use INPUTS;
use dynamics;
use diffusion;
use domains;
use tracers;
use params;
use hadv_3o_upwind;
use remap;
use NetCDF_IO;
//use utils;

use Math;
use AllLocalesBarriers;
use Time;


////////////////////////////////////////////////////////
//                                                    //
//               Runge-Kutta 3th order                //
//                                                    //
//         y_(n+1)=y_n+(1/6)*h*(k_1+4*k_2+k_3)        //
//                      where                         //
//                  k_1=f(x_n,y_n)                    //
//          k2=f(x_n+(1/2)*h,y_n+(1/2)*h*k_1)         //
//            k3=f(x_n+h,y_n+2*h*k_2−h*k_1)           //
//                                                    //
////////////////////////////////////////////////////////


proc TimeStep(ref Dyn: Dynamics, ref Diff: Diffusion, D: Domains, P: Params, step : int) {

  // Load velocity fields for the next timestep
    update_thickness(zeta_np1, H_np1, H0, h, D, P, step+1);
    update_dynamics(Dyn.u_np1, Dyn.v_np1, Dyn.U_np1, Dyn.V_np1, H_np1, D, P, step+1);

  // Calculate horizontal velocities at the (n+1/2) time step.
    calc_half_step_dyn(Dyn, D, P);

  // Calculate thickness at the (n+1/2) time step. This is only used for the diffusion module.
    calc_half_step_tr(D, P);

  /////////////////////////////////////
  //       Update H using RK3        //
  /////////////////////////////////////

    WriteOutput(H_n, "H", "stuff", step);

    RHS_H(k1, Dyn.U_n, Dyn.V_n, D, P);

    forall (t,k,j,i) in D.rho_3D {
      ktmp[t,k,j,i] =  H_n[t,k,j,i] + 0.5*P.dt*k1[t,k,j,i];
    }
    update_halos(ktmp);

    calc_volumetric_fluxes(Dyn.u_np1h, Dyn.v_np1h, Dyn.U_np1h, Dyn.V_np1h, ktmp, D, P);

    RHS_H(k2, Dyn.U_np1h, Dyn.V_np1h, D, P);

    forall (t,k,j,i) in D.rho_3D {
      ktmp[t,k,j,i] =  H_n[t,k,j,i] - P.dt*k1[t,k,j,i] + 2*P.dt*k2[t,k,j,i];
    }
    update_halos(ktmp);

    calc_volumetric_fluxes(Dyn.u_np1, Dyn.v_np1, Dyn.U_np1, Dyn.V_np1, ktmp, D, P);

    RHS_H(k3, Dyn.U_np1, Dyn.V_np1, D, P);

    forall (t,k,j,i) in D.rho_3D {
      H_dagger[t,k,j,i] = H_n[t,k,j,i] + P.one_sixth*P.dt*
                           (k1[t,k,j,i] + 4*k2[t,k,j,i] + k3[t,k,j,i]);
    }

    allLocalesBarrier.barrier();

  //////////////////////////////////////////
  //       Update tracer using RK3        //
  //////////////////////////////////////////

    calc_horizontal_fluxes(Dyn.U_n, Dyn.V_n, Dyn.tmp_U, Dyn.tmp_V, D, P, tracer_n);
    calc_diffusive_fluxes(Diff.tmp_U, Diff.tmp_V, D, P, tracer_n, H_n);

    RHS_tr(k1, Dyn.tmp_U, Dyn.tmp_V, Diff.tmp_U, Diff.tmp_V, H_n, D, P);

    forall (t,k,j,i) in D.rho_3D {
      ktmp[t,k,j,i] =  tracer_n[t,k,j,i] + 0.5*P.dt*k1[t,k,j,i];
    }
    update_halos(ktmp);

    calc_horizontal_fluxes(Dyn.U_np1h, Dyn.V_np1h, Dyn.tmp_U, Dyn.tmp_V, D, P, ktmp);
    calc_diffusive_fluxes(Diff.tmp_U, Diff.tmp_V, D, P, ktmp, H_np1h);

    RHS_tr(k2, Dyn.tmp_U, Dyn.tmp_V, Diff.tmp_U, Diff.tmp_V, H_np1h, D, P);

    forall (t,k,j,i) in D.rho_3D {
      ktmp[t,k,j,i] =  tracer_n[t,k,j,i] - P.dt*k1[t,k,j,i] + 2*P.dt*k2[t,k,j,i];
    }
    update_halos(ktmp);

    calc_horizontal_fluxes(Dyn.U_np1, Dyn.V_np1, Dyn.tmp_U, Dyn.tmp_V, D, P, ktmp);
    calc_diffusive_fluxes(Diff.tmp_U, Diff.tmp_V, D, P, ktmp, H_np1);

    RHS_tr(k3, Dyn.tmp_U, Dyn.tmp_V, Diff.tmp_U, Diff.tmp_V, H_np1, D, P);

    forall (t,k,j,i) in D.rho_3D {
      tracer_dagger[t,k,j,i] = tracer_n[t,k,j,i] + P.one_sixth*P.dt*
                           (k1[t,k,j,i] + 4*k2[t,k,j,i] + k3[t,k,j,i]);
    }

    allLocalesBarrier.barrier();

}


proc RHS_H(ref tmp, ref U, ref V, D: Domains, P: Params) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]}  {
      tmp[t,k,j,i] = - P.iarea * (  (U[t,k,j,i] - U[t,k,j,i-1])
                        + (V[t,k,j,i] - V[t,k,j-1,i])  );
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)}  {
      tmp[t,k,j,i] = - P.iarea * (  (U[t,k,j,i] - U[t,k,j,i-1])
                        + (V[t,k,j,i] - V[t,k,j-1,i])  );
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]}  {
      tmp[t,k,j,i] = - P.iarea *(  (U[t,k,j,i] - U[t,k,j,i-1])
                        + (V[t,k,j,i] - V[t,k,j-1,i])  );
    }
  }

  allLocalesBarrier.barrier();

}


proc RHS_tr(ref tmp, ref U, ref V, ref diff_U, ref diff_V, ref H, D: Domains, P: Params) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]}  {
      tmp[t,k,j,i] = - P.iarea * (  (U[t,k,j,i] - U[t,k,j,i-1])
                                  + (V[t,k,j,i] - V[t,k,j-1,i])
                                  - (diff_U[t,k,j,i] - diff_U[t,k,j,i-1])
                                  - (diff_V[t,k,j,i] - diff_V[t,k,j-1,i]) ) / H[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)}  {
      tmp[t,k,j,i] = - P.iarea * (  (U[t,k,j,i] - U[t,k,j,i-1])
                                  + (V[t,k,j,i] - V[t,k,j-1,i])
                                  - (diff_U[t,k,j,i] - diff_U[t,k,j,i-1])
                                  - (diff_V[t,k,j,i] - diff_V[t,k,j-1,i]) ) / H[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]}  {
      tmp[t,k,j,i] = - P.iarea * (  (U[t,k,j,i] - U[t,k,j,i-1])
                                  + (V[t,k,j,i] - V[t,k,j-1,i])
                                  - (diff_U[t,k,j,i] - diff_U[t,k,j,i-1])
                                  - (diff_V[t,k,j,i] - diff_V[t,k,j-1,i]) ) / H[t,k,j,i];
    }
  }

  allLocalesBarrier.barrier();

}

