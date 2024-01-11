use INPUTS;
use dynamics;
use domains;
use tracers;
use params;
use hadv_3o_upwind;
use vadv_4o_implicit;
use NetCDF_IO;
//use utils;

use Math;
use LinearAlgebra;
use AllLocalesBarriers;
use Time;

////////////////////////////////////////////////////////
//                                                    //
//               Runge-Kutta 4th order                //
//                                                    //
////////////////////////////////////////////////////////

proc TimeStep(ref Dyn: Dynamics, D: Domains, P: Params, step : int) {

  // Load velocity fields for the next timestep
    update_thickness(zeta_np1, H_np1, H0, h, D, P, step+1);
    update_dynamics(Dyn.U_np1, Dyn.V_np1, H_np1, D, P, step+1);

//  // Predict tracer concentration at (n+1/2) using AB3 extrapolation
  // Calculate zeta and thickness at (n-1/2) and (n+1/2)
    calc_half_step_tr(D, P);

  // Calculate horizontal mass fluxes at the (n+1/2) time step.
    calc_half_step_dyn(Dyn, D, P);

//  // Calculate W at n using the continuity equation
    calc_W(Dyn.U_n, Dyn.V_n, Dyn.W_n, D, P, H_nm1h, H_np1h);

  // Calculate W at (n+1/2) using the continuity equation
    calc_W(Dyn.U_np1h, Dyn.V_np1h, Dyn.W_np1h, D, P, H_n, H_np1);

  /////////////////////////////////////
  //       Update using LF-AM3       //
  /////////////////////////////////////

    calc_horizontal_fluxes(Dyn.U_n, Dyn.V_n, Dyn, D, P, tracer, H_n);
    calc_vertical_flux(Dyn.W_n, Dyn, D, P, tracer, H_n);

    calc_auxiliary_thicknesses(Dyn, D, P);

    calc_predictor(Dyn, D, P);
    set_bry(P, P.bryfiles[step], "temp", tracer_np1h, D.rho_3D);



    calc_horizontal_fluxes(Dyn.U_np1h, Dyn.V_np1h, Dyn, D, P, tracer_np1h, H_np1h);
    calc_vertical_flux(Dyn.W_np1h, Dyn, D, P, tracer_np1h, H_np1h);

    calc_courant(Dyn, D, P);
    WriteOutput(courant, "cour", "stuff", step);

    calc_corrector(Dyn, D, P);

    allLocalesBarrier.barrier();

    forall (t,k,j,i) in D.rho_3D {
//      H_nm1[t,k,j,i] = H_n[t,k,j,i];
      H_n[t,k,j,i] = H_np1[t,k,j,i];
      tracer_nm2[t,k,j,i] = tracer_nm1[t,k,j,i];
      tracer_nm1[t,k,j,i] = tracer[t,k,j,i];
      tracer[t,k,j,i] = corrector[t,k,j,i];
    }
    set_bry(P, P.bryfiles[step+1], "temp", tracer, D.rho_3D);

//    update_halos(H_nm1);
    update_halos(H_n);
    update_halos(tracer);
    update_halos(tracer_nm1);
    update_halos(tracer_nm2);

    Dyn.U_nm2 = Dyn.U_nm1;
    Dyn.V_nm2 = Dyn.V_nm1;
    Dyn.U_nm1 = Dyn.U_n;
    Dyn.V_nm1 = Dyn.V_n;
    Dyn.U_n = Dyn.U_np1;
    Dyn.V_n = Dyn.V_np1;
    allLocalesBarrier.barrier();

    WriteOutput(tracer, "temp", "stuff", step);

}

proc calc_predictor(ref Dyn: Dynamics, D: Domains, P: Params) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Dyn) {
      tracer_np1h[t,k,j,i] = ( ( (0.5 + 2*P.gamma) * tracer[t,k,j,i] + (0.5 - 2*P.gamma) * tracer_nm1[t,k,j,i]) * H_minus[t,k,j,i]
              - (1 - 2*P.gamma) * P.dt * P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                                                + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                                                + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) ) ) / H_plus[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Dyn) {
      tracer_np1h[t,k,j,i] = ( ( (0.5 + 2*P.gamma) * tracer[t,k,j,i] + (0.5 - 2*P.gamma) * tracer_nm1[t,k,j,i]) * H_minus[t,k,j,i]
              - (1 - 2*P.gamma) * P.dt * P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                                                + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                                                + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) ) ) / H_plus[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Dyn) {
      tracer_np1h[t,k,j,i] = ( ( (0.5 + 2*P.gamma) * tracer[t,k,j,i] + (0.5 - 2*P.gamma) * tracer_nm1[t,k,j,i]) * H_minus[t,k,j,i]
              - (1 - 2*P.gamma) * P.dt * P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                                                + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                                                + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) ) ) / H_plus[t,k,j,i];
    }
  }

}

proc calc_corrector(ref Dyn: Dynamics, D: Domains, P: Params) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Dyn) {
      corrector[t,k,j,i] = (tracer[t,k,j,i]*H_n[t,k,j,i]
                         - P.dt*P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                                       + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                                       + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) ) ) / H_np1[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Dyn) {
      corrector[t,k,j,i] = (tracer[t,k,j,i]*H_n[t,k,j,i]
                         - P.dt*P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                                       + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                                       + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) ) ) / H_np1[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Dyn) {
      corrector[t,k,j,i] = (tracer[t,k,j,i]*H_n[t,k,j,i]
                         - P.dt*P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                                       + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                                       + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) ) ) / H_np1[t,k,j,i];
    }
  }

}

proc calc_courant(ref Dyn: Dynamics, D: Domains, P: Params) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Dyn) {
      courant[t,k,j,i] = P.dt * P.iarea *
                                        ((max(Dyn.U_np1h[t,k,j,i],0.0) - min(Dyn.U_np1h[t,k,j,i-1],0.0))
                                       + (max(Dyn.V_np1h[t,k,j,i],0.0) - min(Dyn.V_np1h[t,k,j-1,i],0.0))
                                       + (max(Dyn.W_np1h[t,k,j,i],0.0) - min(Dyn.W_np1h[t,k-1,j,i],0.0))) / H_n[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Dyn) {
      courant[t,k,j,i] = P.dt * P.iarea *
                                        ((max(Dyn.U_np1h[t,k,j,i],0.0) - min(Dyn.U_np1h[t,k,j,i-1],0.0))
                                       + (max(Dyn.V_np1h[t,k,j,i],0.0) - min(Dyn.V_np1h[t,k,j-1,i],0.0))
                                       + (max(Dyn.W_np1h[t,k,j,i],0.0) - min(Dyn.W_np1h[t,k-1,j,i],0.0))) / H_n[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Dyn) {
      courant[t,k,j,i] = P.dt * P.iarea *
                                        ((max(Dyn.U_np1h[t,k,j,i],0.0) - min(Dyn.U_np1h[t,k,j,i-1],0.0))
                                       + (max(Dyn.V_np1h[t,k,j,i],0.0) - min(Dyn.V_np1h[t,k,j-1,i],0.0))
                                       + (max(Dyn.W_np1h[t,k,j,i],0.0) - min(Dyn.W_np1h[t,k-1,j,i],0.0))) / H_n[t,k,j,i];
    }
  }

}
