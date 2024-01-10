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

  // Set the boundary conditions for the tracer by reading in the
  // values at the boundaries
    set_bry(P, P.bryfiles[step], "temp", tracer, D.rho_3D);

  // Update halos for current tracer field
    update_halos(tracer);

  // Calculate zeta and thicknesses at (n-1/2) and (n+1/2) and (n+3/2)
    calc_half_step_tr(D, P);

  // Calculate horizontal mass fluxes at the (n+1/2) time step.
    calc_half_step_dyn(Dyn, D);

  // Calculate W at n using the continuity equation
    calc_W(Dyn.U_n, Dyn.V_n, Dyn.W_n, D, P, H_nm1h, H_np1h);

  // Calculate W at (n+1/2) using the continuity equation
    calc_W(Dyn.U_np1h, Dyn.V_np1h, Dyn.W_np1h, D, P, H_n, H_np1);

  // Calculate W at (n+1) using the continuity equation
    calc_W(Dyn.U_np1, Dyn.V_np1, Dyn.W_np1, D, P, H_np1h, H_np3h);

  //////////////////////////////////
  //       Begin RK4 steps        //
  //////////////////////////////////

  // Calculate tracer fluxes at the n time step.
    calc_horizontal_fluxes(Dyn.U_n, Dyn.V_n, Dyn, D, P, tracer, H_n);
    calc_vertical_flux(Dyn.W_n, Dyn, D, P, tracer, H_n);

    RHS(k1, Dyn, D, P, H_n);

    forall (t,k,j,i) in D.rho_3D {
      ktmp[t,k,j,i] =  tracer[t,k,j,i] + 0.5*P.dt*k1[t,k,j,i];
    }
    update_halos(ktmp);

    calc_horizontal_fluxes(Dyn.U_np1h, Dyn.V_np1h, Dyn, D, P, ktmp, H_np1h);
    calc_vertical_flux(Dyn.W_np1h, Dyn, D, P, ktmp, H_np1h);

    RHS(k2, Dyn, D, P, H_np1h);

    forall (t,k,j,i) in D.rho_3D {
      ktmp[t,k,j,i] =  tracer[t,k,j,i] + 0.5*P.dt*k2[t,k,j,i];
    }
    update_halos(ktmp);

    calc_horizontal_fluxes(Dyn.U_np1h, Dyn.V_np1h, Dyn, D, P, ktmp, H_np1h);
    calc_vertical_flux(Dyn.W_np1h, Dyn, D, P, ktmp, H_np1h);

    RHS(k3, Dyn, D, P, H_np1h);

    forall (t,k,j,i) in D.rho_3D {
      ktmp[t,k,j,i] =  tracer[t,k,j,i] + P.dt*k3[t,k,j,i];
    }
    update_halos(ktmp);

    calc_horizontal_fluxes(Dyn.U_np1, Dyn.V_np1, Dyn, D, P, ktmp, H_np1);
    calc_vertical_flux(Dyn.W_np1, Dyn, D, P, ktmp, H_np1);

    RHS(k4, Dyn, D, P, H_np1);

    forall (t,k,j,i) in D.rho_3D {
      tracer[t,k,j,i] = (tracer[t,k,j,i]*H_n[t,k,j,i] + P.one_sixth*P.dt*
                           (k1[t,k,j,i] + 2*k2[t,k,j,i] + 2*k3[t,k,j,i] + k4[t,k,j,i])) / H_np1[t,k,j,i];
    }
    allLocalesBarrier.barrier();

    forall (t,k,j,i) in D.rho_3D {
      H_nm1[t,k,j,i] = H_n[t,k,j,i];
      H_n[t,k,j,i] = H_np1[t,k,j,i];
    }
    update_halos(H_nm1);
    update_halos(H_n);

    Dyn.U_n = Dyn.U_np1;
    Dyn.V_n = Dyn.V_np1;
    allLocalesBarrier.barrier();

    WriteOutput(tracer, "temp", "stuff", step);

}


proc RHS(ref tmp, ref Dyn: Dynamics, D: Domains, P: Params, ref H) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Dyn) {
      tmp[t,k,j,i] = - P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                            + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                            + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) )  / H[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Dyn) {
      tmp[t,k,j,i] = - P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                            + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                            + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) )  / H[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Dyn) {
      tmp[t,k,j,i] = - P.iarea*((Dyn.tmp_U[t,k,j,i] - Dyn.tmp_U[t,k,j,i-1])
                            + (Dyn.tmp_V[t,k,j,i] - Dyn.tmp_V[t,k,j-1,i])
                            + (Dyn.tmp_W[t,k,j,i] - Dyn.tmp_W[t,k-1,j,i]) )  / H[t,k,j,i];
    }
  }

  allLocalesBarrier.barrier();

}

