use INPUTS;
use dynamics;
use domains;
use tracers;
use params;
use adv_4o_cntr;
use NetCDF_IO;
//use utils;

use LinearAlgebra;
use AllLocalesBarriers;

////////////////////////////////////////////////////////
//                                                    //
//               Runge-Kutta 4th order                //
//                                                    //
////////////////////////////////////////////////////////

proc TimeStep(ref Dyn: Dynamics, D: Domains, ref Trc: Tracers, P: Params, step : int) {


  // Load velocity fields for the next timestep
    update_thickness(Trc.zeta_np1, Trc.H_np1, Trc.H0, Trc.h, D, P, step+1);
    update_dynamics(Dyn.U_np1, Dyn.V_np1, Trc.H_np1, D, P, step+1);

WriteOutput(Trc.tracer, "init", "stuff", step);

  // Set the boundary conditions for the tracer by reading in the
  // values at the boundaries
    set_bry(P, P.bryfiles[step], "temp", Trc.tracer, D.rho_3D);

  // Update halos for current tracer field
    update_halos(Trc.tracer);

  // Calculate zeta and thicknesses at (n-1/2) and (n+1/2) and (n+3/2)
    calc_half_step_tr(Trc, D, P);

  // Calculate horizontal mass fluxes at the (n+1/2) time step.
    calc_half_step_dyn(Dyn);

  // Calculate W at n using the continuity equation
    calc_W(Dyn.U_n, Dyn.V_n, Dyn.W_n, D, P, Trc.H_nm1h, Trc.H_np1h);

  // Calculate W at (n+1/2) using the continuity equation
    calc_W(Dyn.U_np1h, Dyn.V_np1h, Dyn.W_np1h, D, P, Trc.H_n, Trc.H_np1);

  // Calculate W at (n+1) using the continuity equation
    calc_W(Dyn.U_np1, Dyn.V_np1, Dyn.W_np1, D, P, Trc.H_np1h, Trc.H_np3h);

  // Calculate tracer fluxes at the n time step.
    calc_fluxes(Dyn.U_n, Dyn.V_n, Dyn.W_n, Dyn, D, P, Trc.tracer, Trc.H_n);

    RHS(Trc.k1, Dyn, D, P, Trc.H_n);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + 0.5*P.dt*Trc.k1[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(Dyn.U_np1h, Dyn.V_np1h, Dyn.W_np1h, Dyn, D, P, Trc.tmp, Trc.H_np1h);

    RHS(Trc.k2, Dyn, D, P, Trc.H_np1h);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + 0.5*P.dt*Trc.k2[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(Dyn.U_np1h, Dyn.V_np1h, Dyn.W_np1h, Dyn, D, P, Trc.tmp, Trc.H_np1h);

    RHS(Trc.k3, Dyn, D, P, Trc.H_np1h);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + P.dt*Trc.k3[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(Dyn.U_np1, Dyn.V_np1, Dyn.W_np1, Dyn, D, P, Trc.tmp, Trc.H_np1);

    Trc.tracer[D.rho_3D] = (Trc.tracer[D.rho_3D]*Trc.H_n[D.rho_3D] + P.one_sixth*P.dt*
                           (Trc.k1[D.rho_3D] + 2*Trc.k2[D.rho_3D] + 2*Trc.k3[D.rho_3D] + Trc.k4[D.rho_3D])) / Trc.H_np1[D.rho_3D];

 WriteOutput(Trc.tracer, "temp", "stuff", step);

    Trc.H_nm1[D.rho_3D] = Trc.H_n[D.rho_3D];
    Trc.H_n[D.rho_3D] = Trc.H_np1[D.rho_3D];

    Dyn.U_n = Dyn.U_np1;
    Dyn.V_n = Dyn.V_np1;

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

}

