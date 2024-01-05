	use INPUTS;
use arrays;
use domains;
use tracers;
use files;
use NetCDF_IO;

use LinearAlgebra;
use AllLocalesBarriers;

////////////////////////////////////////////////////////
//                                                    //
//          Leapfrog - 3rd order Adams-Moulton        //
//                                                    //
//          Predictor step:  SM05, 4.8                //
//          Corrector step:  SM05, 4.2                //
//                                                    //
//        Piecewise Parabolic Splines:  H07, 63       //
//                                                    //
////////////////////////////////////////////////////////

proc TimeStep(ref Ac: Arrays, D: Domains, ref Trc: Tracers, F: Files, step : int) {


  // Load velocity fields for the next timestep into field "An"
    update_thickness(Trc.zeta_np1, Trc.H_np1, Trc.H0, Trc.h, D, F, step+1);
    update_dynamic_arrays(Ac.U_np1, Ac.V_np1, Trc.H_np1, D, F, step+1);

WriteOutput(Trc.tracer, "init", "stuff", step);

  // Set the boundary conditions for the tracer by reading in the
  // values at the boundaries
    set_bry(F.bry[step], "temp", Trc.tracer, D.rho_3D);

  // Update halos for current tracer field
    update_halos(Trc.tracer);

  // Calculate zeta and thicknesses at (n-1/2) and (n+1/2) and (n+3/2)
    calc_half_step_tr(Trc, D);

  // Calculate W at n using the continuity equation
    calc_W(Ac.U_n, Ac.V_n, Ac.W_n, D, Trc.H_nm1h, Trc.H_np1h);

  // Calculate W at (n+1) using the continuity equation
    calc_W(Ac.U_np1, Ac.V_np1, Ac.W_np1, D, Trc.H_np1h, Trc.H_np3h);

  // Calculate horizontal mass fluxes at the (n+1/2) time step.
    calc_half_step_dyn(Ac);

  // Calculate W at (n+1/2) using the continuity equation
    calc_W(Ac.U_np1h, Ac.V_np1h, Ac.W_np1h, D, Trc.H_n, Trc.H_np1);

  // Calculate tracer fluxes at the n time step.
    calc_fluxes(Ac.U_n, Ac.V_n, Ac.W_n, Ac, D, Trc.tracer, Trc.H_n);

    RHS(Trc.k1, Ac, D, Trc.H_n);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + 0.5*dt*Trc.k1[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(Ac.U_np1h, Ac.V_np1h, Ac.W_np1h, Ac, D, Trc.tmp, Trc.H_np1h);

    RHS(Trc.k2, Ac, D, Trc.H_np1h);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + 0.5*dt*Trc.k2[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(Ac.U_np1h, Ac.V_np1h, Ac.W_np1h, Ac, D, Trc.tmp, Trc.H_np1h);

    RHS(Trc.k3, Ac, D, Trc.H_np1h);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + dt*Trc.k3[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(Ac.U_np1, Ac.V_np1, Ac.W_np1, Ac, D, Trc.tmp, Trc.H_np1);

    Trc.tracer[D.rho_3D] = (Trc.tracer[D.rho_3D]*Trc.H_n[D.rho_3D] + 0.166666*dt*(Trc.k1[D.rho_3D] + 2*Trc.k2[D.rho_3D] + 2*Trc.k3[D.rho_3D] + Trc.k4[D.rho_3D])) / Trc.H_np1[D.rho_3D];

 WriteOutput(Trc.tracer, "temp", "stuff", step);

    Trc.H_nm1[D.rho_3D] = Trc.H_n[D.rho_3D];
    Trc.H_n[D.rho_3D] = Trc.H_np1[D.rho_3D];

    Ac.U_n = Ac.U_np1;
    Ac.V_n = Ac.V_np1;

}


proc calc_fluxes(ref U, ref V, ref W, ref A: Arrays, D: Domains, ref arr, ref H) {

  // Horizontal: upstream-biased parabolic interpolation: SM05, after 4.13

  /////////////////////////////////////////
  //              U-fluxes               //
  /////////////////////////////////////////

  if (here.id == 0) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..D.u_3D.first[3]} with (ref A) {
       A.tmp_U[t,k,j,i] = ((5*arr[t,k,j,i] + 8*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], (D.u_3D.first[3]+1)..D.u_3D.last[3]} with (ref A) {
        A.tmp_U[t,k,j,i] = ((-arr[t,k,j,i-1] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * U[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.last[3]..D.u_3D.last[3]} with (ref A) {
      A.tmp_U[t,k,j,i] = ((5*arr[t,k,j,i+1] + 8*arr[t,k,j,i] - arr[t,k,j,i-1])/12) * U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..(D.u_3D.last[3]-1)} with (ref A) {
        A.tmp_U[t,k,j,i] = ((-arr[t,k,j,i-1] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * U[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in D.u_3D with (ref A) {
      A.tmp_U[t,k,j,i] = ((-arr[t,k,j,i-1] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * U[t,k,j,i];
    }
  }

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.first[2]..D.v_3D.first[2], D.v_3D.dim[3]} with (ref A) {
    A.tmp_V[t,k,j,i] = ((5*arr[t,k,j,i] + 8*arr[t,k,j+1,i] - arr[t,k,j+2,i])/12) * V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.last[2]..D.v_3D.last[2], D.v_3D.dim[3]} with (ref A) {
    A.tmp_V[t,k,j,i] = ((5*arr[t,k,j+1,i] + 8*arr[t,k,j,i] - arr[t,k,j-1,i])/12) * V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], (D.v_3D.first[2]+1)..(D.v_3D.last[2]-1), D.v_3D.dim[3]} with (ref A) {
      A.tmp_V[t,k,j,i] = ((-arr[t,k,j-1,i] + 7*arr[t,k,j,i] + 7*arr[t,k,j+1,i] - arr[t,k,j+2,i])/12) * V[t,k,j,i];
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //                                     W-fluxes                                     //
  //  Calculated with Implicit Fourth-order scheme:  White and Adcroft, 2008, Eq. 46  //
  //////////////////////////////////////////////////////////////////////////////////////

  var Dp : domain(1) = {0..Nz};
  var DpDp : domain(2) = {0..Nz, 0..Nz};

  forall (t,j,i) in {0..0,D.rho_3D.dim[2], D.rho_3D.dim[3]} with (ref A) {

    // Bottom boundary extrapolation
    // Calculated using a cubic polynomial over the bottom four cells

    var B : [0..3] real;
    var M : [0..3,0..3] real;

    var h_b : real = 0;
    var h_t : real = H[t,0,j,i];
    var iH = 1.0 / (h_t - h_b);

    for k in 0..3 {
      M[k,0] = iH*(h_t - h_b);
      M[k,1] = 0.5*iH*(h_t**2 - h_b**2);
      M[k,2] = (1.0/3.0)*iH*(h_t**3 - h_b**3);
      M[k,3] = 0.25 * iH*(h_t**4 - h_b**4);

      B[k] = arr[t,k,j,i];

      h_b = h_b + H[t,k,j,i];
      h_t = h_t + H[t,k+1,j,i];
      iH = 1.0/(h_t - h_b);
    }

    var coeffs = solve(M,B);
    var Ts_bot = coeffs[0];

    // Top boundary extrapolation
    h_b = 0;
    h_t = H[t,Nz-1,j,i];
    iH = 1.0/(h_t - h_b);
    for k in 0..3 {
      M[k,0] = iH*(h_t - h_b);
      M[k,1] = 0.5*iH*(h_t**2 - h_b**2);
      M[k,2] = (1.0/3.0)*iH*(h_t**3 - h_b**3);
      M[k,3] = 0.25 * iH*(h_t**4 - h_b**4);

      B[k] = arr[t,Nz-1-k,j,i];

      h_b = h_b + H[t,Nz-1-k,j,i];
      h_t = h_t + H[t,Nz-2-k,j,i];
      iH = 1.0/(h_t - h_b);
    }

    coeffs = solve(M,B) ;
    var Ts_top = coeffs[0];

    /////////////////////////////

    var Bf : [Dp] real;
    var Mf : [DpDp] real;

    Mf[0,0] = 1.0;
    Mf[Nz,Nz] = 1.0;
    Bf[0] = Ts_bot;
    Bf[Nz] = Ts_top;

    for k in 0..(Nz-2) {
      var h0 = H[t,k,j,i];
      var h1 = H[t,k+1,j,i];

      var alpha = (h1**2) / ((h0 + h1)**2);
      var beta = (h0**2) / ((h0 + h1)**2);
      var b = 2*(h1**2)*(h1**2 + 2*h0**2 + 3*h0*h1) / ((h0+h1)**4);
      var c = 2*(h0**2)*(h0**2 + 2*h1**2 + 3*h0*h1) / ((h0+h1)**4);

      Mf[k+1,k] = alpha;
      Mf[k+1,k+1] = 1.0;
      Mf[k+1,k+2] = beta;

      Bf[k+1] = b*arr[t,k,j,i] + c*arr[t,k+1,j,i];
    }

    A.tmp_W[t,..,j,i] = solve(Mf,Bf) * W[t,..,j,i];

  }

}


proc calc_W(ref U, ref V, ref W, D: Domains, ref H_n, ref H_np1) {

  if (here.id == 0) {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} {
//    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (H_np1[t,k,j,i]-H_n[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]-1} {
//    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]-1} with (ref Ac) {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (H_np1[t,k,j,i]-H_n[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }

    }
  }
  else {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]} {
//    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]} with (ref Ac) {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (H_np1[t,k,j,i]-H_n[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }
    }
  }

}


proc calc_half_step_tr(ref Trc: Tracers, D: Domains) {

    forall (t,j,i) in D.rho_2D {
      Trc.zeta_nm1h[t,j,i] = 0.5 * (Trc.zeta_nm1[t,j,i] + Trc.zeta_n[t,j,i]);
      Trc.zeta_np1h[t,j,i] = 0.5 * (Trc.zeta_n[t,j,i] + Trc.zeta_np1[t,j,i]);
      Trc.zeta_np3h[t,j,i] = (1.5 + beta) * Trc.zeta_np1[t,j,i] - (0.5 + 2*beta) * Trc.zeta_n[t,j,i] + beta * Trc.zeta_nm1[t,j,i];
    }

  // From SM09, Eq. 2.13
    forall (t,k,j,i) in D.rho_3D {
      Trc.H_nm1h[t,k,j,i] = Trc.H0[t,k,j,i] * (1 + Trc.zeta_nm1h[t,j,i] / Trc.h[j,i]);
      Trc.H_np1h[t,k,j,i] = Trc.H0[t,k,j,i] * (1 + Trc.zeta_np1h[t,j,i] / Trc.h[j,i]);
      Trc.H_np3h[t,k,j,i] = Trc.H0[t,k,j,i] * (1 + Trc.zeta_np3h[t,j,i] / Trc.h[j,i]);
    }

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      Trc.H_nm1h.updateFluff();
      Trc.H_np1h.updateFluff();
      Trc.H_np3h.updateFluff();
    }
    allLocalesBarrier.barrier();
}

proc calc_half_step_dyn(ref Ac: Arrays) {

  Ac.U_np1h        = 0.5 * (Ac.U_n + Ac.U_np1);
  Ac.V_np1h        = 0.5 * (Ac.V_n + Ac.V_np1);

}

proc RHS(ref tmp, ref Ac: Arrays, D: Domains, ref H) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {
      tmp[t,k,j,i] =                         - iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) )  / H[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Ac) {
      tmp[t,k,j,i] =                       - iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) )  / H[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Ac) {
      tmp[t,k,j,i] =                         - iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) )  / H[t,k,j,i];
    }
  }

}

/*
proc continuity_compare(ref Ac: Arrays, D: Domains, ref thickness_c, ref thickness_n, Trc: Tracers) {

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {
      Trc.div[t,k,j,i] =  -dt*((Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1])
                                                 + (Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i])
                                                 + (Ac.W[t,k,j,i] - Ac.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (thickness_n[t,k,j,i] - thickness_c[t,k,j,i]) * dx * dy;
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Ac) {
      Trc.div[t,k,j,i] =  -dt*((Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1])
                                                 + (Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i])
                                                 + (Ac.W[t,k,j,i] - Ac.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (thickness_n[t,k,j,i] - thickness_c[t,k,j,i]) * dx * dy;
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Ac) {

      Trc.div[t,k,j,i] = -dt*((Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1])
                                                 + (Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i])
                                                 + (Ac.W[t,k,j,i] - Ac.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (thickness_n[t,k,j,i] - thickness_c[t,k,j,i]) * dx * dy;
    }
  }
}
*/
