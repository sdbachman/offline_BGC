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

proc TimeStep(ref Ap2: Arrays, ref Ap: Arrays, ref Ac: Arrays, ref An: Arrays, ref An2: Arrays, D: Domains, ref Trp2: Tracers, ref Trp: Tracers, ref Trc: Tracers, ref Trn: Tracers,ref Trn2: Tracers, F: Files, step : int) {

/*
  // Load velocity fields for the next timestep into field "An"
    update_thickness(Trn, D, F, (step+1));
    update_dynamic_arrays(An, Trn, D, F, (step+1));

WriteOutput(Trc.tracer, "init", "stuff", step);

  // Set the boundary conditions for the tracer by reading in the
  // values at the boundaries
    set_bry(F.bry[step], "temp", Trc.tracer, D.rho_3D);
//    WriteOutput(Trc.tracer, "initial", "stuff", step);

  // Update halos for current tracer field
    update_halos(Trc.tracer);

  // Calculate zeta and thicknesses at n-1/2 and n+1/2
    calc_half_step_tr(Trp2, Trp, Trc, Trn, D);

//
//  // Calculate W at n using the continuity equation
//    calc_W(Ac, D, Trc.thickness_half, Trn.thickness_half);
//
//  // Calculate tracer fluxes based on current velocities
//    calc_fluxes(Ac, Trc, D, Trc.tracer, Trc.thickness);
//
  // Calculate predictor step for LF-AM3. This estimates the tracer
  // states at the (n+1/2) time step.
//    predictor(Ac, D, Trp, Trc);
    predictor(Ac, D, Trp2, Trp, Trc);

  // Set the boundary conditions for the predictor by reading in the
  // tracer values at the boundaries.  Currently only using boundary
  // values at time n, not (n+1/2)
//    set_bry(F.bry[step], "temp", Trc.predictor, D.rho_3D);

WriteOutput(Trc.predictor, "pred", "stuff", step);

  // Update halos for the predictor
 //   update_halos(Trc.predictor);

  // Calculate horizontal mass fluxes at the (n+1/2) time step. Store them in "Ac".
    calc_half_step_dyn(Ap2, Ap, Ac, An);

  // Calculate W at (n+1/2) using the continuity equation
    calc_W(Ac, D, Trc.thickness, Trn.thickness);

  // Calculate tracer fluxes at the (n+1/2) time step.
    calc_fluxes(Ac, Trc, D, Trc.predictor, Trc.thickness_half);

  // Calculate corrector step for LF-AM3. This calculates the tracer
  // states at the (n+1) time step.
    corrector(Ac, D, Trc, Trn);

  // Change the "old" tracer variables to be the "current" tracer variables.
    Trp2.tracer[D.rho_3D] = Trp.tracer[D.rho_3D];
    Trp2.zeta[D.rho_2D] = Trp.zeta[D.rho_2D];
    Trp2.thickness[D.rho_3D] = Trp.thickness[D.rho_3D];

    Trp.tracer[D.rho_3D] = Trc.tracer[D.rho_3D];
    Trp.zeta[D.rho_2D] = Trc.zeta[D.rho_2D];
    Trp.thickness[D.rho_3D] = Trc.thickness[D.rho_3D];

  // Set the "current" tracer variables to be at time (n+1).
    Trc.tracer[D.rho_3D] = Trc.corrector[D.rho_3D];
    Trc.zeta[D.rho_2D] = Trn.zeta[D.rho_2D];
    Trc.thickness[D.rho_3D] = Trn.thickness[D.rho_3D];

  // Set the boundary conditions for tracer_new by reading in the
  // tracer values at the boundaries at time (n+1).
    set_bry(F.bry[step+1], "temp", Trc.tracer, D.rho_3D);

    WriteOutput(Trc.tracer, "temp", "stuff", step);

//allLocalesBarrier.barrier();
//exit();
  // Set "current" dynamical fields to be those from the next timestep
    Ap2.U = Ap.U;
    Ap2.V = Ap2.V;

    Ap.U = Ac.U;
    Ap.V = Ac.V;

    Ac.U = An.U;
    Ac.V = An.V;


//continuity_compare(Ac, D, Trc.thickness, Trn.thickness, Trc);
//WriteOutput(Trc.div, "div", "stuff", step);
//WriteOutput(Trc.dV, "dV", "stuff", step);

*/


  // Load velocity fields for the next timestep into field "An"
    update_thickness(Trn, D, F, (step+1));
    update_thickness(Trn2, D, F, (step+2));
    update_dynamic_arrays(An, Trn, D, F, (step+1));
    update_dynamic_arrays(An2, Trn2, D, F, (step+2));

WriteOutput(Trc.tracer, "init", "stuff", step);

  // Set the boundary conditions for the tracer by reading in the
  // values at the boundaries
    set_bry(F.bry[step], "temp", Trc.tracer, D.rho_3D);

  // Update halos for current tracer field
    update_halos(Trc.tracer);

  // Calculate zeta and thicknesses at (n-1/2) and (n+1/2) and (n+3/2)
    calc_half_step_tr(Trp2, Trp, Trc, Trn, Trn2, D);

  // Calculate W at n using the continuity equation
    calc_W(Ac.U, Ac.V, Ac.W, D, Trc.thickness_half, Trn.thickness_half);

  // Calculate W at (n+1) using the continuity equation
    calc_W(An.U, An.V, An.W, D, Trn.thickness_half, Trn2.thickness_half);

  // Calculate horizontal mass fluxes at the (n+1/2) time step.
    calc_half_step_dyn(Ap2, Ap, Ac, An);

  // Calculate W at (n+1/2) using the continuity equation
    calc_W(Ac.U_half, Ac.V_half, Ac.W_half, D, Trc.thickness, Trn.thickness);

  // Calculate tracer fluxes at the n time step.
    calc_fluxes(Ac.U, Ac.V, Ac.W, Ac, Trc, D, Trc.tracer, Trc.thickness);

    RHS(Trc.k1, Ac, D, Trc.thickness);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + 0.5*dt*Trc.k1[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(Ac.U_half, Ac.V_half, Ac.W_half, Ac, Trc, D, Trc.tmp, Trn.thickness_half);

    RHS(Trc.k2, Ac, D, Trn.thickness_half);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + 0.5*dt*Trc.k2[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(Ac.U_half, Ac.V_half, Ac.W_half, Ac, Trc, D, Trc.tmp, Trn.thickness_half);

    RHS(Trc.k3, Ac, D, Trn.thickness_half);

    Trc.tmp[D.rho_3D] =  Trc.tracer[D.rho_3D] + dt*Trc.k3[D.rho_3D];
    update_halos(Trc.tmp);
    calc_fluxes(An.U, An.V, An.W, Ac, Trc, D, Trc.tmp, Trn.thickness);

    Trc.tracer[D.rho_3D] = (Trc.tracer[D.rho_3D]*Trc.thickness[D.rho_3D] + 0.166666*dt*(Trc.k1[D.rho_3D] + 2*Trc.k2[D.rho_3D] + 2*Trc.k3[D.rho_3D] + Trc.k4[D.rho_3D])) / Trn.thickness[D.rho_3D];

 WriteOutput(Trc.tracer, "temp", "stuff", step);

    Trp.thickness[D.rho_3D] = Trc.thickness[D.rho_3D];
    Trc.thickness[D.rho_3D] = Trn.thickness[D.rho_3D];
    Trn.thickness[D.rho_3D] = Trn2.thickness[D.rho_3D];

    Ap.U = Ac.U;
    Ac.U = An.U;
    An.U = An2.U;

    Ap.V = Ac.V;
    Ac.V = An.V;
    An.V = An2.V;

}


proc calc_fluxes(ref U, ref V, ref W, ref A: Arrays, ref Tr: Tracers, D: Domains, ref arr, ref thickness_c) {

  // Temporary arrays
  //var tmp_U : [D.u_3D] real;
  //var tmp_V : [D.v_3D] real;
  //var tmp_W : [D.w_3D] real;

//  forall (t,k,j,i) in D.rho_3D with (ref A) {
//    var tmp = (0.5 - gamma) * dt * iarea * (A.U[t,k,j,i] - A.U[t,k,j,i-1] + A.V[t,k,j,i] - A.V[t,k,j-1,i] + A.W[t,k,j,i] - A.W[t,k-1,j,i]);
//    Tr.H_plus[t,k,j,i] = thickness_c[t,k,j,i] - tmp;
//    Tr.H_minus[t,k,j,i] = thickness_c[t,k,j,i] + tmp;
//  }

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
    var h_t : real = thickness_c[t,0,j,i];
    var ithickness = 1.0 / (h_t - h_b);

    for k in 0..3 {
      M[k,0] = ithickness*(h_t - h_b);
      M[k,1] = 0.5*ithickness*(h_t**2 - h_b**2);
      M[k,2] = (1.0/3.0)*ithickness*(h_t**3 - h_b**3);
      M[k,3] = 0.25 * ithickness*(h_t**4 - h_b**4);

      B[k] = arr[t,k,j,i];

      h_b = h_b + thickness_c[t,k,j,i];
      h_t = h_t + thickness_c[t,k+1,j,i];
      ithickness = 1.0/(h_t - h_b);
    }

    var coeffs = solve(M,B);
    var Ts_bot = coeffs[0];

    // Top boundary extrapolation
    h_b = 0;
    h_t = thickness_c[t,Nz-1,j,i];
    ithickness = 1.0/(h_t - h_b);
    for k in 0..3 {
      M[k,0] = ithickness*(h_t - h_b);
      M[k,1] = 0.5*ithickness*(h_t**2 - h_b**2);
      M[k,2] = (1.0/3.0)*ithickness*(h_t**3 - h_b**3);
      M[k,3] = 0.25 * ithickness*(h_t**4 - h_b**4);

      B[k] = arr[t,Nz-1-k,j,i];

      h_b = h_b + thickness_c[t,Nz-1-k,j,i];
      h_t = h_t + thickness_c[t,Nz-2-k,j,i];
      ithickness = 1.0/(h_t - h_b);
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
      var h0 = thickness_c[t,k,j,i];
      var h1 = thickness_c[t,k+1,j,i];

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


proc calc_W(ref U, ref V, ref W, D: Domains, ref thickness_c, ref thickness_n) {

  if (here.id == 0) {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} {
//    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (thickness_n[t,k,j,i]-thickness_c[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]-1} {
//    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]-1} with (ref Ac) {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (thickness_n[t,k,j,i]-thickness_c[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }

    }
  }
  else {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]} {
//    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]} with (ref Ac) {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (thickness_n[t,k,j,i]-thickness_c[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }
    }
  }

}


proc calc_half_step_tr(ref Trp2: Tracers, ref Trp: Tracers, ref Trc: Tracers, ref Trn: Tracers, ref Trn2: Tracers, D: Domains) {

  //Trc.zeta_half        = 0.5 * (Trp.zeta + Trc.zeta);
  //Trn.zeta_half        = 0.5 * (Trc.zeta + Trn.zeta);

    forall (t,j,i) in D.rho_2D {
      //Trc.zeta_half[t,j,i] = (1.5 + beta) * Trc.zeta[t,j,i] - (0.5 + 2*beta) * Trp.zeta[t,j,i] + beta * Trp2.zeta[t,j,i];
      Trc.zeta_half[t,j,i] = 0.5 * (Trp.zeta[t,j,i] + Trc.zeta[t,j,i]);
      Trn.zeta_half[t,j,i] = 0.5 * (Trc.zeta[t,j,i] + Trn.zeta[t,j,i]);
      Trn2.zeta_half[t,j,i] = 0.5 * (Trn.zeta[t,j,i] + Trn2.zeta[t,j,i]);
    }

  // From SM09, Eq. 2.13
    forall (t,k,j,i) in D.rho_3D {
      Trc.thickness_half[t,k,j,i] = Trc.thickness0[t,k,j,i] * (1 + Trc.zeta_half[t,j,i] / Trc.h[j,i]);
      Trn.thickness_half[t,k,j,i] = Trn.thickness0[t,k,j,i] * (1 + Trn.zeta_half[t,j,i] / Trn.h[j,i]);
      Trn2.thickness_half[t,k,j,i] = Trn2.thickness0[t,k,j,i] * (1 + Trn2.zeta_half[t,j,i] / Trn2.h[j,i]);
    }

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      Trc.thickness_half.updateFluff();
      Trn.thickness_half.updateFluff();
      Trn2.thickness_half.updateFluff();
    }
    allLocalesBarrier.barrier();
}

proc calc_half_step_dyn(ref Ap2: Arrays, ref Ap: Arrays, ref Ac: Arrays, ref An: Arrays) {

  Ac.U_half        = 0.5 * (Ac.U + An.U);
  Ac.V_half        = 0.5 * (Ac.V + An.V);

  // Ac.U = (1.5 + beta) * Ac.U - (0.5 + 2*beta) * Ap.U + beta * Ap2.U;
  // Ac.V = (1.5 + beta) * Ac.V - (0.5 + 2*beta) * Ap.V + beta * Ap2.V;

}

//proc predictor(ref A: Arrays, D: Domains, Trp: Tracers, Trc: Tracers) {
proc predictor(ref A: Arrays, D: Domains, Trp2: Tracers, Trp: Tracers, Trc: Tracers) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1/2) timestep  //
  /////////////////////////////////////////
/*
  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref A) {
      Trc.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Trc.tracer[t,k,j,i] + (0.5 - 2*gamma) * Trp.tracer[t,k,j,i])*Trc.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k,j,i] - A.tmp_W[t,k-1,j,i]) ) ) / Trc.H_plus[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref A) {
      Trc.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Trc.tracer[t,k,j,i] + (0.5 - 2*gamma) * Trp.tracer[t,k,j,i])*Trc.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k,j,i] - A.tmp_W[t,k-1,j,i]) ) ) / Trc.H_plus[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref A) {
      Trc.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Trc.tracer[t,k,j,i] + (0.5 - 2*gamma) * Trp.tracer[t,k,j,i])*Trc.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k,j,i] - A.tmp_W[t,k-1,j,i]) ) ) / Trc.H_plus[t,k,j,i];
    }
  }
*/
    forall (t,k,j,i) in D.rho_3D with (ref A) {
      Trc.predictor[t,k,j,i] = (1.5 + beta) * Trc.tracer[t,k,j,i] - (0.5 + 2*beta) * Trp.tracer[t,k,j,i] + beta * Trp2.tracer[t,k,j,i];
    }

}


proc corrector(ref Ac: Arrays, D: Domains, ref Trc: Tracers, ref Trn: Tracers) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {
      Trc.corrector[t,k,j,i] = (Trc.tracer[t,k,j,i]*Trc.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) ) ) / Trn.thickness[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Ac) {
      Trc.corrector[t,k,j,i] = (Trc.tracer[t,k,j,i]*Trc.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) ) ) / Trn.thickness[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Ac) {
      Trc.corrector[t,k,j,i] = (Trc.tracer[t,k,j,i]*Trc.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) ) ) / Trn.thickness[t,k,j,i];
    }
  }

}

proc RHS(ref tmp, ref Ac: Arrays, D: Domains, ref thickness) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {
      tmp[t,k,j,i] =                         - iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) )  / thickness[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Ac) {
      tmp[t,k,j,i] =                       - iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) )  / thickness[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Ac) {
      tmp[t,k,j,i] =                         - iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                   + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                   + (Ac.tmp_W[t,k,j,i] - Ac.tmp_W[t,k-1,j,i]) )  / thickness[t,k,j,i];
    }
  }

}

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

