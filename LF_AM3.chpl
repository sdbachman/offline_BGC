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

proc TimeStep(ref Ac: Arrays, ref An: Arrays, D: Domains, ref Trp: Tracers, ref Trc: Tracers, ref Trn: Tracers, F: Files, step : int) {

  // Load velocity fields for the next timestep into field "An"
    update_thickness(Trn, D, F, (step+1));
    update_dynamic_arrays(An, Trn, D, F, (step+1));

  // Set the boundary conditions for the tracer by reading in the
  // values at the boundaries
    set_bry(F.bry[step], "temp", Trc.tracer, D.rho_3D);
//WriteOutput(Trc.tracer, "tr", "stuff", step);

  // Update halos for current tracer field
    update_halos(Trc.tracer);
//WriteOutput(Trc.tracer, "b4p", "stuff", step);

  // Calculate zeta and thicknesses at n-1/2 and n+1/2
    calc_half_step_tr(Trp, Trc, Trn, D);

//WriteOutput(Trp.thickness, "thp", "stuff", step);
//WriteOutput(Trc.thickness, "thc", "stuff", step);
//WriteOutput(Trn.thickness, "thn", "stuff", step);
//WriteOutput(Trc.thickness_half, "thc", "stuff", step);
//WriteOutput(Trn.thickness_half, "thn", "stuff", step);

//    allLocalesBarrier.barrier();
//    exit();

  // Calculate W at n using the continuity equation
    calc_W(Ac, D, Trc.thickness_half, Trn.thickness_half);

  // Calculate tracer fluxes based on current velocities
    calc_fluxes(Ac, Trc, D, Trc.tracer);

  // Calculate predictor step for LF-AM3. This estimates the tracer
  // states at the (n+1/2) time step.
    predictor(Ac, D, Trp, Trc);
//WriteOutput(Trc.predictor, "tracer", "stuff", step);
//allLocalesBarrier.barrier();
//exit();

  // Set the boundary conditions for the predictor by reading in the
  // tracer values at the boundaries.  Currently only using boundary
  // values at time n, not (n+1/2)
    set_bry(F.bry[step], "temp", Trc.predictor, D.rho_3D);
//WriteOutput(Trc.predictor, "pred", "stuff", step);

  // Update halos for the predictor
    update_halos(Trc.predictor);

  // Calculate horizontal mass fluxes at the (n+1/2) time step. Store them in "Ac".
    calc_half_step_dyn(Ac, An);

  // Calculate W at (n+1/2) using the continuity equation
    calc_W(Ac, D, Trc.thickness, Trn.thickness);

  // Calculate tracer fluxes at the (n+1/2) time step.
    calc_fluxes(Ac, Trc, D, Trc.predictor);


//// TESTING
     continuity_compare(Ac, D, Trc, Trn);
WriteOutput(Trc.div, "div", "stuff", step);
//WriteOutput(Trc.div2, "div2", "stuff", step);
WriteOutput(Trc.dV, "dV", "stuff", step);
//WriteOutput(Trc.thickness, "trc", "stuff", step);
//WriteOutput(Trn.thickness, "trn", "stuff", step);
exit();

  // Calculate corrector step for LF-AM3. This calculates the tracer
  // states at the (n+1) time step.
    corrector(Ac, D, Trc, Trn);
WriteOutput(Trc.corrector, "corr", "stuff", step);

    allLocalesBarrier.barrier();
    exit();

  // Change the "old" tracer variables to be the "current" tracer variables.
    Trp.tracer[D.rho_3D.localSubdomain()] = Trc.tracer[D.rho_3D.localSubdomain()];
    Trp.zeta[D.rho_2D.localSubdomain()] = Trc.zeta[D.rho_2D.localSubdomain()];
    Trp.thickness[D.rho_3D.localSubdomain()] = Trc.thickness[D.rho_3D.localSubdomain()];

  // Set the "current" tracer variables to be at time (n+1).
    Trc.tracer[D.rho_3D.localSubdomain()] = Trc.corrector[D.rho_3D.localSubdomain()];
    Trc.zeta[D.rho_2D.localSubdomain()] = Trn.zeta[D.rho_2D.localSubdomain()];
    Trc.thickness[D.rho_3D.localSubdomain()] = Trn.thickness[D.rho_3D.localSubdomain()];

//WriteOutput(Tr.tracer_old, "tracer_old", "stuff", step);

  // Set the boundary conditions for tracer_new by reading in the
  // tracer values at the boundaries at time (n+1).
    set_bry(F.bry[step+1], "temp", Trc.tracer, D.rho_3D);

//WriteOutput(Tr.tracer_new, "tracer", "stuff", step);

  // Set "current" dynamical fields to be those from the next timestep
    Ac = An;

}


proc calc_fluxes(ref A: Arrays, ref Tr: Tracers, D: Domains, ref arr) {

  // Temporary arrays
  //var tmp_U : [D.u_3D] real;
  //var tmp_V : [D.v_3D] real;
  //var tmp_W : [D.w_3D] real;

  forall (t,k,j,i) in D.rho_3D with (ref A) {
    var tmp = (0.5 - gamma) * dt * iarea * (A.U[t,k,j,i] - A.U[t,k,j,i-1] + A.V[t,k,j,i] - A.V[t,k,j-1,i] + A.W[t,k+1,j,i] - A.W[t,k,j,i]);
    Tr.H_plus[t,k,j,i] = Tr.thickness[t,k,j,i] - tmp;
    Tr.H_minus[t,k,j,i] = Tr.thickness[t,k,j,i] + tmp;
  }

  // Horizontal: upstream-biased parabolic interpolation: SM05, after 4.13

  /////////////////////////////////////////
  //              U-fluxes               //
  /////////////////////////////////////////

  if (here.id == 0) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..D.u_3D.first[3]} with (ref A) {
//        A.tmp_U[t,k,j,i] = 0.5*(arr[t,k,j,i+1] + arr[t,k,j,i]) * A.U[t,k,j,i];
       A.tmp_U[t,k,j,i] = ((5*arr[t,k,j,i] + 8*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * A.U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], (D.u_3D.first[3]+1)..D.u_3D.last[3]} with (ref A) {
        var tmp : real;
       // if (A.U[t,k,j,i] > 0) {
       //   tmp = arr[t,k,j,i+1] - 2*arr[t,k,j,i] + arr[t,k,j,i-1];
       // }
       // else {
       //   tmp = arr[t,k,j,i+2] - 2*arr[t,k,j,i+1] + arr[t,k,j,i];
       // }
       // A.tmp_U[t,k,j,i] = (0.5*(arr[t,k,j,i+1] + arr[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
        A.tmp_U[t,k,j,i] = ((-arr[t,k,j,i] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * A.U[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.last[3]..D.u_3D.last[3]} with (ref A) {
      //A.tmp_U[t,k,j,i] = 0.5*(arr[t,k,j,i+1] + arr[t,k,j,i]) * A.U[t,k,j,i];
      A.tmp_U[t,k,j,i] = ((5*arr[t,k,j,i+1] + 8*arr[t,k,j,i] - arr[t,k,j,i-1])/12) * A.U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..(D.u_3D.last[3]-1)} with (ref A) {
        var tmp : real;
       // if (A.U[t,k,j,i] > 0) {
       //   tmp = arr[t,k,j,i+1] - 2*arr[t,k,j,i] + arr[t,k,j,i-1];
       // }
       // else {
       //   tmp = arr[t,k,j,i+2] - 2*arr[t,k,j,i+1] + arr[t,k,j,i];
       // }
       // A.tmp_U[t,k,j,i] = (0.5*(arr[t,k,j,i+1] + arr[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
        A.tmp_U[t,k,j,i] = ((-arr[t,k,j,i] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * A.U[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in D.u_3D with (ref A) {
        var tmp : real;
       // if (A.U[t,k,j,i] > 0) {
       //   tmp = arr[t,k,j,i+1] - 2*arr[t,k,j,i] + arr[t,k,j,i-1];
       // }
       // else {
       //   tmp = arr[t,k,j,i+2] - 2*arr[t,k,j,i+1] + arr[t,k,j,i];
       // }
       // A.tmp_U[t,k,j,i] = (0.5*(arr[t,k,j,i+1] + arr[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
      A.tmp_U[t,k,j,i] = ((-arr[t,k,j,i] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * A.U[t,k,j,i];
    }
  }

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.first[2]..D.v_3D.first[2], D.v_3D.dim[3]} with (ref A) {
    A.tmp_V[t,k,j,i] = 0.5*(arr[t,k,j+1,i] + arr[t,k,j,i]) * A.V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.last[2]..D.v_3D.last[2], D.v_3D.dim[3]} with (ref A) {
    A.tmp_V[t,k,j,i] = 0.5*(arr[t,k,j+1,i] + arr[t,k,j,i]) * A.V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], (D.v_3D.first[2]+1)..(D.v_3D.last[2]-1), D.v_3D.dim[3]} with (ref A) {
      var tmp : real;
      if (A.V[t,k,j,i] > 0) {
        tmp = arr[t,k,j+1,i] - 2*arr[t,k,j,i] + arr[t,k,j-1,i];
      }
      else {
        tmp = arr[t,k,j+2,i] - 2*arr[t,k,j+1,i] + arr[t,k,j,i];
      }
      A.tmp_V[t,k,j,i] = (0.5*(arr[t,k,j+1,i] + arr[t,k,j,i]) - us*tmp) * A.V[t,k,j,i];
  }

  /////////////////////////////////////////
  //              W-fluxes               //
  //  Calculated with Splines:  H07, 63  //
  //   Use upwind advection instead of splines, as recommended
  //  at the bottom of
  // https://www.myroms.org/projects/src/ticket/839
  // Centered fourth-order scheme just below Eq. 23:
  // https://www.myroms.org/wiki/Numerical_Solution_Technique
  // or in SM05, 4.10
  /////////////////////////////////////////

  var Dp : domain(1) = {0..Nz};
  var DpDp : domain(2) = {0..Nz, 0..Nz};

  forall (t,j,i) in {0..0,D.rho_3D.dim[2], D.rho_3D.dim[3]} with (ref A) {
    var b : [Dp] real;
    var M : [DpDp] real;

    M[0,0] = 1.0;
    M[Nz,Nz] = 1.0;

    for k in 1..<Nz {
      M[k,k] = 2*Tr.thickness[t,k,j,i] + Tr.thickness[t,k-1,j,i];
      M[k,k-1] = 2*Tr.thickness[t,k,j,i];
      M[k,k+1] = Tr.thickness[t,k-1,j,i];

      b[k] = 3*(Tr.thickness[t,k,j,i]*arr[t,k-1,j,i] + Tr.thickness[t,k-1,j,i]*arr[t,k,j,i]);
    }

    A.tmp_W[t,..,j,i] = solve(M, b) * A.W[t,..,j,i];

  }

}


proc calc_W(ref Ac: Arrays, D: Domains, ref thickness_c, ref thickness_n) {

  if (here.id == 0) {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {

      Ac.W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        Ac.W[t,k,j,i] = Ac.W[t,k-1,j,i] - (thickness_n[t,k,j,i]-thickness_c[t,k,j,i])*area/dt -
                   (Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1] + Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i]);
      }
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]-1} with (ref Ac) {

      Ac.W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        Ac.W[t,k,j,i] = Ac.W[t,k-1,j,i] - (thickness_n[t,k,j,i]-thickness_c[t,k,j,i])*area/dt -
                   (Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1] + Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i]);
      }

    }
  }
  else {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]} with (ref Ac) {

      Ac.W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        Ac.W[t,k,j,i] = Ac.W[t,k-1,j,i] - (thickness_n[t,k,j,i]-thickness_c[t,k,j,i])*area/dt -
                   (Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1] + Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i]);
      }
    }
  }

}


proc calc_half_step_tr(ref Trp: Tracers, ref Trc: Tracers, ref Trn: Tracers, D: Domains) {

  Trc.zeta_half        = 0.5 * (Trp.zeta + Trc.zeta);
  Trn.zeta_half        = 0.5 * (Trc.zeta + Trn.zeta);

  // From SM09, Eq. 2.13
    forall (t,k,j,i) in D.rho_3D {
      Trc.thickness_half[t,k,j,i] = Trc.thickness0[t,k,j,i] * (1 + Trc.zeta_half[t,j,i] / Trc.h[j,i]);
      Trn.thickness_half[t,k,j,i] = Trn.thickness0[t,k,j,i] * (1 + Trn.zeta_half[t,j,i] / Trn.h[j,i]);
    }

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      Trc.thickness_half.updateFluff();
      Trn.thickness_half.updateFluff();
    }
    allLocalesBarrier.barrier();
}

proc calc_half_step_dyn(ref Ac: Arrays, ref An: Arrays) {

  //Ac.H_plus  = 0.5 * (Ac.H_plus + An.H_plus);
  //Ac.H_minus = 0.5 * (Ac.H_minus + Ac.H_minus);
  Ac.U        = 0.5 * (Ac.U + An.U);
  Ac.V        = 0.5 * (Ac.V + An.V);
  //Ac.W        = 0.5 * (Ac.W + An.W);

}

proc predictor(ref A: Arrays, D: Domains, Trp: Tracers, Trc: Tracers) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1/2) timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref A) {
      Trc.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Trc.tracer[t,k,j,i] + (0.5 - 2*gamma) * Trp.tracer[t,k,j,i])*Trc.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k+1,j,i] - A.tmp_W[t,k,j,i]) ) ) / Trc.H_plus[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref A) {
      Trc.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Trc.tracer[t,k,j,i] + (0.5 - 2*gamma) * Trp.tracer[t,k,j,i])*Trc.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k+1,j,i] - A.tmp_W[t,k,j,i]) ) ) / Trc.H_plus[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref A) {
      Trc.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Trc.tracer[t,k,j,i] + (0.5 - 2*gamma) * Trp.tracer[t,k,j,i])*Trc.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k+1,j,i] - A.tmp_W[t,k,j,i]) ) ) / Trc.H_plus[t,k,j,i];
    }
  }

}


proc corrector(ref Ac: Arrays, D: Domains, ref Trc: Tracers, ref Trn: Tracers) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////
/*
  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {
      Tr.corrector[t,k,j,i] = (Tr.tracer_old[t,k,j,i]*Tr.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                                 + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                                 + (Ac.tmp_W[t,k+1,j,i] - Ac.tmp_W[t,k,j,i]) ) ) / An.thickness[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Ac) {
      Tr.corrector[t,k,j,i] = (Tr.tracer_old[t,k,j,i]*Ac.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                                 + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                                 + (Ac.tmp_W[t,k+1,j,i] - Ac.tmp_W[t,k,j,i]) ) ) / An.thickness[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Ac) {
      Tr.corrector[t,k,j,i] = (Tr.tracer_old[t,k,j,i]*Ac.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                                 + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                                 + (Ac.tmp_W[t,k+1,j,i] - Ac.tmp_W[t,k,j,i]) ) ) / An.thickness[t,k,j,i];
    }
  }
*/

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {
      Trc.corrector[t,k,j,i] = (Trc.tracer[t,k,j,i]*Trc.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                                 + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                                 + (Ac.tmp_W[t,k+1,j,i] - Ac.tmp_W[t,k,j,i]) ) ) / Trn.thickness[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Ac) {
      Trc.corrector[t,k,j,i] = (Trc.tracer[t,k,j,i]*Trc.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                                 + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                                 + (Ac.tmp_W[t,k+1,j,i] - Ac.tmp_W[t,k,j,i]) ) ) / Trn.thickness[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Ac) {
      Trc.corrector[t,k,j,i] = (Trc.tracer[t,k,j,i]*Trc.thickness[t,k,j,i]
                         - dt*iarea*((Ac.tmp_U[t,k,j,i] - Ac.tmp_U[t,k,j,i-1])
                                                 + (Ac.tmp_V[t,k,j,i] - Ac.tmp_V[t,k,j-1,i])
                                                 + (Ac.tmp_W[t,k+1,j,i] - Ac.tmp_W[t,k,j,i]) ) ) / Trn.thickness[t,k,j,i];
    }
  }

}


proc continuity_compare(ref Ac: Arrays, D: Domains, ref Trc: Tracers, ref Trn: Tracers) {

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Ac) {
      Trc.div[t,k,j,i] =  - dt*((Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1])
                                                 + (Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i])
                                                 + (Ac.W[t,k,j,i] - Ac.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (Trn.thickness[t,k,j,i] - Trc.thickness[t,k,j,i]) * dx * dy;
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Ac) {
      Trc.div[t,k,j,i] =  - dt*((Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1])
                                                 + (Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i])
                                                 + (Ac.W[t,k,j,i] - Ac.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (Trn.thickness[t,k,j,i] - Trc.thickness[t,k,j,i]) * dx * dy;
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Ac) {

      Trc.div[t,k,j,i] = - dt*((Ac.U[t,k,j,i] - Ac.U[t,k,j,i-1])
                                                 + (Ac.V[t,k,j,i] - Ac.V[t,k,j-1,i])
                                                 + (Ac.W[t,k,j,i] - Ac.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (Trn.thickness[t,k,j,i] - Trc.thickness[t,k,j,i]) * dx * dy;
    }
  }
}

