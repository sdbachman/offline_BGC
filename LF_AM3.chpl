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

proc TimeStep(ref Ac: Arrays, ref An: Arrays, D: Domains, Tr: Tracers, F: Files, step : int) {

  // Load velocity fields for the current timestep
    if (step == Nt_start) {
      update_dynamic_arrays(Ac, D, F, step);
    }
    else {
      Ac = An;
    }

  // Load velocity fields for the next timestep
    update_dynamic_arrays(An, D, F, step);

////////////
  // Update boundary conditions for the tracers
////////////

  calc_fluxes(Ac, D, Tr);

  predictor(Ac, D, Tr);

  get_bry(F.bry[step], "temp", Tr.predictor, D.rho_3D);

  WriteOutput("predictory.nc", Tr.predictor, "tracer", "stuff", 10);

  exit();

  calc_half_step(Ac, An);

  calc_fluxes(Ac, D, Tr);

////////////
  // corrector
////////////

}


proc calc_fluxes(ref A: Arrays, D: Domains, Tr: Tracers) {

  // Temporary arrays
  //var tmp_U : [D.u_3D] real;
  //var tmp_V : [D.v_3D] real;
  //var tmp_W : [D.w_3D] real;

  forall (t,k,j,i) in D.rho_3D with (ref A) {
    var tmp = (0.5 - gamma) * dt * iarea * (A.U[t,k,j,i] - A.U[t,k,j,i-1] + A.V[t,k,j,i] - A.V[t,k,j-1,i] + A.W[t,k+1,j,i] - A.W[t,k,j,i]);
    A.H_plus[t,k,j,i] = A.thickness[t,k,j,i] - tmp;
    A.H_minus[t,k,j,i] = A.thickness[t,k,j,i] + tmp;
  }

  // Horizontal: upstream-biased parabolic interpolation: SM05, after 4.13

  /////////////////////////////////////////
  //              U-fluxes               //
  /////////////////////////////////////////

  if (here.id == 0) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..D.u_3D.first[3]} with (ref A) {
        A.tmp_U[t,k,j,i] = 0.5*(Tr.tracer_new[t,k,j,i+1] + Tr.tracer_new[t,k,j,i]) * A.U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], (D.u_3D.first[3]+1)..D.u_3D.last[3]} with (ref A) {
        var tmp : real;
        if (A.U[t,k,j,i] > 0) {
          tmp = Tr.tracer_new[t,k,j,i+1] - 2*Tr.tracer_new[t,k,j,i] + Tr.tracer_new[t,k,j,i-1];
        }
        else {
          tmp = Tr.tracer_new[t,k,j,i+2] - 2*Tr.tracer_new[t,k,j,i+1] + Tr.tracer_new[t,k,j,i];
        }
        A.tmp_U[t,k,j,i] = (0.5*(Tr.tracer_new[t,k,j,i+1] + Tr.tracer_new[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.last[3]..D.u_3D.last[3]} with (ref A) {
      A.tmp_U[t,k,j,i] = 0.5*(Tr.tracer_new[t,k,j,i+1] + Tr.tracer_new[t,k,j,i]) * A.U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..(D.u_3D.last[3]-1)} with (ref A) {
        var tmp : real;
        if (A.U[t,k,j,i] > 0) {
          tmp = Tr.tracer_new[t,k,j,i+1] - 2*Tr.tracer_new[t,k,j,i] + Tr.tracer_new[t,k,j,i-1];
        }
        else {
          tmp = Tr.tracer_new[t,k,j,i+2] - 2*Tr.tracer_new[t,k,j,i+1] + Tr.tracer_new[t,k,j,i];
        }
        A.tmp_U[t,k,j,i] = (0.5*(Tr.tracer_new[t,k,j,i+1] + Tr.tracer_new[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in D.u_3D with (ref A) {
        var tmp : real;
        if (A.U[t,k,j,i] > 0) {
          tmp = Tr.tracer_new[t,k,j,i+1] - 2*Tr.tracer_new[t,k,j,i] + Tr.tracer_new[t,k,j,i-1];
        }
        else {
          tmp = Tr.tracer_new[t,k,j,i+2] - 2*Tr.tracer_new[t,k,j,i+1] + Tr.tracer_new[t,k,j,i];
        }
        A.tmp_U[t,k,j,i] = (0.5*(Tr.tracer_new[t,k,j,i+1] + Tr.tracer_new[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
    }
  }

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.first[2]..D.v_3D.first[2], D.v_3D.dim[3]} with (ref A) {
    A.tmp_V[t,k,j,i] = 0.5*(Tr.tracer_new[t,k,j+1,i] + Tr.tracer_new[t,k,j,i]) * A.V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.last[2]..D.v_3D.last[2], D.v_3D.dim[3]} with (ref A) {
    A.tmp_V[t,k,j,i] = 0.5*(Tr.tracer_new[t,k,j+1,i] + Tr.tracer_new[t,k,j,i]) * A.V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], (D.v_3D.first[2]+1)..(D.v_3D.last[2]-1), D.v_3D.dim[3]} with (ref A) {
      var tmp : real;
      if (A.V[t,k,j,i] > 0) {
        tmp = Tr.tracer_new[t,k,j+1,i] - 2*Tr.tracer_new[t,k,j,i] + Tr.tracer_new[t,k,j-1,i];
      }
      else {
        tmp = Tr.tracer_new[t,k,j+2,i] - 2*Tr.tracer_new[t,k,j+1,i] + Tr.tracer_new[t,k,j,i];
      }
      A.tmp_V[t,k,j,i] = (0.5*(Tr.tracer_new[t,k,j+1,i] + Tr.tracer_new[t,k,j,i]) - us*tmp) * A.V[t,k,j,i];
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
      M[k,k] = 2*A.thickness[t,k,j,i] + A.thickness[t,k-1,j,i];
      M[k,k-1] = 2*A.thickness[t,k,j,i];
      M[k,k+1] = A.thickness[t,k-1,j,i];

      b[k] = 3*(A.thickness[t,k,j,i]*Tr.tracer_new[t,k-1,j,i] + A.thickness[t,k-1,j,i]*Tr.tracer_new[t,k,j,i]);
    }

    A.tmp_W[t,0..Nz,j,i] = solve(M, b) * A.W[t,0..Nz,j,i];

  }

}


proc calc_half_step(ref Ac: Arrays, ref An: Arrays) {

  Ac.H_plus  = 0.5 * (Ac.H_plus + An.H_plus);
  Ac.H_minus = 0.5 * (Ac.H_minus + Ac.H_minus);
  Ac.U        = 0.5 * (Ac.U + An.U);
  Ac.V        = 0.5 * (Ac.V + An.V);
  Ac.W        = 0.5 * (Ac.W + An.W);

}

proc predictor(ref A: Arrays, D: Domains, Tr: Tracers) {

  /////////////////////////////////////////
  //  Calculate tracer field at 1/2 timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref A) {
      Tr.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Tr.tracer_new[t,k,j,i] + (0.5 - 2*gamma) * Tr.tracer_old[t,k,j,i])*A.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k+1,j,i] - A.tmp_W[t,k,j,i]) ) ) / A.H_plus[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref A) {
      Tr.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Tr.tracer_new[t,k,j,i] + (0.5 - 2*gamma) * Tr.tracer_old[t,k,j,i])*A.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k+1,j,i] - A.tmp_W[t,k,j,i]) ) ) / A.H_plus[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref A) {
      Tr.predictor[t,k,j,i] = ( ((0.5 + 2*gamma) * Tr.tracer_new[t,k,j,i] + (0.5 - 2*gamma) * Tr.tracer_old[t,k,j,i])*A.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*((A.tmp_U[t,k,j,i] - A.tmp_U[t,k,j,i-1])
                                                 + (A.tmp_V[t,k,j,i] - A.tmp_V[t,k,j-1,i])
                                                 + (A.tmp_W[t,k+1,j,i] - A.tmp_W[t,k,j,i]) ) ) / A.H_plus[t,k,j,i];
    }
  }

}

