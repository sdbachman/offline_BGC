use INPUTS;
use arrays;
use domains;
use NetCDF_IO;

use LinearAlgebra;

//      Predictor step:  SM05, 4.8
//      Corrector step:  SM05, 4.2

//      Splines:  H07, 63



proc TimeStep(ref A: Arrays, D: Domains, tr: ?, tr_old: ?, predictor: ?) {

  var tmp_U : [D.u_3D] real;
  var tmp_V : [D.v_3D] real;
  var tmp_W : [D.w_3D] real;

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
        tmp_U[t,k,j,i] = 0.5*(tr[t,k,j,i+1] + tr[t,k,j,i]) * A.U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], (D.u_3D.first[3]+1)..D.u_3D.last[3]} with (ref A) {
        var tmp : real;
        if (A.U[t,k,j,i] > 0) {
          tmp = tr[t,k,j,i+1] - 2*tr[t,k,j,i] + tr[t,k,j,i-1];
        }
        else {
          tmp = tr[t,k,j,i+2] - 2*tr[t,k,j,i+1] + tr[t,k,j,i];
        }
        tmp_U[t,k,j,i] = (0.5*(tr[t,k,j,i+1] + tr[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.last[3]..D.u_3D.last[3]} with (ref A) {
      tmp_U[t,k,j,i] = 0.5*(tr[t,k,j,i+1] + tr[t,k,j,i]) * A.U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..(D.u_3D.last[3]-1)} with (ref A) {
        var tmp : real;
        if (A.U[t,k,j,i] > 0) {
          tmp = tr[t,k,j,i+1] - 2*tr[t,k,j,i] + tr[t,k,j,i-1];
        }
        else {
          tmp = tr[t,k,j,i+2] - 2*tr[t,k,j,i+1] + tr[t,k,j,i];
        }
        tmp_U[t,k,j,i] = (0.5*(tr[t,k,j,i+1] + tr[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in D.u_3D with (ref A) {
        var tmp : real;
        if (A.U[t,k,j,i] > 0) {
          tmp = tr[t,k,j,i+1] - 2*tr[t,k,j,i] + tr[t,k,j,i-1];
        }
        else {
          tmp = tr[t,k,j,i+2] - 2*tr[t,k,j,i+1] + tr[t,k,j,i];
        }
        tmp_U[t,k,j,i] = (0.5*(tr[t,k,j,i+1] + tr[t,k,j,i]) - us*tmp) * A.U[t,k,j,i];
    }
  }

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.first[2]..D.v_3D.first[2], D.v_3D.dim[3]} with (ref A) {
    tmp_V[t,k,j,i] = 0.5*(tr[t,k,j+1,i] + tr[t,k,j,i]) * A.V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.last[2]..D.v_3D.last[2], D.v_3D.dim[3]} with (ref A) {
    tmp_V[t,k,j,i] = 0.5*(tr[t,k,j+1,i] + tr[t,k,j,i]) * A.V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], (D.v_3D.first[2]+1)..(D.v_3D.last[2]-1), D.v_3D.dim[3]} with (ref A) {
      var tmp : real;
      if (A.V[t,k,j,i] > 0) {
        tmp = tr[t,k,j+1,i] - 2*tr[t,k,j,i] + tr[t,k,j-1,i];
      }
      else {
        tmp = tr[t,k,j+2,i] - 2*tr[t,k,j+1,i] + tr[t,k,j,i];
      }
      tmp_V[t,k,j,i] = (0.5*(tr[t,k,j+1,i] + tr[t,k,j,i]) - us*tmp) * A.V[t,k,j,i];
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

      b[k] = 3*(A.thickness[t,k,j,i]*tr[t,k-1,j,i] + A.thickness[t,k-1,j,i]*tr[t,k,j,i]);
    }

    tmp_W[t,0..Nz,j,i] = solve(M, b) * A.W[t,0..Nz,j,i];

  }


  /////////////////////////////////////////
  //  Calculate tracer field at 1/2 timestep  //
  /////////////////////////////////////////

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref A) {
      predictor[t,k,j,i] = ( 0*((0.5 + 2*gamma) * tr[t,k,j,i] + (0.5 - 2*gamma) * tr_old[t,k,j,i])*A.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*(0*(tmp_U[t,k,j,i] - tmp_U[t,k,j,i-1])
                                                 + (tmp_V[t,k,j,i] - tmp_V[t,k,j-1,i])
                                                 + 0*(tmp_W[t,k+1,j,i] - tmp_W[t,k,j,i]) ) ) / A.H_plus[t,k,j,i];

    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref A) {
      predictor[t,k,j,i] = ( 0*((0.5 + 2*gamma) * tr[t,k,j,i] + (0.5 - 2*gamma) * tr_old[t,k,j,i])*A.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*(0*(tmp_U[t,k,j,i] - tmp_U[t,k,j,i-1])
                                                 + (tmp_V[t,k,j,i] - tmp_V[t,k,j-1,i])
                                                 + 0*(tmp_W[t,k+1,j,i] - tmp_W[t,k,j,i]) ) ) / A.H_plus[t,k,j,i];

    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref A) {
      predictor[t,k,j,i] = ( 0*((0.5 + 2*gamma) * tr[t,k,j,i] + (0.5 - 2*gamma) * tr_old[t,k,j,i])*A.H_minus[t,k,j,i]
                         - (1 - 2*gamma)*dt*iarea*(0*(tmp_U[t,k,j,i] - tmp_U[t,k,j,i-1])
                                                 + (tmp_V[t,k,j,i] - tmp_V[t,k,j-1,i])
                                                 + 0*(tmp_W[t,k+1,j,i] - tmp_W[t,k,j,i]) ) ) / A.H_plus[t,k,j,i];
      //predictor[t,k,j,i] = A.H_minus[t,k,j,i] / A.H_plus[t,k,j,i];

    }
  }

   WriteOutput("pred.nc", predictor, "predictor", "stuff", 10);


}
