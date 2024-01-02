use files;
use tracers;
use domains;
use INPUTS;
use sigma_coordinate;

record Arrays {

  var r0, r1, r2, r3 = 1..0;
  var u0, u1, u2, u3 = 1..0;
  var v0, v1, v2, v3 = 1..0;
  var w0, w1, w2, w3 = 1..0;

/*
  var h : [r2, r3] real;
  var thickness0 : [r1, r2, r3] real;
  var thickness : [r0, r1, r2, r3] real;
  var H_plus : [r0, r1, r2, r3] real;
  var H_minus : [r0, r1, r2, r3] real;

  var zeta_old : [r0, r2, r3] real;
  var zeta_new : [r0, r2, r3] real;
*/

  var u : [u0, u1, u2, u3] real;
  var U : [u0, u1, u2, u3] real;
  var U_half : [u0, u1, u2, u3] real;
  var tmp_U : [u0, u1, u2, u3] real;

  var v : [v0, v1, v2, v3] real;
  var V : [v0, v1, v2, v3] real;
  var V_half : [v0, v1, v2, v3] real;
  var tmp_V : [v0, v1, v2, v3] real;

  var w : [w0, w1, w2, w3] real;
  var W : [w0, w1, w2, w3] real;
  var W_half : [w0, w1, w2, w3] real;
  var tmp_W : [w0, w1, w2, w3] real;


  proc init(arg: Domains) {
    this.r0 = arg.rho_3D.dim[0];
    this.r1 = arg.rho_3D.dim[1];
    this.r2 = arg.rho_3D.dim[2];
    this.r3 = arg.rho_3D.dim[3];

    this.u0 = arg.u_3D.dim[0];
    this.u1 = arg.u_3D.dim[1];
    this.u2 = arg.u_3D.dim[2];
    this.u3 = arg.u_3D.dim[3];

    this.v0 = arg.v_3D.dim[0];
    this.v1 = arg.v_3D.dim[1];
    this.v2 = arg.v_3D.dim[2];
    this.v3 = arg.v_3D.dim[3];

    this.w0 = arg.w_3D.dim[0];
    this.w1 = arg.w_3D.dim[1];
    this.w2 = arg.w_3D.dim[2];
    this.w3 = arg.w_3D.dim[3];
  }

}

/*
proc set_static_arrays(ref A: Arrays, D: Domains, F: Files) {

  A.h = get_var(F.grd, "h", D.grid);

  A.thickness0 = get_thickness0(A.h);

}
*/

proc update_dynamic_arrays(ref A: Arrays, ref Tr: Tracers, D: Domains, F: Files, step : int) {

  // Read in u and v
    A.u = get_var(F.vel[step], "u", D.u_3D);
    A.v = get_var(F.vel[step], "v", D.v_3D);
//    A.w = get_var(F.vel[step], "omega", D.w_3D);

  // Update volumetric fluxes
    forall (t,k,j,i) in D.u_3D with (ref A) {
      A.U[t,k,j,i] = A.u[t,k,j,i] * 0.5 * (Tr.thickness[t,k,j,i] + Tr.thickness[t,k,j,i+1]) * dy;
    }

    forall (t,k,j,i) in D.v_3D with (ref A) {
      A.V[t,k,j,i] = A.v[t,k,j,i] * 0.5 * (Tr.thickness[t,k,j,i] + Tr.thickness[t,k,j+1,i]) * dx;
    }

//    // From SH05, Eq. 1.18
//    forall (t,k,j,i) in D.w_3D with (ref A) {
//      A.W[t,k,j,i] = A.w[t,k,j,i] * dx * dy;
//    }

}

