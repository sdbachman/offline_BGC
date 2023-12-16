use domains;
use INPUTS;
use sigma_coordinate;

record Arrays {

  const r0, r1, r2, r3 = 1..0;
  const u0, u1, u2, u3 = 1..0;
  const v0, v1, v2, v3 = 1..0;
  const w0, w1, w2, w3 = 1..0;

  var h : [r2, r3] real;
  var thickness0 : [r1, r2, r3] real;
  //var volume0 : [r1, r2, r3] real;
  var thickness : [r0, r1, r2, r3] real;
  //var volume : [r0, r1, r2, r3] real;
  //var dHdt : [r0, r1, r2, r3] real;
  //var dVdt : [r0, r1, r2, r3] real;
  var H_plus : [r0, r1, r2, r3] real;
  var H_minus : [r0, r1, r2, r3] real;

  var zeta_old : [r0, r2, r3] real;
  var zeta_new : [r0, r2, r3] real;

  var u : [u0, u1, u2, u3] real;
  var U : [u0, u1, u2, u3] real;

  var v : [v0, v1, v2, v3] real;
  var V : [v0, v1, v2, v3] real;

  var w : [w0, w1, w2, w3] real;
  var W : [w0, w1, w2, w3] real;


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

proc set_static_arrays(ref A: Arrays, D: Domains) {

  A.h = get_var(gridfile, "h", D.grid);

  A.thickness0 = get_thickness0(A.h);

  //A.volume0 = A.thickness0 * dx * dy;

}

proc update_dynamic_arrays(ref A: Arrays, D: Domains, infiles, step : int) {

  // Read in SSH, update thicknesses
    A.zeta_new = get_var(infiles[step], "zeta", D.rho_2D);

    if (step == Nt_start) {
      A.zeta_old = A.zeta_new;
    }

  // From SM09, Eq. 2.13
    forall (t,k,j,i) in D.rho_3D with (ref A) {
      //A.dHdt[t,k,j,i] = A.thickness0[k,j,i] * (A.zeta_new[t,j,i] - A.zeta_old[t,j,i]) / (A.h[j,i] * dt);
      //A.dVdt[t,k,j,i] = A.dHdt[t,k,j,i] * dx * dy;

      A.thickness[t,k,j,i] = A.thickness0[k,j,i] * (1 + A.zeta_new[t,j,i] / A.h[j,i]);
    }

  // Read in u and v
    A.u = get_var(infiles[step], "u", D.u_3D);
    A.v = get_var(infiles[step], "v", D.v_3D);
    A.w = get_var(infiles[step], "omega", D.w_3D);

  // Update volumetric fluxes
    forall (t,k,j,i) in D.u_3D with (ref A) {
      A.U[t,k,j,i] = A.u[t,k,j,i] * 0.5 * (A.thickness[t,k,j,i] + A.thickness[t,k,j,i+1]) * dy;
    }

    forall (t,k,j,i) in D.v_3D with (ref A) {
      A.V[t,k,j,i] = A.v[t,k,j,i] * 0.5 * (A.thickness[t,k,j,i] + A.thickness[t,k,j+1,i]) * dx;
    }

    // From SH05, Eq. 1.18
    //forall (t,k,j,i) in {0..0, 0..<Nz, D.w_3D.dim[1], D.w_3D.dim[2]} with (ref A) {
    forall (t,k,j,i) in D.w_3D with (ref A) {
      A.W[t,k,j,i] = A.w[t,k,j,i] * dx * dy;
    }

}

