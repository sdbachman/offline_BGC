use files;
use tracers;
use domains;
use INPUTS;
use sigma_coordinate;

record Dynamics {

  var r0, r1, r2, r3 = 1..0;
  var u0, u1, u2, u3 = 1..0;
  var v0, v1, v2, v3 = 1..0;
  var w0, w1, w2, w3 = 1..0;

  var U_n : [u0, u1, u2, u3] real;
  var U_np1h : [u0, u1, u2, u3] real;
  var U_np1 : [u0, u1, u2, u3] real;
  var tmp_U : [u0, u1, u2, u3] real;

  var V_n : [v0, v1, v2, v3] real;
  var V_np1h : [v0, v1, v2, v3] real;
  var V_np1 : [v0, v1, v2, v3] real;
  var tmp_V : [v0, v1, v2, v3] real;

  var W_n : [w0, w1, w2, w3] real;
  var W_np1h : [w0, w1, w2, w3] real;
  var W_np1 : [w0, w1, w2, w3] real;
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

proc update_dynamics(ref U, ref V, ref H, D: Domains, F: Files, step : int) {

  // Read in u and v

    var u = get_var(F.vel[step], "u", D.u_3D);
    var v = get_var(F.vel[step], "v", D.v_3D);

  // Get volumetric fluxes for (n+1) timestep
    forall (t,k,j,i) in D.u_3D {
      U[t,k,j,i] = u[t,k,j,i] * 0.5 * (H[t,k,j,i] + H[t,k,j,i+1]) * dy;
    }

    forall (t,k,j,i) in D.v_3D {
      V[t,k,j,i] = v[t,k,j,i] * 0.5 * (H[t,k,j,i] + H[t,k,j+1,i]) * dx;
    }

}

proc calc_W(ref U, ref V, ref W, D: Domains, ref H_n, ref H_np1) {

  if (here.id == 0) {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (H_np1[t,k,j,i]-H_n[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]-1} {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (H_np1[t,k,j,i]-H_n[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }

    }
  }
  else {

    forall (t,j,i) in {0..0, (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3])..D.rho_3D.last[3]} {

      W[t,-1,j,i] = 0;
      for k in 0..<Nz {
        W[t,k,j,i] = W[t,k-1,j,i] - (H_np1[t,k,j,i]-H_n[t,k,j,i])*area/dt -
                   (U[t,k,j,i] - U[t,k,j,i-1] + V[t,k,j,i] - V[t,k,j-1,i]);
      }
    }
  }

}

proc calc_half_step_dyn(ref Dyn: Dynamics) {

  Dyn.U_np1h        = 0.5 * (Dyn.U_n + Dyn.U_np1);
  Dyn.V_np1h        = 0.5 * (Dyn.V_n + Dyn.V_np1);

}
