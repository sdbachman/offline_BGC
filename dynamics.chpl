use params;
use tracers;
use domains;
use INPUTS;
use sigma_coordinate;

use AllLocalesBarriers;

record Dynamics {

  var r0, r1, r2 = 1..0;
  var u0, u1, u2 = 1..0;
  var v0, v1, v2 = 1..0;
  var w0, w1, w2 = 1..0;

  var U_n : [u0, u1, u2] real;
  var U_np1h : [u0, u1, u2] real;
  var U_np1 : [u0, u1, u2] real;
  var tmp_U : [u0, u1, u2] real;

  var V_n : [v0, v1, v2] real;
  var V_np1h : [v0, v1, v2] real;
  var V_np1 : [v0, v1, v2] real;
  var tmp_V : [v0, v1, v2] real;

  var u_n : [u0, u1, u2] real;
  var u_np1h : [u0, u1, u2] real;
  var u_np1 : [u0, u1, u2] real;

  var v_n : [v0, v1, v2] real;
  var v_np1h : [v0, v1, v2] real;
  var v_np1 : [v0, v1, v2] real;


  proc init(arg: Domains) {

    this.r0 = arg.rho_3D.dim[0];
    this.r1 = arg.rho_3D.dim[1];
    this.r2 = arg.rho_3D.dim[2];

    this.u0 = arg.u_3D.dim[0];
    this.u1 = arg.u_3D.dim[1];
    this.u2 = arg.u_3D.dim[2];

    this.v0 = arg.v_3D.dim[0];
    this.v1 = arg.v_3D.dim[1];
    this.v2 = arg.v_3D.dim[2];

  }

}

proc update_dynamics(ref u, ref v, ref U, ref V, ref H, D: Domains, P: Params, step : int) {

  // Read in u and v

//    u = get_var(P.velfiles[step], "u", D.u_3D);
//    v = get_var(P.velfiles[step], "v", D.v_3D);
    u = readZarrArrayLocal(P.velfiles[step] + (here.id : string) + '/u/', real(32), 3);
    v = readZarrArrayLocal(P.velfiles[step] + (here.id : string) + '/v/', real(32), 3);


    calc_volumetric_fluxes(u, v, U, V, H, D, P);

}


proc calc_volumetric_fluxes(ref u, ref v, ref U, ref V, ref H, D: Domains, P: Params) {

  // Get volumetric fluxes for (n+1) timestep
    forall (k,j,i) in D.u_3D {
      U[k,j,i] = u[k,j,i] * 0.5 * (H[k,j,i] + H[k,j,i+1]) * P.dy;
    }

    forall (k,j,i) in D.v_3D {
      V[k,j,i] = v[k,j,i] * 0.5 * (H[k,j,i] + H[k,j+1,i]) * P.dx;
    }

  allLocalesBarrier.barrier();

}


proc calc_half_step_dyn(ref Dyn: Dynamics, D: Domains, P: Params) {

  Dyn.u_np1h = 0.5 * (Dyn.u_n + Dyn.u_np1);
  Dyn.v_np1h = 0.5 * (Dyn.v_n + Dyn.v_np1);

  allLocalesBarrier.barrier();

}

