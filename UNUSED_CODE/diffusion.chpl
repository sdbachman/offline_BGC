use params;
use tracers;
use domains;
use INPUTS;
use sigma_coordinate;

use AllLocalesBarriers;

record Diffusion {

  var u0, u1, u2, u3 = 1..0;
  var v0, v1, v2, v3 = 1..0;

  var U_nm2 : [u0, u1, u2, u3] real;
  var U_nm1 : [u0, u1, u2, u3] real;
  var U_n : [u0, u1, u2, u3] real;
  var U_np1h : [u0, u1, u2, u3] real;
  var U_np1 : [u0, u1, u2, u3] real;
  var tmp_U : [u0, u1, u2, u3] real;

  var V_nm1 : [v0, v1, v2, v3] real;
  var V_nm2 : [v0, v1, v2, v3] real;
  var V_n : [v0, v1, v2, v3] real;
  var V_np1h : [v0, v1, v2, v3] real;
  var V_np1 : [v0, v1, v2, v3] real;
  var tmp_V : [v0, v1, v2, v3] real;

  proc init(arg: Domains) {
    this.u0 = arg.u_3D.dim[0];
    this.u1 = arg.u_3D.dim[1];
    this.u2 = arg.u_3D.dim[2];
    this.u3 = arg.u_3D.dim[3];

    this.v0 = arg.v_3D.dim[0];
    this.v1 = arg.v_3D.dim[1];
    this.v2 = arg.v_3D.dim[2];
    this.v3 = arg.v_3D.dim[3];

  }

}

proc calc_diffusive_fluxes(ref U, ref V, D: Domains, P: Params, ref arr, ref H) {

  /////////////////////////////////////////
  //             Zonal fluxes            //
  /////////////////////////////////////////

  forall (t,k,j,i) in D.u_3D {

    U[t,k,j,i] = 0.5*(sponge[t,k,j,i] + sponge[t,k,j,i+1]) * (arr[t,k,j,i+1] - arr[t,k,j,i])
                         * P.dy * 0.5 * (H[t,k,j,i] + H[t,k,j,i+1]) / P.dx;
  }

  /////////////////////////////////////////
  //         Meridional fluxes           //
  /////////////////////////////////////////

  forall (t,k,j,i) in D.v_3D {

    V[t,k,j,i] = 0.5*(sponge[t,k,j,i] + sponge[t,k,j+1,i]) * (arr[t,k,j+1,i] - arr[t,k,j,i])
                         * P.dx * 0.5 * (H[t,k,j,i] + H[t,k,j+1,i]) / P.dy;
  }

  allLocalesBarrier.barrier();

}

