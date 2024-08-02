use INPUTS;
use dynamics;
use params;
use domains;
use tracers;

use Time;
use AllLocalesBarriers;

proc calc_horizontal_fluxes(ref U, ref V, ref tmp_U, ref tmp_V, D: Domains, P: Params, ref arr) {

  // Horizontal: upstream-biased parabolic interpolation: SM05, after 4.13
  // Near open boundaries this method will assume the point "beyond" the edge has the same value as the edge

  /////////////////////////////////////////
  //              U-fluxes               //
  /////////////////////////////////////////

  if (here.id == 0) {

    forall (k,j,i) in {D.u_3D.dim[0], D.u_3D.dim[1], D.u_3D.first[2]..D.u_3D.first[2]} {
      tmp_U[k,j,i] = 0.5*(arr[k,j,i] + arr[k,j,i+1]) * U[k,j,i] - mask_rho[j,i+2]*
                             (P.one_sixth * max(U[k,j,i], 0.0) * (arr[k,j,i+1] - arr[k,j,i])
                           +  P.one_sixth * min(U[k,j,i], 0.0) * (arr[k,j,i+2] - 2*arr[k,j,i+1] + arr[k,j,i]));
   }

    forall (k,j,i) in {D.u_3D.dim[0], D.u_3D.dim[1], (D.u_3D.first[2]+1)..D.u_3D.last[2]} {
      tmp_U[k,j,i] = 0.5*(arr[k,j,i] + arr[k,j,i+1]) * U[k,j,i] - mask_rho[j,i-1]*mask_rho[j,i+2]*
                             (P.one_sixth * max(U[k,j,i], 0.0) * (arr[k,j,i+1] - 2*arr[k,j,i] + arr[k,j,i-1])
                           +  P.one_sixth * min(U[k,j,i], 0.0) * (arr[k,j,i+2] - 2*arr[k,j,i+1] + arr[k,j,i]));
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (k,j,i) in {D.u_3D.dim[0], D.u_3D.dim[1], D.u_3D.last[2]..D.u_3D.last[2]} {
      tmp_U[k,j,i] = 0.5*(arr[k,j,i] + arr[k,j,i+1]) * U[k,j,i] - mask_rho[j,i-1]*
                             (P.one_sixth * max(U[k,j,i], 0.0) * (arr[k,j,i+1] - 2*arr[k,j,i] + arr[k,j,i-1])
                          +  P.one_sixth * min(U[k,j,i], 0.0) * ( -arr[k,j,i+1] + arr[k,j,i]));
    }

    forall (k,j,i) in {D.u_3D.dim[0], D.u_3D.dim[1], D.u_3D.first[2]..(D.u_3D.last[2]-1)} {
      tmp_U[k,j,i] = 0.5*(arr[k,j,i] + arr[k,j,i+1]) * U[k,j,i] - mask_rho[j,i-1]*mask_rho[j,i+2]*
                             (P.one_sixth * max(U[k,j,i], 0.0) * (arr[k,j,i+1] - 2*arr[k,j,i] + arr[k,j,i-1])
                           +  P.one_sixth * min(U[k,j,i], 0.0) * (arr[k,j,i+2] - 2*arr[k,j,i+1] + arr[k,j,i]));
    }
  }
  else {
    forall (k,j,i) in D.u_3D {
      tmp_U[k,j,i] = 0.5*(arr[k,j,i] + arr[k,j,i+1]) * U[k,j,i] - mask_rho[j,i-1]*mask_rho[j,i+2]*
                             (P.one_sixth * max(U[k,j,i], 0.0) * (arr[k,j,i+1] - 2*arr[k,j,i] + arr[k,j,i-1])
                           +  P.one_sixth * min(U[k,j,i], 0.0) * (arr[k,j,i+2] - 2*arr[k,j,i+1] + arr[k,j,i]));
    }
  }

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (k,j,i) in {D.v_3D.dim[0], D.v_3D.first[1]..D.v_3D.first[1], D.v_3D.dim[2]} {
    tmp_V[k,j,i] = 0.5*(arr[k,j,i] + arr[k,j+1,i]) * V[k,j,i] - mask_rho[j+2,i]*
                             (P.one_sixth * max(V[k,j,i], 0.0) * (arr[k,j+1,i] - arr[k,j,i])
                           +  P.one_sixth * min(V[k,j,i], 0.0) * (arr[k,j+2,i] - 2*arr[k,j+1,i] + arr[k,j,i]));
  }

  forall (k,j,i) in {D.v_3D.dim[0], D.v_3D.last[1]..D.v_3D.last[1], D.v_3D.dim[2]} {
    tmp_V[k,j,i] = 0.5*(arr[k,j,i] + arr[k,j+1,i]) * V[k,j,i] - mask_rho[j-1,i]*
                             (P.one_sixth * max(V[k,j,i], 0.0) * (arr[k,j+1,i] - 2*arr[k,j,i] + arr[k,j-1,i])
                           +  P.one_sixth * min(V[k,j,i], 0.0) * (-arr[k,j+1,i] + arr[k,j,i]));
  }

  forall (k,j,i) in {D.v_3D.dim[0], (D.v_3D.first[1]+1)..(D.v_3D.last[1]-1), D.v_3D.dim[2]} {
    tmp_V[k,j,i] = 0.5*(arr[k,j,i] + arr[k,j+1,i]) * V[k,j,i] - mask_rho[j-1,i]*mask_rho[j+2,i]*
                             (P.one_sixth * max(V[k,j,i], 0.0) * (arr[k,j+1,i] - 2*arr[k,j,i] + arr[k,j-1,i])
                           +  P.one_sixth * min(V[k,j,i], 0.0) * (arr[k,j+2,i] - 2*arr[k,j+1,i] + arr[k,j,i]));
  }

  allLocalesBarrier.barrier();

}

