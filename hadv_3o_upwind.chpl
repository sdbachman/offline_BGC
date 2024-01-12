use INPUTS;
use dynamics;
use params;
use domains;
use tracers;

use Time;
use AllLocalesBarriers;

proc calc_horizontal_fluxes(ref U, ref V, ref Dyn: Dynamics, D: Domains, P: Params, ref arr, ref H) {

  // Horizontal: upstream-biased parabolic interpolation: SM05, after 4.13

  // Near open boundaries this method will assume the point "beyond" the edge has the same value as the edge

  /////////////////////////////////////////
  //              U-fluxes               //
  /////////////////////////////////////////

  if (here.id == 0) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..D.u_3D.first[3]} with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = 0.5*(arr[t,k,j,i] + arr[t,k,j,i+1]) * U[t,k,j,i] - mask_rho[j,i+2]*
                             (P.one_sixth * max(U[t,k,j,i], 0.0) * (arr[t,k,j,i+1] - arr[t,k,j,i])
                           +  P.one_sixth * min(U[t,k,j,i], 0.0) * (arr[t,k,j,i+2] - 2*arr[t,k,j,i+1] + arr[t,k,j,i]));
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], (D.u_3D.first[3]+1)..D.u_3D.last[3]} with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = 0.5*(arr[t,k,j,i] + arr[t,k,j,i+1]) * U[t,k,j,i] - mask_rho[j,i-1]*mask_rho[j,i+2]*
                             (P.one_sixth * max(U[t,k,j,i], 0.0) * (arr[t,k,j,i+1] - 2*arr[t,k,j,i] + arr[t,k,j,i-1])
                           +  P.one_sixth * min(U[t,k,j,i], 0.0) * (arr[t,k,j,i+2] - 2*arr[t,k,j,i+1] + arr[t,k,j,i]));
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.last[3]..D.u_3D.last[3]} with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = 0.5*(arr[t,k,j,i] + arr[t,k,j,i+1]) * U[t,k,j,i] - mask_rho[j,i-1]*
                             (P.one_sixth * max(U[t,k,j,i], 0.0) * (arr[t,k,j,i+1] - 2*arr[t,k,j,i] + arr[t,k,j,i-1])
                           +  P.one_sixth * min(U[t,k,j,i], 0.0) * ( -arr[t,k,j,i+1] + arr[t,k,j,i]));
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..(D.u_3D.last[3]-1)} with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = 0.5*(arr[t,k,j,i] + arr[t,k,j,i+1]) * U[t,k,j,i] - mask_rho[j,i-1]*mask_rho[j,i+2]*
                             (P.one_sixth * max(U[t,k,j,i], 0.0) * (arr[t,k,j,i+1] - 2*arr[t,k,j,i] + arr[t,k,j,i-1])
                           +  P.one_sixth * min(U[t,k,j,i], 0.0) * (arr[t,k,j,i+2] - 2*arr[t,k,j,i+1] + arr[t,k,j,i]));
    }
  }
  else {
    forall (t,k,j,i) in D.u_3D with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = 0.5*(arr[t,k,j,i] + arr[t,k,j,i+1]) * U[t,k,j,i] - mask_rho[j,i-1]*mask_rho[j,i+2]*
                             (P.one_sixth * max(U[t,k,j,i], 0.0) * (arr[t,k,j,i+1] - 2*arr[t,k,j,i] + arr[t,k,j,i-1])
                           +  P.one_sixth * min(U[t,k,j,i], 0.0) * (arr[t,k,j,i+2] - 2*arr[t,k,j,i+1] + arr[t,k,j,i]));
    }
  }

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.first[2]..D.v_3D.first[2], D.v_3D.dim[3]} with (ref Dyn) {
    Dyn.tmp_V[t,k,j,i] = 0.5*(arr[t,k,j,i] + arr[t,k,j,i+1]) * V[t,k,j,i] - mask_rho[j+2,i]*
                             (P.one_sixth * max(V[t,k,j,i], 0.0) * (arr[t,k,j+1,i] - arr[t,k,j,i])
                           +  P.one_sixth * min(V[t,k,j,i], 0.0) * (arr[t,k,j+2,i] - 2*arr[t,k,j+1,i] + arr[t,k,j,i]));
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.last[2]..D.v_3D.last[2], D.v_3D.dim[3]} with (ref Dyn) {
    Dyn.tmp_V[t,k,j,i] = 0.5*(arr[t,k,j,i] + arr[t,k,j,i+1]) * V[t,k,j,i] - mask_rho[j-1,i]*
                             (P.one_sixth * max(V[t,k,j,i], 0.0) * (arr[t,k,j+1,i] - 2*arr[t,k,j,i] + arr[t,k,j-1,i])
                           +  P.one_sixth * min(V[t,k,j,i], 0.0) * (-arr[t,k,j+1,i] + arr[t,k,j,i]));
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], (D.v_3D.first[2]+1)..(D.v_3D.last[2]-1), D.v_3D.dim[3]} with (ref Dyn) {
    Dyn.tmp_V[t,k,j,i] = 0.5*(arr[t,k,j,i] + arr[t,k,j,i+1]) * V[t,k,j,i] - mask_rho[j-1,i]*mask_rho[j+2,i]*
                             (P.one_sixth * max(V[t,k,j,i], 0.0) * (arr[t,k,j+1,i] - 2*arr[t,k,j,i] + arr[t,k,j-1,i])
                           +  P.one_sixth * min(V[t,k,j,i], 0.0) * (arr[t,k,j+2,i] - 2*arr[t,k,j+1,i] + arr[t,k,j,i]));
  }

  allLocalesBarrier.barrier();

}

