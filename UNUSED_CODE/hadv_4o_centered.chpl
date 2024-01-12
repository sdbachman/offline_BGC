use INPUTS;
use dynamics;
use params;
use domains;
use tracers;

use Time;
use AllLocalesBarriers;

proc calc_horizontal_fluxes(ref U, ref V, ref Dyn: Dynamics, D: Domains, P: Params, ref arr, ref H) {

  /////////////////////////////////////////
  //              U-fluxes               //
  /////////////////////////////////////////

  if (here.id == 0) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..D.u_3D.first[3]} with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = ((5*arr[t,k,j,i] + 8*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], (D.u_3D.first[3]+1)..D.u_3D.last[3]} with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = ((-arr[t,k,j,i-1] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * U[t,k,j,i];
    }
  }
  else if (here.id == (Locales.size-1)) {

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.last[3]..D.u_3D.last[3]} with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = ((5*arr[t,k,j,i+1] + 8*arr[t,k,j,i] - arr[t,k,j,i-1])/12) * U[t,k,j,i];
    }

    forall (t,k,j,i) in {0..0, D.u_3D.dim[1], D.u_3D.dim[2], D.u_3D.first[3]..(D.u_3D.last[3]-1)} with (ref Dyn) {
        Dyn.tmp_U[t,k,j,i] = ((-arr[t,k,j,i-1] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * U[t,k,j,i];
    }
  }
  else {
    forall (t,k,j,i) in D.u_3D with (ref Dyn) {
      Dyn.tmp_U[t,k,j,i] = ((-arr[t,k,j,i-1] + 7*arr[t,k,j,i] + 7*arr[t,k,j,i+1] - arr[t,k,j,i+2])/12) * U[t,k,j,i];
    }
  }

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.first[2]..D.v_3D.first[2], D.v_3D.dim[3]} with (ref Dyn) {
     Dyn.tmp_V[t,k,j,i] = ((5*arr[t,k,j,i] + 8*arr[t,k,j+1,i] - arr[t,k,j+2,i])/12) * V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], D.v_3D.last[2]..D.v_3D.last[2], D.v_3D.dim[3]} with (ref Dyn) {
     Dyn.tmp_V[t,k,j,i] = ((5*arr[t,k,j+1,i] + 8*arr[t,k,j,i] - arr[t,k,j-1,i])/12) * V[t,k,j,i];
  }

  forall (t,k,j,i) in {0..0, D.v_3D.dim[1], (D.v_3D.first[2]+1)..(D.v_3D.last[2]-1), D.v_3D.dim[3]} with (ref Dyn) {
     Dyn.tmp_V[t,k,j,i] = ((-arr[t,k,j-1,i] + 7*arr[t,k,j,i] + 7*arr[t,k,j+1,i] - arr[t,k,j+2,i])/12) * V[t,k,j,i];
  }

  allLocalesBarrier.barrier();
}
