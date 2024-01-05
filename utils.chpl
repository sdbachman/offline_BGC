use INPUTS;
use dynamics;
use domains;
use tracers;
use files;


proc continuity_compare(ref Dyn: Dynamics, D: Domains, ref thickness_c, ref thickness_n, Trc: Tracers) {

  if (here.id == 0) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), (D.rho_3D.first[3]+1)..D.rho_3D.last[3]} with (ref Dyn) {
      Trc.div[t,k,j,i] =  -dt*((Dyn.U[t,k,j,i] - Dyn.U[t,k,j,i-1])
                             + (Dyn.V[t,k,j,i] - Dyn.V[t,k,j-1,i])
                             + (Dyn.W[t,k,j,i] - Dyn.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (thickness_n[t,k,j,i] - thickness_c[t,k,j,i]) * dx * dy;
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..(D.rho_3D.last[3]-1)} with (ref Dyn) {
      Trc.div[t,k,j,i] =  -dt*((Dyn.U[t,k,j,i] - Dyn.U[t,k,j,i-1])
                             + (Dyn.V[t,k,j,i] - Dyn.V[t,k,j-1,i])
                             + (Dyn.W[t,k,j,i] - Dyn.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (thickness_n[t,k,j,i] - thickness_c[t,k,j,i]) * dx * dy;
    }
  }
  else {
    forall (t,k,j,i) in {0..0, D.rho_3D.dim[1], (D.rho_3D.first[2]+1)..(D.rho_3D.last[2]-1), D.rho_3D.first[3]..D.rho_3D.last[3]} with (ref Dyn) {

      Trc.div[t,k,j,i] = -dt*((Dyn.U[t,k,j,i] - Dyn.U[t,k,j,i-1])
                                                 + (Dyn.V[t,k,j,i] - Dyn.V[t,k,j-1,i])
                                                 + (Dyn.W[t,k,j,i] - Dyn.W[t,k-1,j,i]) );
      Trc.dV[t,k,j,i] = (thickness_n[t,k,j,i] - thickness_c[t,k,j,i]) * dx * dy;
    }
  }
}
