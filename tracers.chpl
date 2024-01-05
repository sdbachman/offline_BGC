use INPUTS;
use domains;
use files;
use sigma_coordinate;
use NetCDF_IO;
use StencilDist;
use AllLocalesBarriers;

class Tracers {

  // Creating a singleton first dimension to store the time
  const FullDomain_grid = {0..<Ny, 0..<Nx};
  const FullDomain2D = {0..0, 0..<Ny, 0..<Nx};
  const FullDomain3D = {0..0, 0..<Nz, 0..<Ny, 0..<Nx};
  var myTargetLocales_grid = reshape(Locales, {1..1, 1..Locales.size});
  var myTargetLocales2D = reshape(Locales, {1..1, 1..1, 1..Locales.size});
  var myTargetLocales3D = reshape(Locales, {1..1, 1..1, 1..1, 1..Locales.size});

  // Create stencilDist arrays for the tracers. This is necessary to do the halo updates efficiently.
    const stencil_grid = new stencilDist(boundingBox=FullDomain_grid, targetLocales=myTargetLocales_grid, fluff=(0,2));
    const stencil2D = new stencilDist(boundingBox=FullDomain2D, targetLocales=myTargetLocales2D, fluff=(0,0,2));
    const stencil3D = new stencilDist(boundingBox=FullDomain3D, targetLocales=myTargetLocales3D, fluff=(0,0,0,2));
    const D_grid = stencil_grid.createDomain(FullDomain_grid);
    const D2 = stencil2D.createDomain(FullDomain2D);
    const D3 = stencil3D.createDomain(FullDomain3D);

    var tracer : [D3] real;

    var k1 : [D3] real;
    var k2 : [D3] real;
    var k3 : [D3] real;
    var k4 : [D3] real;
    var tmp : [D3] real;

    var h : [D_grid] real;
    var H0 : [D3] real;
    var H_nm1 : [D3] real;
    var H_nm1h : [D3] real;
    var H_n : [D3] real;
    var H_np1h : [D3] real;
    var H_np1 : [D3] real;
    var H_np3h : [D3] real;

    var zeta_nm1 : [D2] real;
    var zeta_nm1h : [D2] real;
    var zeta_n : [D2] real;
    var zeta_np1h : [D2] real;
    var zeta_np1 : [D2] real;
    var zeta_np3h : [D2] real;

    var div : [D3] real;
    var div2 : [D3] real;
    var dV  : [D3] real;
}

proc initialize(ref Tr: Tracers, F: Files, grdfile, D) {

    Tr.h[D.grid] = get_var(grdfile, "h", D.grid);

    Tr.H0[D.rho_3D] = get_H0(Tr.h[D.grid]);

    // Initialize tracer fields
      Tr.tracer[D.rho_3D] = get_var(F.vel[Nt_start+1], "temp", D.rho_3D);

    // Initialize zeta and thicknesses
      update_thickness(Tr.zeta_nm1, Tr.H_nm1, Tr.H0, Tr.h, D, F, Nt_start);
      update_thickness(Tr.zeta_n, Tr.H_n, Tr.H0, Tr.h, D, F, Nt_start+1);

    // Update the halos for the static arrays
      allLocalesBarrier.barrier();
      if (here.id == 0) {
        Tr.h.updateFluff();
        Tr.H0.updateFluff();
      }
      allLocalesBarrier.barrier();

}

proc update_halos(ref arr) {

    // Update the halos. Only need to issue the updateFluff
    // method from one Locale.

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      arr.updateFluff();
    }
    allLocalesBarrier.barrier();

}

proc update_thickness(ref zeta, ref H, ref H0, ref h, D: Domains, F: Files, step : int) {

  // Read in SSH, update thicknesses for (n+1)
    zeta[D.rho_2D] = get_var(F.vel[step], "zeta", D.rho_2D);

  // From SM09, Eq. 2.13
    forall (t,k,j,i) in D.rho_3D {
      H[t,k,j,i] = H0[t,k,j,i] * (1 + zeta[t,j,i] / h[j,i]);
    }

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      H.updateFluff();
    }
    allLocalesBarrier.barrier();

}

proc calc_half_step_tr(ref Tr: Tracers, D: Domains) {

    forall (t,j,i) in D.rho_2D {
      Tr.zeta_nm1h[t,j,i] = 0.5 * (Tr.zeta_nm1[t,j,i] + Tr.zeta_n[t,j,i]);
      Tr.zeta_np1h[t,j,i] = 0.5 * (Tr.zeta_n[t,j,i] + Tr.zeta_np1[t,j,i]);
      Tr.zeta_np3h[t,j,i] = (1.5 + beta) * Tr.zeta_np1[t,j,i] - (0.5 + 2*beta) * Tr.zeta_n[t,j,i] + beta * Tr.zeta_nm1[t,j,i];
    }

  // From SM09, Eq. 2.13
    forall (t,k,j,i) in D.rho_3D {
      Tr.H_nm1h[t,k,j,i] = Tr.H0[t,k,j,i] * (1 + Tr.zeta_nm1h[t,j,i] / Tr.h[j,i]);
      Tr.H_np1h[t,k,j,i] = Tr.H0[t,k,j,i] * (1 + Tr.zeta_np1h[t,j,i] / Tr.h[j,i]);
      Tr.H_np3h[t,k,j,i] = Tr.H0[t,k,j,i] * (1 + Tr.zeta_np3h[t,j,i] / Tr.h[j,i]);
    }

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      Tr.H_nm1h.updateFluff();
      Tr.H_np1h.updateFluff();
      Tr.H_np3h.updateFluff();
    }
    allLocalesBarrier.barrier();
}

