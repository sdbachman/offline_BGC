use INPUTS;
use domains;
use params;
use dynamics;
use sigma_coordinate;
use NetCDF_IO;
use StencilDist;
use AllLocalesBarriers;

// Notation here will be consistent with the description in
// https://adcroft.github.io/assets/pdf/ALE_workshop_NCWCP_2016.pdf


// Creating a singleton first dimension to store the time
  const FullDomain_grid = {0..<Ny_, 0..<Nx_};
  const FullDomain2D = {0..<Ny_, 0..<Nx_};
  const FullDomain3D = {0..<Nz_, 0..<Ny_, 0..<Nx_};
  var myTargetLocales_grid = reshape(Locales, {1..1, 1..Locales.size});
  var myTargetLocales2D = reshape(Locales, {1..1, 1..Locales.size});
  var myTargetLocales3D = reshape(Locales, {1..1, 1..1, 1..Locales.size});

// Create stencilDist arrays for the tracers. This is necessary to do the halo updates efficiently.
  const stencil_grid = new stencilDist(boundingBox=FullDomain_grid, targetLocales=myTargetLocales_grid, fluff=(0,2));
  const stencil2D = new stencilDist(boundingBox=FullDomain2D, targetLocales=myTargetLocales2D, fluff=(0,2));
  const stencil3D = new stencilDist(boundingBox=FullDomain3D, targetLocales=myTargetLocales3D, fluff=(0,0,2));
  const D_grid = stencil_grid.createDomain(FullDomain_grid);
  const D2 = stencil2D.createDomain(FullDomain2D);
  const D3 = stencil3D.createDomain(FullDomain3D);

// For RK3
  var k1 : [D3] real;
  var k2 : [D3] real;
  var k3 : [D3] real;
  var ktmp : [D3] real;

  var tracer_n : [D3] real;
  var tracer_np1h : [D3] real;
  var tracer_np1 : [D3] real;
  var tracer_dagger : [D3] real;

  var mask_rho : [D_grid] real;
  var h : [D_grid] real;
  var H0 : [D3] real;
  var H_n : [D3] real;
  var H_np1h : [D3] real;
  var H_np1 : [D3] real;
  var H_dagger : [D3] real;

  var zeta_n : [D2] real;
  var zeta_np1h : [D2] real;
  var zeta_np1 : [D2] real;

  var sponge : [D3] real;

  var kappa_v : [D3] real;

proc initialize_tr(D, P: Params) {

    mask_rho[D.grid] = get_var(P.grdfile, "mask_rho", D.grid);
    h[D.grid] = get_var(P.grdfile, "h", D.grid);


    H0[D.rho_3D] = get_H0(h[D.grid], P);

//    kappa_v[D.rho_3D] = get_var(P.velfiles[P.Nt_start], "Akt", D.rho_3D);
    kappa_v[D.rho_3D] = readZarrArrayLocal(P.velfiles[P.Nt_start] + (here.id : string) + '/Akt/', real(32), 3);

    // Initialize tracer fields
//    tracer_n[D.rho_3D] = get_var(P.velfiles[P.Nt_start], "temp", D.rho_3D);
    tracer_n[D.rho_3D] = readZarrArrayLocal(P.velfiles[P.Nt_start] + (here.id : string) + '/temp/', real(32), 3);

    // Initialize zeta and thicknesses
      update_thickness(zeta_n, H_n, H0, h, D, P, P.Nt_start);

    // Update the halos for the static arrays
      allLocalesBarrier.barrier();
      if (here.id == 0) {
        mask_rho.updateFluff();
        h.updateFluff();
        H0.updateFluff();
      }
      allLocalesBarrier.barrier();

    // Update the halos for the tracer arrays
      allLocalesBarrier.barrier();
      if (here.id == 0) {
        tracer_n.updateFluff();
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

proc update_thickness(ref zeta, ref H, ref H0, ref h, D: Domains, P: Params, step : int) {

  // Read in SSH, update thicknesses for (n+1)
//    zeta[D.rho_2D] = get_var(P.velfiles[step], "zeta", D.rho_2D);
    zeta[D.rho_2D] = readZarrArrayLocal(P.velfiles[step] + (here.id : string) + '/zeta/', real(32), 2);

  // From SM09, Eq. 2.13
    forall (k,j,i) in D.rho_3D {
      H[k,j,i] = H0[k,j,i] * (1 + zeta[j,i] / h[j,i]);
    }

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      H.updateFluff();
    }
    allLocalesBarrier.barrier();

}

proc calc_half_step_tr(D: Domains, P: Params) {

    forall (j,i) in D.rho_2D {
      zeta_np1h[j,i] = 0.5 * (zeta_n[j,i] + zeta_np1[j,i]);
    }

  // From SM09, Eq. 2.13
    forall (k,j,i) in D.rho_3D {
      H_np1h[k,j,i] = H0[k,j,i] * (1 + zeta_np1h[j,i] / h[j,i]);
    }

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      H_np1h.updateFluff();
    }
    allLocalesBarrier.barrier();
}
