use INPUTS;
use dynamics;
use params;
use domains;


proc calc_fluxes(ref U, ref V, ref W, ref Dyn: Dynamics, D: Domains, P: Params, ref arr, ref H) {

  // Horizontal: upstream-biased parabolic interpolation: SM05, after 4.13

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

  //////////////////////////////////////////////////////////////////////////////////////
  //                                     W-fluxes                                     //
  //  Calculated with Implicit Fourth-order scheme:  White and Adcroft, 2008, Eq. 46  //
  //////////////////////////////////////////////////////////////////////////////////////

  var Dp : domain(1) = {0..P.Nz};
  var DpDp : domain(2) = {0..P.Nz, 0..P.Nz};

  forall (t,j,i) in {0..0,D.rho_3D.dim[2], D.rho_3D.dim[3]} with (ref Dyn) {

    // Bottom boundary extrapolation
    // Calculated using a cubic polynomial over the bottom four cells

    var B : [0..3] real;
    var M : [0..3,0..3] real;

    var h_b : real = 0;
    var h_t : real = H[t,0,j,i];
    var iH = 1.0 / (h_t - h_b);

    for k in 0..3 {
      M[k,0] = iH*(h_t - h_b);
      M[k,1] = 0.5*iH*(h_t**2 - h_b**2);
      M[k,2] = (1.0/3.0)*iH*(h_t**3 - h_b**3);
      M[k,3] = 0.25 * iH*(h_t**4 - h_b**4);

      B[k] = arr[t,k,j,i];

      h_b = h_b + H[t,k,j,i];
      h_t = h_t + H[t,k+1,j,i];
      iH = 1.0/(h_t - h_b);
    }

    var coeffs = solve(M,B);
    var Ts_bot = coeffs[0];

    // Top boundary extrapolation
    h_b = 0;
    h_t = H[t,P.Nz-1,j,i];
    iH = 1.0/(h_t - h_b);
    for k in 0..3 {
      M[k,0] = iH*(h_t - h_b);
      M[k,1] = 0.5*iH*(h_t**2 - h_b**2);
      M[k,2] = (1.0/3.0)*iH*(h_t**3 - h_b**3);
      M[k,3] = 0.25 * iH*(h_t**4 - h_b**4);

      B[k] = arr[t,P.Nz-1-k,j,i];

      h_b = h_b + H[t,P.Nz-1-k,j,i];
      h_t = h_t + H[t,P.Nz-2-k,j,i];
      iH = 1.0/(h_t - h_b);
    }

    coeffs = solve(M,B) ;
    var Ts_top = coeffs[0];

    /////////////////////////////

    var Bf : [Dp] real;
    var Mf : [DpDp] real;

    Mf[0,0] = 1.0;
    Mf[P.Nz,P.Nz] = 1.0;
    Bf[0] = Ts_bot;
    Bf[P.Nz] = Ts_top;

    for k in 0..(P.Nz-2) {
      var h0 = H[t,k,j,i];
      var h1 = H[t,k+1,j,i];

      var alpha = (h1**2) / ((h0 + h1)**2);
      var beta = (h0**2) / ((h0 + h1)**2);
      var b = 2*(h1**2)*(h1**2 + 2*h0**2 + 3*h0*h1) / ((h0+h1)**4);
      var c = 2*(h0**2)*(h0**2 + 2*h1**2 + 3*h0*h1) / ((h0+h1)**4);

      Mf[k+1,k] = alpha;
      Mf[k+1,k+1] = 1.0;
      Mf[k+1,k+2] = beta;

      Bf[k+1] = b*arr[t,k,j,i] + c*arr[t,k+1,j,i];
    }

    Dyn.tmp_W[t,..,j,i] = solve(Mf,Bf) * W[t,..,j,i];

  }

}

