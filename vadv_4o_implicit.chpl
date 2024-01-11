use INPUTS;
use dynamics;
use params;
use domains;

use Math;
use AllLocalesBarriers;
use Time;


proc calc_vertical_flux(ref W, ref Dyn: Dynamics, D: Domains, P: Params, ref arr, ref H) {

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

    var tmmp = thomas(t,j,i,Ts_bot, Ts_top, arr, H, P);
//  var tmmp = gaussian(t,j,i,Ts_bot, Ts_top, arr, H, P, Dp, DpDp);

    Dyn.tmp_W[t,..,j,i] = tmmp * W[t,..,j,i];
//    for k in -1..<P.Nz {
//      Dyn.tmp_W[t,k,j,i] = max(W[t,k,j,i], 0.0) * tmmp[k+1] + min(W[t,k,j,i],0.0) * tmmp[k+2];
//    }

  }

  allLocalesBarrier.barrier();

}

proc thomas(t, j, i, Ts_bot, Ts_top, ref arr, ref H, P: Params) {

    var n = P.Nz+1;
    var Dp : domain(1) = {1..n};

    var a : [Dp] real;
    var b : [Dp] real;
    var c : [Dp] real;
    var d : [Dp] real;

    var cp : [Dp] real;
    var dp : [Dp] real;
    var x  : [Dp] real;

    b[1] = 1.0;
    b[n] = 1.0;
    d[1] = Ts_bot;
    d[n] = Ts_top;

    for k in 1..(n-2) {
      var h0 = H[t,k-1,j,i];
      var h1 = H[t,k,j,i];

      var alpha = (h1**2) / ((h0 + h1)**2);
      var beta = (h0**2) / ((h0 + h1)**2);
      var d1 = 2*(h1**2)*(h1**2 + 2*h0**2 + 3*h0*h1) / ((h0+h1)**4);
      var d2 = 2*(h0**2)*(h0**2 + 2*h1**2 + 3*h0*h1) / ((h0+h1)**4);

      a[k+1] = alpha;
      b[k+1] = 1.0;
      c[k+1] = beta;

      d[k+1] = d1*arr[t,k-1,j,i] + d2*arr[t,k,j,i];
    }

    cp[1] = c[1] / b[1];
    dp[1] = d[1] / b[1];

    for k in 2..(n-1) {
      cp[k] = c[k] / (b[k] - a[k]*cp[k-1]);
      dp[k] = (d[k] - a[k]*dp[k-1]) / (b[k] - a[k]*cp[k-1]);
    }

    dp[n] = (d[n] - a[n]*dp[n-1]) / (b[n] - a[n]*cp[n-1]);

    x[n] = dp[n];
    for k in 1..(n-1) by -1 {
      x[k] = dp[k] - cp[k]*x[k+1];
    }

    return x;
}

proc gaussian(t, j, i, Ts_bot, Ts_top, ref arr, ref H, P: Params, Dp, DpDp) {

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

    var tmmp = solve(Mf,Bf);

    return tmmp;
}

