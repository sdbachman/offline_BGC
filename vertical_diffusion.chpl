use params;
use tracers;
use domains;
use INPUTS;
use sigma_coordinate;

use AllLocalesBarriers;

proc calc_vertical_diffusion(ref arr, ref H, D: Domains, P: Params) {

  // This will update tracer_dagger with an implicit timestep for the vertical diffusion

  var Dp : domain(1) = {0..P.Nz};
  var DpDp : domain(2) = {0..P.Nz, 0..P.Nz};

  forall (j,i) in {D.rho_3D.dim[1], D.rho_3D.dim[2]} {

    var tmmp = thomas_diff(j,i, arr, H, P);

    for kk in 0..<P.Nz {
      arr[kk,j,i] = tmmp[kk+1];
    }
  }

  allLocalesBarrier.barrier();

}

proc thomas_diff(j, i, ref arr, ref H, P: Params) {

    var n = P.Nz;
    var Dp : domain(1) = {1..n};

    var a : [Dp] real;
    var b : [Dp] real;
    var c : [Dp] real;
    var d : [Dp] real;

    var cp : [Dp] real;
    var dp : [Dp] real;
    var x  : [Dp] real;

//    b[1] = 1.0;
//    b[n] = 1.0;
//    d[1] = Ts_bot;
//    d[n] = Ts_top;

    a[1] = 0;
    c[1] = -(2*P.dt*kappa_v[1,j,i]) / (H[0,j,i] * (H[1,j,i]+H[0,j,i]));
    b[1] = 1 - a[1] - c[1];
    d[1] = arr[0,j,i];

    a[n] = -(2*P.dt*kappa_v[n-1,j,i]) / (H[n-1,j,i] * (H[n-1,j,i]+H[n-2,j,i]));
    c[n] = 0;
    b[n] = 1 - a[n] - c[n];
    d[n] = arr[n-1,j,i];

    for k in 2..<n {
//      var h0 = H[k-1,j,i];
//      var h1 = H[k,j,i];

//      var alpha = (h1**2) / ((h0 + h1)**2);
//      var beta = (h0**2) / ((h0 + h1)**2);
//      var d1 = 2*(h1**2)*(h1**2 + 2*h0**2 + 3*h0*h1) / ((h0+h1)**4);
//      var d2 = 2*(h0**2)*(h0**2 + 2*h1**2 + 3*h0*h1) / ((h0+h1)**4);

      a[k] = -(2*P.dt*kappa_v[k-1,j,i]) / (H[k-1,j,i] * (H[k-1,j,i]+H[k-2,j,i]));
      c[k] = -(2*P.dt*kappa_v[k,j,i]) / (H[k-1,j,i] * (H[k,j,i]+H[k-1,j,i]));
      b[k] = 1 - a[k] - c[k];

      d[k] = arr[k-1,j,i];
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
