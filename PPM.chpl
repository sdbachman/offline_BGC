use INPUTS;
use dynamics;
use domains;
use tracers;
use params;
use NetCDF_IO;
//use utils;

use Math;
use AllLocalesBarriers;
use Time;


// This function will apply the piecewise parabolic method (PPM) as described in
// White and Adcroft (2008).

proc Polyfit(D: Domains, P: Params) {

  forall (t,j,i) in {0..0,D.rho_3D.dim[2], D.rho_3D.dim[3]} {

  if (mask_rho[j,i] == 1) {

    //writeln("On point (", j, ", ", i, ")");

    var interface_values = calc_interface_values(t, j, i, D, P, tracer_dagger, H_dagger);

    // Create the coefficients for the polynomial in each cell
    // Going to treat "left" as equivalent to "bottom", and "right" as equivalent to "top"

    var a0 : [0..<P.Nz] real;
    var a1 : [0..<P.Nz] real;
    var a2 : [0..<P.Nz] real;
    for k in 0..<P.Nz {
      a0[k] = interface_values[k];
      a1[k] = 6*tracer_dagger[t,k,j,i] - 4*interface_values[k] - 2*interface_values[k+1];
      a2[k] = 3*(interface_values[k] + interface_values[k+1] - 2*tracer_dagger[t,k,j,i]);
    }

    // These are vector copies of each column.  H_orig is going to have one extra layer
    // to allow the reconstruction loop to exit gracefully.

      var H_orig : [0..<P.Nz] real;
      var H_new  : [0..<P.Nz] real;
      for kk in 0..<P.Nz {
        H_orig[kk] = H_dagger[t,kk,j,i];
        H_new[kk]  = H_np1[t,kk,j,i];
      }

    // Going to normalize the thicknesses to make the forthcoming loop logic and
    // tracer conservation a bit easier. This is just for bookkeeping here and will not
    // affect the actual thicknesses.
      H_new = H_new * (+ reduce H_orig) / (+ reduce H_new);

    // These represent the indices in each column
      var k_orig = 0;
      var curr_k_new  = 0;

    // These represent how much thickness we have left before we exit the current cell
      var curr_H_orig : real = H_orig[0];
      var curr_H_new  : real = H_new[0];

    // These identify our position within the current cell
      var z0 : real = 0;
      var z1 : real = 1;

    // This holds the reconstruction
    var reconstruction : [0..<P.Nz] real;

//    if (i,j) == (24,17) {

    var tol = 1e-5;

    // Going to overwrite tracer_n with the interpolated values
    label ko for k_new in 0..<P.Nz {
//writeln("Starting on new layer ", k_new);
      curr_H_new = H_new[k_new];

      while (k_orig < P.Nz) {

        while (curr_H_new >= tol) {

          curr_H_orig = H_orig[k_orig];

//writeln("Starting on orig layer ", k_orig);

          z1 = min(1, z0 + curr_H_new / curr_H_orig);
          var tmp = integrate(a0[k_orig], a1[k_orig], a2[k_orig], z0, z1, P);
          var frac = (z1 - z0) * curr_H_orig / H_new[k_new];

          reconstruction[k_new] = reconstruction[k_new] + frac*tmp;
//writeln("For reconstruction, z0, z1 = ", z0, " ", z1);
//writeln("a0, a1, a2: ", a0[k_orig], " ", a1[k_orig], " ", a2[k_orig]);
//writeln("curr_H_orig was ", curr_H_orig);
//writeln("Added ", (frac*tmp), " to reconstruction");
          curr_H_new = max(0, curr_H_new - (z1-z0)*curr_H_orig);
          z0 = 0;
          // Increment k_orig only if z1 = 1
          k_orig += (floor(z1) : int);

        }
        // This will ensure z0 = 0 if z1 = 1; otherwise if z1<1 then z0=z1
        z0 = z1 - floor(z1);
//writeln("Finished new layer ", k_new);
//writeln("Final reconstructed value is ", reconstruction[k_new]);
//writeln("Original value in this layer was ", tracer_n[t,k_new,j,i]);
//writeln();

        continue ko;

      }
    }

/*
    // Adjust values to satisfy conservation.
    var integrated_orig = (+ reduce (H_orig * tracer_dagger[t,..,j,i]) );
    var integrated_new  = (+ reduce (H_new * reconstruction) );
    var int_ratio = integrated_orig / integrated_new;
    reconstruction = reconstruction * int_ratio;
*/

    for kk in 0..<P.Nz {
      tracer_n[t,kk,j,i] = reconstruction[kk];
    }

// } // (i,j)

} // mask_rho

  } // end of forall

} // end of subroutine


proc integrate(a0, a1, a2, z0, z1, P: Params) {

  // definite_integral = a0 * (z1 - z0) + 0.5*a1*(z1**2 - z0**2) + P.one_third*a2*(z1**3 - z0**3);
  // Mean = definite integral / interval = definite integral / (z1 - z0);

  var mean_of_integral = a0 + 0.5*a1*(z1 + z0) + P.one_third*a2*(z1**2 + z1*z0 + z0**2);

  return mean_of_integral;
}

proc calc_interface_values(t, j, i, D: Domains, P: Params, ref arr, ref H) {

  //////////////////////////////////////////////////////////////////////////////////////
  //                   Get tracer values at layer interfaces                          //
  //  Calculated with Implicit Fourth-order scheme:  White and Adcroft, 2008, Eq. 46  //
  //////////////////////////////////////////////////////////////////////////////////////

  var Dp : domain(1) = {0..P.Nz};
  var DpDp : domain(2) = {0..P.Nz, 0..P.Nz};

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

    var interface_values : [0..P.Nz] real = thomas(t,j,i,Ts_bot, Ts_top, arr, H, P);

  return interface_values;

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

