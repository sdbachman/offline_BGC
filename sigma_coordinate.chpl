use Math;
use INPUTS;

proc Cs(ref sigma: [?D] real) {

  var C = (1 - cosh(theta_s * sigma)) / (cosh(theta_s) - 1);
  var C2 = (exp(theta_b * C) - 1) / (1 - exp(-theta_b));

  return C2;
}

proc get_thickness0(ref h: [?D] real) {

  var k_w : [0..Nz] int;

  for idx in k_w.domain {
    k_w[idx] = idx;
  }

  var sigma_w = (k_w - Nz) / (1.0*Nz);

  var Cs_w = Cs(sigma_w);

  var z_w : [0..Nz, D.dim[0], D.dim[1]] real;
  var thickness0 : [0..<Nz, D.dim[0], D.dim[1]] real;

  forall (k,j,i) in z_w.domain {
    z_w[k,j,i] = h[j,i] * (hc*sigma_w[k] + h[j,i]*Cs_w[k]) / (hc + h[j,i]);
  }

  forall (k,j,i) in thickness0.domain {
    thickness0[k,j,i] = z_w[k+1,j,i] - z_w[k,j,i];
  }

  return thickness0;
}

