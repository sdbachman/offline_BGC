use Math;
use params;
use INPUTS;

proc Cs(ref sigma: [?D] real, P: Params) {

  var C = (1 - cosh(P.theta_s * sigma)) / (cosh(P.theta_s) - 1);
  var C2 = (exp(P.theta_b * C) - 1) / (1 - exp(-P.theta_b));

  return C2;
}

proc get_H0(ref h: [?D] real, P: Params) {

  var k_w : [0..P.Nz] int;

  for idx in k_w.domain {
    k_w[idx] = idx;
  }

  var sigma_w = (k_w - P.Nz) / (1.0*P.Nz);

  var Cs_w = Cs(sigma_w, P);

  var z_w : [0..P.Nz, D.dim[0], D.dim[1]] real;
  var H0 : [0..<P.Nz, D.dim[0], D.dim[1]] real;

  forall (k,j,i) in z_w.domain {
    z_w[k,j,i] = h[j,i] * (P.hc*sigma_w[k] + h[j,i]*Cs_w[k]) / (P.hc + h[j,i]);
  }

  forall (k,j,i) in H0.domain {
    H0[k,j,i] = z_w[k+1,j,i] - z_w[k,j,i];
  }

  return H0;
}

