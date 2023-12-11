use StencilDist;

proc set_domains(const D_rho: ?) {

  var locD = D_rho.localSubdomain();

  var locD_static = {locD.dim[2], locD.dim[3]};

  var locD_2D = {locD.dim[0], locD.dim[2], locD.dim[3]};

  var locD_u = get_locD_u(locD);

  var locD_v = get_locD_v(locD);

  var locD_w = get_locD_w(locD);

  return (locD, locD_static, locD_2D, locD_u, locD_v, locD_w);
}

proc get_locD_u(ref D: domain(4)) {

  var D_u : domain(4);

  if (here.id == 0) {
    D_u = {D.dim[0], D.dim[1], D.dim[2], 0..D.last[3]};
  }
  else if (here.id == (Locales.size-1)) {
    D_u = {D.dim[0], D.dim[1], D.dim[2], (D.first[3]-1)..(D.last[3]-1)};
  }
  else {
    D_u = {D.dim[0], D.dim[1], D.dim[2], (D.first[3]-1)..(D.last[3])};
  }

  return D_u;
}


proc get_locD_v(ref D: domain(4)) {

  var D_v = {D.dim[0], D.dim[1], 0..(D.last[2]-1), D.dim[3]};

  return D_v;
}


proc get_locD_w(ref D: domain(4)) {

  var D_w = {D.dim[0], 0..(D.last[1]+1), D.dim[2], D.dim[3]};

  return D_w;
}

