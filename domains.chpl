use StencilDist;


class Domains {
  var grid : domain(2);
  var rho_2D : domain(3);
  var rho_3D : domain(4);
  var u_3D : domain(4);
  var v_3D : domain(4);
  var w_3D : domain(4);

}

proc set_domains(arg: Domains, const full_domain: ?, const grid_domain: ?) {

  arg.rho_3D = full_domain.localSubdomain();

  arg.grid = {arg.rho_3D.dim[2], arg.rho_3D.dim[3]};

  arg.rho_2D = {arg.rho_3D.dim[0], arg.rho_3D.dim[2], arg.rho_3D.dim[3]};

  arg.u_3D = get_u(arg.rho_3D);

  arg.v_3D = get_v(arg.rho_3D);

  arg.w_3D = get_w(arg.rho_3D);

}

proc get_u(ref D: domain(4)) {

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


proc get_v(ref D: domain(4)) {

  var D_v = {D.dim[0], D.dim[1], 0..(D.last[2]-1), D.dim[3]};

  return D_v;
}


proc get_w(ref D: domain(4)) {

  var D_w = {D.dim[0], -1..D.last[1], D.dim[2], D.dim[3]};

  return D_w;
}

