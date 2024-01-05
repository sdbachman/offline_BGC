use INPUTS;
use FileSystem;

class Params {

  var grdfile : string;

  var vf : domain(1);
  var bf : domain(1);
  var ff : domain(1);

  var velfiles : [vf] string;
  var bryfiles : [bf] string;
  var frcfiles : [ff] string;

  var Nx : int;
  var Ny : int;
  var Nz : int;

  /* Sigma coordinate parameters */
  var theta_s : real;
  var theta_b : real;
  var hc      : real;

  var dx : real;
  var dy : real;
  var area : real;
  var iarea : real;

  var dt : real;

  /* Timestepping */
  var Nt_start : int;
  var Nt : int;

  // For LF-AM3 scheme
  var gamma : real;
  var us : real;

  // For AB3 scheme
  var beta : real;

  // For RK4 scheme
  var one_sixth : real;



  proc init() {

    this.grdfile = gridfile;

    var tmp = glob(velocity_files);
    var tmp2 = glob(boundary_files);
    var tmp3 = glob(forcing_files);

    this.vf = tmp.domain;
    this.bf = tmp2.domain;
    this.ff = tmp3.domain;

    this.velfiles = tmp;
    this.bryfiles = tmp2;
    this.frcfiles = tmp3;

    this.Nx = Nx_;
    this.Ny = Ny_;
    this.Nz = Nz_;

    /* Sigma coordinate parameters */
    this.theta_s = theta_s_;
    this.theta_b = theta_b_;
    this.hc      = hc_;

    this.dx = dx_;
    this.dy = dy_;
    this.area = area_;
    this.iarea = iarea_;

    this.dt = dt_;

    this.Nt_start = Nt_start_;
    this.Nt = Nt_;

    // For LF-AM3 scheme
    this.gamma = gamma_;
    this.us = us_;

    // For AB3 scheme
    this.beta = beta_;

    // For RK4 scheme
    this.one_sixth = one_sixth_;

  }

}
