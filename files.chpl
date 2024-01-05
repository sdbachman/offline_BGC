use INPUTS;
use FileSystem;

class Files {

  var grd : string;

  var vf : domain(1);
  var vel : [vf] string;

  var bf : domain(1);
  var bry : [bf] string;

  var ff : domain(1);
  var frc : [ff] string;

  proc init() {

    this.grd = gridfile;

    var tmp = glob(velfiles);
    this.vf = tmp.domain;
    this.vel = tmp;

    var tmp2 = glob(bryfiles);
    this.bf = tmp2.domain;
    this.bry = tmp2;

//    tmp = glob(frcfiles);
//    this.ff = tmp.domain;
//    this.frc = tmp;

  }

}
