use IO;
use BlockDist;
use StencilDist;
use Time;
use AutoMath;
use LinearAlgebra;
use IO.FormattedIO;
use Math;
use AllLocalesBarriers;

use NetCDF_IO;
use LF_AM3;
use INPUTS;
use files;
use domains;
use arrays;
use tracers;

proc main() {

  var t : stopwatch;
  t.start();

  var Tr_next2 = new owned Tracers();
  var Tr_next = new owned Tracers();
  var Tr_curr = new owned Tracers();
  var Tr_prev = new owned Tracers();
  var Tr_prev2 = new owned Tracers();

  coforall loc in Locales with (ref Tr_next2, ref Tr_next, ref Tr_curr, ref Tr_prev, ref Tr_prev2) do on loc {

    var locF = new owned Files();
    var locD = new owned Domains();

    set_domains(locD, Tr_curr.D3, Tr_curr.D_grid);

    var locA_prev2 = new Arrays(locD);
    var locA_prev = new Arrays(locD);
    var locA_curr = new Arrays(locD);
    var locA_next = new Arrays(locD);
    var locA_next2 = new Arrays(locD);

    initialize(Tr_next2, locF.vel[Nt_start+4], locF.grd, locD);
    initialize(Tr_next, locF.vel[Nt_start+3], locF.grd, locD);
    initialize(Tr_curr, locF.vel[Nt_start+2], locF.grd, locD);
    initialize(Tr_prev, locF.vel[Nt_start+1], locF.grd, locD);
    initialize(Tr_prev2, locF.vel[Nt_start], locF.grd, locD);

    update_thickness(Tr_next2, locD, locF, Nt_start+4);
    update_thickness(Tr_next, locD, locF, Nt_start+3);
    update_thickness(Tr_curr, locD, locF, Nt_start+2);
    update_thickness(Tr_prev, locD, locF, Nt_start+1);
    update_thickness(Tr_prev2, locD, locF, Nt_start);

//    update_dynamic_arrays(locA_curr, Tr_curr, locD, locF, Nt_start);
    update_dynamic_arrays(locA_prev2, Tr_prev2, locD, locF, Nt_start);
    update_dynamic_arrays(locA_prev, Tr_prev, locD, locF, Nt_start+1);
    update_dynamic_arrays(locA_curr, Tr_curr, locD, locF, Nt_start+2);
    update_dynamic_arrays(locA_next, Tr_next, locD, locF, Nt_start+3);
    update_dynamic_arrays(locA_next2, Tr_next2, locD, locF, Nt_start+4);

    // timestepping loop
      for step in (Nt_start+2)..(Nt_start+Nt) {

        // LF-AM3 timestep
          TimeStep(locA_prev2, locA_prev, locA_curr, locA_next, locA_next2, locD, Tr_prev2, Tr_prev, Tr_curr, Tr_next, Tr_next2, locF, step);

      } // timestepping loop

  } // coforall loop

} // end program
