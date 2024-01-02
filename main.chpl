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

  var Tr_next = new owned Tracers();
  var Tr_curr = new owned Tracers();
  var Tr_prev = new owned Tracers();

  coforall loc in Locales with (ref Tr_next, ref Tr_curr, ref Tr_prev) do on loc {

    var locF = new owned Files();
    var locD = new owned Domains();

    set_domains(locD, Tr_curr.D3, Tr_curr.D_grid);

    var locA_curr = new Arrays(locD);
    var locA_next = new Arrays(locD);

    initialize(Tr_next, locF.vel[Nt_start+2], locF.grd, locD);
    initialize(Tr_curr, locF.vel[Nt_start+1], locF.grd, locD);
    initialize(Tr_prev, locF.vel[Nt_start], locF.grd, locD);

    update_thickness(Tr_next, locD, locF, Nt_start+2);
    update_thickness(Tr_curr, locD, locF, Nt_start+1);
    update_thickness(Tr_prev, locD, locF, Nt_start);

    update_dynamic_arrays(locA_curr, Tr_curr, locD, locF, Nt_start);
//    update_dynamic_arrays(locA_prev, Tr_prev, locD, locF, Nt_start);

    // timestepping loop
      for step in (Nt_start+1)..(Nt_start+Nt) {

        // LF-AM3 timestep
          TimeStep(locA_curr, locA_next, locD, Tr_prev, Tr_curr, Tr_next, locF, step);

      } // timestepping loop

  } // coforall loop

} // end program
