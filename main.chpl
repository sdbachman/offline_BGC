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

    var Tr = new owned Tracers();

  coforall loc in Locales with (ref Tr) do on loc {

    var locF = new owned Files();

    var locD = new owned Domains();
    set_domains(locD, Tr.D);

    var locA_curr = new Arrays(locD);
    var locA_next = new Arrays(locD);
    set_static_arrays(locA_curr, locD, locF);
    set_static_arrays(locA_next, locD, locF);

    initialize(Tr, locF.vel[Nt_start], locD);

    update_halos(Tr);

    WriteOutput("initial.nc", Tr.tracer_old, "tracer", "stuff", 10);

    // timestepping loop
    for step in (Nt_start)..(Nt_start+Nt) {

      // LF-AM3 timestep
      TimeStep(locA_curr, locA_next, locD, Tr, locF, step);

    } // timestepping loop

  } // coforall loop

} // end program
