use IO;
use BlockDist;
use StencilDist;
use Time;
use AutoMath;
use LinearAlgebra;
use IO.FormattedIO;
use Math;
use AllLocalesBarriers;
use FileSystem;

use NetCDF_IO;
use sigma_coordinate;
use LF_AM3;
use INPUTS;
use domains;
use arrays;
use tracers;

proc main(args: [] string) {

  var t : stopwatch;
  t.start();

    var Tr = new owned Tracers();

  coforall loc in Locales with (ref Tr) do on loc {

    var infiles = glob('/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL/SMALL_avg.??????????????.nc');

    var locD = new owned Domains();
    set_domains(locD, Tr.D);

    var locA = new Arrays(locD);
    set_static_arrays(locA, locD);

    initialize(Tr, infiles[Nt_start], locD);

    update_halos(Tr);

    WriteOutput("initial.nc", Tr.tracer_old, "tracer", "stuff", 10);

    // timestepping loop
    for step in (Nt_start)..(Nt_start+Nt) {

      update_dynamic_arrays(locA, locD, infiles, step);

      TimeStep(locA, locD, Tr);

    } // timestepping loop

  } // coforall loop

} // end program
