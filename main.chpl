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


proc main(args: [] string) {

  var t : stopwatch;
  t.start();

  // Creating a singleton first dimension to store the time
  const Full = {0..0, 0..<Nz, 0..<Ny, 0..<Nx};
  var myTargetLocales2D = reshape(Locales, {1..1, 1..1, 1..1, 1..Locales.size});

  // Create stencilDist arrays for the tracers. This is necessary to do the halo updates efficiently.
    const sD = new stencilDist(boundingBox=Full, targetLocales=myTargetLocales2D, fluff=(0,0,0,2));
    const D_stencil = sD.createDomain(Full);
    var tracer_old : [D_stencil] real;
    var tracer_new : [D_stencil] real;
    var predictor  : [D_stencil] real;

  coforall loc in Locales with (ref tracer_old) do on loc {

    var infiles = glob('/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL/SMALL_avg.??????????????.nc');

    var locD = new owned Domains();
    set_domains(locD, D_stencil);

    //writeln(here.id, ": ", locD.rho_3D, " ", locD.u_3D);
    //allLocalesBarrier.barrier();
    //exit();


    var locA = new Arrays(locD);
    set_static_arrays(locA, locD);

    // Initialize tracer fields
    tracer_old[locD.rho_3D] = get_var(infiles[Nt_start], "temp", locD.rho_3D);
    tracer_new[locD.rho_3D] = get_var(infiles[Nt_start], "temp", locD.rho_3D);

    // Update the halos
    allLocalesBarrier.barrier();
    if (here.id == 0) {
      tracer_old.updateFluff();
      tracer_new.updateFluff();
    }
    allLocalesBarrier.barrier();

    // Copy the part of the tracer array that is on the computational stencil
    // to do the NetCDF output
 //   WriteOutput("tmp3.nc", tracer_old, "tracer", "stuff", 10);

    // timestepping loop
    for step in (Nt_start)..(Nt_start+Nt) {

      update_dynamic_arrays(locA, locD, infiles, step);

      TimeStep(locA, locD, tracer_new, tracer_old, predictor);

     } // Timestepping loop

  } // coforall loop

} // end program
