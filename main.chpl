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
use RK3;
use INPUTS;
use params;
use domains;
use dynamics;
use horizontal_diffusion;
use tracers;
use PPM;
use Zarr;

proc main() {

  var t : stopwatch;
  t.start();

  coforall loc in Locales with (ref H_n) do on loc {

    // Create and initialize variables and classes
    var P = new owned Params();
    var D = new owned Domains();

    set_domains(D, D3, D_grid);

    initialize_tr(D, P);
    initialize_sponge(D, P);

    var Dyn = new Dynamics(D);
    var Diff = new Diffusion(D);

    // Load fields for the first timestep
      update_dynamics(Dyn.u_n, Dyn.v_n, Dyn.U_n, Dyn.V_n, H_n, D, P, P.Nt_start);

    // Load fields for the next timestep
      update_thickness(zeta_np1, H_np1, H0, h, D, P, P.Nt_start+1);
      update_dynamics(Dyn.u_np1, Dyn.v_np1, Dyn.U_np1, Dyn.V_np1, H_np1, D, P, P.Nt_start+1);

    // timestepping loop
      for step in (P.Nt_start)..(P.Nt_start+P.Nt) {

        var t1 : stopwatch;
        var t2 : stopwatch;
        var t3 : stopwatch;
        var t4 : stopwatch;

        // Step forward

          t1.start();
          Explicit_TimeStep(Dyn, Diff, D, P, step);
          t1.stop();

          t2.start();
          Implicit_TimeStep(Dyn, Diff, D, P, step);
          t2.stop();

        // Create polynomial fit to current grid
          t3.start();
          Polyfit(D, P);
          t3.stop();

        WriteOutput(tracer_n, "after", "stuff", step);
        allLocalesBarrier.barrier();


        // Update fields to prepare for next time step
          t4.start();
          update_fields(Dyn, Diff, D, P, step);
          t4.stop();

        writeln("Locale ", here.id, " time for explicit step: ", t1.elapsed());
        writeln("Locale ", here.id, " time for implicit step: ", t2.elapsed());
        writeln("Locale ", here.id, " time for polyfit: ", t3.elapsed());
        writeln("Locale ", here.id, " time for reading: ", t4.elapsed());

      } // timestepping loop

  } // coforall loop

  t.stop();
  writeln("Program finished in ", t.elapsed(), " seconds.");

} // end program
