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
use diffusion;
use tracers;
use PPM;

proc main() {

  var t : stopwatch;
  t.start();

  coforall loc in Locales with (ref H_n) do on loc {

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

        // Step forward
          TimeStep(Dyn, Diff, D, P, step);

        // Create polynomial fit to current grid
          Polyfit(D, P);

        // WriteOutput(tracer_n, "after", "stuff", step);
        // allLocalesBarrier.barrier();

        // Update fields to prepare for next time step
          update_fields(Dyn, Diff, D, P, step);

      } // timestepping loop

  } // coforall loop

  t.stop();
  writeln("Program finished in ", t.elapsed(), " seconds.");

} // end program
