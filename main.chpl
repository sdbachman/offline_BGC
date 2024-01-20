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

    update_dynamics(Dyn.u_n, Dyn.v_n, Dyn.U_n, Dyn.V_n, H_n, D, P, P.Nt_start);

    // timestepping loop
      for step in (P.Nt_start)..(P.Nt_start+P.Nt) {

        // Step forward
          TimeStep(Dyn, Diff, D, P, step);

        // Create polynomial fit to current grid
//          Polyfit(D, P, tracer_dagger);

      } // timestepping loop

  } // coforall loop

  t.stop();
  writeln("Program finished in ", t.elapsed(), " seconds.");

} // end program
