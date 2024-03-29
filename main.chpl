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
use RK4;
use INPUTS;
use params;
use domains;
use dynamics;
use tracers;

proc main() {

  var t : stopwatch;
  t.start();

  coforall loc in Locales with (ref H_n) do on loc {

    var P = new owned Params();
    var D = new owned Domains();

    set_domains(D, D3, D_grid);

    initialize_tr(P, D);

    var Dyn = new Dynamics(D);

    update_dynamics(Dyn.U_n, Dyn.V_n, H_n, D, P, P.Nt_start+1);

    // timestepping loop
      for step in (P.Nt_start+1)..(P.Nt_start+P.Nt) {

        // LF-AM3 timestep
          TimeStep(Dyn, D, P, step);

      } // timestepping loop

  } // coforall loop

  t.stop();
  writeln("Program finished in ", t.elapsed(), " seconds.");

} // end program
