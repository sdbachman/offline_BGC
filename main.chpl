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
use params;
use domains;
use dynamics;
use tracers;

proc main() {

  var t : stopwatch;
  t.start();

  var Tr = new owned Tracers();

  coforall loc in Locales with (ref Tr) do on loc {

    var P = new owned Params();
    var D = new owned Domains();

    set_domains(D, Tr.D3, Tr.D_grid);

    var Dyn = new Dynamics(D);

    initialize(Tr, P, D);

    update_dynamics(Dyn.U_n, Dyn.V_n, Tr.H_n, D, P, P.Nt_start+1);

    // timestepping loop
      for step in (P.Nt_start+1)..(P.Nt_start+P.Nt) {

        // LF-AM3 timestep
          TimeStep(Dyn, D, Tr, P, step);

      } // timestepping loop

  } // coforall loop

} // end program
