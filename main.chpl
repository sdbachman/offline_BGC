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
use dynamics;
use tracers;

proc main() {

  var t : stopwatch;
  t.start();

  var Tr = new owned Tracers();

  coforall loc in Locales with (ref Tr) do on loc {

    var locF = new owned Files();
    var locD = new owned Domains();

    set_domains(locD, Tr.D3, Tr.D_grid);

    var locA = new Dynamics(locD);

    initialize(Tr, locF, locF.grd, locD);

    update_dynamics(locA.U_n, locA.V_n, Tr.H_n, locD, locF, Nt_start+1);

    // timestepping loop
      for step in (Nt_start+1)..(Nt_start+Nt) {

        // LF-AM3 timestep
          TimeStep(locA, locD, Tr, locF, step);

      } // timestepping loop

  } // coforall loop

} // end program
