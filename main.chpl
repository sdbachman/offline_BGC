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
//use LF-AM3;
use INPUTS;
use domains;


class local_domains {
  var locD : domain(4);
  var locD_static : domain(2);
  var locD_2D : domain(3);
  var locD_u : domain(4);
  var locD_v : domain(4);
  var locD_w : domain(4);
}

proc main(args: [] string) {

  var t : stopwatch;
  t.start();

  // Creating a singleton first dimension to store the time
  const Full = {0..0, 0..<Nz, 0..<Ny, 0..<Nx};
  var myTargetLocales2D = reshape(Locales, {1..1, 1..1, 1..1, 1..Locales.size});

  // Create stencilDist arrays for the tracers. This is necessary to do the halo updates efficiently.
    const sD = new stencilDist(boundingBox=Full, targetLocales=myTargetLocales2D, fluff=(0,0,0,1));
    const D_stencil = sD.createDomain(Full);
    var tracer_old : [D_stencil] real;
    var tracer_new : [D_stencil] real;

  coforall loc in Locales do on loc {

    var infiles = glob('/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL/SMALL_his.??????????????.nc');
    writeln(infiles);
    var (locD, locD_static, locD_2D, locD_u, locD_v, locD_w) = set_domains(D_stencil);


    var h = get_var(gridfile, "h", locD_static);

    // These are the cell thicknesses and volumes for SSH = 0
    var thickness0 = get_thickness0(h);
    var volume0 = thickness0 * dx * dy;

    // These will be updated as the SSH changes
    var thickness : [locD] real;
    var volume : [locD] real;

    // These will track the change in thickness and volume
    var dHdt : [locD] real;
    var dVdt : [locD] real;

    // Read in initial values for zeta_old and zeta_new
    var zeta_old : [locD_2D] real;
    var zeta_new : [locD_2D] real;

    // Create arrays for the volumetric fluxes
      var U : [locD_u] real;
      var V : [locD_v] real;
      var W : [locD_w] real;

    // Create arrays for auxiliary thicknesses (for LF-AM3 predictor step)
      var H_plus : [locD] real;
      var H_minus : [locD] real;

    // Initialize tracer fields
      tracer_old[locD] = get_var(infiles[Nt_start], "temp", locD);
      tracer_new[locD] = get_var(infiles[Nt_start], "temp", locD);

    // Update the halos
      allLocalesBarrier.barrier();
      if (here.id == 0) {
        tracer_old.updateFluff();
        tracer_new.updateFluff();
      }
      allLocalesBarrier.barrier();

      // Copy the part of the tracer array that is on the computational stencil
      // to do the NetCDF output
      WriteOutput("tmp3.nc", tracer_old, "tracer", "stuff", 10);


    // timestepping loop
    for step in (Nt_start)..(Nt_start+Nt) {

      // Read in SSH, update thicknesses
        zeta_new = get_var(infiles[step], "zeta", locD_2D);
        if (step == Nt_start) {
          zeta_old = zeta_new;
        }

      // From SM09, Eq. 2.13
        forall (t,k,j,i) in locD {
          dHdt[t,k,j,i] = thickness0[k,j,i] * (zeta_new[t,j,i] - zeta_old[t,j,i]) / (h[j,i] / dt);
          dVdt[t,k,j,i] = dHdt[t,k,j,i] * dx * dy;

          thickness[t,k,j,i] = thickness0[k,j,i] * (1 + zeta_new[t,j,i] / h[j,i]);
        }

      // Read in u and v
        var u = get_var(infiles[step], "u", locD_u);
        var v = get_var(infiles[step], "v", locD_v);

      // Update volumetric fluxes
        forall (t,k,j,i) in locD_u {
          U[t,k,j,i] = u[t,k,j,i] * 0.5 * (thickness[t,k,j,i] + thickness[t,k,j,i+1]) * dy;
        }

        forall (t,k,j,i) in locD_v {
          V[t,k,j,i] = v[t,k,j,i] * 0.5 * (thickness[t,k,j,i] + thickness[t,k,j+1,i]) * dx;
        }

        // From SH05, Eq. 1.18
        forall (t,k,j,i) in {0..0, 1..Nz, locD_w.dim[1], locD_w.dim[2]} {
          W[t,k,j,i] = -W[t,k-1,j,i] - (dVdt[t,k-1,j,i] + (U[t,k-1,j,i+1] - U[t,k-i,k,i-1]) + (V[t,k-1,j+1,i] - V[t,k-1,j,i]));
        }

        //TimeStep();
     } // Timestepping loop

  } // coforall loop

} // end program
