use FileSystem;

/* INPUT PARAMETERS */

config const gridfile = '/glade/derecho/scratch/bachman/UCLA-ROMS/Work/SMALL/INPUT/roms_grd.nc';

config const Nx = 66;
config const Ny = 34;
config const Nz = 100;

config const theta_s = 5.0;
config const theta_b = 2.0;
config const hc = 300;

config const dx = 1000;
config const dy = 1000;

config const dt = 3600;

/* Timestepping */
config var Nt_start : int = 0;
config var Nt : int = 1;

// For LF-AM3 scheme
  config const gamma = 0.0833333333333;
