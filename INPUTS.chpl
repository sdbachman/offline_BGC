use FileSystem;

/* INPUT PARAMETERS */

config const gridfile = '/glade/derecho/scratch/bachman/UCLA-ROMS/Work/SMALL/INPUT/roms_grd.nc';
config var velfiles = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL/SMALL_avg.??????????????.nc';
config var bryfiles = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL/SMALL_bry.??????????????.nc';
config const frcfiles = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL/SMALL_frc.??????????????.nc';

config const Nx = 66;
config const Ny = 34;
config const Nz = 100;

config const theta_s = 5.0;
config const theta_b = 2.0;
config const hc = 300;

config const dx = 4000;
config const dy = 4000;
const area = dx * dy;
const iarea = 1.0 / area;

config const dt = 60;

/* Timestepping */
config var Nt_start : int = 1;
config var Nt : int = 0;

// For LF-AM3 scheme
  config const gamma = 0.0833333333333;
  config const us = 0.16666666666666;
