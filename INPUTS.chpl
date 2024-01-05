use FileSystem;

/* INPUT PARAMETERS */

config const gridfile = '/glade/derecho/scratch/bachman/UCLA-ROMS/Work/SMALL/INPUT/roms_grd.nc';
config var velfiles = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL2/SMALL_his.??????????????.nc';
config var bryfiles = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL2/SMALL_bry.??????????????.nc';
config const frcfiles = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL2/SMALL_frc.??????????????.nc';

config const Nx = 66;
config const Ny = 34;
config const Nz = 100;

config const theta_s : real = 5.0;
config const theta_b : real = 2.0;
config const hc      : real = 300.0;

config const dx : real = 4000;
config const dy : real = 4000;
const area = dx * dy;
const iarea = 1.0 / area;

config const dt : real = 60.0;

/* Timestepping */
config var Nt_start : int = 0;
config var Nt : int = 100;

// For LF-AM3 scheme
config const gamma = 0.0833333333333;
config const us = 0.16666666666666;

// For AB3 scheme
config const beta = 5.0/12.0;

// For RK4 scheme
config const one_sixth = 1.0 / 6.0;
