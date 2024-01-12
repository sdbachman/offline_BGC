use FileSystem;

/* INPUT PARAMETERS */

config const gridfile       = '/glade/derecho/scratch/bachman/UCLA-ROMS/Work/SMALL/INPUT/roms_grd.nc';
config var   velocity_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL2/SMALL_his.??????????????.nc';
config var   boundary_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL2/SMALL_bry.??????????????.nc';
config const forcing_files  = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/SMALL2/SMALL_frc.??????????????.nc';

config const Nx_ = 66;
config const Ny_ = 34;
config const Nz_ = 100;

/* Sigma coordinate parameters */
config const theta_s_ : real = 5.0;
config const theta_b_ : real = 2.0;
config const hc_      : real = 300.0;

config const dx_ : real = 4000;
config const dy_ : real = 4000;
const area_ = dx_ * dy_;
const iarea_ = 1.0 / area_;

config const dt_ : real = 60.0;

/* Timestepping */
config var Nt_start_ : int = 0;
config var Nt_ : int = 100;

// For LF-AM3 scheme
config const gamma_ = 0.0833333333333;
config const us_ = 0.16666666666666;

// For AB3 scheme
config const beta_ = 5.0/12.0;

// For RK4 scheme
config const one_sixth_ = 1.0 / 6.0;

// For sponge
config const v_sponge_ : real = 300;
config const sponge_width_ : real = 15;
