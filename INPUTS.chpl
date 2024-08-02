use FileSystem;

/* INPUT PARAMETERS */

//config const gridfile       = '/glade/derecho/scratch/bachman/UCLA-ROMS/Work/Iceland1/INPUT/Iceland1_grd.nc';
//config const velocity_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_avg.??????????????.nc';
//config const boundary_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';
//config const forcing_files  = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';

config const gridfile       = '/glade/derecho/scratch/bachman/UCLA-ROMS/Work/Iceland1/INPUT/Iceland1_grd.nc';
config const velocity_files = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/remove_time/INPUT/Iceland1_avg.??????????????.zarr';
//config const velocity_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_avg.??????????????.nc';
config const boundary_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';
config const forcing_files  = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';

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

// For PPM scheme
config const one_third_ = 1.0 / 3.0;

// For sponge
config const v_sponge_ : real = 300;
config const sponge_width_ : real = 15;

// Order of polynomial for boundary value extrapolation
config const ord_ : int = 3;
