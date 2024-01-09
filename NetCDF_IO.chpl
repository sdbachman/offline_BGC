use INPUTS;
use params;

use NetCDF.C_NetCDF;
use CTypes;
use AllLocalesBarriers;

proc get_shape (filename : string, varName : string) {

  //var filename = args[1];
  //var varName = args[2];

  var ncid : c_int;
  var varid : c_int;
  var ndims : c_int;
  var dimid: c_int;

  // Open the file
  // (1)  int nc_open(const char* path, int mode,     int* ncidp)
  extern proc nc_open(path : c_ptrConst(c_char), mode : c_int, ncidp : c_ptr(c_int)) : c_int;
  nc_open(filename.c_str(), NC_NOWRITE, c_ptrTo(ncid));

  // Get the variable ID
  //
  //      int nc_inq_varid(int ncid,    const char* name,      int* varidp)
  extern proc nc_inq_varid(ncid: c_int, varName: c_ptrConst(c_char), varid: c_ptr(c_int));
  nc_inq_varid(ncid, varName.c_str(), c_ptrTo(varid));

  //writeln("varid: ", varid);


  // Get the number of dimensions for this variable
  //
  //      int nc_inq_varndims(int ncid,    int varid,    int* ndimsp)
  extern proc nc_inq_varndims(ncid: c_int, varid: c_int, ndims: c_ptr(c_int));
  nc_inq_varndims(ncid, varid, c_ptrTo(ndims));

  //writeln("ndims: ", ndims);

  var dimids : [0..#ndims] c_int;

  // Get the IDs of each dimension
  //
  //      int nc_inq_vardimid(int ncid,     int varid,     int* dimidsp)
  extern proc nc_inq_vardimid(ncid : c_int, varid : c_int, dimidsp : c_ptr(c_int)) : c_int;

  nc_inq_vardimid(ncid, varid, c_ptrTo(dimids));

  //writeln("dimids: ", dimids);

  var dimlens : [0..#ndims] c_size_t;

  // Get the size of each dimension
  //
  //      int nc_inq_dimlen(int ncid,     int dimid,     size_t* lenp)
  extern proc nc_inq_dimlen(ncid : c_int, dimid : c_int, lenp : c_ptr(c_size_t)) : c_int;
  for i in 0..#ndims do {
    nc_inq_dimlen(ncid, dimids[i], c_ptrTo(dimlens[i]));
  }

  //writeln("dimlens: ", dimlens);
  //writeln("dimlens size: ", dimlens.size);

  return dimlens;
}


proc get_var(filename : string, varName : string, dom_in) {

      var in_array : [dom_in] real(64);

        /* Some external procedure declarations */
          extern proc nc_get_vara_double(ncid : c_int, varid : c_int, startp : c_ptr(c_size_t), countp : c_ptr(c_size_t), ip : c_ptr(real(64))) : c_int;
          extern proc nc_open(path : c_ptrConst(c_char), mode : c_int, ncidp : c_ptr(c_int)) : c_int;
          extern proc nc_inq_varid(ncid: c_int, varName: c_ptrConst(c_char), varid: c_ptr(c_int));

        /* Determine where to start reading file, and how many elements to read */
          // Start specifies a hyperslab.  It expects an array of dimension sizes
          var start = tuplify(dom_in.localSubdomain().first);
          // Count specifies a hyperslab.  It expects an array of dimension sizes
          var count = tuplify(dom_in.localSubdomain().shape);

        /* Create arrays of c_size_t for compatibility with NetCDF-C functions. */
          var start_c = [i in 0..#start.size] start[i] : c_size_t;
          var count_c = [i in 0..#count.size] count[i] : c_size_t;

          var ncid : c_int;
          var varid : c_int;

        /* Open the file */
          nc_open(filename.c_str(), NC_NOWRITE, ncid);

        /* Get the variable ID */
          nc_inq_varid(ncid, varName.c_str(), c_ptrTo(varid));


          nc_get_vara_double(ncid, varid, c_ptrTo(start_c), c_ptrTo(count_c), c_ptrTo(in_array[start]));

          nc_close(ncid);

      return in_array;

}


proc set_bry(P: Params, filename : string, varName : string, ref arr, dom_in) {

        /* Some external procedure declarations */
          extern proc nc_get_vara_double(ncid : c_int, varid : c_int, startp : c_ptr(c_size_t), countp : c_ptr(c_size_t), ip : c_ptr(real(64))) : c_int;
          extern proc nc_open(path : c_ptrConst(c_char), mode : c_int, ncidp : c_ptr(c_int)) : c_int;
          extern proc nc_inq_varid(ncid: c_int, varName: c_ptrConst(c_char), varid: c_ptr(c_int));

        /* Determine where to start reading file, and how many elements to read */

        if (here.id == 0) {

          // Start specifies a hyperslab.  It expects an array of dimension sizes
          var first = dom_in.localSubdomain().first;
          var last  = dom_in.localSubdomain().last;

          // Count specifies a hyperslab.  It expects an array of dimension sizes
          var shape = dom_in.localSubdomain().shape;

          var tmp_west : [first[1]..last[1], first[2]..last[2]] real;
          var tmp_north : [first[1]..last[1], first[3]..last[3]] real;
          var tmp_south : [first[1]..last[1], first[3]..last[3]] real;

          var start_w = (first[1], first[2]);
          var count_w = (shape[1], shape[2]);

          var start_n = (first[1], first[3]);
          var count_n = (shape[1], shape[3]);

          var start_s = (first[1], first[3]);
          var count_s = (shape[1], shape[3]);

        /* Create arrays of c_size_t for compatibility with NetCDF-C functions. */
          var start_w_c = [i in 0..#start_w.size] start_w[i] : c_size_t;
          var count_w_c = [i in 0..#count_w.size] count_w[i] : c_size_t;

          var start_n_c = [i in 0..#start_n.size] start_n[i] : c_size_t;
          var count_n_c = [i in 0..#count_n.size] count_n[i] : c_size_t;

          var start_s_c = [i in 0..#start_s.size] start_s[i] : c_size_t;
          var count_s_c = [i in 0..#count_s.size] count_s[i] : c_size_t;

          var ncid : c_int;
          var varid_w : c_int;
          var varid_n : c_int;
          var varid_s : c_int;

        /* Open the file */
          nc_open(filename.c_str(), NC_NOWRITE, ncid);

        /* Get the variable ID */
          var vw = varName + "_west";
          nc_inq_varid(ncid, vw.c_str(), c_ptrTo(varid_w));

          var vn = varName + "_north";
          nc_inq_varid(ncid, vn.c_str(), c_ptrTo(varid_n));

          var vs = varName + "_south";
          nc_inq_varid(ncid, vs.c_str(), c_ptrTo(varid_s));

          nc_get_vara_double(ncid, varid_w, c_ptrTo(start_w_c), c_ptrTo(count_w_c), c_ptrTo(tmp_west[start_w]));
          arr[0,0..<P.Nz,..,0] = tmp_west;

          nc_get_vara_double(ncid, varid_n, c_ptrTo(start_n_c), c_ptrTo(count_n_c), c_ptrTo(tmp_north[start_n]));
          arr[0,0..<P.Nz,last[2],first[3]..last[3]] = tmp_north;

          nc_get_vara_double(ncid, varid_s, c_ptrTo(start_s_c), c_ptrTo(count_s_c), c_ptrTo(tmp_south[start_s]));
          arr[0,0..<P.Nz,0,first[3]..last[3]] = tmp_south;

          nc_close(ncid);

         }
         else if (here.id == (Locales.size-1)) {

          // Start specifies a hyperslab.  It expects an array of dimension sizes
          var first = dom_in.localSubdomain().first;
          var last  = dom_in.localSubdomain().last;

          // Count specifies a hyperslab.  It expects an array of dimension sizes
          var shape = dom_in.localSubdomain().shape;

          var tmp_east : [first[1]..last[1], first[2]..last[2]] real;
          var tmp_north : [first[1]..last[1], first[3]..last[3]] real;
          var tmp_south : [first[1]..last[1], first[3]..last[3]] real;

          var start_e = (first[1], first[2]);
          var count_e = (shape[1], shape[2]);

          var start_n = (first[1], first[3]);
          var count_n = (shape[1], shape[3]);

          var start_s = (first[1], first[3]);
          var count_s = (shape[1], shape[3]);

        /* Create arrays of c_size_t for compatibility with NetCDF-C functions. */
          var start_e_c = [i in 0..#start_e.size] start_e[i] : c_size_t;
          var count_e_c = [i in 0..#count_e.size] count_e[i] : c_size_t;

          var start_n_c = [i in 0..#start_n.size] start_n[i] : c_size_t;
          var count_n_c = [i in 0..#count_n.size] count_n[i] : c_size_t;

          var start_s_c = [i in 0..#start_s.size] start_s[i] : c_size_t;
          var count_s_c = [i in 0..#count_s.size] count_s[i] : c_size_t;

          var ncid : c_int;
          var varid_e : c_int;
          var varid_n : c_int;
          var varid_s : c_int;

        /* Open the file */
          nc_open(filename.c_str(), NC_NOWRITE, ncid);

        /* Get the variable ID */
          var ve = varName + "_east";
          nc_inq_varid(ncid, ve.c_str(), c_ptrTo(varid_e));

          var vn = varName + "_north";
          nc_inq_varid(ncid, vn.c_str(), c_ptrTo(varid_n));

          var vs = varName + "_south";
          nc_inq_varid(ncid, vs.c_str(), c_ptrTo(varid_s));

          nc_get_vara_double(ncid, varid_e, c_ptrTo(start_e_c), c_ptrTo(count_e_c), c_ptrTo(tmp_east[start_e]));
          arr[0,0..<P.Nz,..,last[3]] = tmp_east;

          nc_get_vara_double(ncid, varid_n, c_ptrTo(start_n_c), c_ptrTo(count_n_c), c_ptrTo(tmp_north[start_n]));
          arr[0,0..<P.Nz,last[2],first[3]..last[3]] = tmp_north;

          nc_get_vara_double(ncid, varid_s, c_ptrTo(start_s_c), c_ptrTo(count_s_c), c_ptrTo(tmp_south[start_s]));
          arr[0,0..<P.Nz,0,first[3]..last[3]] = tmp_south;

          nc_close(ncid);

         }
         else {
          // Start specifies a hyperslab.  It expects an array of dimension sizes
          var first = dom_in.localSubdomain().first;
          var last  = dom_in.localSubdomain().last;

          // Count specifies a hyperslab.  It expects an array of dimension sizes
          var shape = dom_in.localSubdomain().shape;

          var tmp_north : [first[1]..last[1], first[3]..last[3]] real;
          var tmp_south : [first[1]..last[1], first[3]..last[3]] real;

          var start_n = (first[1], first[3]);
          var count_n = (shape[1], shape[3]);

          var start_s = (first[1], first[3]);
          var count_s = (shape[1], shape[3]);

        /* Create arrays of c_size_t for compatibility with NetCDF-C functions. */

          var start_n_c = [i in 0..#start_n.size] start_n[i] : c_size_t;
          var count_n_c = [i in 0..#count_n.size] count_n[i] : c_size_t;

          var start_s_c = [i in 0..#start_s.size] start_s[i] : c_size_t;
          var count_s_c = [i in 0..#count_s.size] count_s[i] : c_size_t;

          var ncid : c_int;
          var varid_n : c_int;
          var varid_s : c_int;

        /* Open the file */
          nc_open(filename.c_str(), NC_NOWRITE, ncid);

        /* Get the variable ID */

          var vn = varName + "_north";
          nc_inq_varid(ncid, vn.c_str(), c_ptrTo(varid_n));

          var vs = varName + "_south";
          nc_inq_varid(ncid, vs.c_str(), c_ptrTo(varid_s));

          nc_get_vara_double(ncid, varid_n, c_ptrTo(start_n_c), c_ptrTo(count_n_c), c_ptrTo(tmp_north[start_n]));
          arr[0,0..<P.Nz,last[2],first[3]..last[3]] = tmp_north;

          nc_get_vara_double(ncid, varid_s, c_ptrTo(start_s_c), c_ptrTo(count_s_c), c_ptrTo(tmp_south[start_s]));
          arr[0,0..<P.Nz,0,first[3]..last[3]] = tmp_south;

          nc_close(ncid);

         }

  allLocalesBarrier.barrier();

}


var y: atomic int;

proc WriteOutput(ref arr_in: [?D] real, varName : string, units : string, i : int) {

  allLocalesBarrier.barrier();

  /* The timestamp for the filename */
    var currentIter = i : string;
    const maxLen = 10;
    const zero_len = maxLen - currentIter.size;
    const paddedStr = (zero_len * "0") + currentIter;
    var filename = (varName + "." + paddedStr + ".nc");
    //var filename = "tmp.nc";

    if (here.id == 0) {
  /* IDs for the netCDF file, dimensions, and variables. */
    var ncid, varid : c_int;
    var x_dimid, y_dimid, z_dimid, time_dimid : c_int;
    var x_varid, y_varid, z_varid, time_varid : c_int;

    var ndims : int = 4;
    var dimids: [0..#ndims] c_int;

    var shape = arr_in.shape;

    var att_text : string;

    //var current_time = t;

    //writeln("Current time: ", t);
    //var zo = z;
    //var yo = y[..,0];
    //var xo = x[0,..];

    var to : [0..<shape[0]] real;
    var zo : [0..<shape[1]] real;
    var yo : [0..<shape[2]] real;
    var xo : [0..<shape[3]] real;

    //for i in 0..<shape[0] {
      to[0] = dt_ * i;
    //}
    for ii in 0..<shape[1] {
      zo[ii] = ii;
    }
    for ii in 0..<shape[2] {
      yo[ii] = ii;
    }
    for ii in 0..<shape[3] {
      xo[ii] = ii;
    }

    var timeName = "time";
    var zName = "z";
    var yName = "y";
    var xName = "x";

  /* Create the file. */
    extern proc nc_create(path : c_ptrConst(c_char), cmode : c_int, ncidp : c_ptr(c_int)) : c_int;
    nc_create( filename.c_str(), NC_CLOBBER, c_ptrTo(ncid));

  /* Define the dimensions. The record dimension is defined to have
     unlimited length - it can grow as needed. In this example it is
     the time dimension.*/
    extern proc nc_def_dim(ncid : c_int, name : c_ptrConst(c_char), len : c_size_t, idp : c_ptr(c_int)) : c_int;
//    nc_def_dim(ncid, timeName.c_str(), shape[0] : c_size_t, time_dimid);
    nc_def_dim(ncid, timeName.c_str(), NC_UNLIMITED, time_dimid);
    nc_def_dim(ncid, zName.c_str(), shape[1] : c_size_t, z_dimid);
    nc_def_dim(ncid, yName.c_str(), shape[2] : c_size_t, y_dimid);
    nc_def_dim(ncid, xName.c_str(), shape[3] : c_size_t, x_dimid);

  /* Define the coordinate variables. */
    extern proc nc_def_var(ncid : c_int, name : c_ptrConst(c_char), xtype : nc_type, ndims : c_int, dimidsp : c_ptr(c_int), varidp : c_ptr(c_int)) : c_int;
    nc_def_var(ncid, timeName.c_str(), NC_DOUBLE, 1 : c_int, c_ptrTo(time_dimid), c_ptrTo(time_varid));
    nc_def_var(ncid, zName.c_str(), NC_DOUBLE, 1 : c_int, c_ptrTo(z_dimid), c_ptrTo(z_varid));
    nc_def_var(ncid, yName.c_str(), NC_DOUBLE, 1 : c_int, c_ptrTo(y_dimid), c_ptrTo(y_varid));
    nc_def_var(ncid, xName.c_str(), NC_DOUBLE, 1 : c_int, c_ptrTo(x_dimid), c_ptrTo(x_varid));

  /* Assign units attributes to coordinate variables. */
    extern proc nc_put_att_text(ncid : c_int, varid : c_int, name : c_ptrConst(c_char), len : c_size_t, op : c_ptrConst(c_char)) : c_int;
    att_text = "timesteps";
    nc_put_att_text(ncid, time_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());
    att_text = "meters";
    nc_put_att_text(ncid, z_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());
    att_text = "meters";
    nc_put_att_text(ncid, y_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());
    att_text = "meters";
    nc_put_att_text(ncid, x_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());

  /* The dimids array is used to pass the dimids of the dimensions of
     the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
    dimids[0] = time_dimid;
    dimids[1] = z_dimid;
    dimids[2] = y_dimid;
    dimids[3] = x_dimid;

  /* Define the netCDF variable. */
    nc_def_var(ncid, varName.c_str(), NC_DOUBLE, ndims : c_int, c_ptrTo(dimids[0]), c_ptrTo(varid));

  /* Assign units attributes to the netCDF variables. */
    nc_put_att_text(ncid, varid, "units".c_str(), units.numBytes : c_size_t, units.c_str());

  /* End define mode. */
    nc_enddef(ncid);

  /* Write the coordinate variable data. */
    extern proc nc_put_var_double(ncid : c_int, varid : c_int, op : c_ptr(c_double)) : c_int;
    nc_put_var_double(ncid, time_varid, c_ptrTo(to[0]));
    nc_put_var_double(ncid, z_varid, c_ptrTo(zo[0]));
    nc_put_var_double(ncid, y_varid, c_ptrTo(yo[0]));
    nc_put_var_double(ncid, x_varid, c_ptrTo(xo[0]));

    nc_close(ncid); // comment this and nc_open out to work on 1 node
  }
  /////////////////////////////////

  allLocalesBarrier.barrier();

  y.waitFor(here.id%numLocales);

  var ncid, varid : c_int;

  /* Declaration of nc_open */
    extern proc nc_open(path : c_ptrConst(c_char), mode : c_int, ncidp : c_ptr(c_int)) : c_int;
  // Open the file
    nc_open( filename.c_str() , NC_WRITE, c_ptrTo(ncid));

  /* Create arrays of c_size_t for compatibility with NetCDF-C functions. */
  /* Determine where to start reading file, and how many elements to read */
    // Start specifies a hyperslab.  It expects an array of dimension sizes
      var start = tuplify(D.localSubdomain().first);
    // Count specifies a hyperslab.  It expects an array of dimension sizes
      var count = tuplify(D.localSubdomain().shape);

  /* Adding an extra first element to account for the "time" dimension. */
    var start_c : [0..<start.size] c_size_t;
    var count_c : [0..<count.size] c_size_t;

    for i in 0..<start.size {
      start_c[i] = start[i] : c_size_t;
      count_c[i] = count[i] : c_size_t;
    }

    // Need to create a local copy of the input array here that *only* has the part
    // of the array from the localSubdomain held in memory. The halo row of the
    // StencilDist array in contigous in memory so it will not write correctly to the NetCDF file
    var tmp_arr = arr_in[D.localSubdomain()];

    extern proc nc_inq_varid(ncid: c_int, varName: c_ptrConst(c_char), varid: c_ptr(c_int));
    nc_inq_varid(ncid, varName.c_str(), c_ptrTo(varid));

    extern proc nc_put_vara_double(ncid : c_int, varid : c_int, startp : c_ptr(c_size_t), countp : c_ptr(c_size_t), op : c_ptr(c_double)) : c_int;
    nc_put_vara_double(ncid, varid, c_ptrTo(start_c), c_ptrTo(count_c), c_ptrTo(tmp_arr[start]));

    nc_close(ncid);

    const inc = (y.read() + 1) % numLocales;
    y.write(inc);

//    allLocalesBarrier.barrier();
}

inline proc tuplify(x) {
  if isTuple(x) then return x; else return (x,);
}
