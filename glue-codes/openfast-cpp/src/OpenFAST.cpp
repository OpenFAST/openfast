#include "OpenFAST.H"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <array>

inline void check_nc_error(int code, std::string msg) {
    if (code != 0)
        throw std::runtime_error("OpenFAST C++ API:: NetCDF error: " + msg);
}

int fast::OpenFAST::AbortErrLev = ErrID_Fatal; // abort error level; compare with NWTC Library

int time_step_ratio(double fastDt, double driverDt, double epsFactor=1e-6)
{
  // ensure that the ratio is robust to integer conversion by making sure it will always truncate down
  // provide an epsilon that is small relative to dtFast to help with integer conversion
  const double eps = fastDt*epsFactor;
  return static_cast<int>((driverDt+eps)/fastDt);
}

//Constructor
fast::fastInputs::fastInputs():
    nTurbinesGlob(0),
    dryRun(false),
    debug(false),
    tStart(-1.0),
    restartFreq(-1),
    tMax(0.0),
    dtDriver(0.0)
{
    //Nothing to do here
}

//Constructor
fast::OpenFAST::OpenFAST()
{

    ncRstVarNames_ = {"time", "rst_filename", "twr_ref_pos", "bld_ref_pos", "nac_ref_pos", "hub_ref_pos", "twr_def", "twr_vel", "twr_ld", "bld_def", "bld_vel", "bld_ld", "hub_def", "hub_vel", "nac_def", "nac_vel", "bld_root_def", "bld_pitch", "x_vel", "vel_vel", "x_force", "xdot_force", "orient_force", "vel_force", "force"};
    ncRstDimNames_ = {"n_tsteps", "n_states", "n_twr_data", "n_bld_data", "n_pt_data", "n_bld_root_data", "n_bld_pitch_data", "n_vel_pts_data", "n_force_pts_data", "n_force_pts_orient_data"};

    ncOutVarNames_ = {"time", "twr_ref_pos", "twr_ref_orient", "bld_chord", "bld_rloc", "bld_ref_pos", "bld_ref_orient", "hub_ref_pos", "hub_ref_orient", "nac_ref_pos", "nac_ref_orient", "twr_disp", "twr_orient", "twr_vel", "twr_rotvel", "twr_ld", "twr_moment", "bld_disp", "bld_orient", "bld_vel", "bld_rotvel", "bld_ld", "bld_ld_loc", "bld_moment", "hub_disp", "hub_orient", "hub_vel", "hub_rotvel", "nac_disp", "nac_orient", "nac_vel", "nac_rotvel", "bld_root_ref_pos", "bld_root_ref_orient", "bld_root_disp", "bld_root_orient"};
    ncOutDimNames_ = {"n_tsteps", "n_dim", "n_twr_nds", "n_blds", "n_bld_nds"};
}

inline bool fast::OpenFAST::checkFileExists(const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

void fast::OpenFAST::findRestartFile(int iTurbLoc) {

    int ncid;
    size_t n_tsteps;
    size_t count1 = 1;
    double latest_time;

    //Find the file and open it in read only mode
    std::stringstream rstfile_ss;
    rstfile_ss << "turb_" ;
    rstfile_ss << std::setfill('0') << std::setw(2) << turbineData[iTurbLoc].TurbID;
    rstfile_ss << "_rst.nc";
    std::string rst_filename = rstfile_ss.str();
    int ierr = nc_open(rst_filename.c_str(), NC_NOWRITE, &ncid);
    check_nc_error(ierr, "nc_open");


    for (auto const& dim_name: ncRstDimNames_) {
        int tmpDimID;
        ierr = nc_inq_dimid(ncid, dim_name.data(), &tmpDimID);
        if (ierr == NC_NOERR)
            ncRstDimIDs_[dim_name] = tmpDimID;
    }

    for (auto const& var_name: ncRstVarNames_) {
        int tmpVarID;
        ierr = nc_inq_varid(ncid, var_name.data(), &tmpVarID);
        if (ierr == NC_NOERR)
            ncRstVarIDs_[var_name] = tmpVarID;
    }

    ierr = nc_inq_dimlen(ncid, ncRstDimIDs_["n_tsteps"], &n_tsteps);
    check_nc_error(ierr, "nc_inq_dimlen");
    n_tsteps -= 1; //To account for 0 based indexing
    ierr = nc_get_vara_double(ncid, ncRstVarIDs_["time"], &n_tsteps, &count1, &latest_time);
    check_nc_error(ierr, "nc_get_vara_double - getting latest time");
    tStart = latest_time;

    char *tmpOutFileRoot;
    size_t len;
    ierr = nc_inq_attlen(ncid, NC_GLOBAL, "out_file_root", &len);
    check_nc_error(ierr, "nc_inq_attlen - getting out_file_root length");

    tmpOutFileRoot = (char*) malloc(len + 1);
    ierr = nc_get_att_text(ncid, NC_GLOBAL, "out_file_root", tmpOutFileRoot);
    check_nc_error(ierr, "nc_get_att_text - getting out_file_root");
    tmpOutFileRoot[len] = '\0';
    turbineData[iTurbLoc].outFileRoot.assign(tmpOutFileRoot);

    ierr = nc_get_att_double(ncid, NC_GLOBAL, "dt_fast", &dtFAST);
    check_nc_error(ierr, "nc_get_att_double");

    ierr = nc_get_att_double(ncid, NC_GLOBAL, "dt_driver", &dtDriver);
    check_nc_error(ierr, "nc_get_att_double");

    ierr = nc_get_att_int(ncid, NC_GLOBAL, "output_freq", &outputFreq_);
    check_nc_error(ierr, "nc_get_att_int");

    ierr = nc_get_att_int(ncid, NC_GLOBAL, "restart_freq", &restartFreq_);
    check_nc_error(ierr, "nc_get_att_int");

    int tstep = std::round(latest_time/dtFAST);

    std::stringstream rstfilename;
    rstfilename << turbineData[iTurbLoc].outFileRoot <<  "." << tstep ;
    turbineData[iTurbLoc].FASTRestartFileName = rstfilename.str();

    std::cout << "Restarting from time " << latest_time << " at time step " << tstep << " from file name " << turbineData[iTurbLoc].FASTRestartFileName << std::endl ;

    nc_close(ncid);
    free(tmpOutFileRoot);
}

void fast::OpenFAST::prepareRestartFile(int iTurbLoc) {

    int ncid;
    //This will destroy any existing file
    std::stringstream rstfile_ss;
    rstfile_ss << "turb_" ;
    rstfile_ss << std::setfill('0') << std::setw(2) << turbineData[iTurbLoc].TurbID;
    rstfile_ss << "_rst.nc";
    std::string rst_filename = rstfile_ss.str();
    int ierr = nc_create(rst_filename.c_str(), NC_CLOBBER, &ncid);
    check_nc_error(ierr, "nc_create");

    nc_put_att_text(ncid, NC_GLOBAL, "out_file_root", turbineData[iTurbLoc].outFileRoot.size()+1, turbineData[iTurbLoc].outFileRoot.c_str());
    nc_put_att_double(ncid, NC_GLOBAL, "dt_fast", NC_DOUBLE, 1, &dtFAST);
    nc_put_att_double(ncid, NC_GLOBAL, "dt_driver", NC_DOUBLE, 1, &dtDriver);
    nc_put_att_int(ncid,NC_GLOBAL,"output_freq", NC_INT, 1, &outputFreq_);
    nc_put_att_int(ncid,NC_GLOBAL,"restart_freq", NC_INT, 1, &restartFreq_);

    //Define dimensions
    int tmpDimID;
    ierr = nc_def_dim(ncid, "n_tsteps", NC_UNLIMITED, &tmpDimID);
    ncRstDimIDs_["n_tsteps"] = tmpDimID;
    ierr = nc_def_dim(ncid, "n_states", 4, &tmpDimID);
    ncRstDimIDs_["n_states"] = tmpDimID;

    //Define variables
    int tmpVarID;
    ierr = nc_def_var(ncid, "time", NC_DOUBLE, 1, &ncRstDimIDs_["n_tsteps"], &tmpVarID);
    ncRstVarIDs_["time"] = tmpVarID;

    if (turbineData[iTurbLoc].sType == EXTLOADS) {

        ierr = nc_def_dim(ncid, "n_twr_data", turbineData[iTurbLoc].nBRfsiPtsTwr*6, &tmpDimID);
        ncRstDimIDs_["n_twr_data"] = tmpDimID;
        ierr = nc_def_dim(ncid,"n_bld_data", turbineData[iTurbLoc].nTotBRfsiPtsBlade*6, &tmpDimID);
        ncRstDimIDs_["n_bld_data"] = tmpDimID;
        ierr = nc_def_dim(ncid,"n_bld_root_data", turbineData[iTurbLoc].numBlades*6, &tmpDimID);
        ncRstDimIDs_["n_bld_root_data"] = tmpDimID;
        ierr = nc_def_dim(ncid,"n_bld_pitch_data", turbineData[iTurbLoc].numBlades, &tmpDimID);
        ncRstDimIDs_["n_bld_pitch_data"] = tmpDimID;
        ierr = nc_def_dim(ncid,"n_pt_data", 6, &tmpDimID);
        ncRstDimIDs_["n_pt_data"] = tmpDimID;

        const std::vector<int> twrDefLoadsDims{ncRstDimIDs_["n_tsteps"], ncRstDimIDs_["n_states"], ncRstDimIDs_["n_twr_data"]};
        const std::vector<int> bldDefLoadsDims{ncRstDimIDs_["n_tsteps"], ncRstDimIDs_["n_states"], ncRstDimIDs_["n_bld_data"]};
        const std::vector<int> bldRootDefsDims{ncRstDimIDs_["n_tsteps"], ncRstDimIDs_["n_states"], ncRstDimIDs_["n_bld_root_data"]};
        const std::vector<int> bldPitchDims{ncRstDimIDs_["n_tsteps"], ncRstDimIDs_["n_states"], ncRstDimIDs_["n_bld_pitch_data"]};
        const std::vector<int> ptDefLoadsDims{ncRstDimIDs_["n_tsteps"], ncRstDimIDs_["n_states"], ncRstDimIDs_["n_pt_data"],};

        ierr = nc_def_var(ncid, "twr_def", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["twr_def"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_vel", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["twr_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_ld", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["twr_ld"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_def", NC_DOUBLE, 3, bldDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["bld_def"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_vel", NC_DOUBLE, 3, bldDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["bld_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_ld", NC_DOUBLE, 3, bldDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["bld_ld"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_def", NC_DOUBLE, 3, ptDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["hub_def"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_vel", NC_DOUBLE, 3, ptDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["hub_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_def", NC_DOUBLE, 3, ptDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["nac_def"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_vel", NC_DOUBLE, 3, ptDefLoadsDims.data(), &tmpVarID);
        ncRstVarIDs_["nac_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_root_def", NC_DOUBLE, 3, bldRootDefsDims.data(), &tmpVarID);
        ncRstVarIDs_["bld_root_def"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_pitch", NC_DOUBLE, 3, bldPitchDims.data(), &tmpVarID);
        ncRstVarIDs_["bld_pitch"] = tmpVarID;

    } else if (turbineData[iTurbLoc].sType == EXTINFLOW) {

        ierr = nc_def_dim(ncid, "n_vel_pts_data", turbineData[iTurbLoc].numVelPts*3, &tmpDimID);
        ncRstDimIDs_["n_vel_pts_data"] = tmpDimID;
        ierr = nc_def_dim(ncid, "n_force_pts_data", turbineData[iTurbLoc].numForcePts*3, &tmpDimID);
        ncRstDimIDs_["n_force_pts_data"] = tmpDimID;
        ierr = nc_def_dim(ncid, "n_force_pts_orient_data", turbineData[iTurbLoc].numForcePts*9, &tmpDimID);
        ncRstDimIDs_["n_force_pts_orient_data"] = tmpDimID;

        const std::vector<int> velPtsDataDims{ncRstDimIDs_["n_tsteps"], ncRstDimIDs_["n_states"], ncRstDimIDs_["n_vel_pts_data"]};
        const std::vector<int> forcePtsDataDims{ncRstDimIDs_["n_tsteps"], ncRstDimIDs_["n_states"], ncRstDimIDs_["n_force_pts_data"],};
        const std::vector<int> forcePtsOrientDataDims{ncRstDimIDs_["n_tsteps"], ncRstDimIDs_["n_states"], ncRstDimIDs_["n_force_pts_orient_data"],};

        ierr = nc_def_var(ncid, "x_vel", NC_DOUBLE, 3, velPtsDataDims.data(), &tmpVarID);
        ncRstVarIDs_["x_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "vel_vel", NC_DOUBLE, 3, velPtsDataDims.data(), &tmpVarID);
        ncRstVarIDs_["vel_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "xref_force", NC_DOUBLE, 1, &ncRstDimIDs_["n_force_pts_data"], &tmpVarID);
        ncRstVarIDs_["xref_force"] = tmpVarID;
        ierr = nc_def_var(ncid, "x_force", NC_DOUBLE, 3, forcePtsDataDims.data(), &tmpVarID);
        ncRstVarIDs_["x_force"] = tmpVarID;
        ierr = nc_def_var(ncid, "xdot_force", NC_DOUBLE, 3, forcePtsDataDims.data(), &tmpVarID);
        ncRstVarIDs_["xdot_force"] = tmpVarID;
        ierr = nc_def_var(ncid, "vel_force", NC_DOUBLE, 3, forcePtsDataDims.data(), &tmpVarID);
        ncRstVarIDs_["vel_force"] = tmpVarID;
        ierr = nc_def_var(ncid, "force", NC_DOUBLE, 3, forcePtsDataDims.data(), &tmpVarID);
        ncRstVarIDs_["force"] = tmpVarID;
        ierr = nc_def_var(ncid, "orient_force", NC_DOUBLE, 3, forcePtsOrientDataDims.data(), &tmpVarID);
        ncRstVarIDs_["orient_force"] = tmpVarID;

    }

    //! Indicate that we are done defining variables, ready to write data
    ierr = nc_enddef(ncid);
    check_nc_error(ierr, "nc_enddef");

    if ( (turbineData[iTurbLoc].sType == EXTINFLOW) && (turbineData[iTurbLoc].inflowType == 2) ) {

        int nfpts_data = 3*get_numForcePtsLoc(iTurbLoc);
        int ierr = nc_put_var_double(ncid, ncRstVarIDs_["xref_force"], velForceNodeData[iTurbLoc][fast::STATE_NP1].xref_force.data());
    }

    ierr = nc_close(ncid);
    check_nc_error(ierr, "nc_close");


}

void fast::OpenFAST::findOutputFile(int iTurbLoc) {

    int ncid;
    size_t n_tsteps;
    size_t count1 = 1;
    double latest_time;

    //Find the file and open it in read only mode
    std::stringstream outfile_ss;
    outfile_ss << "turb_" ;
    outfile_ss << std::setfill('0') << std::setw(2) << turbineData[iTurbLoc].TurbID;
    outfile_ss << "_output.nc";
    std::string out_filename = outfile_ss.str();
    int ierr = nc_open(out_filename.c_str(), NC_NOWRITE, &ncid);
    check_nc_error(ierr, "nc_open");


    for (auto const& dim_name: ncOutDimNames_) {
        int tmpDimID;
        ierr = nc_inq_dimid(ncid, dim_name.data(), &tmpDimID);
        if (ierr == NC_NOERR)
            ncOutDimIDs_[dim_name] = tmpDimID;
    }

    for (auto const& var_name: ncOutVarNames_) {
        int tmpVarID;
        ierr = nc_inq_varid(ncid, var_name.data(), &tmpVarID);
        if (ierr == NC_NOERR)
            ncOutVarIDs_[var_name] = tmpVarID;
    }

    ierr = nc_inq_dimlen(ncid, ncOutDimIDs_["n_tsteps"], &n_tsteps);
    check_nc_error(ierr, "nc_inq_dimlen");
    n_tsteps -= 1; //To account for 0 based indexing
    ierr = nc_get_vara_double(ncid, ncOutVarIDs_["time"], &n_tsteps, &count1, &latest_time);
    check_nc_error(ierr, "nc_get_vara_double - getting latest time");
    nc_close(ncid);

}

void fast::OpenFAST::prepareOutputFile(int iTurbLoc) {

    int ncid;
    //Create the file - this will destory any file
    std::stringstream defloads_fstream;
    defloads_fstream << "turb_" ;
    defloads_fstream << std::setfill('0') << std::setw(2) << turbineData[iTurbLoc].TurbID;
    defloads_fstream << "_output.nc";
    std::string defloads_filename = defloads_fstream.str();
    int ierr = nc_create(defloads_filename.c_str(), NC_CLOBBER, &ncid);
    check_nc_error(ierr, "nc_create");

    //Define dimensions
    int tmpDimID;
    ierr = nc_def_dim(ncid, "n_dim", 3, &tmpDimID);
    ncOutDimIDs_["n_dim"] = tmpDimID;
    ierr = nc_def_dim(ncid, "n_tsteps", NC_UNLIMITED, &tmpDimID);
    ncOutDimIDs_["n_tsteps"] = tmpDimID;

    //Now define variables
    int tmpVarID;
    ierr = nc_def_var(ncid, "time", NC_DOUBLE, 1, &ncOutDimIDs_["n_tsteps"], &tmpVarID);
    ncOutVarIDs_["time"] = tmpVarID;

    if (turbineData[iTurbLoc].sType == EXTLOADS) {

        int nBlades = turbineData[iTurbLoc].numBlades;
        int nTwrPts = turbineData[iTurbLoc].nBRfsiPtsTwr;
        int nTotBldPts = turbineData[iTurbLoc].nTotBRfsiPtsBlade;
        int nBldPts = nTotBldPts/nBlades;

        ierr = nc_def_dim(ncid, "n_twr_nds", nTwrPts, &tmpDimID);
        ncOutDimIDs_["n_twr_nds"] = tmpDimID;
        ierr = nc_def_dim(ncid,"n_blds", nBlades, &tmpDimID);
        ncOutDimIDs_["n_blds"] = tmpDimID;
        ierr = nc_def_dim(ncid, "n_bld_nds", nBldPts, &tmpDimID);
        ncOutDimIDs_["n_bld_nds"] = tmpDimID;

        const std::vector<int> twrRefDims{ncOutDimIDs_["n_dim"], ncOutDimIDs_["n_twr_nds"]};
        const std::vector<int> twrDefLoadsDims{ncOutDimIDs_["n_tsteps"], ncOutDimIDs_["n_dim"], ncOutDimIDs_["n_twr_nds"]};
        const std::vector<int> bldParamDims{ncOutDimIDs_["n_blds"], ncOutDimIDs_["n_bld_nds"]};
        const std::vector<int> bldRefDims{ncOutDimIDs_["n_blds"], ncOutDimIDs_["n_dim"], ncOutDimIDs_["n_bld_nds"]};
        const std::vector<int> bldRootRefDims{ncOutDimIDs_["n_blds"], ncOutDimIDs_["n_dim"]};
        const std::vector<int> bldDefLoadsDims{ncOutDimIDs_["n_tsteps"], ncOutDimIDs_["n_blds"], ncOutDimIDs_["n_dim"], ncOutDimIDs_["n_bld_nds"]};
        const std::vector<int> bldRootDefDims{ncOutDimIDs_["n_tsteps"], ncOutDimIDs_["n_blds"], ncOutDimIDs_["n_dim"]};
        const std::vector<int> ptRefDims{ncOutDimIDs_["n_dim"]};
        const std::vector<int> ptDefLoadsDims{ncOutDimIDs_["n_tsteps"], ncOutDimIDs_["n_dim"]};

        ierr = nc_def_var(ncid, "twr_ref_pos", NC_DOUBLE, 2, twrRefDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_ref_orient", NC_DOUBLE, 2, twrRefDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_ref_orient"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_chord", NC_DOUBLE, 2, bldParamDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_chord"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_rloc", NC_DOUBLE, 2, bldParamDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_rloc"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_ref_pos", NC_DOUBLE, 3, bldRefDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_ref_orient", NC_DOUBLE, 3, bldRefDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_ref_orient"] = tmpVarID;

        ierr = nc_def_var(ncid, "bld_root_ref_pos", NC_DOUBLE, 2, bldRootRefDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_root_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_root_ref_orient", NC_DOUBLE, 2, bldRootRefDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_root_ref_orient"] = tmpVarID;

        ierr = nc_def_var(ncid, "hub_ref_pos", NC_DOUBLE, 1, ptRefDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_ref_orient", NC_DOUBLE, 1, ptRefDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_ref_orient"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_ref_pos", NC_DOUBLE, 1, ptRefDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_ref_orient", NC_DOUBLE, 1, ptRefDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_ref_orient"] = tmpVarID;

        ierr = nc_def_var(ncid, "twr_disp", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_orient", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_orient"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_vel", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_rotvel", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_rotvel"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_ld", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_ld"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_moment", NC_DOUBLE, 3, twrDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_moment"] = tmpVarID;

        ierr = nc_def_var(ncid, "bld_disp", NC_DOUBLE, 4, bldDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_orient", NC_DOUBLE, 4, bldDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_orient"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_vel", NC_DOUBLE, 4, bldDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_rotvel", NC_DOUBLE, 4, bldDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_rotvel"] = tmpVarID;

        ierr = nc_def_var(ncid, "bld_root_disp", NC_DOUBLE, 3, bldRootDefDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_root_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_root_orient", NC_DOUBLE, 3, bldRootDefDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_root_orient"] = tmpVarID;

        ierr = nc_def_var(ncid, "bld_ld", NC_DOUBLE, 4, bldDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_ld"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_ld_loc", NC_DOUBLE, 4, bldDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_ld_loc"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_moment", NC_DOUBLE, 4, bldDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_moment"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_disp", NC_DOUBLE, 2, ptDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_orient", NC_DOUBLE, 2, ptDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_orient"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_vel", NC_DOUBLE, 2, ptDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_rotvel", NC_DOUBLE, 2, ptDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_rotvel"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_disp", NC_DOUBLE, 2, ptDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_orient", NC_DOUBLE, 2, ptDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_orient"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_vel", NC_DOUBLE, 2, ptDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_rotvel", NC_DOUBLE, 2, ptDefLoadsDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_rotvel"] = tmpVarID;

    } else if (turbineData[iTurbLoc].sType == EXTINFLOW) {

        int nBlades = get_numBladesLoc(iTurbLoc);
        int nBldPts = get_numForcePtsBladeLoc(iTurbLoc);
        int nTwrPts = get_numForcePtsTwrLoc(iTurbLoc);

        ierr = nc_def_dim(ncid, "n_twr_nds", nTwrPts, &tmpDimID);
        ncOutDimIDs_["n_twr_nds"] = tmpDimID;
        ierr = nc_def_dim(ncid,"n_blds", nBlades, &tmpDimID);
        ncOutDimIDs_["n_blds"] = tmpDimID;
        ierr = nc_def_dim(ncid, "n_bld_nds", nBldPts, &tmpDimID);
        ncOutDimIDs_["n_bld_nds"] = tmpDimID;

        const std::vector<int> twrRefDataDims{ncOutDimIDs_["n_dim"], ncOutDimIDs_["n_twr_nds"]};
        const std::vector<int> twrDataDims{ncOutDimIDs_["n_tsteps"], ncOutDimIDs_["n_dim"], ncOutDimIDs_["n_twr_nds"]};
        const std::vector<int> bldParamDims{ncOutDimIDs_["n_blds"], ncOutDimIDs_["n_bld_nds"]};
        const std::vector<int> bldRefDataDims{ncOutDimIDs_["n_blds"], ncOutDimIDs_["n_dim"], ncOutDimIDs_["n_bld_nds"]};
        const std::vector<int> bldDataDims{ncOutDimIDs_["n_tsteps"], ncOutDimIDs_["n_blds"], ncOutDimIDs_["n_dim"], ncOutDimIDs_["n_bld_nds"]};
        const std::vector<int> ptRefDataDims{ncOutDimIDs_["n_dim"]};
        const std::vector<int> ptDataDims{ncOutDimIDs_["n_tsteps"], ncOutDimIDs_["n_dim"]};

        ierr = nc_def_var(ncid, "bld_chord", NC_DOUBLE, 2, bldParamDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_chord"] = tmpVarID;

        ierr = nc_def_var(ncid, "twr_ref_pos", NC_DOUBLE, 2, twrRefDataDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_disp", NC_DOUBLE, 3, twrDataDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_vel", NC_DOUBLE, 3, twrDataDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "twr_ld", NC_DOUBLE, 3, twrDataDims.data(), &tmpVarID);
        ncOutVarIDs_["twr_ld"] = tmpVarID;

        ierr = nc_def_var(ncid, "bld_ref_pos", NC_DOUBLE, 3, bldRefDataDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_disp", NC_DOUBLE, 4, bldDataDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_vel", NC_DOUBLE, 4, bldDataDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_ld", NC_DOUBLE, 4, bldDataDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_ld"] = tmpVarID;
        ierr = nc_def_var(ncid, "bld_ld_loc", NC_DOUBLE, 4, bldDataDims.data(), &tmpVarID);
        ncOutVarIDs_["bld_ld_loc"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_ref_pos", NC_DOUBLE, 1, ptRefDataDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_disp", NC_DOUBLE, 2, ptDataDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_vel", NC_DOUBLE, 2, ptDataDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "hub_rotvel", NC_DOUBLE, 2, ptDataDims.data(), &tmpVarID);
        ncOutVarIDs_["hub_rotvel"] = tmpVarID;

        ierr = nc_def_var(ncid, "nac_ref_pos", NC_DOUBLE, 1, ptRefDataDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_ref_pos"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_disp", NC_DOUBLE, 2, ptDataDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_disp"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_vel", NC_DOUBLE, 2, ptDataDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_vel"] = tmpVarID;
        ierr = nc_def_var(ncid, "nac_rotvel", NC_DOUBLE, 2, ptDataDims.data(), &tmpVarID);
        ncOutVarIDs_["nac_rotvel"] = tmpVarID;

    }

    //! Indicate that we are done defining variables, ready to write data
    ierr = nc_enddef(ncid);
    check_nc_error(ierr, "nc_enddef");

    if (turbineData[iTurbLoc].sType == EXTLOADS) {

        int nBlades = turbineData[iTurbLoc].numBlades;
        int nTwrPts = turbineData[iTurbLoc].nBRfsiPtsTwr;
        int nTotBldPts = turbineData[iTurbLoc].nTotBRfsiPtsBlade;
        int nBldPts = nTotBldPts/nBlades;

        std::vector<double> tmpArray;

        tmpArray.resize(nTwrPts);
        {
            std::vector<size_t> count_dim{1,static_cast<size_t>(nTwrPts)};
            for (size_t idim=0;idim < 3; idim++) {
                for (size_t i=0; i < nTwrPts; i++)
                    tmpArray[i] = brFSIData[iTurbLoc][3].twr_ref_pos[i*6+idim];
                std::vector<size_t> start_dim{idim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_ref_pos"], start_dim.data(),
                                          count_dim.data(), tmpArray.data());
            }
            for (size_t idim=0;idim < 3; idim++) {
                for (size_t i=0; i < nTwrPts; i++)
                    tmpArray[i] = brFSIData[iTurbLoc][3].twr_ref_pos[i*6+3+idim];
                std::vector<size_t> start_dim{idim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_ref_orient"], start_dim.data(),
                                          count_dim.data(), tmpArray.data());
            }
        }

        tmpArray.resize(nBldPts);
        {
            std::vector<size_t> count_dim{1,1,static_cast<size_t>(nBldPts)};
            for (size_t iDim=0;iDim < 3; iDim++) {
                int iStart = 0 ;
                for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                    for (auto i=0; i < nBldPts; i++) {
                        tmpArray[i] = brFSIData[iTurbLoc][3].bld_ref_pos[(iStart*6)+iDim];
                        iStart++;
                    }
                    std::vector<size_t> start_dim{iBlade,iDim,0};
                    ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_ref_pos"], start_dim.data(),
                                              count_dim.data(), tmpArray.data());
                }
            }
            for (size_t iDim=0;iDim < 3; iDim++) {
                int iStart = 0 ;
                for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                    for (auto i=0; i < nBldPts; i++) {
                        tmpArray[i] = brFSIData[iTurbLoc][3].bld_ref_pos[(iStart*6)+iDim+3];
                        iStart++;
                    }
                    std::vector<size_t> start_dim{iBlade,iDim,0};
                    ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_ref_orient"], start_dim.data(),
                                              count_dim.data(), tmpArray.data());
                }
            }

            std::vector<size_t> param_count_dim{1,static_cast<size_t>(nBldPts)};
            int iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (size_t i=0; i < nBldPts; i++) {
                    tmpArray[i] = brFSIData[iTurbLoc][3].bld_chord[iStart];
                    iStart++;
                }
                std::vector<size_t> start_dim{iBlade,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_chord"], start_dim.data(),
                                          param_count_dim.data(), tmpArray.data());
            }
            iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (size_t i=0; i < nBldPts; i++) {
                    tmpArray[i] = brFSIData[iTurbLoc][3].bld_rloc[iStart];
                    iStart++;
                }
                std::vector<size_t> start_dim{iBlade,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_rloc"], start_dim.data(),
                                          param_count_dim.data(), tmpArray.data());
            }
        }

        for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
            std::vector<size_t> start_dim{iBlade,0};
            std::vector<size_t> count_dim{1,3};
            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_root_ref_pos"],
                                     start_dim.data(),
                                     count_dim.data(),
                                     &brFSIData[iTurbLoc][3].bld_root_ref_pos[iBlade*6+0]);

            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_root_ref_orient"],
                                     start_dim.data(),
                                     count_dim.data(),
                                     &brFSIData[iTurbLoc][3].bld_root_ref_pos[iBlade*6+3]);
        }

        ierr = nc_put_var_double(ncid, ncOutVarIDs_["nac_ref_pos"],
                                 &brFSIData[iTurbLoc][3].nac_ref_pos[0]);
        ierr = nc_put_var_double(ncid, ncOutVarIDs_["nac_ref_orient"],
                                 &brFSIData[iTurbLoc][3].nac_ref_pos[3]);

        ierr = nc_put_var_double(ncid, ncOutVarIDs_["hub_ref_pos"],
                                 &brFSIData[iTurbLoc][3].hub_ref_pos[0]);
        ierr = nc_put_var_double(ncid, ncOutVarIDs_["hub_ref_orient"],
                                 &brFSIData[iTurbLoc][3].hub_ref_pos[3]);

    } else if (turbineData[iTurbLoc].sType == EXTINFLOW) {

        int nBlades = get_numBladesLoc(iTurbLoc);
        int nBldPts = get_numForcePtsBladeLoc(iTurbLoc);
        int nTwrPts = get_numForcePtsTwrLoc(iTurbLoc);

        std::vector<double> tmpArray;

        {

            tmpArray.resize(nBldPts);
            std::vector<size_t> count_dim{1,1,static_cast<size_t>(nBldPts)};
            for (size_t iDim=0;iDim < 3; iDim++) {
                for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                    int node_bld_start = (1 + iBlade * nBldPts);
                    for (auto i=0; i < nBldPts; i++)
                        tmpArray[i] = velForceNodeData[iTurbLoc][3].x_force[(node_bld_start+i)*3+iDim] ;
                    std::vector<size_t> start_dim{iBlade,iDim,0};
                    ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_ref_pos"], start_dim.data(), count_dim.data(), tmpArray.data());
                }
            }

            std::vector<size_t> param_count_dim{1,static_cast<size_t>(nBldPts)};
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                int iStart = 1 + iBlade*nBldPts;
                for (size_t i=0; i < nBldPts; i++)
                    tmpArray[i] = extinfw_i_f_FAST[iTurbLoc].forceNodesChord[iStart+i];
                std::vector<size_t> start_dim{iBlade,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_chord"], start_dim.data(),
                                          param_count_dim.data(), tmpArray.data());
            }
        }
    }

    ierr = nc_close(ncid);
    check_nc_error(ierr, "nc_close");

}


void fast::OpenFAST::init() {

    allocateMemory_preInit();

    if (!dryRun) {
        switch (simStart) {

        case fast::trueRestart:

            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

                findRestartFile(iTurb);
                findOutputFile(iTurb);
                std::string tmpRstFileRoot{turbineData[iTurb].FASTRestartFileName};
                tmpRstFileRoot.resize(INTERFACE_STRING_LENGTH, ' ');
                if (turbineData[iTurb].sType == EXTINFLOW) {
                    /* note that this will set nt_global inside the FAST library */
                    FAST_ExtInfw_Restart(
                        &iTurb,
                        tmpRstFileRoot.c_str(),
                        &AbortErrLev,
                        &turbineData[iTurb].dt,
                        &turbineData[iTurb].numBlades,
                        &turbineData[iTurb].numVelPtsBlade,
                        &turbineData[iTurb].numVelPtsTwr,
                        &ntStart,
                        &extinfw_i_f_FAST[iTurb],
                        &extinfw_o_t_FAST[iTurb],
                        &ErrStat,
                        ErrMsg);
                    checkError(ErrStat, ErrMsg);
                } else if(turbineData[iTurb].sType == EXTLOADS) {
                    FAST_ExtLoads_Restart(
                        &iTurb,
                        tmpRstFileRoot.c_str(),
                        &AbortErrLev,
                        &turbineData[iTurb].dt,
                        &turbineData[iTurb].numBlades,
                        &ntStart,
                        &extld_i_f_FAST[iTurb],
                        &extld_p_f_FAST[iTurb],
                        &extld_o_t_FAST[iTurb],
                        &ErrStat,
                        ErrMsg);
                    checkError(ErrStat, ErrMsg);
                    turbineData[iTurb].inflowType = 0;
                }

                nt_global = ntStart;
                allocateMemory_postInit(iTurb);

                get_ref_positions_from_openfast(iTurb);

                readRestartFile(iTurb, nt_global);

            }
            checkAndSetSubsteps();

            break ;

        case fast::init:

            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

                char tmpOutFileRoot[INTERFACE_STRING_LENGTH];
                char inputFileName[INTERFACE_STRING_LENGTH];
                if (turbineData[iTurb].sType == EXTINFLOW) {

                    std::copy(
                        turbineData[iTurb].FASTInputFileName.data(),
                        turbineData[iTurb].FASTInputFileName.data() + (turbineData[iTurb].FASTInputFileName.size() + 1),
                        inputFileName
                        );
                    FAST_ExtInfw_Init(
                        &iTurb,
                        &tMax,
                        inputFileName,
                        &turbineData[iTurb].TurbID,
                        tmpOutFileRoot,
                        &turbineData[iTurb].numForcePtsBlade,
                        &turbineData[iTurb].numForcePtsTwr,
                        turbineData[iTurb].TurbineBasePos.data(),
                        &AbortErrLev,
                        &dtDriver,
                        &turbineData[iTurb].dt,
                        &turbineData[iTurb].inflowType,
                        &turbineData[iTurb].numBlades,
                        &turbineData[iTurb].numVelPtsBlade,
                        &turbineData[iTurb].numVelPtsTwr,
                        &turbineData[iTurb].nodeClusterType,
                        &extinfw_i_f_FAST[iTurb],
                        &extinfw_o_t_FAST[iTurb],
                        &ErrStat,
                        ErrMsg);
                    checkError(ErrStat, ErrMsg);

                    std::cerr << "turbineData[iTurb].inflowType = " << turbineData[iTurb].inflowType << std::endl;

                    turbineData[iTurb].numVelPtsTwr = extinfw_o_t_FAST[iTurb].u_Len - turbineData[iTurb].numBlades*turbineData[iTurb].numVelPtsBlade - 1;
                    if(turbineData[iTurb].numVelPtsTwr == 0) {
                        turbineData[iTurb].numForcePtsTwr = 0;
                        std::cout << "Aerodyn doesn't want to calculate forces on the tower. All actuator points on the tower are turned off for turbine " << turbineMapProcToGlob[iTurb] << "." << std::endl ;
                    }

                } else if(turbineData[iTurb].sType == EXTLOADS) {

                    char inputFileName[INTERFACE_STRING_LENGTH];
                    FAST_ExtLoads_Init(
                        &iTurb,
                        &tMax,
                        turbineData[iTurb].FASTInputFileName.data(),
                        &turbineData[iTurb].TurbID,
                        tmpOutFileRoot,
                        turbineData[iTurb].TurbineBasePos.data(),
                        &AbortErrLev,
                        &dtDriver,
                        &turbineData[iTurb].dt,
                        &turbineData[iTurb].numBlades,
                        &turbineData[iTurb].azBlendMean,
                        &turbineData[iTurb].azBlendDelta,
                        &extld_i_f_FAST[iTurb],
                        &extld_p_f_FAST[iTurb],
                        &extld_o_t_FAST[iTurb],
                        &ErrStat,
                        ErrMsg);
                    checkError(ErrStat, ErrMsg);

                    turbineData[iTurb].inflowType = 0;

                }
                timeZero = true;

                turbineData[iTurb].outFileRoot.assign(tmpOutFileRoot, strlen(tmpOutFileRoot));

                allocateMemory_postInit(iTurb);

                get_data_from_openfast(fast::STATE_NM2);
                get_data_from_openfast(fast::STATE_NM1);
                get_data_from_openfast(fast::STATE_N);
                get_data_from_openfast(fast::STATE_NP1);

                get_ref_positions_from_openfast(iTurb);

            }
            timeZero = true;
            checkAndSetSubsteps();

            break ;

        case fast::restartDriverInitFAST:

            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

                findOutputFile(iTurb);
                findRestartFile(iTurb);
                char tmpOutFileRoot[INTERFACE_STRING_LENGTH];
                char inputFileName[INTERFACE_STRING_LENGTH];
                if (turbineData[iTurb].sType == EXTINFLOW) {

                    std::copy(
                        turbineData[iTurb].FASTInputFileName.data(),
                        turbineData[iTurb].FASTInputFileName.data() + (turbineData[iTurb].FASTInputFileName.size() + 1),
                        inputFileName
                    );

                    FAST_ExtInfw_Init(
                        &iTurb,
                        &tMax,
                        inputFileName,
                        &turbineData[iTurb].TurbID,
                        tmpOutFileRoot,
                        &turbineData[iTurb].numForcePtsBlade,
                        &turbineData[iTurb].numForcePtsTwr,
                        turbineData[iTurb].TurbineBasePos.data(),
                        &AbortErrLev,
                        &dtDriver,
                        &turbineData[iTurb].dt,
                        &turbineData[iTurb].inflowType,
                        &turbineData[iTurb].numBlades,
                        &turbineData[iTurb].numVelPtsBlade,
                        &turbineData[iTurb].numVelPtsTwr,
                        &turbineData[iTurb].nodeClusterType,
                        &extinfw_i_f_FAST[iTurb],
                        &extinfw_o_t_FAST[iTurb],
                        &ErrStat,
                        ErrMsg);
                    checkError(ErrStat, ErrMsg);

                    timeZero = true;

                    turbineData[iTurb].numVelPtsTwr = extinfw_o_t_FAST[iTurb].u_Len - turbineData[iTurb].numBlades*turbineData[iTurb].numVelPtsBlade - 1;
                    if(turbineData[iTurb].numVelPtsTwr == 0) {
                        turbineData[iTurb].numForcePtsTwr = 0;
                        std::cout << "Aerodyn doesn't want to calculate forces on the tower. All actuator points on the tower are turned off for turbine " << turbineMapProcToGlob[iTurb] << "." << std::endl ;
                    }

                    allocateMemory_postInit(iTurb);

                    get_data_from_openfast(fast::STATE_NM2);
                    get_data_from_openfast(fast::STATE_NM1);
                    get_data_from_openfast(fast::STATE_N);
                    get_data_from_openfast(fast::STATE_NP1);

                    get_ref_positions_from_openfast(iTurb);

                    checkAndSetSubsteps();

                    ntStart = int(tStart/dtFAST);
                    int ntStartDriver;
                    if( (dtFAST > 0) && (nSubsteps_ > 0))
                        ntStartDriver = int(tStart/dtFAST/nSubsteps_);
                    else
                        ntStartDriver = 0; //Typically for processors that don't contain any turbines

                    std::vector<int> velfile_ncid;
                    velfile_ncid.resize(nTurbinesProc);

                    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                        velfile_ncid[iTurb] = openVelocityDataFile(iTurb);
                        readVelocityData(iTurb, 0, 0, velfile_ncid[iTurb]);
                    }

                    int nVelPts = get_numVelPtsLoc(iTurb);
                    std::cout << std::endl ;
                    std::cout << "nt_global = " << 0 << " nlin_iter = " << 0 << std::endl ;
                    for (size_t k = 0; k < nVelPts; k++)
                         std::cout << k << ", " << velForceNodeData[iTurb][3].vel_vel[k*3 + 0] << " " << velForceNodeData[iTurb][3].vel_vel[k*3 + 1] << " " << velForceNodeData[iTurb][3].vel_vel[k*3 + 2] << " " << std::endl ;

                    init_velForceNodeData();

                    solution0(false) ;

                    for (int iPrestart=0 ; iPrestart < ntStartDriver; iPrestart++) {
                        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                            int nlinIters = read_nlin_iters(iTurb, iPrestart+1, velfile_ncid[iTurb]);
                            for (int iNlin=0; iNlin < nlinIters; iNlin++) {
                                readVelocityData(iTurb, iPrestart+1, iNlin, velfile_ncid[iTurb]);
                                update_states_driver_time_step(false);
                            }
                            advance_to_next_driver_time_step(false);
                        }
                    }

                    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++)
                        nc_close(velfile_ncid[iTurb]);

                    readRestartFile(iTurb, nt_global);

                } else {

                    throw std::runtime_error("RESTARTDRIVERINITFAST option not supported for blade-resolved FSI yet");
                }
            }

            break;

        case fast::simStartType_END:

            break;
        }

    }
}

void fast::OpenFAST::solution0(bool writeFiles) {

    if (!dryRun) {

        if (writeFiles) {
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                prepareRestartFile(iTurb);
                prepareOutputFile(iTurb);
                prepareVelocityDataFile(iTurb);
            }
        }

        // Unfortunately setVelocity only sets the velocity at 'n+1'. Need to copy 'n+1' to 'n'
        init_velForceNodeData() ;
        send_data_to_openfast(fast::STATE_NP1);

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

            FAST_CFD_Solution0(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);

            FAST_CFD_InitIOarrays_SubStep(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
        }

        get_data_from_openfast(fast::STATE_N);
        get_data_from_openfast(fast::STATE_NM1);
        get_data_from_openfast(fast::STATE_NM2);

        if (writeFiles) {
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                if (turbineData[iTurb].inflowType == 2)
                    writeVelocityData(iTurb, -nSubsteps_, 0);
            }
        }

        timeZero = false;

    }

}

void fast::OpenFAST::set_state_from_state(fast::timeStep fromState, fast::timeStep toState) {

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        if (turbineData[iTurb].sType == EXTINFLOW) {
            int nvelpts = get_numVelPtsLoc(iTurb);
            int nfpts = get_numForcePtsLoc(iTurb);
            for (int i=0; i<nvelpts; i++) {
                for (int j=0 ; j < 3; j++) {
                    velForceNodeData[iTurb][toState].x_vel[i*3+j] = velForceNodeData[iTurb][fromState].x_vel[i*3+j];
                    velForceNodeData[iTurb][toState].vel_vel[i*3+j] = velForceNodeData[iTurb][fromState].vel_vel[i*3+j];
                }
            }
            for (int i=0; i<nfpts; i++) {
                for (int j=0 ; j < 3; j++) {
                    velForceNodeData[iTurb][toState].x_force[i*3+j] = velForceNodeData[iTurb][fromState].x_force[i*3+j];
                    velForceNodeData[iTurb][toState].xdot_force[i*3+j] = velForceNodeData[iTurb][fromState].xdot_force[i*3+j];
                    velForceNodeData[iTurb][toState].vel_force[i*3+j] = velForceNodeData[iTurb][fromState].vel_force[i*3+j];
                    velForceNodeData[iTurb][toState].force[i*3+j] = velForceNodeData[iTurb][fromState].force[i*3+j];
                }
                for (int j=0;j<9;j++)
                    velForceNodeData[iTurb][toState].orient_force[i*9+j] = velForceNodeData[iTurb][fromState].orient_force[i*9+j];
            }
        } else if (turbineData[iTurb].sType == EXTLOADS) {

            int numBlades = get_numBladesLoc(iTurb);
            int nBRfsiPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
            int nTotBRfsiPtsBlade = turbineData[iTurb].nTotBRfsiPtsBlade;

            for (int i=0; i < nBRfsiPtsTwr; i++) {
                for (int j=0; j < 6; j++) {
                    brFSIData[iTurb][toState].twr_ref_pos[i*6+j] = brFSIData[iTurb][fromState].twr_ref_pos[i*6+j];
                    brFSIData[iTurb][toState].twr_def[i*6+j] = brFSIData[iTurb][fromState].twr_def[i*6+j];
                    brFSIData[iTurb][toState].twr_vel[i*6+j] = brFSIData[iTurb][fromState].twr_vel[i*6+j];
                    brFSIData[iTurb][toState].twr_ld[i*6+j] = brFSIData[iTurb][fromState].twr_ld[i*6+j];
                }
            }
            for (int i=0; i < nTotBRfsiPtsBlade; i++) {
                for (int j=0; j < 6; j++) {
                    brFSIData[iTurb][toState].bld_ref_pos[i*6+j] = brFSIData[iTurb][fromState].bld_ref_pos[i*6+j];
                    brFSIData[iTurb][toState].bld_def[i*6+j] = brFSIData[iTurb][fromState].bld_def[i*6+j];
                    brFSIData[iTurb][toState].bld_vel[i*6+j] = brFSIData[iTurb][fromState].bld_vel[i*6+j];
                    brFSIData[iTurb][toState].bld_ld[i*6+j] = brFSIData[iTurb][fromState].bld_ld[i*6+j];
                }
            }
            for (int j=0; j < 6; j++) {
                brFSIData[iTurb][toState].hub_ref_pos[j] = brFSIData[iTurb][fromState].hub_ref_pos[j];
                brFSIData[iTurb][toState].hub_def[j] = brFSIData[iTurb][fromState].hub_def[j];
                brFSIData[iTurb][toState].hub_vel[j] = brFSIData[iTurb][fromState].hub_vel[j];
                brFSIData[iTurb][toState].nac_ref_pos[j] = brFSIData[iTurb][fromState].nac_ref_pos[j];
                brFSIData[iTurb][toState].nac_def[j] = brFSIData[iTurb][fromState].nac_def[j];
                brFSIData[iTurb][toState].nac_vel[j] = brFSIData[iTurb][fromState].nac_vel[j];
            }
            for (int i=0; i < numBlades; i++) {
                for (int k=0; k < 6; k++) {
                    brFSIData[iTurb][toState].bld_root_def[i*6+k] = brFSIData[iTurb][fromState].bld_root_def[i*6+k];
                }
                brFSIData[iTurb][toState].bld_pitch[i] = brFSIData[iTurb][fromState].bld_pitch[i];
            }
        }
    }

}

void fast::OpenFAST::init_velForceNodeData() {

    set_state_from_state(fast::STATE_NP1, fast::STATE_N);
    set_state_from_state(fast::STATE_NP1, fast::STATE_NM1);
    set_state_from_state(fast::STATE_NP1, fast::STATE_NM2);

}

//! Dot product of two vectors
double dot(double * a, double * b) {

    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);

}

//! Cross product of two vectors
void cross(double * a, double * b, double * aCrossb) {

    aCrossb[0] = a[1]*b[2] - a[2]*b[1];
    aCrossb[1] = a[2]*b[0] - a[0]*b[2];
    aCrossb[2] = a[0]*b[1] - a[1]*b[0];

}

//! Compose Wiener-Milenkovic parameters 'p' and 'q' into 'pPlusq'. If a transpose of 'p' is required, set tranposeP to '-1', else leave blank or set to '+1'
void composeWM(double * p, double * q, double * pPlusq, double transposeP, double transposeQ) {

    double p0 = 2.0 - 0.125*dot(p,p);
    double q0 = 2.0 - 0.125*dot(q,q);
    std::vector<double> pCrossq(3,0.0);
    cross(p, q, pCrossq.data());

    double delta1 = (4.0-p0)*(4.0-q0);
    double delta2 = p0*q0 - transposeP*dot(p,q);
    double premultFac = 0.0;
    if (delta2 < 0)
        premultFac = -4.0/(delta1 - delta2);
    else
        premultFac = 4.0/(delta1 + delta2);

    for (size_t i=0; i < 3; i++)
        pPlusq[i] = premultFac * (transposeQ * p0 * q[i] + transposeP * q0 * p[i] + transposeP * transposeQ * pCrossq[i] );

}

//! Extrapolate Wiener-Milenkovic parameters from state 'nm2', 'nm1', 'n' to 'np1'
void extrapRotation(double *rnm2, double *rnm1, double *rn, double *rnp1) {

    std::array<double,3> rrnm1{ {0.0,0.0,0.0} };
    std::array<double,3> rrn{ {0.0,0.0,0.0} };
    std::array<double,3> rrnp1{ {0.0,0.0,0.0} };

    composeWM(rnm2, rnm1, rrnm1.data(), -1.0, 1.0); // Remove rigid body rotaiton of rnm2 from rnm1
    composeWM(rnm2, rn, rrn.data(), -1.0, 1.0); // Remove rigid body rotaiton of rnm2 from rnm1
    for(int i=0; i<3; i++) {
        rrnp1[i] = 3.0 * ( rrn[i] - rrnm1[i]) ;
    }
    composeWM(rnm2, rrnp1.data(), rnp1, 1.0, 1.0); //Add rigid body rotation of nm2 back

}


void fast::OpenFAST::predict_states() {

    if (firstPass_) {
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            int nvelpts = get_numVelPtsLoc(iTurb);
            int nfpts = get_numForcePtsLoc(iTurb);
            for (int i=0; i<nvelpts; i++) {
                for (int j=0 ; j < 3; j++) {
                    velForceNodeData[iTurb][fast::STATE_NP1].x_vel[i*3+j] = velForceNodeData[iTurb][fast::STATE_NM2].x_vel[i*3+j] + 3.0*velForceNodeData[iTurb][fast::STATE_N].x_vel[i*3+j] - 3.0*velForceNodeData[iTurb][fast::STATE_NM1].x_vel[i*3+j];
                    velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[i*3+j] = velForceNodeData[iTurb][fast::STATE_NM2].vel_vel[i*3+j] + 3.0*velForceNodeData[iTurb][fast::STATE_N].vel_vel[i*3+j] - 3.0*velForceNodeData[iTurb][fast::STATE_NM1].vel_vel[i*3+j];
                }
            }
            velForceNodeData[iTurb][fast::STATE_NP1].x_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].vel_vel_resid = 0.0;
            for (int i=0; i<nfpts; i++) {
                for (int j=0 ; j < 3; j++) {
                    velForceNodeData[iTurb][fast::STATE_NP1].x_force[i*3+j] = velForceNodeData[iTurb][fast::STATE_NM2].x_force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::STATE_N].x_force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::STATE_NM1].x_force[i*3+j];
                    velForceNodeData[iTurb][fast::STATE_NP1].xdot_force[i*3+j] = velForceNodeData[iTurb][fast::STATE_NM2].xdot_force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::STATE_N].xdot_force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::STATE_NM1].xdot_force[i*3+j];
                    velForceNodeData[iTurb][fast::STATE_NP1].vel_force[i*3+j] =  velForceNodeData[iTurb][fast::STATE_NM2].vel_force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::STATE_N].vel_force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::STATE_NM1].vel_force[i*3+j];
                    velForceNodeData[iTurb][fast::STATE_NP1].force[i*3+j] = velForceNodeData[iTurb][fast::STATE_NM2].force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::STATE_N].force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::STATE_NM1].force[i*3+j];
                }
                for (int j=0;j<9;j++)
                    velForceNodeData[iTurb][fast::STATE_NP1].orient_force[i*9+j] = velForceNodeData[iTurb][fast::STATE_NM2].orient_force[i*9+j] + 3.0*velForceNodeData[iTurb][fast::STATE_N].orient_force[i*9+j] - 3.0*velForceNodeData[iTurb][fast::STATE_NM1].orient_force[i*9+j];
            }
            velForceNodeData[iTurb][fast::STATE_NP1].x_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].xdot_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].orient_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].vel_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].force_resid = 0.0;

            if(turbineData[iTurb].sType == EXTLOADS) {
                int nTotBladeNodes = turbineData[iTurb].nTotBRfsiPtsBlade;
                for (int j=0; j < nTotBladeNodes; j++) {
                    extrapRotation(&brFSIData[iTurb][fast::STATE_NM2].bld_def[j*6+3], &brFSIData[iTurb][fast::STATE_NM1].bld_def[j*6+3], &brFSIData[iTurb][fast::STATE_N].bld_def[j*6+3], &brFSIData[iTurb][fast::STATE_NP1].bld_def[j*6+3]);
                    for (int k=0; k < 3; k++) {
                        brFSIData[iTurb][fast::STATE_NP1].bld_def[j*6+k] = brFSIData[iTurb][fast::STATE_NM2].bld_def[j*6+k] + 3.0*(brFSIData[iTurb][fast::STATE_N].bld_def[j*6+k] - brFSIData[iTurb][fast::STATE_NM1].bld_def[j*6+k]);
                        brFSIData[iTurb][fast::STATE_NP1].bld_vel[j*6+k] = brFSIData[iTurb][fast::STATE_NM2].bld_vel[j*6+k] + 3.0*(brFSIData[iTurb][fast::STATE_N].bld_vel[j*6+k] - brFSIData[iTurb][fast::STATE_NM1].bld_vel[j*6+k]);
                        brFSIData[iTurb][fast::STATE_NP1].bld_vel[j*6+k+3] = brFSIData[iTurb][fast::STATE_NM2].bld_vel[j*6+k+3] + 3.0*(brFSIData[iTurb][fast::STATE_N].bld_vel[j*6+k+3] - brFSIData[iTurb][fast::STATE_NM1].bld_vel[j*6+k+3]);
                    }
                }

                int nBlades = turbineData[iTurb].numBlades;
                for (int j=0; j < nBlades; j++) {
                    brFSIData[iTurb][fast::STATE_NP1].bld_pitch[j] = brFSIData[iTurb][fast::STATE_NM2].bld_pitch[j] + 3.0*(brFSIData[iTurb][fast::STATE_N].bld_pitch[j] - brFSIData[iTurb][fast::STATE_NM1].bld_pitch[j]);
                    extrapRotation(&brFSIData[iTurb][fast::STATE_NM2].bld_root_def[j*6+3], &brFSIData[iTurb][fast::STATE_NM1].bld_root_def[j*6+3], &brFSIData[iTurb][fast::STATE_N].bld_root_def[j*6+3], &brFSIData[iTurb][fast::STATE_NP1].bld_root_def[j*6+3]);
                    for (int k=0; k < 3; k++) {
                        brFSIData[iTurb][fast::STATE_NP1].bld_root_def[j*6+k] = brFSIData[iTurb][fast::STATE_NM2].bld_root_def[j*6+k] + 3.0*(brFSIData[iTurb][fast::STATE_N].bld_root_def[j*6+k] - brFSIData[iTurb][fast::STATE_NM1].bld_root_def[j*6+k]);
                    }
                }

                for (int k=0; k < 3; k++) {
                    brFSIData[iTurb][fast::STATE_NP1].hub_def[k] = brFSIData[iTurb][fast::STATE_NM2].hub_def[k] + 3.0*(brFSIData[iTurb][fast::STATE_N].hub_def[k] - brFSIData[iTurb][fast::STATE_NM1].hub_def[k]);
                    extrapRotation(&brFSIData[iTurb][fast::STATE_NM2].hub_def[3], &brFSIData[iTurb][fast::STATE_NM1].hub_def[3], &brFSIData[iTurb][fast::STATE_N].hub_def[3], &brFSIData[iTurb][fast::STATE_NP1].hub_def[3]);
                    brFSIData[iTurb][fast::STATE_NP1].hub_vel[k] = brFSIData[iTurb][fast::STATE_NM2].hub_vel[k] + 3.0*(brFSIData[iTurb][fast::STATE_N].hub_vel[k] - brFSIData[iTurb][fast::STATE_NM1].hub_vel[k]);
                    brFSIData[iTurb][fast::STATE_NP1].hub_vel[k+3] = brFSIData[iTurb][fast::STATE_NM2].hub_vel[k+3] + 3.0*(brFSIData[iTurb][fast::STATE_N].hub_vel[k+3] - brFSIData[iTurb][fast::STATE_NM1].hub_vel[k+3]);
                    brFSIData[iTurb][fast::STATE_NP1].nac_def[k] = brFSIData[iTurb][fast::STATE_NM2].nac_def[k] + 3.0*(brFSIData[iTurb][fast::STATE_N].nac_def[k] - brFSIData[iTurb][fast::STATE_NM1].nac_def[k]);
                    extrapRotation(&brFSIData[iTurb][fast::STATE_NM2].nac_def[3], &brFSIData[iTurb][fast::STATE_NM1].nac_def[3], &brFSIData[iTurb][fast::STATE_N].nac_def[3], &brFSIData[iTurb][fast::STATE_NP1].nac_def[3]);
                    brFSIData[iTurb][fast::STATE_NP1].nac_vel[k] = brFSIData[iTurb][fast::STATE_NM2].nac_vel[k] + 3.0*(brFSIData[iTurb][fast::STATE_N].nac_vel[k] - brFSIData[iTurb][fast::STATE_NM1].nac_vel[k]);
                    brFSIData[iTurb][fast::STATE_NP1].nac_vel[k+3] = brFSIData[iTurb][fast::STATE_NM2].nac_vel[k+3] + 3.0*(brFSIData[iTurb][fast::STATE_N].nac_vel[k+3] - brFSIData[iTurb][fast::STATE_NM1].nac_vel[k+3]);
                }

                int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
                for (int j=0; j < nPtsTwr; j++) {
                    extrapRotation(&brFSIData[iTurb][fast::STATE_NM2].twr_def[j*6+3],&brFSIData[iTurb][fast::STATE_NM1].twr_def[j*6+3],&brFSIData[iTurb][fast::STATE_N].twr_def[j*6+3], &brFSIData[iTurb][fast::STATE_NP1].twr_def[j*6+3]);
                    for (int k = 0; k < 3; k++) {
                        brFSIData[iTurb][fast::STATE_NP1].twr_def[j*6+k] = brFSIData[iTurb][fast::STATE_NM2].twr_def[j*6+k] + 3.0*(brFSIData[iTurb][fast::STATE_N].twr_def[j*6+k] - brFSIData[iTurb][fast::STATE_NM1].twr_def[j*6+k]);
                        brFSIData[iTurb][fast::STATE_NP1].twr_vel[j*6+k] = brFSIData[iTurb][fast::STATE_NM2].twr_vel[j*6+k] + 3.0*(brFSIData[iTurb][fast::STATE_N].twr_vel[j*6+k] - brFSIData[iTurb][fast::STATE_NM1].twr_vel[j*6+k]);
                        brFSIData[iTurb][fast::STATE_NP1].twr_vel[j*6+k+3] = brFSIData[iTurb][fast::STATE_NM2].twr_vel[j*6+k+3] + 3.0*(brFSIData[iTurb][fast::STATE_N].twr_vel[j*6+k+3] - brFSIData[iTurb][fast::STATE_NM1].twr_vel[j*6+k+3]);
                    }
                }
            }

        }
    }
}

void fast::OpenFAST::prework() {

    if (nSubsteps_ > 1) {

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_CFD_Store_SubStep(&iTurb, &nt_global, &ErrStat, ErrMsg) ;
            checkError(ErrStat, ErrMsg);
        }

    } else {

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_CFD_Prework(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
        }
    }
}

void fast::OpenFAST::update_states_driver_time_step(bool writeFiles) {

    if (firstPass_)
        prework();

    if (nSubsteps_ > 1) {

        if (!firstPass_) {
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                FAST_CFD_Reset_SubStep(&iTurb, &nSubsteps_, &ErrStat, ErrMsg);
                checkError(ErrStat, ErrMsg);
            }
        }

        for (int iSubstep=1; iSubstep < nSubsteps_+1; iSubstep++) {
            double ss_time = double(iSubstep)/double(nSubsteps_);
            step(ss_time);
        }

        get_data_from_openfast(fast::STATE_NP1);

        if (writeFiles) {
            if ( isDebug() ) {
                for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                    std::ofstream fastcpp_velocity_file;
                    fastcpp_velocity_file.open("fastcpp_residual." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv", std::ios_base::app) ;
                    fastcpp_velocity_file << "Time step " << nt_global << " Velocity residual at the force nodes = " << velForceNodeData[iTurb][fast::STATE_NP1].vel_force_resid << std::endl ;
                    fastcpp_velocity_file << "          " << nt_global << " Position residual at the force nodes = " << velForceNodeData[iTurb][fast::STATE_NP1].x_force_resid << std::endl ;
                    fastcpp_velocity_file << "          " << nt_global << " Force residual at the force nodes = " << velForceNodeData[iTurb][fast::STATE_NP1].force_resid << std::endl ;
                    fastcpp_velocity_file.close() ;
                }
            }
        }

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            velForceNodeData[iTurb][fast::STATE_NP1].x_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].vel_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].x_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].xdot_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].orient_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].vel_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].force_resid = 0.0;
        }
    } else {

        send_data_to_openfast(fast::STATE_NP1);

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_CFD_UpdateStates(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);

            // Compute the force from the nacelle only if the drag coefficient is
            //   greater than zero
            if (get_nacelleCdLoc(iTurb) > 0.) {

                calc_nacelle_force (

                                    extinfw_o_t_FAST[iTurb].u[0],
                                    extinfw_o_t_FAST[iTurb].v[0],
                                    extinfw_o_t_FAST[iTurb].w[0],
                                    get_nacelleCdLoc(iTurb),
                                    get_nacelleAreaLoc(iTurb),
                                    get_airDensityLoc(iTurb),
                                    extinfw_i_f_FAST[iTurb].fx[0],
                                    extinfw_i_f_FAST[iTurb].fy[0],
                                    extinfw_i_f_FAST[iTurb].fz[0]

                                    );

            }

        }

        get_data_from_openfast(fast::STATE_NP1);

        if ( writeFiles ) {
            if ( isDebug() ) {
                for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                    std::ofstream fastcpp_velocity_file;
                    fastcpp_velocity_file.open("fastcpp_residual." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv", std::ios_base::app) ;
                    fastcpp_velocity_file << "Time step " << nt_global << " Velocity residual at the force nodes = " << velForceNodeData[iTurb][fast::STATE_NP1].vel_force_resid << std::endl ;
                    fastcpp_velocity_file << "          " << nt_global << " Position residual at the force nodes = " << velForceNodeData[iTurb][fast::STATE_NP1].x_force_resid << std::endl ;
                    fastcpp_velocity_file << "          " << nt_global << " Force residual at the force nodes = " << velForceNodeData[iTurb][fast::STATE_NP1].force_resid << std::endl ;
                    fastcpp_velocity_file.close() ;
                }
            }
        }

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            velForceNodeData[iTurb][fast::STATE_NP1].x_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].vel_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].x_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].xdot_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].orient_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].vel_force_resid = 0.0;
            velForceNodeData[iTurb][fast::STATE_NP1].force_resid = 0.0;
        }

    }

    if (writeFiles) {
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            if (turbineData[iTurb].inflowType == 2)
                writeVelocityData(iTurb, nt_global, nlinIter_);
        }
    }

    firstPass_ = false;
    nlinIter_ +=1 ;
}

void fast::OpenFAST::advance_to_next_driver_time_step(bool writeFiles) {

    if (nSubsteps_ > 1) {
        //Nothing to do here

    } else {

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_CFD_AdvanceToNextTimeStep(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
        }

    }

    nt_global = nt_global + nSubsteps_;


    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        FAST_CFD_WriteOutput(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
    }

    set_state_from_state(fast::STATE_NM1, fast::STATE_NM2);
    set_state_from_state(fast::STATE_N, fast::STATE_NM1);
    set_state_from_state(fast::STATE_NP1, fast::STATE_N);

    if (writeFiles) {
      for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
          int tStepRatio = time_step_ratio(dtFAST, dtDriver);
          if ( (restartFreq_*tStepRatio > 0) && (((nt_global - ntStart) % (restartFreq_*tStepRatio)) == 0 )  && (nt_global != ntStart) ) {
              turbineData[iTurb].FASTRestartFileName = " "; // if blank, it will use FAST convention <RootName>.nt_global
              std::string tmpRstFileRoot{turbineData[iTurb].FASTRestartFileName};
              tmpRstFileRoot.resize(INTERFACE_STRING_LENGTH, ' ');
              FAST_CreateCheckpoint(&iTurb, tmpRstFileRoot.c_str(), &ErrStat, ErrMsg);
              checkError(ErrStat, ErrMsg);
              writeRestartFile(iTurb, nt_global);
          }

          if ( (((nt_global - ntStart) % (outputFreq_ * tStepRatio) ) == 0 )  && (nt_global != ntStart) ) {
              writeOutputFile(iTurb, nt_global);
          }
      }

    }

    nlinIter_ = 0;
    firstPass_ = true ; // Set firstPass_ to true for the next time step
}

void fast::OpenFAST::calc_nacelle_force(const float & u, const float & v, const float & w, const float & cd, const float & area, const float & rho, float & fx, float & fy, float & fz) {
    // Calculate the force on the nacelle (fx,fy,fz) given the
    //   velocity sampled at the nacelle point (u,v,w),
    //   drag coefficient 'cd' and nacelle area 'area'
    // The velocity magnitude
    float Vmag = std::sqrt(u * u + v * v + w * w);

    // Velocity correction based on Martinez-Tossas PhD Thesis 2017
    // The correction samples the velocity at the center of the
    // Gaussian kernel and scales it to obtain the inflow velocity
    float epsilon_d = std::sqrt(2.0 / M_PI * cd * area);
    float correction = 1. / (1.0 - cd * area / (4.0 * M_PI * epsilon_d * epsilon_d));

    // Compute the force for each velocity component
    fx = rho * 1./2. * cd * area * Vmag * u * correction * correction;
    fy = rho * 1./2. * cd * area * Vmag * v * correction * correction;
    fz = rho * 1./2. * cd * area * Vmag * w * correction * correction;

}

/* A version of step allowing for sub-timesteps when the driver program has a larger time step than OpenFAST */
void fast::OpenFAST::step(double ss_time) {

    /* ******************************
       set inputs from this code and call FAST:
       ********************************* */

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        // this advances the states, calls CalcOutput, and solves for next inputs. Predictor-corrector loop is imbeded here:
        FAST_CFD_Prework(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        send_data_to_openfast(ss_time);
        FAST_CFD_UpdateStates(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        FAST_CFD_AdvanceToNextTimeStep(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);

    }

}

void fast::OpenFAST::step(bool writeFiles) {

    /* ******************************
       set inputs from this code and call FAST:
       ********************************* */

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        // this advances the states, calls CalcOutput, and solves for next inputs. Predictor-corrector loop is imbeded here:
        // (note CFD could do subcycling around this step)

        if (turbineData[iTurb].inflowType == 2)

            writeVelocityData(iTurb, nt_global, 0);

            if (writeFiles) {
                if ( isDebug() && (turbineData[iTurb].inflowType == 2) ) {

                    std::ofstream fastcpp_velocity_file;
                    fastcpp_velocity_file.open("fastcpp_velocity." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
                    fastcpp_velocity_file << "# x, y, z, Vx, Vy, Vz" << std::endl ;
                    for (int iNode=0; iNode < get_numVelPtsLoc(iTurb); iNode++) {
                        fastcpp_velocity_file << extinfw_i_f_FAST[iTurb].pxVel[iNode] << ", " << extinfw_i_f_FAST[iTurb].pyVel[iNode] << ", " << extinfw_i_f_FAST[iTurb].pzVel[iNode] << ", " << extinfw_o_t_FAST[iTurb].u[iNode] << ", " << extinfw_o_t_FAST[iTurb].v[iNode] << ", " << extinfw_o_t_FAST[iTurb].w[iNode] << " " << std::endl ;
                    }
                    fastcpp_velocity_file.close() ;
                }
            }

        FAST_CFD_Prework(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        send_data_to_openfast(fast::STATE_NP1);
        FAST_CFD_UpdateStates(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        get_data_from_openfast(fast::STATE_NP1);
        FAST_CFD_AdvanceToNextTimeStep(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);

        // Compute the force from the nacelle only if the drag coefficient is
        //   greater than zero
        if (get_nacelleCdLoc(iTurb) > 0.) {

            calc_nacelle_force (

                extinfw_o_t_FAST[iTurb].u[0],
                extinfw_o_t_FAST[iTurb].v[0],
                extinfw_o_t_FAST[iTurb].w[0],
                get_nacelleCdLoc(iTurb),
                get_nacelleAreaLoc(iTurb),
                get_airDensityLoc(iTurb),
                extinfw_i_f_FAST[iTurb].fx[0],
                extinfw_i_f_FAST[iTurb].fy[0],
                extinfw_i_f_FAST[iTurb].fz[0]

                );

        }

        if (writeFiles) {
            if ( isDebug() && (turbineData[iTurb].inflowType == 2) ) {
                std::ofstream actuatorForcesFile;
                actuatorForcesFile.open("actuator_forces." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
                actuatorForcesFile << "# x, y, z, fx, fy, fz" << std::endl ;
                for (int iNode=0; iNode < get_numForcePtsLoc(iTurb); iNode++) {
                    actuatorForcesFile << extinfw_i_f_FAST[iTurb].pxForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].pyForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].pzForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].fx[iNode] << ", " << extinfw_i_f_FAST[iTurb].fy[iNode] << ", " << extinfw_i_f_FAST[iTurb].fz[iNode] << " " << std::endl ;
                }
                actuatorForcesFile.close() ;
            }
        }

    }

    nt_global = nt_global + 1;

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        FAST_CFD_WriteOutput(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
    }

    if (writeFiles) {
        int tStepRatio = time_step_ratio(dtFAST, dtDriver);
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            if ( (((nt_global - ntStart) % (restartFreq_ * tStepRatio)) == 0 )  && (nt_global != ntStart) ) {
                turbineData[iTurb].FASTRestartFileName = " "; // if blank, it will use FAST convention <RootName>.nt_global
                std::string tmpRstFileRoot{turbineData[iTurb].FASTRestartFileName};
                tmpRstFileRoot.resize(INTERFACE_STRING_LENGTH, ' ');
                FAST_CreateCheckpoint(&iTurb, tmpRstFileRoot.c_str(), &ErrStat, ErrMsg);
                checkError(ErrStat, ErrMsg);
                writeRestartFile(iTurb, nt_global);
            }
            if ( (((nt_global - ntStart) % (outputFreq_ * tStepRatio) ) == 0 )  && (nt_global != ntStart) ) {
                writeOutputFile(iTurb, nt_global);
            }
        }
    }

}

void fast::OpenFAST::setInputs(const fast::fastInputs & fi ) {

    mpiComm = fi.comm;

    MPI_Comm_rank(mpiComm, &worldMPIRank);
    MPI_Comm_group(mpiComm, &worldMPIGroup);

    nTurbinesGlob = fi.nTurbinesGlob;

    if (nTurbinesGlob > 0) {
        dryRun = fi.dryRun;

        debug = fi.debug;

        tStart = fi.tStart;
        simStart = fi.simStart;
        restartFreq_ = fi.restartFreq;
        outputFreq_ = fi.outputFreq;
        tMax = fi.tMax;
        dtDriver = fi.dtDriver;

        ///TODO: Check if this is right and necessary
        // if (simStart == fast::restartDriverInitFAST) {
        //     nt_global = 0;
        // } else {
        //     nt_global = ntStart;
        // }

        globTurbineData.resize(nTurbinesGlob);
        globTurbineData = fi.globTurbineData;

    } else {
        throw std::runtime_error("Number of turbines < 0 ");
    }
}

int fast::OpenFAST::checkAndSetSubsteps() {

    if ( nTurbinesProc > 0) {
        if (dtDriver > 0) {
            dtFAST = turbineData[0].dt;
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                if (dtFAST != turbineData[iTurb].dt) {
                    throw std::runtime_error("All turbines don't have the same time step ");
                }
            }
            if (dtFAST > 0) {
                int tStepRatio = time_step_ratio(dtFAST, dtDriver);
                if (std::abs(dtDriver - tStepRatio * dtFAST) < 0.001) {// TODO: Fix arbitrary number 0.001
                    nSubsteps_ = tStepRatio;
                    return 1;
                } else {
                    return -1;
                }
            } else {
                throw std::runtime_error("FAST time step is zero");
            }

        } else {
            throw std::runtime_error("Driver time step is not set or set to zero");
        }
    } else {
        return 1;
    }

}


void fast::OpenFAST::setDriverTimeStep(double dt_driver) {
    dtDriver = dt_driver;
}


void fast::OpenFAST::setDriverCheckpoint(int nt_checkpoint_driver) {

    if (nTurbinesProc > 0) {
        if (nSubsteps_ > 0) {
            restartFreq_ = nt_checkpoint_driver;
        } else {
            throw std::runtime_error("Trying to set driver checkpoint when nSubsteps_ is zero. Set driver time step first may be?");
        }
    }
}

void fast::OpenFAST::get_turbineParams(int iTurbGlob, turbineDataType & turbData) {

    //TODO: Figure out a better copy operator for the turbineDataType struct
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    turbData.TurbID = turbineData[iTurbLoc].TurbID;
    turbData.FASTInputFileName = turbineData[iTurbLoc].FASTInputFileName;
    turbData.FASTRestartFileName = turbineData[iTurbLoc].FASTRestartFileName;
    if(turbineData[iTurbLoc].TurbineBasePos.size() > 0) {
        turbData.TurbineBasePos.resize(turbineData[iTurbLoc].TurbineBasePos.size());
        for(int i=0; i < turbineData[iTurbLoc].TurbineBasePos.size(); i++)
            turbData.TurbineBasePos[i] = turbineData[iTurbLoc].TurbineBasePos[i];
    }
    if(turbineData[iTurbLoc].TurbineHubPos.size() > 0) {
        turbData.TurbineHubPos.resize(turbineData[iTurbLoc].TurbineHubPos.size());
        for(int i=0; i < turbineData[iTurbLoc].TurbineHubPos.size(); i++)
            turbData.TurbineHubPos[i] = turbineData[iTurbLoc].TurbineHubPos[i];
    }
    turbData.sType = turbineData[iTurbLoc].sType;
    turbData.numBlades = turbineData[iTurbLoc].numBlades;
    turbData.numVelPtsBlade = turbineData[iTurbLoc].numVelPtsBlade;
    turbData.numVelPtsTwr = turbineData[iTurbLoc].numVelPtsTwr;
    turbData.numVelPts = turbineData[iTurbLoc].numVelPts;
    turbData.numForcePtsBlade = turbineData[iTurbLoc].numForcePtsBlade;
    turbData.numForcePtsTwr = turbineData[iTurbLoc].numForcePtsTwr;
    turbData.numForcePts = turbineData[iTurbLoc].numForcePts;
    turbData.inflowType = turbineData[iTurbLoc].inflowType;
    turbData.nacelle_cd = turbineData[iTurbLoc].nacelle_cd;
    turbData.nacelle_area = turbineData[iTurbLoc].nacelle_area;
    turbData.air_density = turbineData[iTurbLoc].air_density;
    turbData.nBRfsiPtsBlade.resize(turbData.numBlades);
    turbData.nTotBRfsiPtsBlade = 0;
    for (int i=0; i < turbData.numBlades; i++) {
        turbData.nBRfsiPtsBlade[i] = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        turbData.nTotBRfsiPtsBlade += turbData.nBRfsiPtsBlade[i];
    }
    turbData.nBRfsiPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
    turbData.azBlendMean = turbineData[iTurbLoc].azBlendMean;
    turbData.azBlendDelta = turbineData[iTurbLoc].azBlendDelta;

}


void fast::OpenFAST::checkError(const int ErrStat, const char * ErrMsg)
{
    if (ErrStat != ErrID_None){

        if (ErrStat >= AbortErrLev){
            throw std::runtime_error(std::string(ErrMsg));
        } else {
            std::cout << "Warning from OpenFAST: " << ErrMsg << std::endl;
        }
    }
}

// Actuator stuff

void fast::OpenFAST::setExpLawWindSpeed(double t){

    double sinOmegat = 0.1 * std::sin(10.0*t);
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        // routine sets the u-v-w wind speeds used in FAST
        int nVelPts = get_numVelPts(iTurb);
        int iTurbGlob = turbineMapProcToGlob[iTurb];
        for (int j = 0; j < nVelPts; j++){
            std::vector<double> coords(3,0.0);
            std::vector<double> tmpVel(3,0.0);
            getVelNodeCoordinates(coords, j, iTurbGlob, fast::STATE_NP1);
            tmpVel[0] = (float) 10.0*pow((coords[2] / 90.0), 0.2) + sinOmegat; // 0.2 power law wind profile using reference 10 m/s at 90 meters + a perturbation
            setVelocity(tmpVel, j, iTurbGlob);
        }
    }
}

void fast::OpenFAST::getApproxHubPos(double* currentCoords, int iTurbGlob, int nSize) {
  assert(nSize==3);
  // Get hub position of Turbine 'iTurbGlob'
  for(int i =0; i<nSize; ++i){
    currentCoords[i] = globTurbineData[iTurbGlob].TurbineHubPos[i];
  }
}

void fast::OpenFAST::getHubPos(double* currentCoords, int iTurbGlob, fast::timeStep t, int nSize) {
    assert(nSize==3);
    // Get hub position of Turbine 'iTurbGlob'
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int i=0; i < nSize; i++)
        currentCoords[i] = velForceNodeData[iTurbLoc][t].x_force[i] + turbineData[iTurbLoc].TurbineBasePos[i] ;
}

void fast::OpenFAST::getHubShftDir(double* hubShftVec, int iTurbGlob, fast::timeStep t, int nSize) {
    assert(nSize==3);
    // Get hub shaft direction of current turbine - pointing downwind
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int i=0; i < nSize; i++)
        hubShftVec[i] = velForceNodeData[iTurbLoc][t].orient_force[i*3] ;
}


void fast::OpenFAST::getVelNodeCoordinates(double* currentCoords, int iNode, int iTurbGlob, fast::timeStep t, int nSize) {
    assert(nSize==3);
    // Set coordinates at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbLoc);
    for (int i=0; i < nSize; i++)
        currentCoords[i] = velForceNodeData[iTurbLoc][t].x_vel[iNode*3+i] + turbineData[iTurbLoc].TurbineBasePos[i] ;

}

void fast::OpenFAST::getForceNodeCoordinates(double* currentCoords, int iNode, int iTurbGlob, fast::timeStep t, int nSize) {
    assert(nSize==3);
    // Set coordinates at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int i=0; i < nSize; i++)
        currentCoords[i] = velForceNodeData[iTurbLoc][t].x_force[iNode*3+i] + turbineData[iTurbLoc].TurbineBasePos[i] ;

}

void fast::OpenFAST::getForceNodeOrientation(double* currentOrientation, int iNode, int iTurbGlob, fast::timeStep t, int nSize) {
    assert(nSize==9);
    // Set orientation at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    for(int i=0; i < nSize; i++)
        currentOrientation[i] = velForceNodeData[iTurbLoc][t].orient_force[iNode*9+i] ;
}

void fast::OpenFAST::getRelativeVelForceNode(double* currentVelocity, int iNode, int iTurbGlob, fast::timeStep t, int nSize) {
    assert(nSize==3);
    // Get relative velocity at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    currentVelocity[0] = velForceNodeData[iTurbLoc][t].vel_force[iNode*3+0] - velForceNodeData[iTurbLoc][t].xdot_force[iNode*3+0];
    currentVelocity[1] = velForceNodeData[iTurbLoc][t].vel_force[iNode*3+1] - velForceNodeData[iTurbLoc][t].xdot_force[iNode*3+1];
    currentVelocity[2] = velForceNodeData[iTurbLoc][t].vel_force[iNode*3+2] - velForceNodeData[iTurbLoc][t].xdot_force[iNode*3+2];
}

void fast::OpenFAST::getForce(double* currentForce, int iNode, int iTurbGlob, fast::timeStep t, int nSize) {
    assert(nSize==3);
    // Set forces at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    for(int i=0; i < nSize; i++)
        currentForce[i] = -velForceNodeData[iTurbLoc][t].force[iNode*3+i] ;
}

double fast::OpenFAST::getRHloc(int iNode, int iTurbGlob) {

    // Return radial location/height along blade/tower at current node of current turbine
    // Inactive for now
    return -1.0;

}

double fast::OpenFAST::getChord(int iNode, int iTurbGlob) {
    // Return blade chord/tower diameter at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    if (turbineData[iTurbLoc].sType == EXTINFLOW) {
        for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
        return extinfw_i_f_FAST[iTurbLoc].forceNodesChord[iNode] ;
    } else {
        return -1.0;
    }

}

void fast::OpenFAST::setVelocity(double* currentVelocity, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Set velocity at current node of current turbine -
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbLoc);
    for(int k=0; k < 3; k++) {
        velForceNodeData[iTurbLoc][fast::STATE_NP1].vel_vel_resid += (velForceNodeData[iTurbLoc][fast::STATE_NP1].vel_vel[iNode*3+k] - currentVelocity[k])*(velForceNodeData[iTurbLoc][fast::STATE_NP1].vel_vel[iNode*3+k] - currentVelocity[k]);
        velForceNodeData[iTurbLoc][fast::STATE_NP1].vel_vel[iNode*3+k] = currentVelocity[k];
    }

    // Put this in send_data_to_openfast
    // extinfw_o_t_FAST[iTurbLoc].u[iNode] = currentVelocity[0];
    // extinfw_o_t_FAST[iTurbLoc].v[iNode] = currentVelocity[1];
    // extinfw_o_t_FAST[iTurbLoc].w[iNode] = currentVelocity[2];
}

void fast::OpenFAST::setVelocityForceNode(double* currentVelocity, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Set velocity at current node of current turbine -
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    for(int k=0; k < nSize; k++) {
        velForceNodeData[iTurbLoc][fast::STATE_NP1].vel_force_resid += (velForceNodeData[iTurbLoc][fast::STATE_NP1].vel_force[iNode*3+k] - currentVelocity[k])*(velForceNodeData[iTurbLoc][fast::STATE_NP1].vel_force[iNode*3+k] - currentVelocity[k]);
        velForceNodeData[iTurbLoc][fast::STATE_NP1].vel_force[iNode*3+k] = currentVelocity[k];
    }
}

void fast::OpenFAST::interpolateVel_ForceToVelNodes() {

    // Interpolates the velocity from the force nodes to the velocity nodes
    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        // Hub location

        for (int k=0; k < 3; k++) {
            double tmp = velForceNodeData[iTurb][fast::STATE_NP1].vel_force[k];
            velForceNodeData[iTurb][fast::STATE_NP1].vel_vel_resid += (velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[k] - tmp)*(velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[k] - tmp);
            velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[k] = tmp;
        }

        // Do the blades first
        int nBlades = get_numBladesLoc(iTurb);
        for(int iBlade=0; iBlade < nBlades; iBlade++) {

            // Create interpolating parameter - Distance from hub
            int nForcePtsBlade = get_numForcePtsBladeLoc(iTurb);
            std::vector<double> rDistForce(nForcePtsBlade) ;
            for(int j=0; j < nForcePtsBlade; j++) {
                int iNodeForce = 1 + iBlade * nForcePtsBlade + j ; //The number of actuator force points is always the same for all blades
                rDistForce[j] = sqrt(
                    (extinfw_i_f_FAST[iTurb].pxForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pxForce[0])*(extinfw_i_f_FAST[iTurb].pxForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pxForce[0])
                    + (extinfw_i_f_FAST[iTurb].pyForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pyForce[0])*(extinfw_i_f_FAST[iTurb].pyForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pyForce[0])
                    + (extinfw_i_f_FAST[iTurb].pzForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pzForce[0])*(extinfw_i_f_FAST[iTurb].pzForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pzForce[0])
                    );
            }

            // Interpolate to the velocity nodes
            int nVelPtsBlade = get_numVelPtsBladeLoc(iTurb);
            for(int j=0; j < nVelPtsBlade; j++) {
                int iNodeVel = 1 + iBlade * nVelPtsBlade + j ; //Assumes the same number of velocity (Aerodyn) nodes for all blades
                double rDistVel = sqrt(
                    (extinfw_i_f_FAST[iTurb].pxVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pxVel[0])*(extinfw_i_f_FAST[iTurb].pxVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pxVel[0])
                    + (extinfw_i_f_FAST[iTurb].pyVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pyVel[0])*(extinfw_i_f_FAST[iTurb].pyVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pyVel[0])
                    + (extinfw_i_f_FAST[iTurb].pzVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pzVel[0])*(extinfw_i_f_FAST[iTurb].pzVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pzVel[0])
                    );
                //Find nearest two force nodes
                int jForceLower = 0;
                while ( (rDistForce[jForceLower+1] < rDistVel) && ( jForceLower < (nForcePtsBlade-2)) )   {
                    jForceLower = jForceLower + 1;
                }
                int iNodeForceLower = 1 + iBlade * nForcePtsBlade + jForceLower ;
                double rInterp = (rDistVel - rDistForce[jForceLower])/(rDistForce[jForceLower+1]-rDistForce[jForceLower]);

                for (int k=0; k < 3; k++) {
                    double tmp = velForceNodeData[iTurb][fast::STATE_NP1].vel_force[iNodeForceLower*3+k] + rInterp * (velForceNodeData[iTurb][fast::STATE_NP1].vel_force[(iNodeForceLower+1)*3+k] - velForceNodeData[iTurb][fast::STATE_NP1].vel_force[iNodeForceLower*3+k]);
                    velForceNodeData[iTurb][fast::STATE_NP1].vel_vel_resid += (velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+k] - tmp)*(velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+k] - tmp);
                    velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+k] = tmp;
                }
            }
        }

        // Now the tower if present and used
        int nVelPtsTower = get_numVelPtsTwrLoc(iTurb);
        if ( nVelPtsTower > 0 ) {

            // Create interpolating parameter - Distance from first node from ground
            int nForcePtsTower = get_numForcePtsTwrLoc(iTurb);
            std::vector<double> hDistForce(nForcePtsTower) ;
            int iNodeBotTowerForce = 1 + nBlades * get_numForcePtsBladeLoc(iTurb); // The number of actuator force points is always the same for all blades
            for(int j=0; j < nForcePtsTower; j++) {
                int iNodeForce = iNodeBotTowerForce + j ;
                hDistForce[j] = sqrt(
                    (extinfw_i_f_FAST[iTurb].pxForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pxForce[iNodeBotTowerForce])*(extinfw_i_f_FAST[iTurb].pxForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pxForce[iNodeBotTowerForce])
                    + (extinfw_i_f_FAST[iTurb].pyForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pyForce[iNodeBotTowerForce])*(extinfw_i_f_FAST[iTurb].pyForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pyForce[iNodeBotTowerForce])
                    + (extinfw_i_f_FAST[iTurb].pzForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pzForce[iNodeBotTowerForce])*(extinfw_i_f_FAST[iTurb].pzForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pzForce[iNodeBotTowerForce])
                    );
            }

            int iNodeBotTowerVel = 1 + nBlades * get_numVelPtsBladeLoc(iTurb); // Assumes the same number of velocity (Aerodyn) nodes for all blades
            for(int j=0; j < nVelPtsTower; j++) {
                int iNodeVel = iNodeBotTowerVel + j ;
                double hDistVel = sqrt(
                    (extinfw_i_f_FAST[iTurb].pxVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pxVel[iNodeBotTowerVel])*(extinfw_i_f_FAST[iTurb].pxVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pxVel[iNodeBotTowerVel])
                    + (extinfw_i_f_FAST[iTurb].pyVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pyVel[iNodeBotTowerVel])*(extinfw_i_f_FAST[iTurb].pyVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pyVel[iNodeBotTowerVel])
                    + (extinfw_i_f_FAST[iTurb].pzVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pzVel[iNodeBotTowerVel])*(extinfw_i_f_FAST[iTurb].pzVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pzVel[iNodeBotTowerVel])
                    );
                //Find nearest two force nodes
                int jForceLower = 0;
                while ( (hDistForce[jForceLower+1] < hDistVel) && ( jForceLower < (nForcePtsTower-2)) )   {
                    jForceLower = jForceLower + 1;
                }
                int iNodeForceLower = iNodeBotTowerForce + jForceLower ;
                double rInterp = (hDistVel - hDistForce[jForceLower])/(hDistForce[jForceLower+1]-hDistForce[jForceLower]);
                for (int k=0; k < 3; k++) {
                    double tmp = velForceNodeData[iTurb][fast::STATE_NP1].vel_force[iNodeForceLower*3+k] + rInterp * (velForceNodeData[iTurb][fast::STATE_NP1].vel_force[(iNodeForceLower+1)*3+k] - velForceNodeData[iTurb][fast::STATE_NP1].vel_force[iNodeForceLower*3+k]);
                    velForceNodeData[iTurb][fast::STATE_NP1].vel_vel_resid += (velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+k] - tmp)*(velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+k] - tmp);
                    velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+k] = tmp;
                }
            }
        }

    } // End loop over turbines

}

void fast::OpenFAST::computeTorqueThrust(int iTurbGlob, double* torque, double* thrust, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob) ;
    if (turbineData[iTurbLoc].sType != EXTINFLOW)
        return;

    //Compute the torque and thrust based on the forces at the actuator nodes
    std::vector<double> relLoc(3,0.0);
    std::vector<double> rPerpShft(3);
    thrust[0] = 0.0; thrust[1] = 0.0; thrust[2] = 0.0;
    torque[0] = 0.0; torque[1] = 0.0; torque[2] = 0.0;

    std::vector<double> hubShftVec(3);
    getHubShftDir(hubShftVec, iTurbGlob, fast::STATE_NP1);

    int nfpts = get_numForcePtsBlade(iTurbLoc);
    for (int k=0; k < get_numBladesLoc(iTurbLoc); k++) {
        for (int j=0; j < nfpts; j++) {
            int iNode = 1 + nfpts*k + j ;

            thrust[0] = thrust[0] + extinfw_i_f_FAST[iTurbLoc].fx[iNode] ;
            thrust[1] = thrust[1] + extinfw_i_f_FAST[iTurbLoc].fy[iNode] ;
            thrust[2] = thrust[2] + extinfw_i_f_FAST[iTurbLoc].fz[iNode] ;

            relLoc[0] = extinfw_i_f_FAST[iTurbLoc].pxForce[iNode] - extinfw_i_f_FAST[iTurbLoc].pxForce[0] ;
            relLoc[1] = extinfw_i_f_FAST[iTurbLoc].pyForce[iNode] - extinfw_i_f_FAST[iTurbLoc].pyForce[0];
            relLoc[2] = extinfw_i_f_FAST[iTurbLoc].pzForce[iNode] - extinfw_i_f_FAST[iTurbLoc].pzForce[0];

            double rDotHubShftVec = relLoc[0]*hubShftVec[0] + relLoc[1]*hubShftVec[1] + relLoc[2]*hubShftVec[2];
            for (int j=0; j < 3; j++)  rPerpShft[j] = relLoc[j] - rDotHubShftVec * hubShftVec[j];

            torque[0] = torque[0] + rPerpShft[1] * extinfw_i_f_FAST[iTurbLoc].fz[iNode] - rPerpShft[2] * extinfw_i_f_FAST[iTurbLoc].fy[iNode] + extinfw_i_f_FAST[iTurbLoc].momentx[iNode] ;
            torque[1] = torque[1] + rPerpShft[2] * extinfw_i_f_FAST[iTurbLoc].fx[iNode] - rPerpShft[0] * extinfw_i_f_FAST[iTurbLoc].fz[iNode] + extinfw_i_f_FAST[iTurbLoc].momenty[iNode] ;
            torque[2] = torque[2] + rPerpShft[0] * extinfw_i_f_FAST[iTurbLoc].fy[iNode] - rPerpShft[1] * extinfw_i_f_FAST[iTurbLoc].fx[iNode] + extinfw_i_f_FAST[iTurbLoc].momentz[iNode] ;

        }
    }
}

fast::ActuatorNodeType fast::OpenFAST::getVelNodeType(int iTurbGlob, int iNode) {
    // Return the type of velocity node for the given node number. The node ordering (from FAST) is
    // Node 0 - Hub node
    // Blade 1 nodes
    // Blade 2 nodes
    // Blade 3 nodes
    // Tower nodes

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbGlob);
    if (iNode) {
        if ( (iNode + 1 - (get_numVelPts(iTurbLoc) - get_numVelPtsTwr(iTurbLoc)) ) > 0) {
            return TOWER;
        }
        else {
            return BLADE;
        }
    }
    else {
        return HUB;
    }

}

fast::ActuatorNodeType fast::OpenFAST::getForceNodeType(int iTurbGlob, int iNode) {
    // Return the type of actuator force node for the given node number. The node ordering (from FAST) is
    // Node 0 - Hub node
    // Blade 1 nodes
    // Blade 2 nodes
    // Blade 3 nodes
    // Tower nodes

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbGlob);
    if (iNode) {
        if ( (iNode + 1 - (get_numForcePts(iTurbLoc) - get_numForcePtsTwr(iTurbLoc)) ) > 0) {
            return TOWER;
        }
        else {
            return BLADE;
        }
    }
    else {
        return HUB;
    }
}

void fast::OpenFAST::allocateMemory_preInit() {

    for (int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
        if (dryRun) {
            if(worldMPIRank == 0) {
                std::cout << "iTurb = " << iTurb << " turbineMapGlobToProc[iTurb] = " << turbineMapGlobToProc[iTurb] << std::endl ;
            }
        }
        if(worldMPIRank == turbineMapGlobToProc[iTurb]) {
            turbineMapProcToGlob[nTurbinesProc] = iTurb;
            reverseTurbineMapProcToGlob[iTurb] = nTurbinesProc;
            nTurbinesProc++ ;
        }
        turbineSetProcs.insert(turbineMapGlobToProc[iTurb]);
    }

    int nProcsWithTurbines=0;
    turbineProcs.resize(turbineSetProcs.size());

    for (std::set<int>::const_iterator p = turbineSetProcs.begin(); p != turbineSetProcs.end(); p++) {
        turbineProcs[nProcsWithTurbines] = *p;
        nProcsWithTurbines++ ;
    }

    if (dryRun) {
        if (nTurbinesProc > 0) {
            std::ofstream turbineAllocFile;
            turbineAllocFile.open("turbineAlloc." + std::to_string(worldMPIRank) + ".txt") ;
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                turbineAllocFile << "Proc " << worldMPIRank << " loc iTurb " << iTurb << " glob iTurb " << turbineMapProcToGlob[iTurb] << std::endl ;
            }
            turbineAllocFile.flush();
            turbineAllocFile.close() ;
        }
    }

    // // Construct a group containing all procs running atleast 1 turbine in FAST
    // MPI_Group_incl(worldMPIGroup, nProcsWithTurbines, &turbineProcs[0], &fastMPIGroup) ;
    // int fastMPIcommTag = MPI_Comm_create(mpiComm, fastMPIGroup, &fastMPIComm);
    // if (MPI_COMM_NULL != fastMPIComm) {
    //     MPI_Comm_rank(fastMPIComm, &fastMPIRank);
    // }

    turbineData.resize(nTurbinesProc);
    velForceNodeData.resize(nTurbinesProc);
    brFSIData.resize(nTurbinesProc);

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        turbineData[iTurb].TurbineBasePos.resize(3);
        turbineData[iTurb].TurbineHubPos.resize(3);

        int iTurbGlob = turbineMapProcToGlob[iTurb];
        turbineData[iTurb].TurbID = globTurbineData[iTurbGlob].TurbID;
        turbineData[iTurb].sType = globTurbineData[iTurbGlob].sType;
        turbineData[iTurb].FASTInputFileName = globTurbineData[iTurbGlob].FASTInputFileName ;
        turbineData[iTurb].FASTRestartFileName = globTurbineData[iTurbGlob].FASTRestartFileName ;
        for(int i=0;i<3;i++) {
            turbineData[iTurb].TurbineBasePos[i] = globTurbineData[iTurbGlob].TurbineBasePos[i];
            turbineData[iTurb].TurbineHubPos[i] = globTurbineData[iTurbGlob].TurbineHubPos[i];
        }
        turbineData[iTurb].numForcePtsBlade = globTurbineData[iTurbGlob].numForcePtsBlade;
        turbineData[iTurb].numForcePtsTwr = globTurbineData[iTurbGlob].numForcePtsTwr;
        turbineData[iTurb].azBlendMean = globTurbineData[iTurbGlob].azBlendMean;
        turbineData[iTurb].azBlendDelta = globTurbineData[iTurbGlob].azBlendDelta;

        velForceNodeData[iTurb].resize(4); // To hold data for 4 time steps
        brFSIData[iTurb].resize(4);

    }

    // Allocate memory for Turbine datastructure for all turbines
    FAST_AllocateTurbines(&nTurbinesProc, &ErrStat, ErrMsg);

    // Allocate memory for ExtInfw Input types in FAST
    extinfw_i_f_FAST.resize(nTurbinesProc) ;
    extinfw_o_t_FAST.resize(nTurbinesProc) ;

    // Allocate memory for ExtLd Input/Parameter/Output types in FAST
    extld_i_f_FAST.resize(nTurbinesProc) ;
    extld_p_f_FAST.resize(nTurbinesProc) ;
    extld_o_t_FAST.resize(nTurbinesProc) ;

}

void fast::OpenFAST::allocateMemory_postInit(int iTurbLoc) {

    if (turbineData[iTurbLoc].sType == EXTINFLOW) {
        turbineData[iTurbLoc].nBRfsiPtsBlade = std::vector<int>(turbineData[iTurbLoc].numBlades,0);
        turbineData[iTurbLoc].nBRfsiPtsTwr = 0;

        if ( turbineData[iTurbLoc].inflowType == 1) {
            // Inflow data is coming from inflow module
            turbineData[iTurbLoc].numForcePtsTwr = 0;
            turbineData[iTurbLoc].numForcePtsBlade = 0;
            turbineData[iTurbLoc].numForcePts = 0;
            turbineData[iTurbLoc].numVelPtsTwr = 0;
            turbineData[iTurbLoc].numVelPtsBlade = 0;
            turbineData[iTurbLoc].numVelPts = 0;
        } else {
            //Inflow data is coming from external program like a CFD solver
            turbineData[iTurbLoc].numForcePts = 1 + turbineData[iTurbLoc].numForcePtsTwr + turbineData[iTurbLoc].numBlades * turbineData[iTurbLoc].numForcePtsBlade ;
            turbineData[iTurbLoc].numVelPts = 1 + turbineData[iTurbLoc].numVelPtsTwr + turbineData[iTurbLoc].numBlades * turbineData[iTurbLoc].numVelPtsBlade ;

            int nfpts = get_numForcePtsLoc(iTurbLoc);
            int nvelpts = get_numVelPtsLoc(iTurbLoc);

            velForceNodeData[iTurbLoc][3].xref_force.resize(3*nfpts);
            for(int k=0; k<4; k++) {
                velForceNodeData[iTurbLoc][k].x_vel.resize(3*nvelpts) ;
                velForceNodeData[iTurbLoc][k].vel_vel.resize(3*nvelpts) ;
                velForceNodeData[iTurbLoc][k].x_force.resize(3*nfpts) ;
                velForceNodeData[iTurbLoc][k].xdot_force.resize(3*nfpts) ;
                velForceNodeData[iTurbLoc][k].orient_force.resize(9*nfpts) ;
                velForceNodeData[iTurbLoc][k].vel_force.resize(3*nfpts) ;
                velForceNodeData[iTurbLoc][k].force.resize(3*nfpts) ;
            }

            if ( isDebug() ) {
                for (int iNode=0; iNode < get_numVelPtsLoc(iTurbLoc); iNode++) {
                    std::cout << "Node " << iNode << " Position = " << extinfw_i_f_FAST[iTurbLoc].pxVel[iNode] << " " << extinfw_i_f_FAST[iTurbLoc].pyVel[iNode] << " " << extinfw_i_f_FAST[iTurbLoc].pzVel[iNode] << " " << std::endl ;
                }
            }
        }
        std::cerr << "turbineData[iTurbLoc].inflowType " << turbineData[iTurbLoc].inflowType << std::endl;
        std::cerr << "turbineData[iTurbLoc].numForcePtsTwr = " << turbineData[iTurbLoc].numForcePtsTwr << std::endl;
        std::cerr << "turbineData[iTurbLoc].numForcePtsBlade = " << turbineData[iTurbLoc].numForcePtsBlade << std::endl;
        std::cerr << "turbineData[iTurbLoc].numForcePts = " << turbineData[iTurbLoc].numForcePts << std::endl;


    } else if (turbineData[iTurbLoc].sType == EXTLOADS) {
        turbineData[iTurbLoc].nBRfsiPtsBlade.resize(turbineData[iTurbLoc].numBlades);
        int nTotBldNds = 0;
        for(int i=0; i < turbineData[iTurbLoc].numBlades; i++) {
            nTotBldNds += extld_p_f_FAST[iTurbLoc].nBladeNodes[i];
            turbineData[iTurbLoc].nBRfsiPtsBlade[i] = extld_p_f_FAST[iTurbLoc].nBladeNodes[i];
            turbineData[iTurbLoc].nTotBRfsiPtsBlade += turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        }
        turbineData[iTurbLoc].nBRfsiPtsTwr = extld_p_f_FAST[iTurbLoc].nTowerNodes[0];

        // Allocate memory for reference position only for 1 time step - np1
        for(int k=0; k<4; k++) {
            brFSIData[iTurbLoc][k].twr_ref_pos.resize(6*turbineData[iTurbLoc].nBRfsiPtsTwr);
            brFSIData[iTurbLoc][k].twr_def.resize(6*turbineData[iTurbLoc].nBRfsiPtsTwr);
            brFSIData[iTurbLoc][k].twr_vel.resize(6*turbineData[iTurbLoc].nBRfsiPtsTwr);
            brFSIData[iTurbLoc][k].bld_rloc.resize(nTotBldNds);
            brFSIData[iTurbLoc][k].bld_chord.resize(nTotBldNds);
            brFSIData[iTurbLoc][k].bld_ref_pos.resize(6*nTotBldNds);
            brFSIData[iTurbLoc][k].bld_def.resize(6*nTotBldNds);
            brFSIData[iTurbLoc][k].bld_vel.resize(6*nTotBldNds);
            brFSIData[iTurbLoc][k].twr_ld.resize(6*turbineData[iTurbLoc].nBRfsiPtsTwr);
            brFSIData[iTurbLoc][k].bld_ld.resize(6*nTotBldNds);
            brFSIData[iTurbLoc][k].hub_ref_pos.resize(6);
            brFSIData[iTurbLoc][k].hub_def.resize(6);
            brFSIData[iTurbLoc][k].hub_vel.resize(6);
            brFSIData[iTurbLoc][k].nac_ref_pos.resize(6);
            brFSIData[iTurbLoc][k].nac_def.resize(6);
            brFSIData[iTurbLoc][k].nac_vel.resize(6);
            brFSIData[iTurbLoc][k].hub_ref_pos.resize(6);
            brFSIData[iTurbLoc][k].bld_pitch.resize(turbineData[iTurbLoc].numBlades);
            brFSIData[iTurbLoc][k].bld_root_ref_pos.resize(6*turbineData[iTurbLoc].numBlades);
            brFSIData[iTurbLoc][k].bld_root_def.resize(6*turbineData[iTurbLoc].numBlades);
        }
    }

}

void fast::OpenFAST::allocateTurbinesToProcsSimple() {
    // Allocate turbines to each processor - round robin fashion
    int nProcs ;
    MPI_Comm_size(mpiComm, &nProcs);
    for(int j = 0; j < nTurbinesGlob; j++)  turbineMapGlobToProc[j] = j % nProcs ;
}

void fast::OpenFAST::end() {

    // Deallocate types we allocated earlier

    if ( !dryRun) {
        bool stopTheProgram = false;
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_End(&iTurb, &stopTheProgram);
        }
    }

    // MPI_Group_free(&fastMPIGroup);
    // if (MPI_COMM_NULL != fastMPIComm) {
    //     MPI_Comm_free(&fastMPIComm);
    // }
    // MPI_Group_free(&worldMPIGroup);

}

int fast::OpenFAST::read_nlin_iters(int iTurb, int n_t_global, int ncid) {

    int nlin_iters = 0;
    size_t count1 = 1;
    size_t n_tsteps = n_t_global;
    int ierr = nc_get_vara_int(ncid, 1, &n_tsteps, &count1, &nlin_iters);

    return nlin_iters;

}


void fast::OpenFAST::readVelocityData(int iTurb, int n_t_global, int nlinIter, int ncid) {

    size_t n_tsteps = n_t_global;
    const std::vector<size_t> start_dim{n_tsteps, static_cast<size_t>(nlinIter), 0};
    int nVelPts = get_numVelPtsLoc(iTurb);
    const std::vector<size_t> velPtsDataDims{1, 1, static_cast<size_t>(3*nVelPts)};
    int ierr = nc_get_vara_double(ncid, 2, start_dim.data(), velPtsDataDims.data(), velForceNodeData[iTurb][fast::STATE_NP1].vel_vel.data());
}

int fast::OpenFAST::openVelocityDataFile(int iTurb) {

    int ncid;
    std::stringstream velfile_fstream;
    velfile_fstream << "turb_" ;
    velfile_fstream << std::setfill('0') << std::setw(2) << turbineData[iTurb].TurbID;
    velfile_fstream << "_veldata.nc";
    std::string velfile_filename = velfile_fstream.str();
    int ierr = nc_open(velfile_filename.c_str(), NC_WRITE, &ncid);
    check_nc_error(ierr, "nc_open");
    return ncid;

}

void fast::OpenFAST::prepareVelocityDataFile(int iTurb) {

    // Open the file in create mode - this will destory any file
    int ncid;
    std::stringstream velfile_fstream;
    velfile_fstream << "turb_" ;
    velfile_fstream << std::setfill('0') << std::setw(2) << turbineData[iTurb].TurbID;
    velfile_fstream << "_veldata.nc";
    std::string velfile_filename = velfile_fstream.str();
    int ierr = nc_create(velfile_filename.c_str(), NC_CLOBBER, &ncid);
    check_nc_error(ierr, "nc_create");

    //Define dimensions
    int tmpDimID;
    ierr = nc_def_dim(ncid, "n_tsteps", NC_UNLIMITED, &tmpDimID);
    ierr = nc_def_dim(ncid, "n_nonlin_iters_max", 2, &tmpDimID);
    ierr = nc_def_dim(ncid, "n_vel_pts_data", turbineData[iTurb].numVelPts*3, &tmpDimID);

    int tmpVarID;
    tmpDimID = 0;
    ierr = nc_def_var(ncid, "time", NC_DOUBLE, 1, &tmpDimID, &tmpVarID);
    ierr = nc_def_var(ncid, "nlin_iters", NC_INT, 1, &tmpDimID, &tmpVarID);
    const std::vector<int> velPtsDataDims{0, 1, 2};
    ierr = nc_def_var(ncid, "vel_vel", NC_DOUBLE, 3, velPtsDataDims.data(), &tmpVarID);

    //! Indicate that we are done defining variables, ready to write data
    ierr = nc_enddef(ncid);
    check_nc_error(ierr, "nc_enddef");
    ierr = nc_close(ncid);
    check_nc_error(ierr, "nc_close");
}

void fast::OpenFAST::writeVelocityData(int iTurb, int n_t_global, int nlinIter) {

    /* // NetCDF stuff to write velocity data to file */
    int ncid;
    //Find the file and open it in append mode
    std::stringstream velfile_ss;
    velfile_ss << "turb_" ;
    velfile_ss << std::setfill('0') << std::setw(2) << turbineData[iTurb].TurbID;
    velfile_ss << "_veldata.nc";
    std::string vel_filename = velfile_ss.str();
    int ierr = nc_open(vel_filename.c_str(), NC_WRITE, &ncid);
    check_nc_error(ierr, "nc_open");

    size_t count1=1;
    size_t n_tsteps = (n_t_global/nSubsteps_)+1;
    double curTime = (n_t_global + nSubsteps_) * dtFAST;
    ierr = nc_put_vara_double(ncid, 0, &n_tsteps, &count1, &curTime);
    int nVelPts = get_numVelPtsLoc(iTurb) ;
    const std::vector<size_t> velPtsDataDims{1, 1, static_cast<size_t>(3*nVelPts)};
    const std::vector<size_t> start_dim{static_cast<size_t>(n_tsteps),static_cast<size_t>(nlinIter),0};

    std::cout << "Writing velocity data at time step " << n_tsteps << ", nonlinear iteration " << nlinIter << std::endl ;
    ierr = nc_put_vara_double(ncid, 2, start_dim.data(), velPtsDataDims.data(), velForceNodeData[iTurb][3].vel_vel.data());
    nlinIter += 1; // To account for 0-based indexing
    ierr = nc_put_vara_int(ncid, 1, &n_tsteps, &count1, &nlinIter);

    nc_close(ncid);

}

void fast::OpenFAST::send_data_to_openfast(fast::timeStep t) {

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        if ( (turbineData[iTurb].sType == EXTINFLOW) && (turbineData[iTurb].inflowType == 2) ) {
            int nvelpts = get_numVelPtsLoc(iTurb);
            for (int iNodeVel=0; iNodeVel < nvelpts; iNodeVel++) {
                extinfw_o_t_FAST[iTurb].u[iNodeVel] = velForceNodeData[iTurb][t].vel_vel[iNodeVel*3+0];
                extinfw_o_t_FAST[iTurb].v[iNodeVel] = velForceNodeData[iTurb][t].vel_vel[iNodeVel*3+1];
                extinfw_o_t_FAST[iTurb].w[iNodeVel] = velForceNodeData[iTurb][t].vel_vel[iNodeVel*3+2];
            }
        } else if(turbineData[iTurb].sType == EXTLOADS) {

            int nBlades = turbineData[iTurb].numBlades;
            int iRunTot = 0;
            for(int i=0; i < nBlades; i++) {
                int nPtsBlade = turbineData[iTurb].nBRfsiPtsBlade[i];
                for (int j=0; j < nPtsBlade; j++) {
                    for (int k=0; k<6; k++) {
                        extld_o_t_FAST[iTurb].bldLd[iRunTot*6+k] = brFSIData[iTurb][t].bld_ld[iRunTot*6+k];
                    }
                    iRunTot++;
                }
            }

            int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
            for (int i=0; i < nPtsTwr*6; i++)
                extld_o_t_FAST[iTurb].twrLd[i] = brFSIData[iTurb][t].twr_ld[i];

        }

    }
}

void fast::OpenFAST::send_data_to_openfast(double ss_time) {

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        if (turbineData[iTurb].inflowType == 2) {
            int nvelpts = get_numVelPtsLoc(iTurb);
            for (int iNodeVel=0; iNodeVel < nvelpts; iNodeVel++) {
                extinfw_o_t_FAST[iTurb].u[iNodeVel] = velForceNodeData[iTurb][fast::STATE_N].vel_vel[iNodeVel*3+0] + ss_time * (velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+0] - velForceNodeData[iTurb][fast::STATE_N].vel_vel[iNodeVel*3+0]);
                extinfw_o_t_FAST[iTurb].v[iNodeVel] = velForceNodeData[iTurb][fast::STATE_N].vel_vel[iNodeVel*3+1] + ss_time * (velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+1] - velForceNodeData[iTurb][fast::STATE_N].vel_vel[iNodeVel*3+1]);
                extinfw_o_t_FAST[iTurb].w[iNodeVel] = velForceNodeData[iTurb][fast::STATE_N].vel_vel[iNodeVel*3+2] + ss_time * (velForceNodeData[iTurb][fast::STATE_NP1].vel_vel[iNodeVel*3+2] - velForceNodeData[iTurb][fast::STATE_N].vel_vel[iNodeVel*3+2]);
            }
        } else if(turbineData[iTurb].sType == EXTLOADS) {

            int nBlades = turbineData[iTurb].numBlades;
            int iRunTot = 0;
            for(int i=0; i < nBlades; i++) {
                int nPtsBlade = turbineData[iTurb].nBRfsiPtsBlade[i];
                for (int j=0; j < nPtsBlade; j++) {
                    for (int k=0; k<6; k++) {
                        extld_o_t_FAST[iTurb].bldLd[iRunTot*6+k] = brFSIData[iTurb][fast::STATE_N].bld_ld[iRunTot*6+k] + ss_time * (brFSIData[iTurb][fast::STATE_NP1].bld_ld[iRunTot*6+k] - brFSIData[iTurb][fast::STATE_N].bld_ld[iRunTot*6+k]);
                    }
                    iRunTot++;
                }
            }

            int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
            for (int i=0; i < nPtsTwr*6; i++)
                extld_o_t_FAST[iTurb].twrLd[i] = brFSIData[iTurb][fast::STATE_N].twr_ld[i] + ss_time * (brFSIData[iTurb][fast::STATE_NP1].twr_ld[i] - brFSIData[iTurb][fast::STATE_N].twr_ld[i]);

        }
    }

}

void fast::OpenFAST::get_data_from_openfast(timeStep t) {


    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        if(turbineData[iTurb].sType == EXTINFLOW) {

            if (turbineData[iTurb].inflowType == 2) {
                int nvelpts = get_numVelPtsLoc(iTurb);
                int nfpts = get_numForcePtsLoc(iTurb);
                // std::cerr << "nvelpts = " << nvelpts << std::endl;
                // std::cerr << "nfpts = " << nfpts << "  " << get_numForcePtsBladeLoc(iTurb) << " " << get_numForcePtsTwrLoc(iTurb) << std::endl;
                for (int i=0; i<nvelpts; i++) {
                    velForceNodeData[iTurb][t].x_vel_resid += (velForceNodeData[iTurb][t].x_vel[i*3+0] - extinfw_i_f_FAST[iTurb].pxVel[i])*(velForceNodeData[iTurb][t].x_vel[i*3+0] - extinfw_i_f_FAST[iTurb].pxVel[i]);
                    velForceNodeData[iTurb][t].x_vel[i*3+0] = extinfw_i_f_FAST[iTurb].pxVel[i];
                    velForceNodeData[iTurb][t].x_vel_resid += (velForceNodeData[iTurb][t].x_vel[i*3+1] - extinfw_i_f_FAST[iTurb].pyVel[i])*(velForceNodeData[iTurb][t].x_vel[i*3+1] - extinfw_i_f_FAST[iTurb].pyVel[i]);
                    velForceNodeData[iTurb][t].x_vel[i*3+1] = extinfw_i_f_FAST[iTurb].pyVel[i];
                    velForceNodeData[iTurb][t].x_vel_resid += (velForceNodeData[iTurb][t].x_vel[i*3+2] - extinfw_i_f_FAST[iTurb].pzVel[i])*(velForceNodeData[iTurb][t].x_vel[i*3+2] - extinfw_i_f_FAST[iTurb].pzVel[i]);
                    velForceNodeData[iTurb][t].x_vel[i*3+2] = extinfw_i_f_FAST[iTurb].pzVel[i];
                }

                for (int i=0; i<nfpts; i++) {
                    velForceNodeData[iTurb][t].x_force_resid += (velForceNodeData[iTurb][t].x_force[i*3+0] - extinfw_i_f_FAST[iTurb].pxForce[i])*(velForceNodeData[iTurb][t].x_force[i*3+0] - extinfw_i_f_FAST[iTurb].pxForce[i]);
                    velForceNodeData[iTurb][t].x_force[i*3+0] = extinfw_i_f_FAST[iTurb].pxForce[i];
                    velForceNodeData[iTurb][t].x_force_resid += (velForceNodeData[iTurb][t].x_force[i*3+1] - extinfw_i_f_FAST[iTurb].pyForce[i])*(velForceNodeData[iTurb][t].x_force[i*3+1] - extinfw_i_f_FAST[iTurb].pyForce[i]);
                    velForceNodeData[iTurb][t].x_force[i*3+1] = extinfw_i_f_FAST[iTurb].pyForce[i];
                    velForceNodeData[iTurb][t].x_force_resid += (velForceNodeData[iTurb][t].x_force[i*3+2] - extinfw_i_f_FAST[iTurb].pzForce[i])*(velForceNodeData[iTurb][t].x_force[i*3+2] - extinfw_i_f_FAST[iTurb].pzForce[i]);
                    velForceNodeData[iTurb][t].x_force[i*3+2] = extinfw_i_f_FAST[iTurb].pzForce[i];
                    velForceNodeData[iTurb][t].xdot_force_resid += (velForceNodeData[iTurb][t].xdot_force[i*3+0] - extinfw_i_f_FAST[iTurb].xdotForce[i])*(velForceNodeData[iTurb][t].xdot_force[i*3+0] - extinfw_i_f_FAST[iTurb].xdotForce[i]);
                    velForceNodeData[iTurb][t].xdot_force[i*3+0] = extinfw_i_f_FAST[iTurb].xdotForce[i];
                    velForceNodeData[iTurb][t].xdot_force_resid += (velForceNodeData[iTurb][t].xdot_force[i*3+1] - extinfw_i_f_FAST[iTurb].ydotForce[i])*(velForceNodeData[iTurb][t].xdot_force[i*3+1] - extinfw_i_f_FAST[iTurb].ydotForce[i]);
                    velForceNodeData[iTurb][t].xdot_force[i*3+1] = extinfw_i_f_FAST[iTurb].ydotForce[i];
                    velForceNodeData[iTurb][t].xdot_force_resid += (velForceNodeData[iTurb][t].xdot_force[i*3+2] - extinfw_i_f_FAST[iTurb].zdotForce[i])*(velForceNodeData[iTurb][t].xdot_force[i*3+2] - extinfw_i_f_FAST[iTurb].zdotForce[i]);
                    velForceNodeData[iTurb][t].xdot_force[i*3+2] = extinfw_i_f_FAST[iTurb].zdotForce[i];
                    for (int j=0;j<9;j++) {
                        velForceNodeData[iTurb][t].orient_force_resid += (velForceNodeData[iTurb][t].orient_force[i*9+j] - extinfw_i_f_FAST[iTurb].pOrientation[i*9+j])*(velForceNodeData[iTurb][t].orient_force[i*9+j] - extinfw_i_f_FAST[iTurb].pOrientation[i*9+j]);
                        velForceNodeData[iTurb][t].orient_force[i*9+j] = extinfw_i_f_FAST[iTurb].pOrientation[i*9+j];
                    }
                    velForceNodeData[iTurb][t].force_resid += (velForceNodeData[iTurb][t].force[i*3+0] - extinfw_i_f_FAST[iTurb].fx[i])*(velForceNodeData[iTurb][t].force[i*3+0] - extinfw_i_f_FAST[iTurb].fx[i]);
                    velForceNodeData[iTurb][t].force[i*3+0] = extinfw_i_f_FAST[iTurb].fx[i];
                    velForceNodeData[iTurb][t].force_resid += (velForceNodeData[iTurb][t].force[i*3+1] - extinfw_i_f_FAST[iTurb].fy[i])*(velForceNodeData[iTurb][t].force[i*3+1] - extinfw_i_f_FAST[iTurb].fy[i]);
                    velForceNodeData[iTurb][t].force[i*3+1] = extinfw_i_f_FAST[iTurb].fy[i];
                    velForceNodeData[iTurb][t].force_resid += (velForceNodeData[iTurb][t].force[i*3+2] - extinfw_i_f_FAST[iTurb].fz[i])*(velForceNodeData[iTurb][t].force[i*3+2] - extinfw_i_f_FAST[iTurb].fz[i]);
                    velForceNodeData[iTurb][t].force[i*3+2] = extinfw_i_f_FAST[iTurb].fz[i];
                }
            }
        } else if(turbineData[iTurb].sType == EXTLOADS) {

            int nBlades = turbineData[iTurb].numBlades;
            int iRunTot = 0;
            for (int i=0; i < nBlades; i++) {
                int nPtsBlade = turbineData[iTurb].nBRfsiPtsBlade[i];
                for (int j=0; j < nPtsBlade; j++) {
                    for (int k=0; k < 3; k++) {
                        brFSIData[iTurb][t].bld_def[iRunTot*6+k] = extld_i_f_FAST[iTurb].bldDef[iRunTot*12+k];
                        brFSIData[iTurb][t].bld_vel[iRunTot*6+k] = extld_i_f_FAST[iTurb].bldDef[iRunTot*12+3+k];
                        brFSIData[iTurb][t].bld_def[iRunTot*6+3+k] = extld_i_f_FAST[iTurb].bldDef[iRunTot*12+6+k];
                        brFSIData[iTurb][t].bld_vel[iRunTot*6+3+k] = extld_i_f_FAST[iTurb].bldDef[iRunTot*12+9+k];
                    }
                    iRunTot++;
                }
                for (int k=0; k < 3; k++) {
                    brFSIData[iTurb][t].bld_root_def[i*6+k] = extld_i_f_FAST[iTurb].bldRootDef[i*12+k];
                    brFSIData[iTurb][t].bld_root_def[i*6+3+k] = extld_i_f_FAST[iTurb].bldRootDef[i*12+6+k];
                }
                brFSIData[iTurb][t].bld_pitch[i] = extld_i_f_FAST[iTurb].bldPitch[i];
            }

            int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
            for (int i=0; i < nPtsTwr; i++) {
                for (int j = 0; j < 3; j++) {
                    brFSIData[iTurb][t].twr_def[i*6+j] = extld_i_f_FAST[iTurb].twrDef[i*12+j];
                    brFSIData[iTurb][t].twr_vel[i*6+j] = extld_i_f_FAST[iTurb].twrDef[i*12+3+j];
                    brFSIData[iTurb][t].twr_def[i*6+3+j] = extld_i_f_FAST[iTurb].twrDef[i*12+6+j];
                    brFSIData[iTurb][t].twr_vel[i*6+3+j] = extld_i_f_FAST[iTurb].twrDef[i*12+9+j];
                }
            }

            for (int j = 0; j < 3; j++) {
                brFSIData[iTurb][t].hub_def[j] = extld_i_f_FAST[iTurb].hubDef[j];
                brFSIData[iTurb][t].hub_vel[j] = extld_i_f_FAST[iTurb].hubDef[3+j];
                brFSIData[iTurb][t].hub_def[3+j] = extld_i_f_FAST[iTurb].hubDef[6+j];
                brFSIData[iTurb][t].hub_vel[3+j] = extld_i_f_FAST[iTurb].hubDef[9+j];
                brFSIData[iTurb][t].nac_def[j] = extld_i_f_FAST[iTurb].nacDef[j];
                brFSIData[iTurb][t].nac_vel[j] = extld_i_f_FAST[iTurb].nacDef[3+j];
                brFSIData[iTurb][t].nac_def[3+j] = extld_i_f_FAST[iTurb].nacDef[6+j];
                brFSIData[iTurb][t].nac_vel[3+j] = extld_i_f_FAST[iTurb].nacDef[9+j];
            }
            //TODO: May be calculate the residual here as well
        }
    }
}

void fast::OpenFAST::readRestartFile(int iTurbLoc, int n_t_global) {

    int ncid;
    //Find the file and open it in append mode
    std::stringstream rstfile_ss;
    rstfile_ss << "turb_" ;
    rstfile_ss << std::setfill('0') << std::setw(2) << turbineData[iTurbLoc].TurbID;
    rstfile_ss << "_rst.nc";
    std::string rst_filename = rstfile_ss.str();
    int ierr = nc_open(rst_filename.c_str(), NC_NOWRITE, &ncid);
    check_nc_error(ierr, "nc_open");

    size_t count1 = 1;
    size_t n_tsteps;
    ierr = nc_inq_dimlen(ncid, ncRstDimIDs_["n_tsteps"], &n_tsteps);
    check_nc_error(ierr, "nc_inq_dimlen");
    n_tsteps -= 1; // To account for 0 based indexing

    if (turbineData[iTurbLoc].sType == EXTINFLOW) {

        int nvelpts = get_numVelPtsLoc(iTurbLoc);
        int nfpts = get_numForcePtsLoc(iTurbLoc);

        const std::vector<size_t> velPtsDataDims{1, 1, static_cast<size_t>(3*nvelpts)};
        const std::vector<size_t> forcePtsDataDims{1, 1, static_cast<size_t>(3*nfpts)};
        const std::vector<size_t> forcePtsOrientDataDims{1, 1, static_cast<size_t>(9*nfpts)};

        ierr = nc_get_var_double(ncid, ncRstVarIDs_["xref_force"], velForceNodeData[iTurbLoc][fast::STATE_NP1].xref_force.data());

        for (size_t j=0; j < 4; j++) {  // Loop over states - NM2, STATE_NM1, N, NP1

            const std::vector<size_t> start_dim{n_tsteps,j,0};

            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["x_vel"], start_dim.data(), velPtsDataDims.data(), velForceNodeData[iTurbLoc][j].x_vel.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["vel_vel"], start_dim.data(), velPtsDataDims.data(), velForceNodeData[iTurbLoc][j].vel_vel.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["x_force"], start_dim.data(), forcePtsDataDims.data(), velForceNodeData[iTurbLoc][j].x_force.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["xdot_force"], start_dim.data(), forcePtsDataDims.data(), velForceNodeData[iTurbLoc][j].xdot_force.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["vel_force"], start_dim.data(), forcePtsDataDims.data(), velForceNodeData[iTurbLoc][j].vel_force.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["force"], start_dim.data(), forcePtsDataDims.data(), velForceNodeData[iTurbLoc][j].force.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["orient_force"], start_dim.data(), forcePtsOrientDataDims.data(), velForceNodeData[iTurbLoc][j].orient_force.data());

        }

    } else if (turbineData[iTurbLoc].sType == EXTLOADS) {

        int nBRfsiPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
        int nTotBRfsiPtsBlade = turbineData[iTurbLoc].nTotBRfsiPtsBlade;
        int nBlades = turbineData[iTurbLoc].numBlades;
        const std::vector<size_t> twrDataDims{1, 1, static_cast<size_t>(6*nBRfsiPtsTwr)};
        const std::vector<size_t> bldDataDims{1, 1, static_cast<size_t>(6*nTotBRfsiPtsBlade)};
        const std::vector<size_t> bldRootDataDims{1, 1, static_cast<size_t>(6*nBlades)};
        const std::vector<size_t> bldPitchDataDims{1, 1, static_cast<size_t>(nBlades)};
        const std::vector<size_t> ptDataDims{1, 1, 6};

        for (size_t j=0; j < 4; j++) {  // Loop over states - NM2, STATE_NM1, N, NP1

            const std::vector<size_t> start_dim{n_tsteps, j, 0};

            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["twr_def"], start_dim.data(), twrDataDims.data(), brFSIData[iTurbLoc][j].twr_def.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["twr_vel"], start_dim.data(), twrDataDims.data(), brFSIData[iTurbLoc][j].twr_vel.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["twr_ld"], start_dim.data(), twrDataDims.data(), brFSIData[iTurbLoc][j].twr_ld.data());

            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["bld_def"], start_dim.data(), bldDataDims.data(), brFSIData[iTurbLoc][j].bld_def.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["bld_vel"], start_dim.data(), bldDataDims.data(), brFSIData[iTurbLoc][j].bld_vel.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["bld_ld"], start_dim.data(), bldDataDims.data(), brFSIData[iTurbLoc][j].bld_ld.data());

            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["hub_def"], start_dim.data(), ptDataDims.data(), brFSIData[iTurbLoc][j].hub_def.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["hub_vel"], start_dim.data(), ptDataDims.data(), brFSIData[iTurbLoc][j].hub_vel.data());

            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["nac_def"], start_dim.data(), ptDataDims.data(), brFSIData[iTurbLoc][j].nac_def.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["nac_vel"], start_dim.data(), ptDataDims.data(), brFSIData[iTurbLoc][j].nac_vel.data());

            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["bld_root_def"], start_dim.data(), bldRootDataDims.data(), brFSIData[iTurbLoc][j].bld_root_def.data());
            ierr = nc_get_vara_double(ncid, ncRstVarIDs_["bld_pitch"], start_dim.data(), bldPitchDataDims.data(), brFSIData[iTurbLoc][j].bld_pitch.data());

        }



    }

    nc_close(ncid);

}

void fast::OpenFAST::cross(double * a, double * b, double * aCrossb) {

    aCrossb[0] = a[1]*b[2] - a[2]*b[1];
    aCrossb[1] = a[2]*b[0] - a[0]*b[2];
    aCrossb[2] = a[0]*b[1] - a[1]*b[0];

}


//! Apply a DCM rotation 'dcm' to a vector 'r' into 'rRot'. To optionally transpose the rotation, set 'tranpose=-1.0'.
void fast::OpenFAST::applyDCMrotation(double * dcm, double * r, double *rRot, double transpose) {

    if (transpose > 0) {
        for(size_t i=0; i < 3; i++) {
            rRot[i] = 0.0;
            for(size_t j=0; j < 3; j++)
                rRot[i] += dcm[i*3+j] * r[j];
        }
    } else {
        for(size_t i=0; i < 3; i++) {
            rRot[i] = 0.0;
            for(size_t j=0; j < 3; j++)
                rRot[i] += dcm[j*3+i] * r[j];
        }
    }
}

//! Apply a Wiener-Milenkovic rotation 'wm' to a vector 'r' into 'rRot'. To optionally transpose the rotation, set 'tranpose=-1.0'.
void fast::OpenFAST::applyWMrotation(double * wm, double * r, double *rRot, double transpose) {

    double wm0 = 2.0-0.125*dot(wm, wm);
    double nu = 2.0/(4.0-wm0);
    double cosPhiO2 = 0.5*wm0*nu;
    std::vector<double> wmCrossR(3,0.0);
    cross(wm, r, wmCrossR.data());
    std::vector<double> wmCrosswmCrossR(3,0.0);
    cross(wm, wmCrossR.data(), wmCrosswmCrossR.data());

    for(size_t i=0; i < 3; i++)
        rRot[i] = r[i] + transpose * nu * cosPhiO2 * wmCrossR[i] + 0.5 * nu * nu * wmCrosswmCrossR[i];

}


void fast::OpenFAST::writeOutputFile(int iTurbLoc, int n_t_global) {

    int ncid;
    //Open the file in append mode
    std::stringstream outfile_ss;
    outfile_ss << "turb_" ;
    outfile_ss << std::setfill('0') << std::setw(2) << turbineData[iTurbLoc].TurbID;
    outfile_ss << "_output.nc";
    std::string defloads_filename = outfile_ss.str();
    int ierr = nc_open(defloads_filename.c_str(), NC_WRITE, &ncid);
    check_nc_error(ierr, "nc_open");

    size_t count1=1;
    int tStepRatio = time_step_ratio(dtFAST, dtDriver);
    size_t n_tsteps = n_t_global/tStepRatio/outputFreq_ - 1;
    double curTime = n_t_global * dtFAST;
    ierr = nc_put_vara_double(ncid, ncOutVarIDs_["time"], &n_tsteps, &count1, &curTime);

    if ( (turbineData[iTurbLoc].sType == EXTINFLOW) && (turbineData[iTurbLoc].inflowType == 2) ) {

        // Nothing to do here yet
        int nBlades = get_numBladesLoc(iTurbLoc);
        int nBldPts = get_numForcePtsBladeLoc(iTurbLoc);
        int nTwrPts = get_numForcePtsTwrLoc(iTurbLoc);
        std::vector<double> tmpArray;

        tmpArray.resize(nTwrPts);
        {
            int node_twr_start = (1 + nBlades * nBldPts)*3;
            std::vector<size_t> count_dim{1,1,static_cast<size_t>(nTwrPts)};
            for (size_t iDim=0; iDim < 3; iDim++) {
                for (auto i=0; i < nTwrPts; i++)
                    tmpArray[i] = velForceNodeData[iTurbLoc][3].x_force[node_twr_start+i*3+iDim] - velForceNodeData[iTurbLoc][3].xref_force[node_twr_start+i*3+iDim] ;
                std::vector<size_t> start_dim{n_tsteps,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_disp"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
            for (size_t iDim=0; iDim < 3; iDim++) {
                for (auto i=0; i < nTwrPts; i++)
                    tmpArray[i] = velForceNodeData[iTurbLoc][3].xdot_force[node_twr_start+i*3+iDim] ;
                std::vector<size_t> start_dim{n_tsteps,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_vel"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
            for (size_t iDim=0;iDim < 3; iDim++) {
                for (auto i=0; i < nTwrPts; i++)
                    tmpArray[i] = velForceNodeData[iTurbLoc][3].force[node_twr_start+i*3+iDim] ;
                std::vector<size_t> start_dim{n_tsteps,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_ld"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
        }

        tmpArray.resize(nBldPts);
        {
            std::vector<size_t> count_dim{1,1,1,static_cast<size_t>(nBldPts)};
            for (size_t iDim=0;iDim < 3; iDim++) {
                for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                    int node_bld_start = (1 + iBlade * nBldPts);
                    for (auto i=0; i < nBldPts; i++)
                        tmpArray[i] = velForceNodeData[iTurbLoc][3].x_force[(node_bld_start+i)*3+iDim] - velForceNodeData[iTurbLoc][3].xref_force[(node_bld_start+i)*3+iDim] ;
                    std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                    ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_disp"], start_dim.data(), count_dim.data(), tmpArray.data());
                }
            }
            for (size_t iDim=0;iDim < 3; iDim++) {
                int iStart = 0 ;
                for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                    int node_bld_start = (1 + iBlade * nBldPts);
                    for (auto i=0; i < nBldPts; i++)
                        tmpArray[i] = velForceNodeData[iTurbLoc][3].xdot_force[(node_bld_start+i)*3+iDim] ;
                    std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                    ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_vel"], start_dim.data(), count_dim.data(), tmpArray.data());
                }
            }
            for (size_t iDim=0;iDim < 3; iDim++) {
                int iStart = 0 ;
                for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                    int node_bld_start = (1 + iBlade * nBldPts);
                    for (auto i=0; i < nBldPts; i++)
                        tmpArray[i] = velForceNodeData[iTurbLoc][3].force[(node_bld_start+i)*3+iDim] ;
                    std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                    ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_ld"], start_dim.data(), count_dim.data(), tmpArray.data());
                }
            }

            std::vector<double> ld_loc(3*nBlades*nBldPts,0.0);
            for (auto iBlade=0; iBlade < nBlades; iBlade++) {
                int node_bld_start = (1 + iBlade * nBldPts);
                for (auto i=0; i < nBldPts; i++) {
                    applyDCMrotation(&velForceNodeData[iTurbLoc][3].orient_force[(node_bld_start + i)*9], &velForceNodeData[iTurbLoc][3].force[(node_bld_start+i)*3], &ld_loc[(node_bld_start-1)*3]);
                }
            }
            for (size_t iDim=0;iDim < 3; iDim++) {
                int iStart = 0 ;
                for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                    int node_bld_start = (iBlade * nBldPts);
                    for (auto i=0; i < nBldPts; i++)
                        tmpArray[i] = ld_loc[(node_bld_start+i)*3+iDim];
                    std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                    ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_ld_loc"], start_dim.data(), count_dim.data(), tmpArray.data());
                }
            }
        }

        tmpArray.resize(3);
        for (auto i=0; i < 3; i++)
            tmpArray[i] = velForceNodeData[iTurbLoc][3].x_force[i] - velForceNodeData[iTurbLoc][3].xref_force[i];
        std::vector<size_t> start_dim{n_tsteps, 0};
        std::vector<size_t> count_dim{1,3};
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["hub_disp"], start_dim.data(), count_dim.data(), tmpArray.data());
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["hub_vel"], start_dim.data(), count_dim.data(), &velForceNodeData[iTurbLoc][3].xdot_force[0]);
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["hub_ld"], start_dim.data(), count_dim.data(), &velForceNodeData[iTurbLoc][3].force[0]);

    } else if (turbineData[iTurbLoc].sType == EXTLOADS) {

        int nBlades = turbineData[iTurbLoc].numBlades;
        int nTwrPts = turbineData[iTurbLoc].nBRfsiPtsTwr;
        int nTotBldPts = turbineData[iTurbLoc].nTotBRfsiPtsBlade;
        int nBldPts = nTotBldPts/nBlades;

        std::vector<double> tmpArray;
        tmpArray.resize(nTwrPts);
        std::vector<size_t> count_dim{1,1,static_cast<size_t>(nTwrPts)};
        for (size_t iDim=0;iDim < 3; iDim++) {
            for (auto i=0; i < nTwrPts; i++) {
                tmpArray[i] = brFSIData[iTurbLoc][3].twr_def[i*6+iDim] ;
                // std::cerr << "Twr displacement Node " << i << ", dimension " << iDim << " = "
                //           << brFSIData[iTurbLoc][3].twr_ref_pos[i*6+iDim] <<  " "
                //           << brFSIData[iTurbLoc][3].twr_def[i*6+iDim] << std::endl;
            }
            std::vector<size_t> start_dim{n_tsteps,iDim,0};
            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_disp"], start_dim.data(), count_dim.data(), tmpArray.data());
        }

        for (size_t iDim=0;iDim < 3; iDim++) {
            for (auto i=0; i < nTwrPts; i++)
                tmpArray[i] = brFSIData[iTurbLoc][3].twr_def[i*6+3+iDim] ;
            std::vector<size_t> start_dim{n_tsteps,iDim,0};
            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_orient"], start_dim.data(), count_dim.data(), tmpArray.data());
        }
        for (size_t iDim=0;iDim < 3; iDim++) {
            for (auto i=0; i < nTwrPts; i++)
                tmpArray[i] = brFSIData[iTurbLoc][3].twr_vel[i*6+iDim] ;
            std::vector<size_t> start_dim{n_tsteps,iDim,0};
            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_vel"], start_dim.data(), count_dim.data(), tmpArray.data());
        }
        for (size_t iDim=0;iDim < 3; iDim++) {
            for (auto i=0; i < nTwrPts; i++)
                tmpArray[i] = brFSIData[iTurbLoc][3].twr_vel[i*6+3+iDim] ;
            std::vector<size_t> start_dim{n_tsteps,iDim,0};
            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_rotvel"], start_dim.data(), count_dim.data(), tmpArray.data());
        }

        for (size_t iDim=0;iDim < 3; iDim++) {
            for (auto i=0; i < nTwrPts; i++)
                tmpArray[i] = brFSIData[iTurbLoc][3].twr_ld[i*6+iDim] ;
            std::vector<size_t> start_dim{n_tsteps,iDim,0};
            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_ld"], start_dim.data(), count_dim.data(), tmpArray.data());
        }
        for (size_t iDim=0;iDim < 3; iDim++) {
            for (auto i=0; i < nTwrPts; i++)
                tmpArray[i] = brFSIData[iTurbLoc][3].twr_ld[i*6+3+iDim];
            std::vector<size_t> start_dim{n_tsteps,iDim,0};
            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["twr_moment"], start_dim.data(), count_dim.data(), tmpArray.data());
        }

    tmpArray.resize(nBldPts);
    {
        std::vector<size_t> count_dim{1,1,1,static_cast<size_t>(nBldPts)};
        for (size_t iDim=0;iDim < 3; iDim++) {
            int iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (auto i=0; i < nBldPts; i++) {
                    tmpArray[i] = brFSIData[iTurbLoc][3].bld_def[(iStart*6)+iDim];
                    iStart++;
                }
                std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_disp"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
        }
        for (size_t iDim=0;iDim < 3; iDim++) {
            int iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (auto i=0; i < nBldPts; i++) {
                    tmpArray[i] = brFSIData[iTurbLoc][3].bld_def[(iStart*6)+3+iDim];
                    iStart++;
                }
                std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_orient"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
        }
        for (size_t iDim=0;iDim < 3; iDim++) {
            int iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (auto i=0; i < nBldPts; i++) {
                    tmpArray[i] = brFSIData[iTurbLoc][3].bld_vel[(iStart*6)+iDim];
                    iStart++;
                }
                std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_vel"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
        }
        for (size_t iDim=0; iDim < 3; iDim++) {
            int iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (auto i=0; i < nBldPts; i++) {
                    tmpArray[i] = brFSIData[iTurbLoc][3].bld_vel[(iStart*6)+3+iDim];
                    iStart++;
                }
                std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_rotvel"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
        }
        for (size_t iDim=0;iDim < 3; iDim++) {
            int iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (auto i=0; i < nBldPts; i++) {
                    tmpArray[i] = brFSIData[iTurbLoc][3].bld_ld[(iStart*6)+iDim];
                    iStart++;
                }
                std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_ld"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
        }

        std::vector<double> ld_loc(3*nTotBldPts,0.0);
        for (auto i=0; i < nTotBldPts; i++) {
            applyWMrotation(&brFSIData[iTurbLoc][3].bld_def[i*6+3], &brFSIData[iTurbLoc][3].bld_ld[i*6], &ld_loc[i*3]);
        }
        for (size_t iDim=0;iDim < 3; iDim++) {
            int iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (auto i=0; i < nBldPts; i++) {
                    tmpArray[i] = ld_loc[iStart*3+iDim];
                    iStart++;
                }
                std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_ld_loc"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
        }

        for (size_t iDim=0; iDim < 3; iDim++) {
            int iStart = 0 ;
            for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
                for (auto i=0; i < nBldPts; i++) {
                    tmpArray[i] = brFSIData[iTurbLoc][3].bld_ld[(iStart*6)+3+iDim];
                    iStart++;
                }
                std::vector<size_t> start_dim{n_tsteps,iBlade,iDim,0};
                ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_moment"], start_dim.data(), count_dim.data(), tmpArray.data());
            }
        }

    }

    {
        for (size_t iBlade=0; iBlade < nBlades; iBlade++) {

            std::vector<size_t> start_dim{n_tsteps, iBlade, 0};
            std::vector<size_t> count_dim{1,1,3};

            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_root_disp"],
                                      start_dim.data(),
                                      count_dim.data(),
                                      &brFSIData[iTurbLoc][3].bld_root_def[iBlade*6+0]);

            ierr = nc_put_vara_double(ncid, ncOutVarIDs_["bld_root_orient"],
                                      start_dim.data(),
                                      count_dim.data(),
                                      &brFSIData[iTurbLoc][3].bld_root_def[iBlade*6+3]);
        }
    }

    {
        std::vector<size_t> start_dim{n_tsteps, 0};
        std::vector<size_t> count_dim{1,3};
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["hub_disp"], start_dim.data(), count_dim.data(), &brFSIData[iTurbLoc][3].hub_def[0]);
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["hub_orient"], start_dim.data(), count_dim.data(), &brFSIData[iTurbLoc][3].hub_def[3]);
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["hub_vel"], start_dim.data(), count_dim.data(), &brFSIData[iTurbLoc][3].hub_vel[0]);
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["hub_rotvel"], start_dim.data(), count_dim.data(), &brFSIData[iTurbLoc][3].hub_vel[3]);

        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["nac_disp"], start_dim.data(), count_dim.data(), &brFSIData[iTurbLoc][3].nac_def[0]);
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["nac_orient"], start_dim.data(), count_dim.data(), &brFSIData[iTurbLoc][3].nac_def[3]);
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["nac_vel"], start_dim.data(), count_dim.data(), &brFSIData[iTurbLoc][3].nac_vel[0]);
        ierr = nc_put_vara_double(ncid, ncOutVarIDs_["nac_rotvel"], start_dim.data(), count_dim.data(), &brFSIData[iTurbLoc][3].nac_vel[3]);
    }

    }

    nc_close(ncid);



}

void fast::OpenFAST::writeRestartFile(int iTurbLoc, int n_t_global) {

    /* // NetCDF stuff to write states to restart file or read back from it */

    int ncid;
    //Find the file and open it in append mode
    std::stringstream rstfile_ss;
    rstfile_ss << "turb_" ;
    rstfile_ss << std::setfill('0') << std::setw(2) << turbineData[iTurbLoc].TurbID;
    rstfile_ss << "_rst.nc";
    std::string rst_filename = rstfile_ss.str();
    int ierr = nc_open(rst_filename.c_str(), NC_WRITE, &ncid);
    check_nc_error(ierr, "nc_open");

    size_t count1=1;
    int tStepRatio = time_step_ratio(dtFAST, dtDriver);
    size_t n_tsteps = n_t_global/tStepRatio/restartFreq_ - 1;
    double curTime = n_t_global * dtFAST;
    ierr = nc_put_vara_double(ncid, ncRstVarIDs_["time"], &n_tsteps, &count1, &curTime);

    if ( (turbineData[iTurbLoc].sType == EXTINFLOW) && (turbineData[iTurbLoc].inflowType == 2) ){

        int nvelpts = get_numVelPtsLoc(iTurbLoc);
        int nfpts = get_numForcePtsLoc(iTurbLoc);

        const std::vector<size_t> velPtsDataDims{1, 1, static_cast<size_t>(3*nvelpts)};
        const std::vector<size_t> forcePtsDataDims{1, 1, static_cast<size_t>(3*nfpts)};
        const std::vector<size_t> forcePtsOrientDataDims{1, 1, static_cast<size_t>(9*nfpts)};

        for (size_t j=0; j < 4; j++) { // Loop over states - NM2, STATE_NM1, N, NP1

            const std::vector<size_t> start_dim{n_tsteps,j,0};

            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["x_vel"], start_dim.data(), velPtsDataDims.data(), velForceNodeData[iTurbLoc][j].x_vel.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["vel_vel"], start_dim.data(), velPtsDataDims.data(), velForceNodeData[iTurbLoc][j].vel_vel.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["x_force"], start_dim.data(), forcePtsDataDims.data(), velForceNodeData[iTurbLoc][j].x_force.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["xdot_force"], start_dim.data(), forcePtsDataDims.data(), velForceNodeData[iTurbLoc][j].xdot_force.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["vel_force"], start_dim.data(), forcePtsDataDims.data(), velForceNodeData[iTurbLoc][j].vel_force.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["force"], start_dim.data(), forcePtsDataDims.data(), velForceNodeData[iTurbLoc][j].force.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["orient_force"], start_dim.data(), forcePtsOrientDataDims.data(), velForceNodeData[iTurbLoc][j].orient_force.data());
        }

    } else if (turbineData[iTurbLoc].sType == EXTLOADS) {

        int nPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
        int nTotBldPts = turbineData[iTurbLoc].nTotBRfsiPtsBlade;
        int nBlades = turbineData[iTurbLoc].numBlades;
        const std::vector<size_t> twrDataDims{1, 1, static_cast<size_t>(6*nPtsTwr)};
        const std::vector<size_t> bldDataDims{1, 1, static_cast<size_t>(6*nTotBldPts)};
        const std::vector<size_t> bldRootDataDims{1, 1, static_cast<size_t>(6*nBlades)};
        const std::vector<size_t> bldPitchDataDims{1, 1, static_cast<size_t>(nBlades)};
        const std::vector<size_t> ptDataDims{1, 1, 6};

        for (size_t j=0; j < 4; j++) { // Loop over states - STATE_NM2, STATE_NM1, STATE_N, STATE_NP1

            const std::vector<size_t> start_dim{n_tsteps, j, 0};

            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["twr_def"], start_dim.data(), twrDataDims.data(), brFSIData[iTurbLoc][j].twr_def.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["twr_vel"], start_dim.data(), twrDataDims.data(), brFSIData[iTurbLoc][j].twr_vel.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["twr_ld"], start_dim.data(), twrDataDims.data(), brFSIData[iTurbLoc][j].twr_ld.data());

            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["bld_def"], start_dim.data(), bldDataDims.data(), brFSIData[iTurbLoc][j].bld_def.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["bld_vel"], start_dim.data(), bldDataDims.data(), brFSIData[iTurbLoc][j].bld_vel.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["bld_ld"], start_dim.data(), bldDataDims.data(), brFSIData[iTurbLoc][j].bld_ld.data());

            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["hub_def"], start_dim.data(), ptDataDims.data(), brFSIData[iTurbLoc][j].hub_def.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["hub_vel"], start_dim.data(), ptDataDims.data(), brFSIData[iTurbLoc][j].hub_vel.data());

            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["nac_def"], start_dim.data(), ptDataDims.data(), brFSIData[iTurbLoc][j].nac_def.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["nac_vel"], start_dim.data(), ptDataDims.data(), brFSIData[iTurbLoc][j].nac_vel.data());

            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["bld_root_def"], start_dim.data(), bldRootDataDims.data(), brFSIData[iTurbLoc][j].bld_root_def.data());
            ierr = nc_put_vara_double(ncid, ncRstVarIDs_["bld_pitch"], start_dim.data(), bldPitchDataDims.data(), brFSIData[iTurbLoc][j].bld_pitch.data());
        }

    }

    nc_close(ncid);


}

// Mostly Blade-resolved stuff after this

void fast::OpenFAST::get_ref_positions_from_openfast(int iTurb) {

    if(turbineData[iTurb].sType == EXTLOADS) {

        for (int i=0; i < 3; i++) {
            brFSIData[iTurb][fast::STATE_NP1].hub_ref_pos[i] = extld_p_f_FAST[iTurb].hubRefPos[i] + turbineData[iTurb].TurbineBasePos[i];
            brFSIData[iTurb][fast::STATE_NP1].nac_ref_pos[i] = extld_p_f_FAST[iTurb].nacRefPos[i] + turbineData[iTurb].TurbineBasePos[i];
            brFSIData[iTurb][fast::STATE_NP1].hub_ref_pos[i+3] = extld_p_f_FAST[iTurb].hubRefPos[i+3];
            brFSIData[iTurb][fast::STATE_NP1].nac_ref_pos[i+3] = extld_p_f_FAST[iTurb].nacRefPos[i+3];
        }

        int nBlades = turbineData[iTurb].numBlades;
        int iRunTot = 0;
        for (int i=0; i < nBlades; i++) {
            int nPtsBlade = turbineData[iTurb].nBRfsiPtsBlade[i];
            for (int j=0; j < nPtsBlade; j++) {
                for (int k=0; k < 3; k++) {
                    brFSIData[iTurb][fast::STATE_NP1].bld_ref_pos[iRunTot*6+k] = extld_p_f_FAST[iTurb].bldRefPos[iRunTot*6+k] + turbineData[iTurb].TurbineBasePos[k];
                    brFSIData[iTurb][fast::STATE_NP1].bld_ref_pos[iRunTot*6+k+3] = extld_p_f_FAST[iTurb].bldRefPos[iRunTot*6+k+3];
                }
                brFSIData[iTurb][fast::STATE_NP1].bld_chord[iRunTot] = extld_p_f_FAST[iTurb].bldChord[iRunTot];
                brFSIData[iTurb][fast::STATE_NP1].bld_rloc[iRunTot] = extld_p_f_FAST[iTurb].bldRloc[iRunTot];
                iRunTot++;
            }

            for (int k=0; k < 3; k++) {
                brFSIData[iTurb][fast::STATE_NP1].bld_root_ref_pos[i*6+k] = extld_p_f_FAST[iTurb].bldRootRefPos[i*6+k] + turbineData[iTurb].TurbineBasePos[k];
                brFSIData[iTurb][fast::STATE_NP1].bld_root_ref_pos[i*6+k+3] = extld_p_f_FAST[iTurb].bldRootRefPos[i*6+k+3];
            }

        }

        int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
        for (int i=0; i < nPtsTwr; i++) {
            for (int j = 0; j < 3; j++) {
                brFSIData[iTurb][fast::STATE_NP1].twr_ref_pos[i*6+j] = extld_p_f_FAST[iTurb].twrRefPos[i*6+j] + turbineData[iTurb].TurbineBasePos[j];
                brFSIData[iTurb][fast::STATE_NP1].twr_ref_pos[i*6+j+3] = extld_p_f_FAST[iTurb].twrRefPos[i*6+j+3];
            }
        }

    } else if(turbineData[iTurb].sType == EXTINFLOW) {

        if (turbineData[iTurb].inflowType == 2) {
            int nfpts = get_numForcePtsLoc(iTurb);
            for (auto i=0; i<nfpts; i++) {
                for (auto j=0; j < 3; j++)
                    velForceNodeData[iTurb][fast::STATE_NP1].xref_force[i*3+j] = velForceNodeData[iTurb][fast::STATE_NP1].x_force[i*3+j] ;
            }
        }
    }

}

void fast::OpenFAST::getBladeChord(double * bldChord, int iTurbGlob) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    for (int i=0; i < nBlades; i++) {
        int nPtsBlade = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        for(int j=0; j<nPtsBlade; j++) {
            bldChord[iRunTot] = brFSIData[iTurbLoc][fast::STATE_NP1].bld_chord[iRunTot];
            iRunTot++;
        }
    }
}

void fast::OpenFAST::getBladeRloc(double * bldRloc, int iTurbGlob) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    for (int i=0; i < nBlades; i++) {
        int nPtsBlade = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        for(int j=0; j<nPtsBlade; j++) {
            bldRloc[iRunTot] = brFSIData[iTurbLoc][fast::STATE_NP1].bld_rloc[iRunTot];
            iRunTot++;
        }
    }
}

void fast::OpenFAST::getBladeRefPositions(double* bldRefPos, int iTurbGlob, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    for (int i=0; i < nBlades; i++) {
        int nPtsBlade = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        for(int j=0; j<nPtsBlade; j++) {
            for (int k=0; k < nSize; k++) {
                bldRefPos[iRunTot*6+k] = brFSIData[iTurbLoc][fast::STATE_NP1].bld_ref_pos[iRunTot*6+k];
            }
            iRunTot++;
        }
    }

}

void fast::OpenFAST::getBladeRootRefPositions(double* bldRootRefPos, int iTurbGlob, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    for (int i=0; i < nBlades; i++) {
        for (int k=0; k < nSize; k++) {
            bldRootRefPos[i*6+k] = brFSIData[iTurbLoc][fast::STATE_NP1].bld_root_ref_pos[i*6+k];
        }
    }

}

void fast::OpenFAST::getBladeDisplacements(double* bldDefl, double* bldVel, int iTurbGlob, fast::timeStep t, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    std::ofstream bld_bm_mesh;
    bld_bm_mesh.open("blade_beam_mesh.csv", std::ios_base::out) ;

    for (int i=0; i < nBlades; i++) {
        int nPtsBlade = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        for(int j=0; j<nPtsBlade; j++) {

            bld_bm_mesh << brFSIData[iTurbLoc][t].bld_ref_pos[iRunTot*6] << ","
                        << brFSIData[iTurbLoc][t].bld_ref_pos[iRunTot*6+1] << ","
                        << brFSIData[iTurbLoc][t].bld_ref_pos[iRunTot*6+2] << ","
                        << brFSIData[iTurbLoc][t].bld_def[iRunTot*6] << ","
                        << brFSIData[iTurbLoc][t].bld_def[iRunTot*6+1] << ","
                        << brFSIData[iTurbLoc][t].bld_def[iRunTot*6+2] << ","
                        << brFSIData[iTurbLoc][t].bld_def[iRunTot*6+3] << ","
                        << brFSIData[iTurbLoc][t].bld_def[iRunTot*6+4] << ","
                        << brFSIData[iTurbLoc][t].bld_def[iRunTot*6+5] << ","
                        << brFSIData[iTurbLoc][t].bld_vel[iRunTot*6] << ","
                        << brFSIData[iTurbLoc][t].bld_vel[iRunTot*6+1] << ","
                        << brFSIData[iTurbLoc][t].bld_vel[iRunTot*6+2] << ","
                        << brFSIData[iTurbLoc][t].bld_vel[iRunTot*6+3] << ","
                        << brFSIData[iTurbLoc][t].bld_vel[iRunTot*6+4] << ","
                        << brFSIData[iTurbLoc][t].bld_vel[iRunTot*6+5] << std::endl;

            for (int k=0; k < nSize; k++) {
                bldDefl[iRunTot*6+k] = brFSIData[iTurbLoc][t].bld_def[iRunTot*6+k];
                bldVel[iRunTot*6+k] = brFSIData[iTurbLoc][t].bld_vel[iRunTot*6+k];
            }
            iRunTot++;
        }
    }
    bld_bm_mesh.close();

}

void fast::OpenFAST::getBladeRootDisplacements(double* bldRootDefl, int iTurbGlob, fast::timeStep t, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);

    for (int i=0; i < nBlades; i++) {
        for (int k=0; k < nSize; k++)
            bldRootDefl[i*6+k] = brFSIData[iTurbLoc][t].bld_root_def[i*6+k];
    }
}

void fast::OpenFAST::getBladePitch(double* bldPitch, int iTurbGlob, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < nSize; j++)
        bldPitch[j] = brFSIData[iTurbLoc][fast::STATE_NP1].bld_pitch[j];

}

void fast::OpenFAST::getTowerRefPositions(double* twrRefPos, int iTurbGlob, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
    for (int i=0; i < nPtsTwr; i++) {
        for (int j=0; j < nSize; j++) {
            twrRefPos[i*6+j] = brFSIData[iTurbLoc][fast::STATE_NP1].twr_ref_pos[i*6+j];
        }
    }

}

void fast::OpenFAST::getTowerDisplacements(double* twrDefl, double* twrVel, int iTurbGlob, fast::timeStep t, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
    for (int i=0; i < nPtsTwr; i++) {
        for (int j=0; j < nSize; j++) {
            twrDefl[i*6+j] = brFSIData[iTurbLoc][t].twr_def[i*6+j];
            twrVel[i*6+j] = brFSIData[iTurbLoc][t].twr_vel[i*6+j];
        }
    }

}

void fast::OpenFAST::getHubRefPosition(double* hubRefPos, int iTurbGlob, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < nSize; j++)
        hubRefPos[j] = brFSIData[iTurbLoc][fast::STATE_NP1].hub_ref_pos[j];

}

void fast::OpenFAST::getHubDisplacement(double* hubDefl, double* hubVel, int iTurbGlob, fast::timeStep t, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < nSize; j++) {
        hubDefl[j] = brFSIData[iTurbLoc][t].hub_def[j];
        hubVel[j] = brFSIData[iTurbLoc][t].hub_vel[j];
    }

}

void fast::OpenFAST::getNacelleRefPosition(double* nacRefPos, int iTurbGlob, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < nSize; j++)
        nacRefPos[j] = brFSIData[iTurbLoc][fast::STATE_NP1].nac_ref_pos[j];

}


void fast::OpenFAST::getNacelleDisplacement(double* nacDefl, double* nacVel, int iTurbGlob, fast::timeStep t, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < nSize; j++) {
        nacDefl[j] = brFSIData[iTurbLoc][t].nac_def[j];
        nacVel[j] = brFSIData[iTurbLoc][t].nac_vel[j];
    }

}

void fast::OpenFAST::setBladeForces(double* bldForces, int iTurbGlob, fast::timeStep t, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    for (int i=0; i < nBlades; i++) {
        int nPtsBlade = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        for(int j=0; j < nPtsBlade; j++) {
            for(int k=0; k < nSize; k++) {
                brFSIData[iTurbLoc][t].bld_ld[6*iRunTot+k] = bldForces[6*iRunTot+k];
            }
            iRunTot++;
        }
    }

    //TODO: May be calculate the residual as well.
}

void fast::OpenFAST::setTowerForces(double* twrForces, int iTurbGlob, fast::timeStep t, int nSize) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
    for (int i=0; i < nPtsTwr; i++)
        for (int j=0; j < nSize; j++)
            brFSIData[iTurbLoc][t].twr_ld[i*6+j] = twrForces[i*6+j];
    //TODO: May be calculate the residual as well.

}

//! Sets a uniform X force at all blade nodes
void fast::OpenFAST::setUniformXBladeForces(double loadX) {

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        int iTurbGlob = turbineMapProcToGlob[iTurb];
        int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
        std::vector<double> fsiForceTower(6*nPtsTwr,0.0);
        setTowerForces(fsiForceTower, iTurbGlob, fast::STATE_NP1);

        size_t nBlades = get_numBladesLoc(iTurb);
        size_t nTotPtsBlade = 0;
        for(int iBlade=0; iBlade < nBlades; iBlade++)
            nTotPtsBlade +=  turbineData[iTurb].nBRfsiPtsBlade[iBlade];

        std::vector<double> fsiForceBlade(6*nTotPtsBlade, 0.0);
        std::vector<double> dr(nTotPtsBlade, 0.0);

        size_t iNode=0;
        for(int iBlade=0; iBlade < nBlades; iBlade++) {
            int nBldPts = turbineData[iTurb].nBRfsiPtsBlade[iBlade];
            dr[iNode] = 0.5*(brFSIData[iTurb][3].bld_rloc[iNode+1] - brFSIData[iTurb][3].bld_rloc[iNode]);
            iNode++;

            for(int i=1; i < nBldPts-1; i++) {
                dr[iNode] = 0.5*(brFSIData[iTurb][3].bld_rloc[iNode+1] - brFSIData[iTurb][3].bld_rloc[iNode-1]);
                iNode++;
            }
            dr[iNode] = 0.5*(brFSIData[iTurb][3].bld_rloc[iNode] - brFSIData[iTurb][3].bld_rloc[iNode-1]);
            iNode++;
        }

        for(int i=0; i < nTotPtsBlade; i++)
            fsiForceBlade[i*6] = loadX * dr[i]; // X component of force

        setBladeForces(fsiForceBlade, iTurbGlob, fast::STATE_NP1);

    }
}
