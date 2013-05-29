#include "NCHelper.hpp"
#include "moab/ReadUtilIface.hpp"
#include "FileOptions.hpp"

namespace moab {

NCHelper* NCHelper::get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts)
{
  // Not a CAM file, will fill this in later for POP, CICE and CLM
  if (!readNC->isCam)
    return new (std::nothrow) NCHNotCam(readNC, fileId);

  // If a CAM file, which type?
  if (NCHEuler::can_read_file(readNC))
    return new (std::nothrow) NCHEuler(readNC, fileId);
  else if (NCHFV::can_read_file(readNC))
    return new (std::nothrow) NCHFV(readNC, fileId);
  else if (NCHHomme::can_read_file(readNC, opts))
    return new (std::nothrow) NCHHomme(readNC, fileId);

  // Unknown CAM grid
  return new (std::nothrow) NCHUnknownCam(readNC, fileId);
}

bool NCHEuler::can_read_file(ReadNC* readNC)
{
  std::vector<std::string>& dimNames = readNC->dimNames;

  // If dimension names "lon" AND "lat' exist then it's the Eulerian Spectral grid or the FV grid
  if ((std::find(dimNames.begin(), dimNames.end(), std::string("lon")) != dimNames.end()) && (std::find(dimNames.begin(),
    dimNames.end(), std::string("lat")) != dimNames.end()))
  {
    // If dimension names "lon" AND "lat" AND "slon" AND "slat" exist then it's the FV grid
    if ((std::find(dimNames.begin(), dimNames.end(), std::string("slon")) != dimNames.end()) && (std::find(dimNames.begin(),
        dimNames.end(), std::string("slat")) != dimNames.end()))
      return false;
    else
      return true;
  }

  return false;
}

ErrorCode NCHEuler::init_mesh_vals(const FileOptions& opts, EntityHandle file_set)
{
  ErrorCode rval = _readNC->init_EulSpcscd_vals(opts, file_set);
  if (MB_SUCCESS != rval)
    _readNC->readMeshIface->report_error("%s", "Trouble initializing Euler grid.");

  return rval;
}

bool NCHFV::can_read_file(ReadNC* readNC)
{
  std::vector<std::string>& dimNames = readNC->dimNames;

  // If dimension names "lon" AND "lat" AND "slon" AND "slat" exist then it's the FV grid
  if ((std::find(dimNames.begin(), dimNames.end(), std::string("lon")) != dimNames.end()) && (std::find(dimNames.begin(),
      dimNames.end(), std::string("lat")) != dimNames.end()) && (std::find(dimNames.begin(), dimNames.end(), std::string("slon"))
      != dimNames.end()) && (std::find(dimNames.begin(), dimNames.end(), std::string("slat")) != dimNames.end()))
    return true;

  return false;
}

ErrorCode NCHFV::init_mesh_vals(const FileOptions& opts, EntityHandle file_set)
{
  ErrorCode rval = _readNC->init_FVCDscd_vals(opts, file_set);
  if (MB_SUCCESS != rval)
    _readNC->readMeshIface->report_error("%s", "Trouble initializing FV grid.");

  return rval;
}

bool NCHHomme::can_read_file(ReadNC* readNC, const FileOptions& opts)
{
  // If global attribute "np" exists then it's the HOMME grid
  ErrorCode rval = readNC->check_np_attribute();
  if (MB_SUCCESS == rval) {
    if (MB_SUCCESS == opts.match_option("PARTITION_METHOD", "NODAL_PARTITION"))
      readNC->partMethod = -1;

    return true;
  }

  return false;
}

ErrorCode NCHHomme::init_mesh_vals(const FileOptions& opts, EntityHandle file_set)
{
  ErrorCode rval = _readNC->init_HOMMEucd_vals();
  if (MB_SUCCESS != rval)
    _readNC->readMeshIface->report_error("%s", "Trouble initializing Homme grid.");

  return rval;
}

ErrorCode NCHUnknownCam::init_mesh_vals(const FileOptions& opts, EntityHandle file_set)
{
  _readNC->readMeshIface->report_error("%s", "Unknown CAM grid.");

  return MB_FAILURE;
}

ErrorCode NCHNotCam::init_mesh_vals(const FileOptions& opts, EntityHandle file_set)
{
  _readNC->readMeshIface->report_error("%s", "Unknown grid.");

  return MB_FAILURE;
}

} // namespace moab
