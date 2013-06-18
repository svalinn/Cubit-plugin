#include "NCHelperMPAS.hpp"
#include "moab/ReadUtilIface.hpp"
#include "FileOptions.hpp"
#include "moab/SpectralMeshTool.hpp"

#include <cmath>

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {_readNC->readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
    if (err) {_readNC->readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

NCHelperMPAS::NCHelperMPAS(ReadNC* readNC, int fileId, const FileOptions& opts) : UcdNCHelper(readNC, fileId)
{
  if (MB_SUCCESS == opts.match_option("PARTITION_METHOD", "NODAL_PARTITION"))
    readNC->partMethod = -1;
}

bool NCHelperMPAS::can_read_file(ReadNC* readNC, int fileId)
{
  std::vector<std::string>& dimNames = readNC->dimNames;

  // If dimension name "vertexDegree" exists then it should be the MPAS grid
  if (std::find(dimNames.begin(), dimNames.end(), std::string("vertexDegree")) != dimNames.end())
    return true;

  return false;
}

ErrorCode NCHelperMPAS::init_mesh_vals(const FileOptions& opts, EntityHandle file_set)
{
  // TBD
  return MB_SUCCESS;
}

ErrorCode NCHelperMPAS::create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& faces)
{
  // TBD
  return MB_SUCCESS;
}

ErrorCode NCHelperMPAS::read_ucd_variable_setup(std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
                                                 std::vector<ReadNC::VarData>& vdatas, std::vector<ReadNC::VarData>& vsetdatas)
{
  // TBD
  return MB_SUCCESS;
}

ErrorCode NCHelperMPAS::read_ucd_variable_to_nonset_allocate(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  // TBD
  return MB_SUCCESS;
}

#ifdef PNETCDF_FILE
ErrorCode NCHelperMPAS::read_ucd_variable_to_nonset_async(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  // TBD
  return MB_SUCCESS;
}
#else
ErrorCode NCHelperMPAS::read_ucd_variable_to_nonset(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  // TBD
  return MB_SUCCESS;
}
#endif

ErrorCode NCHelperMPAS::convert_ucd_variable(ReadNC::VarData& var_data, int tstep_num)
{
  // TBD
  return MB_SUCCESS;
}

} // namespace moab
