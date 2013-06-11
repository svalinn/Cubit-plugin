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

NCHelperMPAS::NCHelperMPAS(ReadNC* readNC, int fileId, const FileOptions& opts) : NCHelper(readNC, fileId)
{
  if (MB_SUCCESS == opts.match_option("PARTITION_METHOD", "NODAL_PARTITION"))
    readNC->partMethod = -1;
}

bool NCHelperMPAS::can_read_file(ReadNC* readNC, int fileId)
{
  std::vector<std::string>& dimNames = readNC->dimNames;

  // If dimension name "vertexDegree" exists then it should be the MPAS grid
  if (std::find(dimNames.begin(), dimNames.end(), std::string("vertexDegree")) != dimNames.end()) {
    return true;
  }

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

} // namespace moab
