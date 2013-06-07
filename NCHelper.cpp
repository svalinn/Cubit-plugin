#include "NCHelper.hpp"
#include "NCHelperEuler.hpp"
#include "NCHelperFV.hpp"
#include "NCHelperHOMME.hpp"
#include "moab/ReadUtilIface.hpp"

namespace moab {

NCHelper* NCHelper::get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts)
{
  if (NCHelperEuler::can_read_file(readNC, fileId))
    return new (std::nothrow) NCHelperEuler(readNC, fileId);
  else if (NCHelperFV::can_read_file(readNC, fileId))
    return new (std::nothrow) NCHelperFV(readNC, fileId);
  else if (NCHelperHOMME::can_read_file(readNC, fileId))
    return new (std::nothrow) NCHelperHOMME(readNC, fileId, opts);
  else // Unknown NetCDF grid (will fill this in later for POP, CICE and CLM)
    return NULL;
}

} // namespace moab
