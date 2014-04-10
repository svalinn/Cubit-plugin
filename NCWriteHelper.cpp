/*
 * NCWriteHelper.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: iulian
 */

#include "NCWriteHelper.hpp"
#include "NCWriteEuler.hpp"
#include "NCWriteFV.hpp"
#include "NCWriteHOMME.hpp"
#include "NCWriteMPAS.hpp"

namespace moab {

//! Get appropriate helper instance for WriteNC class; based on some info in the file set
NCWriteHelper* NCWriteHelper::get_nc_helper(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
{
  std::string& grid_type = writeNC->grid_type;
  if (grid_type == "CAM_EUL")
    return new (std::nothrow) NCWriteEuler(writeNC, fileId, opts, fileSet);
  else if (grid_type == "CAM_FV")
    return new (std::nothrow) NCWriteFV(writeNC, fileId, opts, fileSet);
  else if (grid_type == "CAM_SE")
    return new (std::nothrow) NCWriteHOMME(writeNC, fileId, opts, fileSet);
  else if (grid_type == "MPAS")
    return new (std::nothrow) NCWriteMPAS(writeNC, fileId, opts, fileSet);

  // Unknown NetCDF grid
  return NULL;
}

} /* namespace moab */
