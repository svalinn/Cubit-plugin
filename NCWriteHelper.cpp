/*
 * NCWriteHelper.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: iulian
 */

#include "NCWriteHelper.hpp"
#include "NCWriteEuler.hpp"

namespace moab {

//! Get appropriate helper instance for WriteNC class; based on some info in the file set
NCWriteHelper* NCWriteHelper::get_nc_helper(WriteNC* writeNC, const FileOptions& opts, EntityHandle fileSet)
{

  // right now, get a cam euler helper
  return new (std::nothrow) NCWriteEuler(writeNC, opts, fileSet);
}


} /* namespace moab */
