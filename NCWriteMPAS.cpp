/*
 * NCWriteMPAS.cpp
 *
 *  Created on: April 9, 2014
 */

#include "NCWriteMPAS.hpp"
#include "moab/WriteUtilIface.hpp"

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { _writeNC->mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { _writeNC->mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

NCWriteMPAS::~NCWriteMPAS()
{
  // TODO Auto-generated destructor stub
}

ErrorCode NCWriteMPAS::collect_mesh_info()
{
  return MB_NOT_IMPLEMENTED;
}

ErrorCode NCWriteMPAS::write_values(std::vector<std::string>& /* var_names */)
{
  return MB_NOT_IMPLEMENTED;
}

} /* namespace moab */
