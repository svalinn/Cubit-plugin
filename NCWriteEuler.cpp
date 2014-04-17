/*
 * NCWriteEuler.cpp
 *
 *  Created on: Mar 28, 2014
 */

#include "NCWriteEuler.hpp"
#include "moab/WriteUtilIface.hpp"

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { _writeNC->mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { _writeNC->mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

NCWriteEuler::~NCWriteEuler()
{
  // TODO Auto-generated destructor stub
}

} /* namespace moab */
