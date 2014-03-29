/*
 * NCWriteEuler.hpp
 *
 *  nc write helper for euler type data (CAM)
 *  Created on: Mar 28, 2014
 *
 */

#ifndef NCWRITEEULER_HPP_
#define NCWRITEEULER_HPP_

#include "NCWriteHelper.hpp"

namespace moab {

class NCWriteEuler: public NCWriteHelper
{
public:
  NCWriteEuler(WriteNC* writeNC,  const FileOptions& opts, EntityHandle fileSet) :
    NCWriteHelper(writeNC, opts, fileSet) {};
  virtual ~NCWriteEuler();
};

} /* namespace moab */
#endif /* NCWRITEEULER_HPP_ */
