/*
 * NCWriteMPAS.hpp
 *
 *  nc write helper for MPAS type data (CAM)
 *  Created on: April 9, 2014
 *
 */

#ifndef NCWRITEMPAS_HPP_
#define NCWRITEMPAS_HPP_

#include "NCWriteHelper.hpp"

namespace moab {

class NCWriteMPAS: public NCWriteHelper
{
public:
  NCWriteMPAS(WriteNC* writeNC,  const FileOptions& opts, EntityHandle fileSet) :
    NCWriteHelper(writeNC, opts, fileSet) {};

  virtual ~NCWriteMPAS();
};

} /* namespace moab */
#endif /* NCWRITEMPAS_HPP_ */
