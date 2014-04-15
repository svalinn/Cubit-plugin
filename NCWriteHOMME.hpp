/*
 * NCWriteHOMME.hpp
 *
 *  nc write helper for HOMME type data (CAM)
 *  Created on: April 9, 2014
 *
 */

#ifndef NCWRITEHOMME_HPP_
#define NCWRITEHOMME_HPP_

#include "NCWriteHelper.hpp"

namespace moab {

class NCWriteHOMME: public NCWriteHelper
{
public:
  NCWriteHOMME(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet) :
    NCWriteHelper(writeNC, fileId, opts, fileSet) {}

  virtual ~NCWriteHOMME();

private:
  ErrorCode write_values(std::vector<std::string>& var_names);
};

} // namespace moab

#endif
