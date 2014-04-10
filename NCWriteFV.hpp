/*
 * NCWriteFV.hpp
 *
 *  nc write helper for FV type data (CAM)
 *  Created on: April 9, 2014
 *
 */

#ifndef NCWRITEFV_HPP_
#define NCWRITEFV_HPP_

#include "NCWriteHelper.hpp"

namespace moab {

class NCWriteFV: public NCWriteHelper
{
public:
  NCWriteFV(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet) :
    NCWriteHelper(writeNC, fileId, opts, fileSet) {}

  virtual ~NCWriteFV();

private:
  ErrorCode write_values(std::vector<std::string>& var_names, EntityHandle fileSet);
};

} // namespace moab

#endif // NCWRITEFV_HPP_
