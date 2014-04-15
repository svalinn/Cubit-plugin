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
  NCWriteEuler(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet) :
    NCWriteHelper(writeNC, fileId, opts, fileSet) {}

  virtual ~NCWriteEuler();

private:
  ErrorCode write_values(std::vector<std::string>& var_names);
};

} // namespace moab

#endif
