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

class NCWriteMPAS: public UcdNCWriteHelper
{
public:
  NCWriteMPAS(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
: UcdNCWriteHelper(writeNC, fileId, opts, fileSet) {}

  virtual ~NCWriteMPAS();

private:
  //! Implementation of NCWriteHelper::collect_mesh_info()
  virtual ErrorCode collect_mesh_info();

  //! Collect data for specified variables
  virtual ErrorCode write_values(std::vector<std::string>& var_names);
};

} // namespace moab

#endif
