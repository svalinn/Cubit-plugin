/*
 * NCWriteHelper.hpp
 *
 * Purpose       : Climate NC writer file helper; abstract, will be implemented for each type
 *
 *  Created on: Mar 28, 2014
 */

#ifndef NCWRITEHELPER_HPP_
#define NCWRITEHELPER_HPP_
#include "WriteNC.hpp"

namespace moab {

class NCWriteHelper
{
public:
  NCWriteHelper(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
  :_writeNC(writeNC), _fileId(fileId), _opts(opts), _fileSet(fileSet) {}

  virtual ~NCWriteHelper() {};

  //! Get appropriate helper instance for WriteNC class; based on some info in the file set
  static NCWriteHelper* get_nc_helper(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet);

  //! Take the info from VarData and write first the coordinates, then the actual variables
  virtual ErrorCode write_values(std::vector<std::string>& var_names) = 0;

protected:
  template <typename T> void jik_to_kji(size_t ni, size_t nj, size_t nk, T* dest, T* source)
  {
    size_t nik = ni * nk, nij = ni * nj;
    for (std::size_t k = 0; k != nk; k++)
      for (std::size_t j = 0; j != nj; j++)
        for (std::size_t i = 0; i != ni; i++)
          dest[k*nij + j*ni + i] = source[j*nik + i*nk + k];
  }

  //! Allow NCWriteHelper to directly access members of WriteNC
  WriteNC* _writeNC;

  //! Cache some information from ReadNC
  int _fileId;
  const FileOptions& _opts;
  EntityHandle _fileSet;
};

} // namespace moab

#endif

