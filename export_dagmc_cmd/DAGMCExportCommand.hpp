#ifndef MYEXPORTCOMMAND_HPP
#define MYEXPORTCOMMAND_HPP

#include "CubitCommandInterface.hpp"

// CGM includes
#include "RefEntity.hpp"

// MOAB includes
#include "moab/Core.hpp"

class MeshExportInterface;

/*!
 * \brief The DAGMCExportCommand class implements all the steps necessary
 * to load faceted data into a MOAB instance and export as a MOAB mesh.
 */
class DAGMCExportCommand: public CubitCommand
{
public:
  DAGMCExportCommand();
  ~DAGMCExportCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

protected:

  ErrorCode create_entity_sets(std::map<RefEntity*, EntityHandle> (&entmap)[5]);
};

#endif // MYEXPORTCOMMAND_HPP
