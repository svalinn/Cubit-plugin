#ifndef MYEXPORTCOMMAND_HPP
#define MYEXPORTCOMMAND_HPP

#include "CubitCommandInterface.hpp"

// CGM includes
#include "RefEntity.hpp"

// MOAB includes
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"

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

  moab::ErrorCode create_entity_sets(std::map<RefEntity*, moab::EntityHandle> (&entmap)[5]);
  moab::ErrorCode create_topology(std::map<RefEntity*, moab::EntityHandle> (&entitymap)[5]);
  moab::ErrorCode store_surface_senses(std::map<RefEntity*, moab::EntityHandle>& surface_map,
                                       std::map<RefEntity*, moab::EntityHandle>& volume_map);

private:

  moab::Interface* mdbImpl;
  moab::GeomTopoTool* myGeomTool;

  moab::Tag geom_tag, id_tag, name_tag, category_tag, faceting_tol_tag, geometry_resabs_tag;
};

#endif // MYEXPORTCOMMAND_HPP
