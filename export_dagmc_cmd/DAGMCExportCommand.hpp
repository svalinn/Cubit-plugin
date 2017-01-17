#ifndef MYEXPORTCOMMAND_HPP
#define MYEXPORTCOMMAND_HPP

#include "CubitCommandInterface.hpp"
#include "CubitMessageHandler.hpp"
#include "TtyProgressTool.hpp"

// CGM includes
#include "RefEntity.hpp"

// MOAB includes
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"

// make_watertight includes
#include "make_watertight/MakeWatertight.hpp"

typedef std::map<RefEntity*, moab::EntityHandle> refentity_handle_map;
typedef std::map<RefEntity*, moab::EntityHandle>::iterator refentity_handle_map_itor;

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
  MakeWatertight* mw;
  
protected:
  TtyProgressTool* progTool;

  moab::ErrorCode create_tags();
  moab::ErrorCode parse_options(CubitCommandData &data, moab::EntityHandle* file_set = 0);
  moab::ErrorCode create_entity_sets(refentity_handle_map (&entmap)[5]);
  moab::ErrorCode create_topology(refentity_handle_map (&entitymap)[5]);
  moab::ErrorCode store_surface_senses(refentity_handle_map& surface_map,
                                       refentity_handle_map& volume_map);
  moab::ErrorCode store_curve_senses(refentity_handle_map& curve_map,
                                     refentity_handle_map& surface_map);
  moab::ErrorCode store_groups(refentity_handle_map (&entitymap)[5]);
  moab::ErrorCode create_group_entsets(refentity_handle_map& group_map);
  moab::ErrorCode store_group_content(refentity_handle_map (&entitymap)[5]);
  moab::ErrorCode create_vertices(refentity_handle_map &vertex_map);
  moab::ErrorCode create_curve_facets(refentity_handle_map& curve_map,
                                      refentity_handle_map& vertex_map);
  moab::ErrorCode create_surface_facets(refentity_handle_map& surface_map,
                                        refentity_handle_map& vertex_map);
  moab::ErrorCode gather_ents(moab::EntityHandle gather_set);  
  moab::ErrorCode teardown();

private:

  moab::Interface* mdbImpl;
  moab::GeomTopoTool* myGeomTool;
  CubitMessageHandler* console;

  std::ostringstream message;

  moab::Tag geom_tag, id_tag, name_tag, category_tag, faceting_tol_tag, geometry_resabs_tag;

  int norm_tol;
  double faceting_tol;
  double len_tol;
  bool verbose_warnings;
  bool fatal_on_curves;
  bool make_watertight;
  
  int failed_curve_count;
  std::vector<int> failed_curves;

  int failed_surface_count;
  std::vector<int> failed_surfaces;


};

#endif // MYEXPORTCOMMAND_HPP
