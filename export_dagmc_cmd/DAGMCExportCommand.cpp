#include "DAGMCExportCommand.hpp"
#include "CubitInterface.hpp"

// CGM includes
#include "GeometryQueryTool.hpp"
#include "Surface.hpp"
#include "RefFace.hpp"
#include "Curve.hpp"
#include "RefEdge.hpp"
#include "SenseEntity.hpp"

// MOAB includes
#include "MBTagConventions.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"

DAGMCExportCommand::DAGMCExportCommand() :
  geom_tag(0), id_tag(0), name_tag(0), category_tag(0), faceting_tol_tag(0), geometry_resabs_tag(0)
{
  moab::ErrorCode rval;

  mdbImpl = new moab::Core();
  myGeomTool = new moab::GeomTopoTool(mdbImpl);
  
  // get some tag handles
  int negone = -1, zero = 0 /*, negonearr[] = {-1, -1, -1, -1}*/;
  rval = mdbImpl->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                                 geom_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT, &negone);
  assert(!rval);
  rval = mdbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                                 id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &zero);
  assert(!rval);
  rval = mdbImpl->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                                 name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  assert(!rval);

  rval = mdbImpl->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                                 category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  assert(!rval);
  rval = mdbImpl->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
                                 moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  assert(!rval);
  rval = mdbImpl->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, 
                                 geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  assert(!rval);

}

DAGMCExportCommand::~DAGMCExportCommand()
{}

std::vector<std::string> DAGMCExportCommand::get_syntax()
{
  // Define the syntax for the command. Note the syntax is a modified BNF
  // format. Full documentation on the command specification syntax can be
  // found in the documentation.
  std::string syntax =
      "export dagmc <string:label='filename',help='<filename>'> "
      "[overwrite]";

  std::vector<std::string> syntax_list;
  syntax_list.push_back(syntax);

  return syntax_list;
}

std::vector<std::string> DAGMCExportCommand::get_syntax_help()
{
  std::vector<std::string> help;
  return help;
}

std::vector<std::string> DAGMCExportCommand::get_help()
{
  std::vector<std::string> help;
  return help;
}

bool DAGMCExportCommand::execute(CubitCommandData &data)
{

  bool result = true;
  moab::ErrorCode rval;

  // Create entity sets for all geometric entities
  std::map<RefEntity*, moab::EntityHandle> entmap[5];

  rval = create_entity_sets(entmap);
  // if (MB_SUCCESS != rval ) // what should error handling look like?

  rval = create_topology(entmap);

  rval = store_surface_senses(entmap[2], entmap[3]);

  rval = store_curve_senses(entmap[1], entmap[2]);

  return result;
}


moab::ErrorCode DAGMCExportCommand::create_entity_sets(std::map<RefEntity*, moab::EntityHandle> (&entmap)[5])
{

  moab::ErrorCode rval;
  const char geom_categories[][CATEGORY_TAG_SIZE] =
              {"Vertex\0", "Curve\0", "Surface\0", "Volume\0", "Group\0"};
  const char* const names[] = {"Vertex", "Curve", "Surface", "Volume"};
  DLIList<RefEntity*> entlist;

  for (int dim = 0; dim < 4; dim++) {
    entlist.clean_out();
    GeometryQueryTool::instance()->ref_entity_list(names[dim], entlist, true);
    entlist.reset();

    std::ostringstream message;
    message << "Found " << entlist.size() << " entities of dimension " << dim << std::endl;

    console = CubitInterface::get_cubit_message_handler();
    console->print_message(message.str().c_str());


    for (int i = entlist.size(); i--; ) {
      RefEntity* ent = entlist.get_and_step();
      moab::EntityHandle handle;
      // Create the new meshset
      rval = mdbImpl->create_meshset(dim == 1 ? moab::MESHSET_ORDERED : moab::MESHSET_SET, handle);
      if (moab::MB_SUCCESS != rval)
        return rval;

      // Map the geom reference entity to the corresponding moab meshset
      entmap[dim][ent] = handle;

      // Create tags for the new meshset
      rval = mdbImpl->tag_set_data(geom_tag, &handle, 1, &dim);
      if (moab::MB_SUCCESS != rval)
        return rval;

      int id = ent->id();
      rval = mdbImpl->tag_set_data(id_tag, &handle, 1, &id);
      if (moab::MB_SUCCESS != rval)
        return rval;

      rval = mdbImpl->tag_set_data(category_tag, &handle, 1, &geom_categories[dim]);
      if (moab::MB_SUCCESS != rval)
        return rval;
    }
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::create_topology(std::map<RefEntity*, moab::EntityHandle> (&entitymap)[5])
{
  moab::ErrorCode rval;
  DLIList<RefEntity*> entitylist;
  std::map<RefEntity*, moab::EntityHandle>::iterator ci;

  for (int dim = 1; dim < 4; ++dim) {
    for (ci = entitymap[dim].begin(); ci != entitymap[dim].end(); ++ci) {
      entitylist.clean_out();
      ci->first->get_child_ref_entities(entitylist);

      entitylist.reset();
      for (int i = entitylist.size(); i--; ) {
        RefEntity* ent = entitylist.get_and_step();
        moab::EntityHandle h = entitymap[dim - 1][ent];
        rval = mdbImpl->add_parent_child(ci->second, h);

        // std::ostringstream message;
        // message << "Created parent-child relationship between " << ci->second << " and " << h << std::endl;
        
        // console = CubitInterface::get_cubit_message_handler();
        // console->print_message(message.str().c_str());

        if (moab::MB_SUCCESS != rval)
          return rval;
      }
    }
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::store_surface_senses(std::map<RefEntity*, moab::EntityHandle>& surface_map,
                                                         std::map<RefEntity*, moab::EntityHandle>& volume_map)
{
  moab::ErrorCode rval;
  std::map<RefEntity*, moab::EntityHandle>::iterator ci;

  for (ci = surface_map.begin(); ci != surface_map.end(); ++ci) {
    RefFace* face = (RefFace*)(ci->first);
    BasicTopologyEntity *forward = 0, *reverse = 0;
    for (SenseEntity* cf = face->get_first_sense_entity_ptr();
         cf; cf = cf->next_on_bte()) {
      BasicTopologyEntity* vol = cf->get_parent_basic_topology_entity_ptr();
      // Allocate vol to the proper topology entity (forward or reverse)
      if (cf->get_sense() == CUBIT_UNKNOWN ||
          cf->get_sense() != face->get_surface_ptr()->bridge_sense()) {
        // Check that each surface has a sense for only one volume
        if (reverse) {
          // std::cout << "Surface " << face->id() << " has reverse sense " <<
          //              "with multiple volume " << reverse->id() << " and " <<
          //              "volume " << vol->id() << std::endl;
          return moab::MB_FAILURE;
        }
        reverse = vol;
      }
      if (cf->get_sense() == CUBIT_UNKNOWN ||
          cf->get_sense() == face->get_surface_ptr()->bridge_sense()) {
        // Check that each surface has a sense for only one volume
        if (forward) {
          // std::cout << "Surface " << face->id() << " has forward sense " <<
          //              "with multiple volume " << forward->id() << " and " <<
          //              "volume " << vol->id() << std::endl;
          return moab::MB_FAILURE;
        }
        forward = vol;
      }
    }

    console = CubitInterface::get_cubit_message_handler();

    if (forward) {
      rval = myGeomTool->set_sense(ci->second, volume_map[forward], moab::SENSE_FORWARD);
      std::ostringstream message;
      message << "Surface " << ci->second << " has forward sense with respect to volume "  << volume_map[forward] << std::endl;
      
      console->print_message(message.str().c_str());

      if (moab::MB_SUCCESS != rval)
        return rval;
    }
    if (reverse) {
      rval = myGeomTool->set_sense(ci->second, volume_map[reverse], moab::SENSE_REVERSE);
      std::ostringstream message;
      message << "Surface " << ci->second << " has reverse sense with respect to volume "  << volume_map[reverse] << std::endl;
      
      console->print_message(message.str().c_str());
      if (moab::MB_SUCCESS != rval)
        return rval;
    }
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::store_curve_senses(std::map<RefEntity*, moab::EntityHandle>& curve_map,
                                      std::map<RefEntity*, moab::EntityHandle>& surface_map)
{
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> ents;
  std::vector<int> senses;
  std::map<RefEntity*, moab::EntityHandle>::iterator ci;
  for (ci = curve_map.begin(); ci != curve_map.end(); ++ci) {
    RefEdge* edge = (RefEdge*)(ci->first);
    ents.clear();
    senses.clear();
    for (SenseEntity* ce = edge->get_first_sense_entity_ptr();
         ce; ce = ce->next_on_bte()) {
      BasicTopologyEntity* fac = ce->get_parent_basic_topology_entity_ptr();
      moab::EntityHandle face = surface_map[fac];
      if (ce->get_sense() == CUBIT_UNKNOWN ||
          ce->get_sense() != edge->get_curve_ptr()->bridge_sense()) {
        ents.push_back(face);
        senses.push_back(moab::SENSE_REVERSE);
      }
      if (ce->get_sense() == CUBIT_UNKNOWN ||
          ce->get_sense() == edge->get_curve_ptr()->bridge_sense()) {
        ents.push_back(face);
        senses.push_back(moab::SENSE_FORWARD);
      }
    }

    rval = myGeomTool->set_senses(ci->second, ents, senses);
    std::ostringstream message;
    message << "Curve " << ci->second << " has "  ;
    for (int ent_num=0;ent_num < ents.size();ent_num++)
      {
        message << std::endl << "\t" << (ent_num>0?"and ":"")
                << (senses[ent_num]==moab::SENSE_FORWARD?"forward":"reverse")
                << " sense with respect to surface " << ents[ent_num];
      }
    message << std::endl;

    console->print_message(message.str().c_str());


    if (moab::MB_SUCCESS != rval)
      return rval;
  }
  return moab::MB_SUCCESS;
}
