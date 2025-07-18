#include "DAGMCExportCommand.hpp"

#include "CubitInterface.hpp"

// CGM includes
#include "GeometryQueryTool.hpp"
#include "ModelQueryEngine.hpp"
#include "GMem.hpp"

#include "RefEntityName.hpp"

#include "Body.hpp"
#include "Surface.hpp"
#include "Curve.hpp"

#include "RefGroup.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "SenseEntity.hpp"

// MOAB includes
#include "MBTagConventions.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"

#define CHK_MB_ERR_RET(A,B)  if (moab::MB_SUCCESS != (B)) { \
  message << (A) << (B) << std::endl;                                   \
  CubitInterface::get_cubit_message_handler() ->print_message(message.str().c_str()); \
  return false;                                                         \
  }

#define CHK_MB_ERR_RET_MB(A,B)  if (moab::MB_SUCCESS != (B)) { \
  message << (A) << (B) << std::endl;                                   \
  return (B);                                                         \
  }

DAGMCExportCommand::DAGMCExportCommand() :
  geom_tag(0), id_tag(0), name_tag(0), category_tag(0), faceting_tol_tag(0), geometry_resabs_tag(0)
{
  // set default values
  norm_tol = 5;
  faceting_tol = 1e-3;
  len_tol = 0.0;
  verbose_warnings = false;
  fatal_on_curves = false;

  CubitMessageHandler* console = CubitInterface::get_cubit_message_handler();
  if (console) {
    std::ostringstream load_message;
    load_message.str("");
    load_message << "-- DAGMC export command available." << std::endl;
    console->print_error(load_message.str().c_str());
  }
}

DAGMCExportCommand::~DAGMCExportCommand()
{}

std::vector<std::string> DAGMCExportCommand::get_syntax()
{
  // Define the syntax for the command. Note the syntax is a modified BNF
  // format. Full documentation on the command specification syntax can be
  // found in the documentation.
  std::string syntax =
      "export dagmc "
      "<string:label='filename',help='<filename>'> "
      "[faceting_tolerance <value:label='faceting_tolerance',help='<faceting tolerance>'>] "
      "[length_tolerance <value:label='length_tolerance',help='<length tolerance>'>] "
      "[normal_tolerance <value:label='normal_tolerance',help='<normal tolerance>'>] "
      "[make_watertight]"
      "[verbose] [fatal_on_curves]";

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

  static moab::Core instance;
  mdbImpl = &instance;
  mw = new MakeWatertight(mdbImpl);
  myGeomTool = new moab::GeomTopoTool(mdbImpl);
  message.str("");

  bool result = true;
  moab::ErrorCode rval;

  // Create entity sets for all geometric entities
  refentity_handle_map entmap[5];

  rval = create_tags();
  CHK_MB_ERR_RET("Error initializing DAGMC export: ",rval);

  // create a file set for storage of tolerance values
  moab::EntityHandle file_set;
  rval = mdbImpl->create_meshset(0, file_set);
  CHK_MB_ERR_RET("Error creating file set.",rval);
  
  rval = parse_options(data, &file_set);
  CHK_MB_ERR_RET("Error parsing options: ",rval);

  rval = create_entity_sets(entmap);
  CHK_MB_ERR_RET("Error creating entity sets: ",rval);

  rval = create_topology(entmap);
  CHK_MB_ERR_RET("Error creating topology: ",rval);
  
  rval = store_surface_senses(entmap[2], entmap[3]);
  CHK_MB_ERR_RET("Error storing surface senses: ",rval);
  
  rval = store_curve_senses(entmap[1], entmap[2]);
  CHK_MB_ERR_RET("Error storing curve senses: ",rval);
    
  rval = store_groups(entmap);
  CHK_MB_ERR_RET("Error storing groups: ",rval);
  
  entmap[3].clear();
  entmap[4].clear();
  
  rval = create_vertices(entmap[0]);
  CHK_MB_ERR_RET("Error creating vertices: ",rval);
  
  rval = create_curve_facets(entmap[1], entmap[0]);
  CHK_MB_ERR_RET("Error faceting curves: ",rval);

  rval = create_surface_facets(entmap[2], entmap[0]);
  CHK_MB_ERR_RET("Error faceting surfaces: ",rval);

  rval = gather_ents(file_set);
  CHK_MB_ERR_RET("Could not gather entities into file set.", rval);

  if (make_watertight) {
    rval = mw->make_mesh_watertight(file_set, faceting_tol, false);
    CHK_MB_ERR_RET("Could not make the model watertight.", rval);
  }
  
  std::string filename;
  data.get_string("filename",filename);
  rval = mdbImpl->write_file(filename.c_str());
  CHK_MB_ERR_RET("Error writing file: ",rval);

  rval = teardown();
  CHK_MB_ERR_RET("Error tearing down export command.",rval);
  
  return result;
}

moab::ErrorCode DAGMCExportCommand::parse_options(CubitCommandData &data, moab::EntityHandle* file_set)
{
  moab::ErrorCode rval;

  // read parsed command for faceting tolerance
  data.get_value("faceting_tolerance",faceting_tol);
  message << "Setting faceting tolerance to " << faceting_tol << std::endl;

  // read parsed command for length tolerance
  data.get_value("length_tolerance",len_tol);
  message << "Setting length tolerance to " << len_tol << std::endl;

  // Always tag with the faceting_tol and geometry absolute resolution
  // If file_set is defined, use that, otherwise (file_set == NULL) tag the interface  
  moab::EntityHandle set = file_set ? *file_set : 0;
  rval = mdbImpl->tag_set_data(faceting_tol_tag, &set, 1, &faceting_tol);
  CHK_MB_ERR_RET_MB("Error setting faceting tolerance tag",rval);

  // read parsed command for normal tolerance
  data.get_value("normal_tolerance",norm_tol);
  message << "Setting normal tolerance to " << norm_tol << std::endl;

  rval = mdbImpl->tag_set_data(geometry_resabs_tag, &set, 1, &GEOMETRY_RESABS);
  CHK_MB_ERR_RET_MB("Error setting geometry_resabs_tag",rval);
  
  // read parsed command for verbosity
  verbose_warnings = data.find_keyword("verbose");
  fatal_on_curves = data.find_keyword("fatal_on_curves");
  make_watertight = data.find_keyword("make_watertight");
  
  if (verbose_warnings && fatal_on_curves)
    message << "This export will fail if curves fail to facet" << std::endl;

  return rval;
}

moab::ErrorCode DAGMCExportCommand::create_tags()
{
  moab::ErrorCode rval;


  // get some tag handles
  int negone = -1, zero = 0 /*, negonearr[] = {-1, -1, -1, -1}*/;
  rval = mdbImpl->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                                 geom_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_ANY, &negone);
  CHK_MB_ERR_RET_MB("Error creating geom_tag",rval);
    
  rval = mdbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                                 id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_ANY, &zero);
  CHK_MB_ERR_RET_MB("Error creating id_tag",rval);
  
  rval = mdbImpl->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                                 name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_ANY);
  CHK_MB_ERR_RET_MB("Error creating name_tag",rval);
  
  rval = mdbImpl->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                                 category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  CHK_MB_ERR_RET_MB("Error creating category_tag",rval);

  rval = mdbImpl->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
                                 moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  CHK_MB_ERR_RET_MB("Error creating faceting_tol_tag",rval);

  rval = mdbImpl->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, 
                                 geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  CHK_MB_ERR_RET_MB("Error creating geometry_resabs_tag",rval);

  return rval;
}

moab::ErrorCode DAGMCExportCommand::teardown()
{
  message  << "***** Faceting Summary Information *****" << std::endl;
  if (0 < failed_curve_count) {
    message << "----- Curve Fail Information -----" << std::endl
            << "There were " << failed_curve_count << " curves that could not be faceted." << std::endl;
  } else {
    message << "----- All curves faceted correctly  -----" << std::endl;
  }
  if (0 < failed_surface_count) {
    message << "----- Facet Fail Information -----" << std::endl
            << "There were " << failed_surface_count << " surfaces that could not be faceted." << std::endl;
  } else {
    message << "----- All surfaces faceted correctly  -----" << std::endl;
  }
  message << "***** End of Faceting Summary Information *****" << std::endl;

  CubitInterface::get_cubit_message_handler()->print_message(message.str().c_str());
  message.str("");

  
  moab::ErrorCode rval = mdbImpl->delete_mesh();
  CHK_MB_ERR_RET_MB("Error cleaning up mesh instance.", rval);
  delete myGeomTool;

  return rval;

}


moab::ErrorCode DAGMCExportCommand::create_entity_sets(refentity_handle_map (&entmap)[5])
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

    message << "Found " << entlist.size() << " entities of dimension " << dim << std::endl;

    for (int i = entlist.size(); i--;) {
      RefEntity* ent = entlist.get_and_step();
      moab::EntityHandle handle;

      // Create the new meshset
      rval = mdbImpl->create_meshset(dim == 1 ? moab::MESHSET_ORDERED : moab::MESHSET_SET, handle);
      if (moab::MB_SUCCESS != rval) return rval;

      // Map the geom reference entity to the corresponding moab meshset
      entmap[dim][ent] = handle;

      // Create tags for the new meshset
      rval = mdbImpl->tag_set_data(geom_tag, &handle, 1, &dim);
      if (moab::MB_SUCCESS != rval) return rval;

      int id = ent->id();
      rval = mdbImpl->tag_set_data(id_tag, &handle, 1, &id);
      if (moab::MB_SUCCESS != rval) return rval;

      rval = mdbImpl->tag_set_data(category_tag, &handle, 1, &geom_categories[dim]);
      if (moab::MB_SUCCESS != rval) return rval;
    }
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::create_topology(refentity_handle_map (&entitymap)[5])
{
  moab::ErrorCode rval;
  DLIList<RefEntity*> entitylist;
  refentity_handle_map_itor ci;

  for (int dim = 1; dim < 4; ++dim) {
    for (ci = entitymap[dim].begin(); ci != entitymap[dim].end(); ++ci) {
      entitylist.clean_out();
      ci->first->get_child_ref_entities(entitylist);

      entitylist.reset();
      for (int i = entitylist.size(); i--; ) {
        RefEntity* ent = entitylist.get_and_step();
        moab::EntityHandle h = entitymap[dim - 1][ent];
        rval = mdbImpl->add_parent_child(ci->second, h);

        if (moab::MB_SUCCESS != rval)
          return rval;
      }
    }
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::store_surface_senses(refentity_handle_map& surface_map,
                                                         refentity_handle_map& volume_map)
{
  moab::ErrorCode rval;
  refentity_handle_map_itor ci;

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
          message << "Surface " << face->id() << " has reverse sense " <<
            "with multiple volume " << reverse->id() << " and " <<
            "volume " << vol->id() << std::endl;
          return moab::MB_FAILURE;
        }
        reverse = vol;
      }
      if (cf->get_sense() == CUBIT_UNKNOWN ||
          cf->get_sense() == face->get_surface_ptr()->bridge_sense()) {
        // Check that each surface has a sense for only one volume
        if (forward) {
          message << "Surface " << face->id() << " has forward sense " <<
            "with multiple volume " << forward->id() << " and " <<
            "volume " << vol->id() << std::endl;
          return moab::MB_FAILURE;
        }
        forward = vol;
      }
    }

    if (forward) {
      rval = myGeomTool->set_sense(ci->second, volume_map[forward], moab::SENSE_FORWARD);
      if (moab::MB_SUCCESS != rval) return rval;
    }
    if (reverse) {
      rval = myGeomTool->set_sense(ci->second, volume_map[reverse], moab::SENSE_REVERSE);
      if (moab::MB_SUCCESS != rval) return rval;
    }
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::store_curve_senses(refentity_handle_map& curve_map,
                                      refentity_handle_map& surface_map)
{
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> ents;
  std::vector<int> senses;
  refentity_handle_map_itor ci;
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
    if (moab::MB_SUCCESS != rval) return rval;
  }
  return moab::MB_SUCCESS;
}


moab::ErrorCode DAGMCExportCommand::store_groups(refentity_handle_map (&entitymap)[5])
{
  moab::ErrorCode rval;

  // Create entity sets for all ref groups
  rval = create_group_entsets(entitymap[4]);
  if (moab::MB_SUCCESS != rval) return rval;

  // Store group names and entities in the mesh
  rval = store_group_content(entitymap);
  if (moab::MB_SUCCESS != rval) return rval;

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::create_group_entsets(refentity_handle_map& group_map)
{
  moab::ErrorCode rval;
  const char geom_categories[][CATEGORY_TAG_SIZE] =
      {"Vertex\0", "Curve\0", "Surface\0", "Volume\0", "Group\0"};
  DLIList<RefEntity*> entitylist;

  // Create entity sets for all ref groups
  std::vector<moab::Tag> extra_name_tags;
  DLIList<CubitString> name_list;
  entitylist.clean_out();

  // Get all entity groups from the CGM model
  GeometryQueryTool::instance()->ref_entity_list("group", entitylist);
  entitylist.reset();

  // Loop over all groups
  for (int i = entitylist.size(); i--; ) {
    // Take the next group
    RefEntity* grp = entitylist.get_and_step();
    name_list.clean_out();
    // Get the names of all entities in this group from the solid model
    RefEntityName::instance()->get_refentity_name(grp, name_list);
    if (name_list.size() == 0)
      continue;
    // Set pointer to first name of the group and set the first name to name1
    name_list.reset();
    CubitString name1 = name_list.get();

    // Create entity handle for the group
    moab::EntityHandle h;
    rval = mdbImpl->create_meshset(moab::MESHSET_SET, h);
    if (moab::MB_SUCCESS != rval) return rval;

    // Set tag data for the group
    char namebuf[NAME_TAG_SIZE];
    memset(namebuf, '\0', NAME_TAG_SIZE);
    strncpy(namebuf, name1.c_str(), NAME_TAG_SIZE - 1);
    if (name1.length() >= (unsigned)NAME_TAG_SIZE)
      {
        message << "WARNING: group name '" << name1.c_str()
                << "' truncated to '" << namebuf << "'" << std::endl;
      }
    rval = mdbImpl->tag_set_data(name_tag, &h, 1, namebuf);
    if (moab::MB_SUCCESS != rval) return rval;


    int id = grp->id();
    rval = mdbImpl->tag_set_data(id_tag, &h, 1, &id);
    if (moab::MB_SUCCESS != rval)
      return moab::MB_FAILURE;

    rval = mdbImpl->tag_set_data(category_tag, &h, 1, &geom_categories[4]);
    if (moab::MB_SUCCESS != rval)
      return moab::MB_FAILURE;
    // Check for extra group names
    if (name_list.size() > 1) {
      for (int j = extra_name_tags.size(); j < name_list.size(); ++j) {
        sprintf(namebuf, "EXTRA_%s%d", NAME_TAG_NAME, j);
        moab::Tag t;
        rval = mdbImpl->tag_get_handle(namebuf, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE, t, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
        assert(!rval);
        extra_name_tags.push_back(t);
      }
      // Add extra group names to the group handle
      for (int j = 0; j < name_list.size(); ++j) {
        name1 = name_list.get_and_step();
        memset(namebuf, '\0', NAME_TAG_SIZE);
        strncpy(namebuf, name1.c_str(), NAME_TAG_SIZE - 1);
        if (name1.length() >= (unsigned)NAME_TAG_SIZE)
          {
            message << "WARNING: group name '" << name1.c_str()
                    << "' truncated to '" << namebuf << "'" << std::endl;
          }
        rval = mdbImpl->tag_set_data(extra_name_tags[j], &h, 1, namebuf);
        if (moab::MB_SUCCESS != rval) return rval;
      }
    }
    // Add the group handle
    group_map[grp] = h;
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::store_group_content(refentity_handle_map (&entitymap)[5])
{
  moab::ErrorCode rval;
  DLIList<RefEntity*> entlist;
  refentity_handle_map_itor ci;
  // Store contents for each group
  entlist.reset();
  for (ci = entitymap[4].begin(); ci != entitymap[4].end(); ++ci) {
    RefGroup* grp = (RefGroup*)(ci->first);
    entlist.clean_out();
    grp->get_child_ref_entities(entlist);

    moab::Range entities;
    while (entlist.size()) {
      RefEntity* ent = entlist.pop();
      int dim = ent->dimension();

      if (dim < 0) {
        Body* body;
        if (entitymap[4].find(ent) != entitymap[4].end()) {
          // Child is another group; examine its contents
          entities.insert(entitymap[4][ent]);
        }
        else if ((body = dynamic_cast<Body*>(ent)) != NULL) {
          // Child is a CGM Body, which presumably comprises some volumes--
          // extract volumes as if they belonged to group.
          DLIList<RefVolume*> vols;
          body->ref_volumes(vols);
          for (int vi = vols.size(); vi--; ) {
            RefVolume* vol = vols.get_and_step();
            if (entitymap[3].find(vol) != entitymap[3].end()) {
              entities.insert(entitymap[3][vol]);
            } else {
              message << "Warning: CGM Body has orphan RefVolume" << std::endl;
            }
          }
        } else {
          // Otherwise, warn user.
          message << "Warning: A dim<0 entity is being ignored by ReadCGM." << std::endl;
        }
      }
      else if (dim < 4) {
        if (entitymap[dim].find(ent) != entitymap[dim].end())
          entities.insert(entitymap[dim][ent]);
      }
    }

    if (!entities.empty()) {
      rval = mdbImpl->add_entities(ci->second, entities);
      if (moab::MB_SUCCESS != rval) return rval;
    }
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::create_vertices(refentity_handle_map &vertex_map)
{
  moab::ErrorCode rval;
  refentity_handle_map_itor ci;

  for (ci = vertex_map.begin(); ci != vertex_map.end(); ++ci) {
    CubitVector pos = dynamic_cast<RefVertex*>(ci->first)->coordinates();
    double coords[3] = {pos.x(), pos.y(), pos.z()};
    moab::EntityHandle vh;

    rval = mdbImpl->create_vertex(coords, vh);
    if (moab::MB_SUCCESS != rval) return rval;

    // Add the vertex to its tagged meshset
    rval = mdbImpl->add_entities(ci->second, &vh, 1);
    if (moab::MB_SUCCESS != rval) return rval;

    // point map entry at vertex handle instead of meshset handle to
    // simplify future operations
    ci->second = vh;
  }

  return moab::MB_SUCCESS;
}


moab::ErrorCode DAGMCExportCommand::create_curve_facets(refentity_handle_map& curve_map,
                                       refentity_handle_map& vertex_map)
{
  moab::ErrorCode rval;
  CubitStatus s;

  // Maximum allowable curve-endpoint proximity warnings
  // If this integer becomes negative, then abs(curve_warnings) is the
  // number of warnings that were suppressed.
  int curve_warnings = 0;
  failed_curve_count = 0;

  // Map iterator
  refentity_handle_map_itor ci;
  
  // Create geometry for all curves
  GMem data;
  for (ci = curve_map.begin(); ci != curve_map.end(); ++ci) {
    // Get the start and end points of the curve in the form of a reference edge
    RefEdge* edge = dynamic_cast<RefEdge*>(ci->first);
    // Get the edge's curve information
    Curve* curve = edge->get_curve_ptr();
    // Clean out previous curve information
    data.clear();
    // Facet curve according to parameters and CGM version
    s = edge->get_graphics(data, norm_tol, faceting_tol);
    
    if( s != CUBIT_SUCCESS )
      {
        // if we fatal on curves
        if(fatal_on_curves)
          {  
             message << "Failed to facet the curve " << edge->id() << std::endl;
             return moab::MB_FAILURE;
           }
        // otherwise record them
        else
          {
            failed_curve_count++;
            failed_curves.push_back(edge->id());
          }
        continue;
      }
    
    std::vector<CubitVector> points = data.point_list();
    
    // Need to reverse data?
    if (curve->bridge_sense() == CUBIT_REVERSED) 
      std::reverse(points.begin(), points.end());
    
    // Check for closed curve
    RefVertex *start_vtx, *end_vtx;
    start_vtx = edge->start_vertex();
    end_vtx = edge->end_vertex();
    
    // Special case for point curve
    if (points.size() < 2) {
      if (start_vtx != end_vtx || curve->measure() > GEOMETRY_RESABS) {
        message << "Warning: No facetting for curve " << edge->id() << std::endl;
        continue;
      }
      moab::EntityHandle h = vertex_map[start_vtx];
      rval = mdbImpl->add_entities(ci->second, &h, 1);
      if (moab::MB_SUCCESS != rval)
        return moab::MB_FAILURE;
      continue;
    }
    // Check to see if the first and last interior vertices are considered to be
    // coincident by CUBIT
    const bool closed = (points.front() - points.back()).length() < GEOMETRY_RESABS;
    if (closed != (start_vtx == end_vtx)) {
      message << "Warning: topology and geometry inconsistant for possibly closed curve "
              << edge->id() << std::endl;
    }
    
    // Check proximity of vertices to end coordinates
    if ((start_vtx->coordinates() - points.front()).length() > GEOMETRY_RESABS ||
        (end_vtx->coordinates() - points.back()).length() > GEOMETRY_RESABS) {
      
      curve_warnings--;
      if (curve_warnings >= 0 || verbose_warnings) {
        message << "Warning: vertices not at ends of curve " << edge->id() << std::endl;
        if (curve_warnings == 0 && !verbose_warnings) {
          message << "         further instances of this warning will be suppressed..." << std::endl;
        }
      }
    }

    // Create interior points
    std::vector<moab::EntityHandle> verts, edges;
    verts.push_back(vertex_map[start_vtx]);
    for (size_t i = 1; i < points.size() - 1; ++i) {
      double coords[] = {points[i].x(), points[i].y(), points[i].z()};
      moab::EntityHandle h;
      // Create vertex entity
      rval = mdbImpl->create_vertex(coords, h);
      if (moab::MB_SUCCESS != rval)
        return moab::MB_FAILURE;
      verts.push_back(h);
    }
    verts.push_back(vertex_map[end_vtx]);

    // Create edges
    for (size_t i = 0; i < verts.size() - 1; ++i) {
      moab::EntityHandle h;
      rval = mdbImpl->create_element(moab::MBEDGE, &verts[i], 2, h);
      if (moab::MB_SUCCESS != rval)
        return moab::MB_FAILURE;
      edges.push_back(h);
    }

    // If closed, remove duplicate
    if (verts.front() == verts.back())
      verts.pop_back();
    // Add entities to the curve meshset from entitymap
    rval = mdbImpl->add_entities(ci->second, &verts[0], verts.size());
    if (moab::MB_SUCCESS != rval)
      return moab::MB_FAILURE;
    rval = mdbImpl->add_entities(ci->second, &edges[0], edges.size());
    if (moab::MB_SUCCESS != rval)
      return moab::MB_FAILURE;
  }

  if (!verbose_warnings && curve_warnings < 0) {
    message << "Suppressed " << -curve_warnings
            << " 'vertices not at ends of curve' warnings." << std::endl;
    //std::cerr << "To see all warnings, use reader param VERBOSE_CGM_WARNINGS." << std::endl;
  }

  return moab::MB_SUCCESS;
}


moab::ErrorCode DAGMCExportCommand::create_surface_facets(refentity_handle_map& surface_map,
                                                          refentity_handle_map& vertex_map)
{
  moab::ErrorCode rval;
  refentity_handle_map_itor ci;
  CubitStatus s;
  failed_surface_count =0;

  DLIList<TopologyEntity*> me_list;

  GMem data;
  // Create geometry for all surfaces
  for (ci = surface_map.begin(); ci != surface_map.end(); ++ci) {
    RefFace* face = dynamic_cast<RefFace*>(ci->first);

    data.clear();
    s = face->get_graphics(data, norm_tol, faceting_tol, len_tol);

    if (CUBIT_SUCCESS != s)
      return moab::MB_FAILURE;

    std::vector<CubitVector> points = data.point_list();

    // Declare array of all vertex handles
    std::vector<moab::EntityHandle> verts(points.size(), 0);

    // Get list of geometric vertices in surface
    me_list.clean_out();
    ModelQueryEngine::instance()->query_model(*face, DagType::ref_vertex_type(), me_list);

    // For each geometric vertex, find a single coincident point in facets
    // Otherwise, print a warning
    for (int i = me_list.size(); i--; ) {
      // Assign geometric vertex
      RefVertex* vtx = dynamic_cast<RefVertex*>(me_list.get_and_step());
      CubitVector pos = vtx->coordinates();

      for (int j = 0; j < points.size(); ++j) {
        // Assign facet vertex
        CubitVector vpos = points[j];

        // Check to see if they are considered coincident
        if ((pos - vpos).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS) {
          // If this facet vertex has already been found coincident, print warning
          if (verts[j]) {
            message << "Warning: Coincident vertices in surface " << face->id() << std::endl;
          }
          // If a coincidence is found, keep track of it in the verts vector
          verts[j] = vertex_map[vtx];
          break;
        }
      }
    }

    // Now create vertices for the remaining points in the facetting
    for (int i = 0; i < points.size(); ++i) {
      if (verts[i]) // If a geometric vertex
        continue;
      double coords[] = {points[i].x(), points[i].y(), points[i].z()};
      // Return vertex handle to verts to fill in all remaining facet
      // vertices
      rval = mdbImpl->create_vertex(coords, verts[i]);
      if (moab::MB_SUCCESS != rval)
        return rval;
    }

    std::vector<int> facet_list = data.facet_list();
    
    // record the failures for information
    if (facet_list.size() == 0)
      {
        failed_surface_count++;
        failed_surfaces.push_back(face->id());
      }

    // Now create facets
    moab::Range facets;
    std::vector<moab::EntityHandle> corners;
    for (int i = 0; i < facet_list.size(); i += facet_list[i] + 1) {
      // Get number of facet verts
      int num_verts = facet_list[i];
      corners.resize(num_verts);
      for (int j = 1; j <= num_verts; ++j) {
        if (facet_list[i+j] >= (int)verts.size()) {
          message << "ERROR: Invalid facet data for surface " << face->id() << std::endl;
          return moab::MB_FAILURE;
        }
        corners[j - 1] = verts[facet_list[i+j]];
      }
      moab::EntityType type;
      if (num_verts == 3)
        type = moab::MBTRI;
      else {
        message << "Warning: non-triangle facet in surface " << face->id() << std::endl;
        message << "  entity has " << num_verts << " edges" << std::endl;
        if (num_verts == 4)
          type = moab::MBQUAD;
        else
          type = moab::MBPOLYGON;
      }

      //if (surf->bridge_sense() == CUBIT_REVERSED)
        //std::reverse(corners.begin(), corners.end());

      moab::EntityHandle h;
      rval = mdbImpl->create_element(type, &corners[0], corners.size(), h);
      if (moab::MB_SUCCESS != rval)
        return moab::MB_FAILURE;

      facets.insert(h);
    }

    // Add vertices and facets to surface set
    rval = mdbImpl->add_entities(ci->second, &verts[0], verts.size());
    if (moab::MB_SUCCESS != rval)
      return moab::MB_FAILURE;
    rval = mdbImpl->add_entities(ci->second, facets);
    if (moab::MB_SUCCESS != rval)
      return moab::MB_FAILURE;
  }

  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCExportCommand::gather_ents(moab::EntityHandle gather_set)
{
  moab::ErrorCode rval;
  moab::Range new_ents;
  rval = mdbImpl->get_entities_by_handle(0,new_ents);
  CHK_MB_ERR_RET_MB("Could not get all entity handles.",rval);

  //make sure there the gather set is empty
  moab::Range gather_ents;
  rval = mdbImpl->get_entities_by_handle(gather_set,gather_ents);
  CHK_MB_ERR_RET_MB("Could not get the gather set entities.",rval);

  if( 0 != gather_ents.size() ){
    CHK_MB_ERR_RET_MB("Unknown entities found in the gather set.",rval);
  }

  rval = mdbImpl->add_entities(gather_set,new_ents);
  CHK_MB_ERR_RET_MB("Could not add newly created entities to the gather set.",rval);

  return rval;
}

