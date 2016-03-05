#include "DAGMCExportCommand.hpp"
#include "CubitInterface.hpp"
#include "CubitMessageHandler.hpp"

// CGM includes
#include "GeometryQueryTool.hpp"

DAGMCExportCommand::DAGMCExportCommand()
{}

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

  ErrorCode rval;

  // Create entity sets for all geometric entities
  std::map<RefEntity*, EntityHandle> entmap[5];

  rval = create_entity_sets(entmap);
  // if (MB_SUCCESS != rval ) // what should error handling look like?


  return result;
}


ErrorCode DAGMCExportCommand::create_entity_sets(std::map<RefEntity*, EntityHandle> (&entmap)[5])
{
  ErrorCode rval;
  const char geom_categories[][CATEGORY_TAG_SIZE] =
              {"Vertex\0", "Curve\0", "Surface\0", "Volume\0", "Group\0"};
  const char* const names[] = {"Vertex", "Curve", "Surface", "Volume"};
  DLIList<RefEntity*> entlist;

  for (int dim = 0; dim < 4; dim++) {
    entlist.clean_out();
    GeometryQueryTool::instance()->ref_entity_list(names[dim], entlist, true);
    entlist.reset();

    for (int i = entlist.size(); i--; ) {
      RefEntity* ent = entlist.get_and_step();
      EntityHandle handle;
      // Create the new meshset
      rval = mdbImpl->create_meshset(dim == 1 ? MESHSET_ORDERED : MESHSET_SET, handle);
      if (MB_SUCCESS != rval)
        return rval;

      // Map the geom reference entity to the corresponding moab meshset
      entmap[dim][ent] = handle;

      // Create tags for the new meshset
      rval = mdbImpl->tag_set_data(geom_tag, &handle, 1, &dim);
      if (MB_SUCCESS != rval)
        return rval;

      int id = ent->id();
      rval = mdbImpl->tag_set_data(id_tag, &handle, 1, &id);
      if (MB_SUCCESS != rval)
        return rval;

      rval = mdbImpl->tag_set_data(category_tag, &handle, 1, &geom_categories[dim]);
      if (MB_SUCCESS != rval)
        return rval;
    }
  }

  return MB_SUCCESS;
}
