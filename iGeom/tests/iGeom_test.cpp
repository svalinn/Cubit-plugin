#include "iGeom_test.hpp"

#include <vector>

#include "CubitVersionCompatibility.hpp"
#include "CubitInterface.hpp"

#include "iGeom.h"

// CGM includes
#include "GeometryModifyTool.hpp"
#include "Body.hpp"
#include "RefEntity.hpp"

iGeom_test::iGeom_test()
{
  //default values
  radius = 2.0;
  radius2 = 1.0;

  CubitMessageHandler* console = CubitInterface::get_cubit_message_handler();
  if (console) {
    std::ostringstream load_message;
    load_message.str("");
    load_message << "-- iGeom_test command available." << std::endl;
    console->print_error(load_message.str().c_str());
  }
}

iGeom_test::~iGeom_test()
{}

std::vector<std::string> iGeom_test::get_syntax()
{
  // Define the syntax for the command. Note the syntax is a modified BNF
  // format. Full documentation on the command specification syntax can be
  // found in the documentation.
  std::string syntax =
      "iGeom_test"
      "[<value:label='radius',help='<radius>'> ]"
      "[<value:label='radius2',help='<radius2>'>]";

  std::vector<std::string> syntax_list;
  syntax_list.push_back(syntax);

  return syntax_list;
}

std::vector<std::string> iGeom_test::get_syntax_help()
{
  std::vector<std::string> help;
  return help;
}

std::vector<std::string> iGeom_test::get_help()
{
  std::vector<std::string> help;
  return help;
}

bool iGeom_test::execute(CubitCommandData &data)
{
  int igm_result;

  iBase_EntitySetHandle set;
  iBase_EntityHandle datum;
  iBase_EntityHandle move[2];
//  iBase_EntityHandle move[1];
  data.get_value("radius",radius);
  data.get_value("radius2",radius2);
  iGeom_createEntSet(NULL, 0, &set, &igm_result);
  iGeom_createSphere(NULL, radius, &datum, &igm_result);
//  iGeom_createBrick(NULL, radius, radius + 1 , radius2, &move[0], &igm_result);
//  iGeom_createCylinder(NULL, 5, radius, radius2, &move[0], &igm_result);
//  iGeom_createCone(NULL, 5, radius, 0.0, radius2, &move[1], &igm_result);
  iGeom_createTorus(NULL, radius, radius2, &move[0], &igm_result);
  iGeom_copyEnt(NULL, move[0], &move[1], &igm_result);
//  iGeom_reflectEnt(NULL, move[1], 0, 0, 0, 0, 0, 1, &igm_result);
  iGeom_rotateEnt(NULL, move[0], 90, 1, 0, 0, &igm_result);
  iGeom_rotateEnt(NULL, move[1], 90, 0, 1, 0, &igm_result);
  iGeom_scaleEnt(NULL, move[1], 0, 0, 0, 1, 0.5, 2, &igm_result);
  iGeom_moveEnt(NULL, move[0], 1, 2, 4, &igm_result);
  iGeom_moveEnt(NULL, move[1], 1, 2, 4, &igm_result);
  iBase_EntityHandle result;
  iGeom_uniteEnts(NULL, move, 2, &result, &igm_result);
//  iGeom_subtractEnts(NULL, move[0], move[1], &result, &igm_result);
//  iGeom_intersectEnts(NULL, move[0], move[1], &result, &igm_result);
  iGeom_sectionEnt(NULL, datum, 1, 0, 0, 0.4, 1, &datum, &igm_result);
  iGeom_addEntToSet(NULL, result, set, &igm_result);
  iGeom_createCone(NULL, 5, radius, 0.0, radius2, &move[1], &igm_result);
  iGeom_deleteEnt(NULL, move[1], &igm_result);

  //Testing how functions handle working with deleted entities
  iGeom_moveEnt(NULL, move[1], 1, 2, 4, &igm_result);
  iGeom_deleteEnt(NULL, move[1], &igm_result);
  iGeom_reflectEnt(NULL, move[1], 0, 0, 0, 0, 0, 1, &igm_result);
  iGeom_sectionEnt(NULL, move[1], 1, 0, 0, 0.4, 0, &datum, &igm_result);
  iGeom_copyEnt(NULL, move[1], &move[0], &igm_result);
  iGeom_reflectEnt(NULL, move[1], 0, 0, 0, 0, 0, 1, &igm_result);

  return true;
}
