#include "iGeom_test.hpp"
#include "iGeom_funcs.hpp"
#include <vector>
#include "CubitInterface.hpp"

// CGM includes
#include "GeometryModifyTool.hpp"
#include "Body.hpp"
#include "RefEntity.hpp"
#include <iostream>

iGeom_test::iGeom_test()
{
  //default values
  radius = 2.0;
  radius2 = 1.0;
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
  iGeom_createEntSet( &set, &igm_result );
  iGeom_createSphere( radius, &datum, &igm_result );
  iGeom_createBrick( radius, radius + 1 , radius2, &move[0], &igm_result );
//  iGeom_createCylinder( 5, radius, radius2, &move[0], &igm_result );
//  iGeom_createCone( 5, radius, 0.0, radius2, &move[1], &igm_result );
//  iGeom_createTorus( radius, radius2, &move[0], &igm_result );
  iGeom_copyEnt( move[0], &move[1], &igm_result );
//  iGeom_reflectEnt( move[1], 0, 0, 0, 0, 0, 1, &igm_result );
  iGeom_rotateEnt( move[0], 90, 1, 0, 0, &igm_result ); 
  iGeom_rotateEnt( move[1], 90, 0, 1, 0, &igm_result );
  iGeom_scaleEnt( move[1], 0, 0, 0, 1, 0.5, 2, &igm_result );
  iGeom_moveEnt( move[0], 1, 2, 4, &igm_result );
  iGeom_moveEnt( move[1], 1, 2, 4, &igm_result );
  iBase_EntityHandle result;
  iGeom_uniteEnts( move, 2, &result, &igm_result );
//  iGeom_subtractEnts( move[0], move[1], &result, &igm_result );
//  iGeom_intersectEnts( move[0], move[1], &result, &igm_result );
  iGeom_sectionEnt( datum, 1, 0, 0, 0.4, 1, &datum, &igm_result );
  iGeom_addEntToSet( result, set, &igm_result );
  iGeom_createCone( 5, radius, 0.0, radius2, &move[1], &igm_result );
  iGeom_deleteEnt( move[1], &igm_result );

  //Testing how functions handle working with deleted entities
  iGeom_moveEnt( move[1], 1, 2, 4, &igm_result );
  iGeom_deleteEnt( move[1], &igm_result );
  iGeom_reflectEnt( move[1], 0, 0, 0, 0, 0, 1, &igm_result );
  iGeom_sectionEnt( move[1], 1, 0, 0, 0.4, 0, &datum, &igm_result );
  iGeom_copyEnt( move[1], &move[0], &igm_result );
  iGeom_reflectEnt( move[1], 0, 0, 0, 0, 0, 1, &igm_result );

  return true;
}


