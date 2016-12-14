#include "iGeom_test.hpp"
#include "iGeom_funcs.hpp"
#include <vector>
#include "CubitInterface.hpp"

// CGM includes
//#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
/*
#include "GeometryModifyEngine.hpp"
#include "ModelQueryEngine.hpp"
#include "GMem.hpp"

#include "RefEntityName.hpp"

*/
#include "Body.hpp"
#include "RefEntity.hpp"
/*
#include "Surface.hpp"
#include "Curve.hpp"

#include "RefGroup.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "SenseEntity.hpp"
*/
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
      /*
     "[faceting_tolerance <value:label='faceting_tolerance',help='<faceting tolerance>'>] "
      "[length_tolerance <value:label='length_tolerance',help='<length tolerance>'>] "
      "[normal_tolerance <value:label='normal_tolerance',help='<normal tolerance>'>] "
      "[verbose] [fatal_on_curves]";
      */

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

//iGeom_getDescription doesn't currently work without CGM.
//Not sure how to test iGeom_getEntBoundBox.

  iBase_EntitySetHandle set;
  iBase_EntityHandle datum;
  iBase_EntityHandle move[2];
//  iBase_EntityHandle move[1];
  data.get_value("radius",radius);
  data.get_value("radius2",radius2);
  iGeom_createEntSet( &set );
  iGeom_createSphere( radius, &datum );
  iGeom_createBrick( radius, radius + 1 , radius2, &move[0] );
//  iGeom_createCylinder( 5, radius, radius2, &move[0] );
//  iGeom_createCone( 5, radius, 0.0, radius2, &move[1] );
//  iGeom_createTorus( radius, radius2, &move[0] );
  iGeom_copyEnt( move[0], &move[1] );
//  iGeom_reflectEnt( move[1], 0, 0, 0, 0, 0, 1 );
  iGeom_rotateEnt( move[0], 90, 1, 0, 0 ); 
  iGeom_rotateEnt( move[1], 90, 0, 1, 0 );
  iGeom_scaleEnt( move[1], 0, 0, 0, 1, 0.5, 2 );
  iGeom_moveEnt( move[0], 1, 2, 4 );
  iGeom_moveEnt( move[1], 1, 2, 4 );
  iBase_EntityHandle result;
  iGeom_uniteEnts( move, 2, &result );
//  iGeom_subtractEnts( move[0], move[1], &result );
//  iGeom_intersectEnts( move[0], move[1], &result );
  iGeom_sectionEnt( datum, 1, 0, 0, 0.4, 1, &datum );
  iGeom_addEntToSet( result, set );
  iGeom_createCone( 5, radius, 0.0, radius2, &move[1] );
  iGeom_deleteEnt( move[1] );

  //Testing how functions handle working with deleted entities
  iGeom_moveEnt( move[1], 1, 2, 4 );
  iGeom_deleteEnt( move[1] );
  iGeom_reflectEnt( move[1], 0, 0, 0, 0, 0, 1 );
  iGeom_sectionEnt( move[1], 1, 0, 0, 0.4, 0, &datum );
  iGeom_copyEnt( move[1], &move[0] );
  iGeom_reflectEnt( move[1], 0, 0, 0, 0, 0, 1 );

  return true;
}


