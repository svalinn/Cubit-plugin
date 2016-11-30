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

iGeom_test::iGeom_test()
{
  //default values
  radius = 2.0;
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
      "[<value:label='radius',help='<radius>'> ]";
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

  data.get_value("radius",radius);
  std::stringstream ss;
  ss << radius;
//  std::string output = "Sphere of radius " + ss.str().c_str() + "coming up!\n"
  iGeom_createSphere( radius );
  return true;
}


