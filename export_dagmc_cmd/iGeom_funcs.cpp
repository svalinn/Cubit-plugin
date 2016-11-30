#include "iGeom_funcs.hpp"
#include "CubitInterface.hpp"
#include "GeometryModifyTool.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"

void
iGeom_createSphere( double radius )
{
  if (radius <= 0.0) {
    std::string output = "Sphere radius must be positive!\n";
  }
  else{
  
    RefEntity* tmp_body = GeometryModifyTool::instance()->sphere( radius );
  }
}

