#include "iGeom_funcs.hpp"
#include "CubitInterface.hpp"
#include "GeometryModifyTool.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"

void
iGeom_createSphere( double radius )
{
/*  if (radius <= 0.0) {
    std::ostringstream error_message;
    error_message.str("");
    error_message << "Sphere radius must be positive!" << std::endl;
    CubitInterface::get_cubit_message_handler()->print_error(error_message.str().c_str());
  }
  else{
  */
  
    RefEntity* tmp_body = GeometryModifyTool::instance()->sphere( radius );
  //}
}

void
iGeom_createBrick( /*in*/ double x,
                   /*in*/ double y,
                   /*in*/ double z )//,
//                   /*out*/ iBase_EntityHandle *geom_entity,
//                   int* err )
{
  double tmp_x = x;
  double tmp_y = y;
  double tmp_z = z;
  
  if (0.0 == y && 0.0 == z) {
    tmp_y = x;
    tmp_z = x;
  }
  
  if (0.0 >= tmp_x || 0.0 >= tmp_y || 0.0 >= tmp_z) {
    std::ostringstream error_message;
    error_message.str("");
    error_message << "Dimensions must be >= 0, or y & z must both be zero." << std::endl;
    CubitInterface::get_cubit_message_handler()->print_error(error_message.str().c_str());
  }
  else{
    
  RefEntity *temp_body = GeometryModifyTool::instance()->brick(tmp_x, tmp_y, tmp_z);

  }
/*  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);

  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
  */
}

void
iGeom_createCylinder( /*in*/ double height,
                      /*in*/ double major_rad,
                      /*in*/ double minor_rad )//,
//                      /*out*/ iBase_EntityHandle *geom_entity,
 //                     int* err )
{
  double tmp_minor = (0.0 == minor_rad ? major_rad : minor_rad);
  RefEntity *temp_body = 
    GeometryModifyTool::instance()->cylinder(height, major_rad, tmp_minor, major_rad);
/*  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);


  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
  */
}

void
iGeom_createCone( /*in*/ double height,
                  /*in*/ double major_rad_base,
                  /*in*/ double minor_rad_base,
                  /*in*/ double rad_top )//,
//                  /*out*/ iBase_EntityHandle *geom_entity,
//                  int* err )
{
  double tmp_minor = (0.0 == minor_rad_base ? major_rad_base : minor_rad_base);
  RefEntity *temp_body = 
    GeometryModifyTool::instance()->cylinder(height, major_rad_base, tmp_minor, rad_top);
/*  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);


  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
  */
}

void
iGeom_createTorus( /*in*/ double major_rad,
                   /*in*/ double minor_rad )//,
//                   /*out*/ iBase_EntityHandle *geom_entity,
//                   int* err )
{
/*  if (minor_rad >= major_rad) {
    ERROR(iBase_INVALID_ARGUMENT, "Major radius must be greater than minor radius for tori.");
  }
  */
  
  RefEntity *temp_body = GeometryModifyTool::instance()->torus(major_rad, minor_rad);
/*  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);
   
  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
  */
}
