#include "iGeom_funcs.hpp"
#include "CubitInterface.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"

void
iGeom_createSphere( double radius,
                   /*out*/ iBase_EntityHandle *geom_entity )
{
/*  if (radius <= 0.0) {
    std::ostringstream error_message;
    error_message.str("");
    error_message << "Sphere radius must be positive!" << std::endl;
    CubitInterface::get_cubit_message_handler()->print_error(error_message.str().c_str());
  }
  else{
  */
  
  RefEntity* temp_body = GeometryModifyTool::instance()->sphere( radius );
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);
  //}
}

void
iGeom_createBrick( /*in*/ double x,
                   /*in*/ double y,
                   /*in*/ double z,
                   /*out*/ iBase_EntityHandle *geom_entity )
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

  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);
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
                      /*in*/ double minor_rad,
                      /*out*/ iBase_EntityHandle *geom_entity )
 //                     int* err )
{
  double tmp_minor = (0.0 == minor_rad ? major_rad : minor_rad);
  RefEntity *temp_body = 
    GeometryModifyTool::instance()->cylinder(height, major_rad, tmp_minor, major_rad);
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);
  /*


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
                  /*in*/ double rad_top,
                  /*out*/ iBase_EntityHandle *geom_entity )
//                  int* err )
{
  double tmp_minor = (0.0 == minor_rad_base ? major_rad_base : minor_rad_base);
  RefEntity *temp_body = 
    GeometryModifyTool::instance()->cylinder(height, major_rad_base, tmp_minor, rad_top);
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);
  /*


  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
  */
}

void
iGeom_createTorus( /*in*/ double major_rad,
                   /*in*/ double minor_rad,
                   /*out*/ iBase_EntityHandle *geom_entity )//,
//                   int* err )
{
/*  if (minor_rad >= major_rad) {
    ERROR(iBase_INVALID_ARGUMENT, "Major radius must be greater than minor radius for tori.");
  }
  */
  
  RefEntity *temp_body = GeometryModifyTool::instance()->torus(major_rad, minor_rad);
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);
  /*
   
  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
  */
}

void
iGeom_moveEnt( /*inout*/ iBase_EntityHandle geom_entity,
               /*in*/ double x,
               /*in*/ double y,
               /*in*/ double z )//,
//               int* err )
{
  CubitVector vec(x, y, z);
  Body *this_bod = dynamic_cast<Body*>(reinterpret_cast<RefEntity*>(geom_entity));
  DLIList<Body*> bodies = DLIList<Body*>();
  bodies.insert( this_bod );
  CubitStatus result;
  if (NULL != this_bod) {
    result = GeometryQueryTool::instance()->translate(bodies, vec);
/*    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to move body.");
      RETURN(iBase_FAILURE);
    }
    
    RETURN(iBase_SUCCESS);
    */
  }
  
  /*
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(reinterpret_cast<RefEntity*>(geom_entity));
  if (NULL != this_bte) {
      // non-body move; check to see if there are any siblings to this entity in the
      // same body; if so, we can't move it; if not, get the body and move that; if
      // there is no body, it's a free entity and we can move it anyway
    Body *this_body = this_bte->body();
    if (NULL == this_body) {
      result = GeometryQueryTool::instance()->translate(this_bte, vec);
    /*  if (CUBIT_SUCCESS != result) {
        ERROR(iBase_FAILURE, "Failed to move entity.");
      }
    }
    else {
      int num_sibs = -1;
      switch (this_bte->dimension()) {
        case 0: num_sibs = this_body->num_ref_vertices(); break;
        case 1: num_sibs = this_body->num_ref_edges(); break;
        case 2: num_sibs = this_body->num_ref_faces(); break;
      }
      if (num_sibs == 1) {
          // ok to move the body instead
        result = GeometryQueryTool::instance()->translate(this_body, vec);
        /*
        if (CUBIT_SUCCESS != result) {
          ERROR(iBase_FAILURE, "Failed to move body even only one entity of that"
                             " dimension in the body.");
        }
      }
      else {
        ERROR(iBase_FAILURE, "Too many siblings for an entity to move it.");
      }
    }

//    RETURN(iBase_SUCCESS);
  }
        */
  
 // ERROR(iBase_INVALID_ENTITY_TYPE, "Wrong type of entity specified for move.");
}
