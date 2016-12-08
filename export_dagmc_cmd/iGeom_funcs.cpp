#include "iGeom_funcs.hpp"
#include "CubitInterface.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "RefEntityFactory.hpp"

#include <iostream>

void
iGeom_copyEnt( /*in*/ iBase_EntityHandle geom_entity,
               /*out*/ iBase_EntityHandle *geom_entity2 )//,
//               int* err )
{
  Body *this_body = dynamic_cast<Body*>(reinterpret_cast<RefEntity*>(geom_entity));
  RefVolume *this_vol = dynamic_cast<RefVolume*>(reinterpret_cast<RefEntity*>(geom_entity));
  if (NULL != this_vol || NULL != this_body) {
      // need to get the associated body, since cgm only supports copying bodies,
      // not volumes
    if (NULL == this_body) {
      this_body = this_vol->get_body_ptr();
      if (NULL == this_body) {
//TODO        ERROR(iBase_FAILURE, "Can't get body from volume.");
        return;
      }
    }

    RefEntity *temp_entity = GeometryModifyTool::instance()->copy_body(this_body);
    *geom_entity2 = reinterpret_cast<iBase_EntityHandle>(temp_entity);
  }
  else {
    RefEntity *this_ent = reinterpret_cast<RefEntity*>(geom_entity);
    RefEntity *temp_entity = GeometryModifyTool::instance()->copy_refentity(this_ent);
    *geom_entity2 = reinterpret_cast<iBase_EntityHandle>(temp_entity);
  }

/*
  if (NULL == *geom_entity2) {
    ERROR(iBase_FAILURE, "NULL returned from CGM copy.");
  }
  
  RETURN(iBase_SUCCESS);
  */
}

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
iGeom_deleteEnt( /*in*/ iBase_EntityHandle geom_entity )//,
//                 int* err )
{
  RefEntity *this_ent = reinterpret_cast<RefEntity*>(geom_entity);

    // special case: if this is a volume, delete the body instead
  CubitStatus result;
  RefVolume *this_vol = dynamic_cast<RefVolume*>(this_ent);
  if (NULL != this_vol){
    result = GeometryQueryTool::instance()->delete_Body(this_vol->body());}
  else if(NULL != this_ent && NULL != dynamic_cast<Body*>(this_ent)) { 
    result = GeometryQueryTool::instance()->delete_RefEntity(this_ent);}

/*
  if (CUBIT_FAILURE == result) {
    ERROR(iBase_FAILURE, "Problems deleting entity.");
  }
  */

    // check to see if this was last thing deleted; if so, reset ids
  RefEntityFactory *rfi = RefEntityFactory::instance();
  if (rfi->num_bodies() == 0 && 
      rfi->num_ref_volumes() == 0 &&
      rfi->num_ref_faces() == 0 &&
      rfi->num_ref_edges() == 0 &&
      rfi->num_ref_vertices() == 0)
    rfi->reset_ids();

//  RETURN(iBase_SUCCESS);
}

void
iGeom_intersectEnts( /*in*/ iBase_EntityHandle ent1,
                     /*in*/ iBase_EntityHandle ent2,
                     /*out*/ iBase_EntityHandle *geom_entity )//,
//                     int* err )
{
  Body *this_ent1 = dynamic_cast<Body*>(reinterpret_cast<RefEntity*>(ent1));
  if( NULL == this_ent1 ){
    return;
  }
  Body *ent1_copy = GeometryModifyTool::instance()->copy_body(this_ent1);
//  if (NULL == ent1_copy) {
//    ERROR(iBase_FAILURE, "Trouble copying blank.");
//  }
  Body *this_ent2 = dynamic_cast<Body*>(reinterpret_cast<RefEntity*>(ent2));
  if( NULL == this_ent2 ){
    return;
  }
  Body *ent2_copy = GeometryModifyTool::instance()->copy_body(this_ent2);
/*  if (NULL == ent2_copy) {
    ERROR(iBase_FAILURE, "Trouble copying tool.");
    GeometryQueryTool::instance()->delete_RefEntity(ent1_copy);
    RETURN(iBase_FAILURE);
  }
  */

  DLIList<Body*> ent1_list, new_body_list;
  ent1_list.append(ent1_copy);
  
  RefEntity *new_body = NULL;
  CubitStatus result = GeometryModifyTool::instance()->intersect(ent2_copy, ent1_list, new_body_list);
/*  if (CUBIT_SUCCESS != result || 0 == new_body_list.size()) {
    ERROR(iBase_FAILURE, "Intersect failed.");
  }
  else {
  */
    new_body = new_body_list.get();
    *geom_entity = reinterpret_cast<iBase_EntityHandle>(new_body);
    GeometryQueryTool::instance()->delete_RefEntity(this_ent2);
    GeometryQueryTool::instance()->delete_RefEntity(this_ent1);
//  }

//  RETURN(iBase_SUCCESS);
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

void
iGeom_reflectEnt( /*inout*/ iBase_EntityHandle geom_entity,
                  /*in*/ double point_x,
                  /*in*/ double point_y,
                  /*in*/ double point_z,
                  /*in*/ double plane_normal_x,
                  /*in*/ double plane_normal_y,
                  /*in*/ double plane_normal_z )//,
  //                int* err )
{
  CubitVector this_plane(plane_normal_x, plane_normal_y, plane_normal_z);
  CubitVector point(point_x, point_y, point_z);
  Body *this_bod = dynamic_cast<Body*>(reinterpret_cast<RefEntity*>(geom_entity));
  DLIList<Body*> bods;
  bods.append(this_bod);
  CubitStatus result;
  if (NULL != this_bod) {
    result = GeometryQueryTool::instance()->reflect(bods, point , this_plane);
    /*
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to reflect body.");
    }
    
    RETURN(iBase_SUCCESS);
    */
    return;
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(reinterpret_cast<RefEntity*>(geom_entity));
  DLIList<BasicTopologyEntity*> btes;
  btes.append(this_bte);
  if (NULL != this_bte) {
    result = GeometryQueryTool::instance()->reflect(btes, point, this_plane);
    /*
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to reflect entity.");
    }
    
    RETURN(iBase_SUCCESS);
    */
    return;
  }
  
  //TODO Output error message here.
//  ERROR(iBase_INVALID_ENTITY_TYPE, "Wrong type of entity specified for reflect.");
}

void
iGeom_rotateEnt( /*inout*/ iBase_EntityHandle geom_entity,
                 /*in*/ double angle,
                 /*in*/ double axis_normal_x,
                 /*in*/ double axis_normal_y,
                 /*in*/ double axis_normal_z )//,
//                 int* err )
{
  CubitVector this_axis(axis_normal_x, axis_normal_y, axis_normal_z);
  Body *this_bod = dynamic_cast<Body*>(reinterpret_cast<RefEntity*>(geom_entity));
  DLIList<Body*> bodies = DLIList<Body*>();
  bodies.insert( this_bod );
  CubitStatus result;
  if (NULL != this_bod) {
    result = GeometryQueryTool::instance()->rotate(bodies, this_axis, angle);
    /*
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to rotate body.");
    }
    
    RETURN(iBase_SUCCESS);
    */
  }
  
  /*
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(reinterpret_cast<RefEntity*>(geom_entity));
  if (NULL != this_bte) {
    result = GeometryQueryTool::instance()->rotate(this_bte, this_axis, angle);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to rotate entity.");
    }
    
    RETURN(iBase_SUCCESS);
  }
  
  ERROR(iBase_INVALID_ENTITY_TYPE, "Wrong type of entity specified for move.");
  */
}

void
iGeom_scaleEnt( /*inout*/ iBase_EntityHandle geom_entity,
                /*in*/ double point_x,
                /*in*/ double point_y,
                /*in*/ double point_z,
                /*in*/ double scale_x,
                /*in*/ double scale_y,
                /*in*/ double scale_z )//,
//                int* err )
{
  CubitVector factor(scale_x, scale_y, scale_z);
  CubitVector point(point_x, point_y, point_z);
  RefEntity *this_ent = reinterpret_cast<RefEntity*>(geom_entity);
  DLIList<RefEntity*> ents = DLIList<RefEntity*>();
  ents.insert( this_ent );
  CubitStatus result;
  if (NULL != this_ent) {
    GeometryQueryTool::instance()->scale(ents, point, scale_x, scale_y, scale_z, true, ents);
    /*
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to scale body.");
    }
    
    RETURN(iBase_SUCCESS);
    */
  }
  
  //I don't think we need to worry about siblings in mcnp2cad


  /*
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(reinterpret_cast<RefEntity*>(geom_entity));
    // non-body move; check to see if there are any siblings to this entity in the
    // same body; if so, we can't move it; if not, get the body and move that; if
    // there is no body, it's a free entity and we can move it anyway
  Body *this_body = this_bte->body();
  if (NULL == this_body && NULL != this_bte) {
    result = GeometryQueryTool::instance()->scale(this_bte, point, factor);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to scale entity.");
    }
  }
  else if (NULL != this_body && NULL != this_bte) {
    int num_sibs = -1;
    switch (this_bte->dimension()) {
      case 0: num_sibs = this_body->num_ref_vertices(); break;
      case 1: num_sibs = this_body->num_ref_edges(); break;
      case 2: num_sibs = this_body->num_ref_faces(); break;
          // for TSTT, always scale volumes
      case 3: num_sibs = 1; break;
    }
    if (num_sibs == 1) {
        // ok to scale the body instead
        */
//    result = GeometryQueryTool::instance()->scale(this_body, point, factor);
/*      if (CUBIT_SUCCESS != result) {
        ERROR(iBase_FAILURE, "Failed to scale body even only one entity of that"
                           " dimension in the body.");
      }
    }
    else {
      ERROR(iBase_FAILURE, "Too many siblings for an entity to scale it.");
    }
  }

  RETURN(iBase_SUCCESS);
  */
}

void
iGeom_subtractEnts( /*in*/ iBase_EntityHandle blank,
                    /*in*/ iBase_EntityHandle tool,
                    /*out*/ iBase_EntityHandle *geom_entity )//,
//                    int* err )
{
  Body *this_blank = dynamic_cast<Body*>(reinterpret_cast<RefEntity*>(blank));
  if( NULL == blank ){
    return;
  }
  Body *blank_copy = GeometryModifyTool::instance()->copy_body(this_blank);
  /*
  if (NULL == blank_copy) {
//    ERROR(iBase_FAILURE, "Trouble copying blank.");
    return;
  }
  */
  Body *this_tool = dynamic_cast<Body*>(reinterpret_cast<RefEntity*>(tool));
  if (NULL == this_tool) {
    return;
  }
  Body *tool_copy = GeometryModifyTool::instance()->copy_body(this_tool);
  /*
  if (NULL == tool_copy) {
  //  ERROR(iBase_FAILURE, "Trouble copying tool.");
    GeometryQueryTool::instance()->delete_RefEntity(blank_copy);
  //  RETURN(iBase_FAILURE);
    return;
  }
  */

  DLIList<Body*> blank_list, new_body_list;
  blank_list.append(blank_copy);
  
  RefEntity *new_body = NULL;
  CubitStatus result = GeometryModifyTool::instance()->subtract(tool_copy, blank_list, new_body_list);
  /*
  if (CUBIT_SUCCESS != result || 0 == new_body_list.size()) {
    ERROR(iBase_FAILURE, "Subtract failed.");
  }
  else {
  */
    new_body = new_body_list.get();
    *geom_entity = reinterpret_cast<iBase_EntityHandle>(new_body);
    GeometryQueryTool::instance()->delete_RefEntity(this_blank);
    GeometryQueryTool::instance()->delete_RefEntity(this_tool);
//  }

//  RETURN(iBase_SUCCESS);
}

void
iGeom_uniteEnts( /*in*/ iBase_EntityHandle const* geom_entities,
                 int geom_entities_size,
                 /*out*/ iBase_EntityHandle *geom_entity )//,
//                 int* err )
{
  DLIList<Body*> bods, orig_bods;
  RefEntity* const* handle_array = reinterpret_cast<RefEntity* const*>(geom_entities);
  for (int i = 0; i < geom_entities_size; i++) {
    Body *this_body = dynamic_cast<Body*>(handle_array[i]);
    if (NULL != this_body) {
//      Body *new_body = GeometryModifyTool::instance()->copy_body(this_body);
  //    if (NULL != new_body) {
        bods.append(this_body);
 //       orig_bods.append(new_body);
  //    }
    }
  }
  /*
  if (bods.size() < geom_entities_size) {
    ERROR(iBase_INVALID_ARGUMENT, "Not all entities input were regions.");
    for (int i = bods.size(); i > 0; i--)
      GeometryQueryTool::instance()->delete_RefEntity(bods.get_and_step());
    
    RETURN(iBase_SUCCESS);
  }
  */
  
  DLIList<Body*> new_bods;
  CubitStatus result = GeometryModifyTool::instance()->unite(bods, new_bods, false);
  /*
  if (CUBIT_SUCCESS != result || 1 != new_bods.size()) {
    ERROR(iBase_FAILURE, "Unite failed.");
  }
  */
    
//  else {
    *geom_entity = reinterpret_cast<iBase_EntityHandle>(dynamic_cast<RefEntity*>(new_bods.get()));
//    for (int i = orig_bods.size(); i > 0; i--)
 //     GeometryQueryTool::instance()->delete_RefEntity(orig_bods.get_and_step());
//  }

// RETURN(iBase_SUCCESS);
}
