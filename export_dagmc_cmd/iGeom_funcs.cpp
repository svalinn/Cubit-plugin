#include "CubitEntity.hpp"
#include "iGeom_funcs.hpp"
#include "CubitInterface.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "RefEntityFactory.hpp"
#include "RefGroup.hpp"
#include <cstring>

//setData
#include "RefEntityName.hpp"

//#include "GeometryModifyTool.hpp"
#include "RefEntity.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include <iostream>

const char *iGeom_entity_type_names[] = {"vertex", "curve", "surface", "body"};

class iGeomArrayManager
{
  void** arrayPtr;

public:


  iGeomArrayManager( void** array_ptr,
                     int& array_allocated_space,
                     int& array_size,
                     int count,
                     int val_size/*,
                     int* err */) : arrayPtr(0)
  {
    if (!*array_ptr) {
      *array_ptr = malloc(val_size * count);
      array_allocated_space = array_size = count;
      /*
      if (!*array_ptr) {
        ERROR(iBase_MEMORY_ALLOCATION_FAILED, "Couldn't allocate array.");
      }
      */
      arrayPtr = array_ptr;
    }
    else {
      array_size = count;
      /*
      if (array_allocated_space < count) {
        ERROR(iBase_BAD_ARRAY_DIMENSION, 
          "Allocated array not large enough to hold returned contents.");
      }
      */
    }
//    RETURN(iBase_SUCCESS);
    return;
  }
  
  ~iGeomArrayManager() 
  {
    if (arrayPtr) {
      free(*arrayPtr);
      *arrayPtr = 0;
    }
  }
  
  void keep_array()
    { arrayPtr = 0; }
};










void
iGeom_addEntToSet( /*in*/ iBase_EntityHandle entity_to_add,
                   /*inout*/ iBase_EntitySetHandle entity_set_handle )//,
//                   int* err )
{
  if (NULL == entity_to_add) {
    //TODO Make error message
    //RETURN(iBase_INVALID_ARGUMENT);
    return;
  }
  
  CubitStatus status = reinterpret_cast<RefGroup*>(entity_set_handle)->
    add_ref_entity(const_cast<RefEntity*>(reinterpret_cast<RefEntity*>(entity_to_add)));
  
  /*
  if (CUBIT_SUCCESS != status) {
    ERROR(iBase_FAILURE, "Problem adding entity to another set.");
  }

  RETURN(iBase_SUCCESS);
  */
}

void
iGeom_copyEnt( /*in*/ iBase_EntityHandle geom_entity,
               /*out*/ iBase_EntityHandle *geom_entity2 )//,
//               int* err )
{
  RefEntity *this_ent = reinterpret_cast<RefEntity*>(geom_entity);
//  if (NULL == dynamic_cast<Body*>(this_ent)) return;
  Body *this_body = dynamic_cast<Body*>(this_ent);
  RefVolume *this_vol = dynamic_cast<RefVolume*>(this_ent);
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
iGeom_createEntSet( /*out*/ iBase_EntitySetHandle *entity_set )//,
//                    int* err )
{
  RefGroup* grp = RefEntityFactory::instance()->construct_RefGroup();
  *entity_set = reinterpret_cast<iBase_EntitySetHandle>(grp);
    // need to set a tag denoting multiset or not...
  /*
  if (*entity_set == NULL) {
    RETURN(iBase_FAILURE);
  }
  
  else {
    RETURN(iBase_SUCCESS);
  }
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

//XXX This uses CGM, so would need a new way to check errors.  But we aren't using the iBase errors at the moment, so not high priority.
void
iGeom_getDescription( char* description_buffer,
                      int description_buffer_length )
{
//    CGM_iGeom_getLastErrorDesc(description_buffer, description_buffer_length);

    std::string lastErrorDesc = "No Error";
    if (description_buffer && description_buffer_length > 0) {
        lastErrorDesc.copy( description_buffer, description_buffer_length );
        if (lastErrorDesc.length() < (unsigned)description_buffer_length)
            description_buffer[lastErrorDesc.length()] = '\0';
    }
}

void
iGeom_getEntBoundBox( /*in*/ iBase_EntityHandle entity_handle,
                      /*out*/ double* min_x,
                      /*out*/ double* min_y,
                      /*out*/ double* min_z,
                      /*out*/ double* max_x,
                      /*out*/ double* max_y,
                      /*out*/ double* max_z )//,
//                      int* err )
{
  RefEntity* entity = (RefEntity*)entity_handle;
  CubitVector minc, maxc;
  //CubitStatus status = iGeom_bounding_box( entity, minc, maxc );
  iGeom_bounding_box( entity, minc, maxc );
  minc.get_xyz( *min_x, *min_y, *min_z );
  maxc.get_xyz( *max_x, *max_y, *max_z );
//  RETURN( (status == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}

void
iGeom_getEntities( /*in*/ iBase_EntitySetHandle set_handle,
                   /*in*/ int gentity_type,
                   /*out*/ iBase_EntityHandle **gentity_handles,
                   int *gentity_handles_allocated,
                   int *gentity_handles_size )//,
//                   int* err )
{
  if (RefGroup *this_set = reinterpret_cast<RefGroup*>(set_handle)) {
    static DLIList<CubitEntity*> centities;
    centities.clean_out();
    this_set->get_child_entities(centities);
    copy_ibase_type( gentity_type, centities, 
                     gentity_handles,
                     gentity_handles_allocated,
                     gentity_handles_size );//,
//                     err );
  }
  else {
    static DLIList<RefEntity*> dim_entities;
    dim_entities.clean_out();
    append_all_ibase_type( gentity_type, dim_entities/*, err*/ );
//    if (iBase_SUCCESS != *err)
//      return;
    
    iGeomArrayManager gentity_handles_manager ( reinterpret_cast<void**>(gentity_handles), *(gentity_handles_allocated), *(gentity_handles_size), dim_entities.size(), sizeof(**gentity_handles) );
    gentity_handles_manager.keep_array();
    dim_entities.copy_to((RefEntity**)*gentity_handles);
//    RETURN(iBase_SUCCESS);
  }  
}

void
iGeom_getNumOfType( /*in*/ iBase_EntitySetHandle set_handle,
                    /*in*/ int gentity_type,
                    int* count )//,
//                    int* err )
{
  const RefGroup *this_set = reinterpret_cast<RefGroup*>(set_handle);
  if (0 == this_set) {
    switch (gentity_type) {
      case iBase_ALL_TYPES:
        *count  = GeometryQueryTool::instance()->num_bodies();
        *count += GeometryQueryTool::instance()->num_ref_faces();
        *count += GeometryQueryTool::instance()->num_ref_edges();
        *count += GeometryQueryTool::instance()->num_ref_vertices();
        break;
      case iBase_REGION:
        *count = GeometryQueryTool::instance()->num_bodies();
        break;
      case iBase_FACE:
        *count = GeometryQueryTool::instance()->num_ref_faces();
        break;
      case iBase_EDGE:
        *count = GeometryQueryTool::instance()->num_ref_edges();
        break;
      case iBase_VERTEX:
        *count = GeometryQueryTool::instance()->num_ref_vertices();
        break;
      default:
//        RETURN(iBase_BAD_TYPE_AND_TOPO);
        break;
    }
//    RETURN (iBase_SUCCESS);
  }
  else {
    static DLIList<CubitEntity*> centities;
    centities.clean_out();
    const_cast<RefGroup*>(this_set)->get_child_entities(centities);
    *count = count_ibase_type( gentity_type, centities/*, err*/ );
  }
}

void
iGeom_getRootSet( //iGeom_Instance,
                  iBase_EntitySetHandle* root )//,
//                  int* err )
{
  *root = NULL;
//  RETURN(iBase_SUCCESS);
}

void
iGeom_getTagHandle( iGeom_Instance instance,
                    /*in*/ const char *tag_name,
                    iBase_TagHandle* tag_handle,
//                    int* err,
                    int tag_name_len )
{
    // make sure string is null-terminated
  std::string tag_name_buf( tag_name, tag_name_len );
  tag_name = tag_name_buf.c_str();
  *tag_handle = reinterpret_cast<iBase_TagHandle>(static_cast<size_t>(reinterpret_cast<CGMTagManager*>(instance)->getTagHandle( tag_name )));

    // XXX: this seems really wrong...
//  iGeom_getErrorType(err);
}

void
iGeom_getTagSizeBytes( iGeom_Instance instance,
                       /*in*/ iBase_TagHandle tag_handle,
                       int* tag_size )//,
//                       int* err )
{
  *tag_size = reinterpret_cast<CGMTagManager*>(instance)->getTagSize(reinterpret_cast<long>(tag_handle));
//  RETURN(iBase_SUCCESS);
}


void
iGeom_imprintEnts( /*in*/ iBase_EntityHandle const* gentity_handles,
                   int gentity_handles_size )//,
//                   int* err )
{
  if (gentity_handles_size < 1) // GMT::imprint segfaults if passed an empty list
//    RETURN(iBase_SUCCESS);
    return;

  DLIList<Body*> bods;
  DLIList<RefVolume*> vols, temp_vols;
  RefEntity* const* handle_array = reinterpret_cast<RefEntity* const*>(gentity_handles);
//  CubitStatus status = CUBIT_SUCCESS;
  bool status = true;
  for (int i = 0; i < gentity_handles_size; i++) {
    Body *temp_bod = dynamic_cast<Body*>(handle_array[i]);
    if (NULL != temp_bod) {
      bods.append_unique(temp_bod);
      continue;
    }
    
    RefVolume *temp_vol = dynamic_cast<RefVolume*>(handle_array[i]);
    if (NULL != temp_vol) {
      TopologyEntity *topo_ent = dynamic_cast<TopologyEntity*>(handle_array[i]);
      if (NULL == topo_ent) {
 //       status = CUBIT_FAILURE;
        status = false;
        continue;
      }
      temp_bod = topo_ent->body();
      if (NULL == temp_bod) {
//        status = CUBIT_FAILURE;
        status = false;
        continue;
      }
      bods.append_unique(temp_bod);
      continue;
    }

      // if we've gotten here, it's an error
 //   status = CUBIT_FAILURE;
    status = false;
  }
  
  if (true != status) return; //RETURN(iBase_FAILURE);

  DLIList<Body*> temp_bods;
  status = GeometryModifyTool::instance()->imprint(bods, temp_bods, false);
  
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
iGeom_newGeom( const char* options,
               iGeom_Instance* instance_out,
//               int* err,
               const int options_size ) 
{
    // scan options for default engine option
  std::string engine;
  if (options && options_size) {
    std::string tmp(options, options_size);
    char f[] = ";engine="; f[0]=tmp[0]; // correct delimiter
    size_t p = tmp.find( f );
    //I'm pretty sure we'll never use engine if we aren't using CGMA.
    //"engine" only occurs in mcnp2cad within ifdef USING_CGMA blocks,
    //and it would require carrying over much of INITCGMA.
    /*
    if (p != std::string::npos) { // if we found engine option
      p += strlen(f); // advance to value (past '=')
      size_t e = tmp.find( tmp[0], p ); // find end delimiter
      if (e == std::string::npos) // if no end delim, then must be last option
        engine = tmp.substr( p, std::string::npos );
      else
        engine = tmp.substr( p, e-p );
    }
    */
  }
  
    // initialize static var with result so that call happens only once
//  static const CubitStatus status = init_cgm( engine );
    // but check the result for every call
  /*
  if (CUBIT_SUCCESS != status)
    RETURN (iBase_FAILURE);
    */

    // return the tagmanager as the instance
#ifdef ITAPS_SHIM
//  CGMTagManager::instance().vtable = &CGM_iGeom_vtable;
#endif
  *instance_out = reinterpret_cast<iGeom_Instance>(&CGMTagManager::instance());
// RETURN(iBase_SUCCESS);
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
iGeom_sectionEnt( /*inout*/ iBase_EntityHandle geom_entity,
                  /*in*/ double plane_normal_x,
                  /*in*/ double plane_normal_y,
                  /*in*/ double plane_normal_z,
                  /*in*/ double offset,
                  /*in*/ int reverse,
                  /*out*/ iBase_EntityHandle *geom_entity2 )//,
//                  int* err )
{
  RefEntity* this_ent = reinterpret_cast<RefEntity*>(geom_entity);

  //For some reason the dynamic_cast can be in the if statement, but if it gets assigned a
  //pointer, Trelis crashes.
  if(NULL == dynamic_cast<Body*>(this_ent)) { 
    return;
  }
  Body *this_body = dynamic_cast<Body*>(this_ent);
/*
  if (NULL == this_body) {
    RefVolume *this_vol = dynamic_cast<RefVolume*>(reinterpret_cast<RefEntity*>(geom_entity));
    if (NULL != this_vol)
      this_body = this_vol->get_body_ptr();
  }
  */
  if (NULL == this_body) {
//    ERROR(iBase_INVALID_ARGUMENT, "Can only section bodies.");
      //TODO add error message
      return;
  }

  CubitVector normal(plane_normal_x, plane_normal_y, plane_normal_z);
  if (normal.length_squared() == 0.0) {
 //   ERROR(iBase_INVALID_ARGUMENT, "Zero-length vector input.");
    //TODO add error message
    return;
  }
  
  CubitVector point1 = normal * CubitVector(1.0, 0.0, 0.0);
  if (point1.length_squared() == 0.0)
    point1 = normal * CubitVector(0.0, 1.0, 0.0);
    
  CubitVector point2 = normal * point1;
  CubitVector point3(0.0, 0.0, 0.0);

  if (0.0 != offset) {
    normal.normalize();
    normal *= offset;
    point1 += normal;
    point2 += normal;
    point3 += normal;
  }

  
  DLIList<Body*> blank_list, new_body_list;
  //Since we have no way to check error, may as well just use the original.
  blank_list.append(this_body);
  //blank_list.append(GeometryModifyTool::instance()->copy_body(this_body));
  CubitStatus result = GeometryModifyTool::instance()->section(blank_list, point1, point2, point3,
                                    new_body_list, !reverse,
                                    false);
  /*
  if (CUBIT_SUCCESS != result || 0 == new_body_list.size()) {
    GeometryQueryTool::instance()->delete_RefEntity(blank_list.get());
    ERROR(iBase_FAILURE, "Section failed.");
  }
    
  else {
  */
      // need to assign it to a RE* first so the void cast gets done right
    RefEntity *new_body = new_body_list.get();
    *geom_entity2 = reinterpret_cast<iBase_EntityHandle>(new_body);
      // also, delete the original body, now that the section worked
//    GeometryQueryTool::instance()->delete_RefEntity(this_body);
//  }

//  RETURN(iBase_SUCCESS);
}


void
iGeom_setData( iGeom_Instance instance,
               /*in*/ iBase_EntityHandle entity_handle,
               /*in*/ iBase_TagHandle tag_handle,
               /*in*/ const void *tag_value_tmp )//,
//               int* err )
{
  const char *tag_value = reinterpret_cast<const char *>(tag_value_tmp);
  RefEntity *tmp_entity = reinterpret_cast<RefEntity*>(entity_handle);
  iBase_ErrorType retval = reinterpret_cast<CGMTagManager*>(instance)->setArrData(&tmp_entity, 1, 
                                          reinterpret_cast<long>(tag_handle), 
                                          tag_value);
//  RETURN(retval);
}


void
iGeom_setEntSetData( iGeom_Instance instance,
                     /*in*/ iBase_EntitySetHandle entity_set,
                     /*in*/ iBase_TagHandle tag_handle,
                     /*in*/ const void *tag_value_tmp )//,
//                     int* err )
{
  const char *tag_value = reinterpret_cast<const char *>(tag_value_tmp);
    // have to go through RefEntity* so that RefEntity** gets set right
  RefEntity *tmp_entity = reinterpret_cast<RefGroup*>(entity_set);
  iBase_ErrorType retval = reinterpret_cast<CGMTagManager*>(instance)->setArrData(&tmp_entity, 1, reinterpret_cast<long>(tag_handle), 
                                          tag_value);
//  RETURN(retval);
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

/*****************
*Helper Functions*
*****************/
//static CubitStatus
static void
iGeom_bounding_box( RefEntity* entity, CubitVector& minc, CubitVector& maxc )
{
  CubitBox box;
  if (BasicTopologyEntity* bte = dynamic_cast<BasicTopologyEntity*>(entity))
    box = bte->bounding_box();
  else if(Body* body = dynamic_cast<Body*>(entity))
    box = body->bounding_box();
  else {
//    CGM_iGeom_setLastError(iBase_INVALID_ENTITY_HANDLE, "Entities passed into gentityBoundingBox must be vertex, edge, face, or region."); 
//    return CUBIT_FAILURE;
    return;
  }
  
  minc = box.minimum();
  maxc = box.maximum();
//  return CUBIT_SUCCESS;
}

static 
void copy_ibase_type( int ibase_type, 
                      const DLIList<CubitEntity*>& list,
                      iBase_EntityHandle** entity_handles,
                      int* entity_handles_alloc,
                      int* entity_handles_size )//,
//                      int* err )
{
  int count;
  if (*entity_handles_alloc == 0) {
    count = count_ibase_type( ibase_type, list/*, err */);
    if (count < 0)
      return;
    *entity_handles = (iBase_EntityHandle*)malloc( count * sizeof(iBase_EntityHandle) );
    if (!*entity_handles) 
//      RETURN(iBase_MEMORY_ALLOCATION_FAILED);
      return;
    *entity_handles_alloc = count;
  }
  
  switch (ibase_type) {
    case iBase_ALL_TYPES:
      count = append_not_type<RefGroup>(list,*entity_handles, *entity_handles_alloc);
      break;
    case iBase_REGION:
      count = append_type<Body>(list,*entity_handles, *entity_handles_alloc);
      break;
    case iBase_FACE:
      count = append_type<RefFace>(list,*entity_handles, *entity_handles_alloc);
      break;
    case iBase_EDGE:
      count = append_type<RefEdge>(list,*entity_handles, *entity_handles_alloc);
      break;
    case iBase_VERTEX:
      count = append_type<RefVertex>(list,*entity_handles, *entity_handles_alloc);
      break;
    default:
   //   RETURN(iBase_INVALID_ENTITY_TYPE);
      return;
      break;
  }
  
  *entity_handles_size = count;
  /*
  if (count > *entity_handles_alloc)
    RETURN(iBase_BAD_ARRAY_DIMENSION);
    */
  
//  RETURN(iBase_SUCCESS);
}

static 
void append_all_ibase_type( int ibase_type, 
                            DLIList<RefEntity*>& target_list )//,
//                            int* err )
{
  RefEntityFactory *const ref = RefEntityFactory::instance();
  if (ibase_type == iBase_ALL_TYPES) {
    for (int i = 0; i < 4; ++i) {
      DLIList<RefEntity*> tmp;
      ref->ref_entity_list( iGeom_entity_type_names[i], tmp );
      target_list += tmp;
    }
  }
  else if (abs(ibase_type) < iBase_ALL_TYPES) {
    ref->ref_entity_list( iGeom_entity_type_names[ibase_type], target_list );
  }
  /*
  else {
    RETURN(iBase_INVALID_ENTITY_TYPE);
  }
  
  RETURN(iBase_SUCCESS);
  */
}

static inline
int count_ibase_type( int ibase_type, const DLIList<CubitEntity*>& list/*, int* err */)
{
//  *err = iBase_SUCCESS;
  switch (ibase_type) {
    case iBase_ALL_TYPES: return list.size() - count_type<RefGroup>(list);
    case iBase_REGION:    return count_type<Body>(list);
    case iBase_FACE:      return count_type<RefFace>(list);
    case iBase_EDGE:      return count_type<RefEdge>(list);
    case iBase_VERTEX:    return count_type<RefVertex>(list);
    default:
//      *err = iBase_INVALID_ENTITY_TYPE;
//      CGM_iGeom_setLastError( *err );
      return -1;
  }
}

template <typename SKIP_TYPE> static inline 
int append_not_type( const DLIList<CubitEntity*>& source_list,
                     iBase_EntityHandle* array, int array_size )
{
  int len = source_list.size();
  int count = 0;
  for (int i = 0; i < len; ++i) {
    if (!dynamic_cast<SKIP_TYPE*>(source_list[i])) {
      if (count == array_size) 
        return -1;
      else if (RefEntity* ent = dynamic_cast<RefEntity*>(source_list[i]))
        array[count++] = reinterpret_cast<iBase_EntityHandle>(ent);
    }
  }
  return count;
}

// Returns count of entities appended.
template <typename TARGET_TYPE> static inline 
int append_type( const DLIList<CubitEntity*>& source_list,
                 iBase_EntityHandle* array, int array_size )
{
  RefEntity* re_ptr;
  int len = source_list.size();
  int count = 0;
  for (int i = 0; i < len; ++i) {
    if (TARGET_TYPE* ent = dynamic_cast<TARGET_TYPE*>(source_list[i])) {
      if (count < array_size)
        array[count] = reinterpret_cast<iBase_EntityHandle>(re_ptr = ent);
      ++count;
    }
  }
  return count;
}

template <typename T> static inline 
int count_type( const DLIList<CubitEntity*>& list )
{
  int count = 0, size = list.size();
  for (int i = 0; i < size; ++i)
    if (dynamic_cast<T*>(list[i]))
      ++count;
  return count;
}

long CGMTagManager::getTagHandle (/*in*/ const char *tag_name)
{
  std::map<std::string,long>::iterator it =
    tagNameMap.find(std::string(tag_name));
  if (it != tagNameMap.end()) {
    bool active = (it->second > 0 ? tagInfo[it->second] :
                   presetTagInfo[-it->second]).isActive;
    if (active) {
//      CGM_iGeom_clearLastError();
      return it->second;
    }
  }

//  CGM_iGeom_setLastError( iBase_TAG_NOT_FOUND );
  return 0;
}

static CGMTagManager::TagInfo preset_tag_list[] = {
   // tag size      tag name           tag data type  default active
 { 0,              "",                 iBase_BYTES,   NULL,  false },
 { 32,             "NAME",             iBase_BYTES,   NULL,   true },
 { sizeof(int),    "GLOBAL_ID",        iBase_INTEGER, NULL,   true },
 { sizeof(int),    "UNIQUE_ID",        iBase_INTEGER, NULL,   true },
 { sizeof(int),    "MESH_INTERVAL",    iBase_INTEGER, NULL,   true },
 { sizeof(double), "MESH_SIZE",        iBase_DOUBLE,  NULL,   true },
 { 4,              "SIZE_FIRMNESS",    iBase_BYTES,   NULL,   true } };
 

CGMTagManager::TagInfo* const CGMTagManager::presetTagInfo = preset_tag_list;

int CGMTagManager::getTagSize (/*in*/ const long tag_handle)
{
//  CGM_iGeom_clearLastError();
  return (tag_handle > 0 ? 
          tagInfo[tag_handle].tagLength : 
          presetTagInfo[-tag_handle].tagLength);
}

iBase_ErrorType CGMTagManager::setArrData (/*in*/ RefEntity* const* entity_handles, 
                                                  const int entity_handles_size,
                                           /*in*/ const long tag_handle,
                                           /*in*/ const char *tag_values/*, const int tag_values_size*/)
{
  TagInfo *tinfo = (tag_handle > 0 ? &tagInfo[tag_handle] : &presetTagInfo[-tag_handle]);
  int tag_size = tinfo->tagLength;
  
  const char *val_ptr = tag_values;

/*
  iBase_ErrorType result = iBase_SUCCESS, tmp_result;
  
  if (tag_handle < 0) {
    for (int i = 0; i < entity_handles_size; i++) {
      if (NULL == entity_handles[i])
        tmp_result = setPresetTagData(interface_group(), tag_handle, val_ptr, tag_size);
      else
        tmp_result = setPresetTagData(entity_handles[i], tag_handle, val_ptr, tag_size);

      val_ptr += tag_size;
//      if (iBase_SUCCESS != tmp_result) result = tmp_result;
//    }
//    RETURN(result);
  }
*/

  for (int i = 0; i < entity_handles_size; i++) {
    RefEntity *this_ent = (NULL == entity_handles[i] ? interface_group() : 
                           entity_handles[i]);
    CATag *catag = get_catag(this_ent, true);
    assert(NULL != catag);
    catag->set_tag_data(tag_handle, val_ptr);
    val_ptr += tag_size;
  }
}


RefGroup *CGMTagManager::interface_group(const bool create_if_missing) 
{
  if (NULL == interfaceGroup) 
    interfaceGroup = 
      dynamic_cast<RefGroup*>(RefEntityName::instance()->get_refentity("interface_group"));
  
  if (NULL == interfaceGroup && create_if_missing)
    interfaceGroup = RefEntityFactory::instance()->construct_RefGroup("interface_group");

  return interfaceGroup;
}

CATag *CGMTagManager::get_catag(RefEntity *ent, 
                                  const bool create_if_missing) 
{
  CubitAttrib *this_attrib = ent->get_cubit_attrib(CATag_att_type, create_if_missing);
  if (NULL != this_attrib)
    return dynamic_cast<CATag*>(this_attrib);
  else
    return NULL;
}
  
iBase_ErrorType CATag::set_tag_data(long tag_handle, const void *tag_data, 
                                     const bool can_shallow_copy)
{
  CGMTagManager::TagInfo *tinfo = (tag_handle > 0 ? 
                                   &(myManager->tagInfo[tag_handle]) : 
                                   &(myManager->presetTagInfo[-tag_handle]));
  
    // check if this attribute has this tag
  std::map<int, void*>::iterator tdpos = tagData.find(tag_handle);
  if (tdpos == tagData.end())
    tdpos = tagData.insert(tagData.end(),
                           std::pair<int,void*>(tag_handle, NULL));
    
  if (!can_shallow_copy) {
      // need to copy the data; might need to allocate first
    if ((*tdpos).second == NULL)
      (*tdpos).second = malloc(tinfo->tagLength);

    memcpy((*tdpos).second, tag_data, tinfo->tagLength);
  }
  else {
      // should shallow copy; might have to delete what's there already
    if ((*tdpos).second != NULL) free((*tdpos).second);
  
      // if shallow copying, caller is saying we can copy, so cast away const
    (*tdpos).second = const_cast<void*>(tag_data);
  }

 // RETURN(iBase_SUCCESS);
}
