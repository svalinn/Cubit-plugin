#ifndef IGEOM_FUNCS_HPP
#define IGEOM_FUNCS_HPP
#include "iBase.h"
#include "GeometryModifyTool.hpp"
#include "RefGroup.hpp"

class CATag;

class CGMTagManager 
{
public:
#ifdef ITAPS_SHIM
//  iGeom_vtable *vtable;
#endif

  friend class CATag;

  
//  ~CGMTagManager();
  
  struct TagInfo
  {
    int tagLength;
    std::string tagName;
    int tagType;
    char *defaultValue;
    bool isActive;
  };

//  static CubitAttrib* CATag_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

 
//  iBase_ErrorType createTag (/*in*/ const char *tag_name,
//                             /*in*/ const int tag_size,
//                             /*in*/ const int tag_type,
//                             /*in*/ char* default_value,
//                             /*out*/ long *tag_handle);

  int getTagSize (/*in*/ const long tag_handle);

  long getTagHandle (/*in*/ const char *tag_name);

  iBase_ErrorType setArrData (/*in*/ RefEntity* const* entity_handles, 
                                     const int entity_handles_size,
                              /*in*/ const long tag_handle,
                              /*in*/ const char *tag_values);

  static inline CGMTagManager& instance()
  {
    static CGMTagManager static_instance;
    return static_instance;
  }

private:
//  CGMTagManager();

  int CATag_att_type;
   
//  long pcTag;
   
  std::vector<TagInfo> tagInfo;
   
  static TagInfo* const presetTagInfo;
//  static const int numPresetTag;
   
  std::map<std::string, long> tagNameMap;
//  static const char *CATag_NAME;
//  static const char *CATag_NAME_INTERNAL;
   
  RefGroup *interfaceGroup;

  iBase_ErrorType setPresetTagData(RefEntity *entity, const long tag_num, 
                                   const char *tag_value, const int tag_size);
  
   
  CATag *get_catag(RefEntity *ent, 
                   const bool create_if_missing = false);

  RefGroup *interface_group(const bool create_if_missing = true);
  
};

class CATag: public CubitAttrib
{
private:
  friend class CGMTagManager;
 
  std::map<int, void*> tagData;

  CGMTagManager *myManager;

  CATag(CGMTagManager *manager, RefEntity *owner, CubitSimpleAttrib *csa_ptr);
    //- create a CATag from a simple attribute

public:

//  virtual ~CATag();

//  CubitStatus actuate() {return CUBIT_SUCCESS;}

//  CubitStatus update();

//  CubitStatus reset();

//  CubitSimpleAttrib* cubit_simple_attrib();
  
//  int int_attrib_type() {return myManager->CATag_att_type;}

  void add_csa_data(CubitSimpleAttrib *csa_ptr);

//  void print();
  
  iBase_ErrorType set_tag_data(long tag_num, const void *tag_data,
                               const bool can_shallow_copy = false);
};


typedef struct iGeom_Instance_Private* iGeom_Instance;

void iGeom_addEntToSet( iBase_EntityHandle entity_to_add,
                        iBase_EntitySetHandle entity_set_handle, 
                        int* err );

void iGeom_createEntSet( iBase_EntitySetHandle *entity_set,
                         int* err );

void iGeom_createSphere( double radius,
                         iBase_EntityHandle *geom_entity,
                         int* err );

void iGeom_createBrick( double x, 
                        double y, 
                        double z, 
                        iBase_EntityHandle *geom_entity,
                        int* err );

void iGeom_createCylinder( double height, 
                           double major_rad, 
                           double minor_rad,
                           iBase_EntityHandle *geom_entity,
                           int* err ); 

void iGeom_createCone( double height,
                       double major_rad_base,
                       double minor_rad_base,
                       double rad_top,
                       iBase_EntityHandle *geom_entity,
                       int* err );
                       
void iGeom_createTorus( double major_rad, 
                        double minor_rad,
                        iBase_EntityHandle *geom_entity,
                        int* err );

void iGeom_getEntBoundBox( iBase_EntityHandle entity_handle,
                           double* min_x,
                           double* min_y,
                           double* min_z,
                           double* max_x,
                           double* max_y,
                           double* max_z,
                           int* err );


void iGeom_moveEnt( iBase_EntityHandle geom_entity, 
                    double x, double y, double z,
                    int* err );

void iGeom_rotateEnt( iBase_EntityHandle geom_entity,
                      double angle,
                      double axis_normal_x,
                      double axis_normal_y,
                      double axis_normal_z,
                      int* err );

void iGeom_copyEnt( iBase_EntityHandle geom_entity,
                    iBase_EntityHandle *geom_entity2,
                    int* err );

void iGeom_scaleEnt( iBase_EntityHandle geom_entity,
                     double point_x,
                     double point_y,
                     double point_z,
                     double scale_x,
                     double scale_y,
                     double scale_z,
                     int* err );

void iGeom_uniteEnts( iBase_EntityHandle const* geom_entities,
                      int geom_entities_size,
                      iBase_EntityHandle *geom_enttiy,
                      int* err );

void iGeom_subtractEnts( iBase_EntityHandle blank,
                         iBase_EntityHandle tool,
                         iBase_EntityHandle *geom_entity,
                         int* err );

void iGeom_deleteEnt( iBase_EntityHandle geom_entity,
                      int* err );

void iGeom_intersectEnts( iBase_EntityHandle ent1,
                          iBase_EntityHandle ent2,
                          iBase_EntityHandle *geom_entity,
                          int* err );

void iGeom_reflectEnt( iBase_EntityHandle geom_entity,
                        double point_x,
                        double point_y,
                        double point_z,
                        double plane_normal_x,
                        double plane_normal_y,
                        double plane_normal_z,
                        int* err );

void iGeom_sectionEnt( iBase_EntityHandle geom_entity,
                       double plane_normal_x,
                       double plane_normal_y,
                       double plane_normal_z,
                       double offset,
                       int reverse,
                       iBase_EntityHandle *geom_entity2,
                       int* err );

void iGeom_imprintEnts( iBase_EntityHandle const* gentity_handles,
                        int gentity_handles_size,
                        int* err );


void iGeom_getDescription( char* description_buffer,
                      int description_buffer_length);//,
//                      int* err );

void iGeom_getEntities( iBase_EntitySetHandle set_handle,
                        int gentity_type,
                        iBase_EntityHandle **gentity_handles,
                        int *gentity_handles_allocated,
                        int *gentity_handles_size,
                        int* err );

void iGeom_getNumOfType( iBase_EntitySetHandle set_handle,
                         int gentity_type,
                         int count,
                         int* err );

void iGeom_getRootSet( //iGeom_Instance instance,
                       iBase_EntitySetHandle* root,
                       int* err );

void iGeom_getTagHandle( iGeom_Instance instance,
                         const char *tag_name,
                         iBase_TagHandle* tag_handle,
                         int* err,
                         int tag_name_len );

void iGeom_getTagSizeBytes( iGeom_Instance instance,
                            iBase_TagHandle tag_handle,
                            int* tag_size,
                            int* err );

void iGeom_newGeom( const  char* options,
                    iGeom_Instance* instance_out,
                    int* err,
                    const int options_size );

void iGeom_setData( iGeom_Instance instance,
                    iBase_EntityHandle entity_handle,
                    iBase_TagHandle tag_handle,
                    const void *tag_value_tmp,
                    int* err );

void iGeom_setEntSetData( iGeom_Instance instance,
                          iBase_EntitySetHandle entity_set,
                          iBase_TagHandle tag_handle,
                          const void *tag_value_tmp,
                          int* err );

void iGeom_mergeEnts( iBase_EntityHandle const* gentity_handles,
                      int gentity_handles_size,
                      double tolerance,
                      int* err );


//Helper Functions
static CubitStatus iGeom_bounding_box( RefEntity* entity,
//static void iGeom_bounding_box( RefEntity* entity,
                                       CubitVector& minc,
                                       CubitVector& maxc );

static void copy_ibase_type( int ibase_type,
                             const DLIList<CubitEntity*>& list,
                             iBase_EntityHandle** entity_handles,
                             int* entith_handles_alloc,
                             int* entity_handles_size,
                             int* err );

static void append_all_ibase_type( int ibase_type,
                                   DLIList<RefEntity*>& target_list,
                                   int* err );

static int count_ibase_type( int ibase_type,
                             const DLIList<CubitEntity*>& list,
                             int* err );

template<typename SKIP_TYPE> static
int append_not_type( const DLIList<CubitEntity*>& source_list,
                     iBase_EntityHandle* array, 
                     int array_size );

template <typename TARGET_TYPE> static
int append_type( const DLIList<CubitEntity*>& source_list,
                 iBase_EntityHandle* array,
                 int array_size );

template <typename T> static
int count_type( const DLIList<CubitEntity*>& list );


#endif
