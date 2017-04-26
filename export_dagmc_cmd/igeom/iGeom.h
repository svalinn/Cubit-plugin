#ifndef _ITAPS_iBase
#define _ITAPS_iBase

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compile time version number digits
 *
 * iBase maintains a major, minor and patch digit in its version number.
 * Technically speaking, there is not much practical value in patch digit
 * for an interface specification. A patch release is typically only used
 * for bug fix releases. Although it is rare, sometimes a bug fix
 * necessitates an API change. So, we define a patch digit for iMesh.
 *
 * Although each interface in ITAPS has been designed to support its own
 * uniqe version numbers, apart from other ITAPS interfaces, as currently
 * used, we require all ITAPS interfaces to use the same ITAPS-wide version
 * number derived from the version number defined by these three digits.
 ******************************************************************************/
#define IBASE_VERSION_MAJOR 1
#define IBASE_VERSION_MINOR 4
#define IBASE_VERSION_PATCH 1

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Version Comparison
 *
 * Evaluates to true at CPP time if the version of iBase currently being
 * compiled is greater than or equal to the version specified.
 ******************************************************************************/
#define IBASE_VERSION_GE(Maj,Min,Pat) \
  (((IBASE_VERSION_MAJOR==(Maj)) && (IBASE_VERSION_MINOR==(Min)) && (IBASE_VERSION_PATCH>=(Pat))) || \
   ((IBASE_VERSION_MAJOR==(Maj)) && (IBASE_VERSION_MINOR>(Min))) || \
   (IBASE_VERSION_MAJOR>(Maj)))

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compose compile-time string represention of the version number
 ******************************************************************************/
#define IBASE_VERSION_STRING___(I,X,Y,Z) #I "_Version_" #X "." #Y "." #Z
#define IBASE_VERSION_STRING__(I,X,Y,Z) IBASE_VERSION_STRING___(I,X,Y,Z)
#define IBASE_VERSION_STRING_(I) IBASE_VERSION_STRING__(I,IBASE_VERSION_MAJOR,IBASE_VERSION_MINOR,IBASE_VERSION_PATCH)
#define IBASE_VERSION_STRING IBASE_VERSION_STRING_(iBase)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compose compile-time symbol name derived from the version number.
 ******************************************************************************/
#define IBASE_VERSION_TAG__(I,X,Y,Z) I##_Version_##X##_##Y##_##Z
#define IBASE_VERSION_TAG_(I,X,Y,Z) IBASE_VERSION_TAG__(I,X,Y,Z)
#define IBASE_VERSION_TAG(I) IBASE_VERSION_TAG_(I,IBASE_VERSION_MAJOR,IBASE_VERSION_MINOR,IBASE_VERSION_PATCH)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief ITAPS-wide (across all ITAPS APIs) version handling
 ******************************************************************************/
#define ITAPS_VERSION_MAJOR IBASE_VERSION_MAJOR
#define ITAPS_VERSION_MINOR IBASE_VERSION_MINOR
#define ITAPS_VERSION_PATCH IBASE_VERSION_PATCH
#define ITAPS_VERSION_GE(Maj,Min,Pat) IBASE_VERSION_GE(Maj,Min,Pat)
#define ITAPS_VERSION_STRING_(I) IBASE_VERSION_STRING_(I)
#define ITAPS_VERSION_STRING ITAPS_VERSION_STRING_(ITAPS)
#define ITAPS_VERSION_TAG_(I) IBASE_VERSION_TAG(I)
#define ITAPS_VERSION_TAG ITAPS_VERSION_TAG_(I)

/***************************************************************************//**
 * \defgroup EnumIterators Enum-Iterators
 * \ingroup iBase
 * \brief Convenience macros for iterating over all possible values in an enum
 *
 * These convenience macros are provided to facilitate iterating over all
 * possible values in an enumerated type. To use these macros, for example...
 * \code 
 * for (iBase_EntityType i  = IBASE_MINENUM(iBase_EntityType);
 *                       i <= IBASE_MAXENUM(iBase_EntityType);
 *                            IBASE_INCENUM(i,iBase_EntityType))
 * {
 * }
 * \endcode
 * Be aware that some enumerated types include a <em>wild card</em> often used
 * in queries to represent all possible values and you may or may not want to
 * include such a value in your iteration.
 ******************************************************************************/

/***************************************************************************//**
 * \ingroup EnumIterators
 * @{
 ******************************************************************************/
#define IBASE_MINENUM(enumName) enumName ## _MIN
#define IBASE_MAXENUM(enumName) enumName ## _MAX
#define IBASE_NUMENUM(enumName) ((int)IBASE_MAXENUM(enumName) - (int)IBASE_MINENUM(enumName) + 1)
#define IBASE_INCENUM(enumName,I) (I = (enum enumName)((int)I+1))
/** @} */

#ifdef __cplusplus
extern "C" {
#endif

  typedef void* iBase_Instance;
  typedef struct iBase_EntityHandle_Private* iBase_EntityHandle;
  typedef struct iBase_EntitySetHandle_Private* iBase_EntitySetHandle;
  typedef struct iBase_TagHandle_Private* iBase_TagHandle;
  typedef struct iBase_EntityIterator_Private* iBase_EntityIterator;
  typedef struct iBase_EntityArrIterator_Private* iBase_EntityArrIterator;

  struct iBase_Error
  {
    int error_type;
    char description[120];
  };

  enum iBase_EntityType {
    iBase_EntityType_MIN = 0,
    /**< facilitates iteration over all values */
    iBase_VERTEX = iBase_EntityType_MIN,
    /**< A topological dimension 0 entity */
    iBase_EDGE,
    /**< A topological dimension 1 entity */
    iBase_FACE,
    /**< A topological dimension 2 entity */
    iBase_REGION,
    /**< A topological dimension 3 entity */
    iBase_ALL_TYPES,
    /**< used only in queires to request information about all types */
    iBase_EntityType_MAX = iBase_ALL_TYPES
      /**< facilitates iteration over all values */
  };

  enum iBase_AdjacencyCost {
    iBase_AdjacencyCost_MIN = 0,
    /**< facilitates iteration over all values */
    iBase_UNAVAILABLE = iBase_AdjacencyCost_MIN,
    /**< Adjacency information not supported */
    iBase_ALL_ORDER_1,
    /**< No more than local mesh traversal required (i!=j) */
    iBase_ALL_ORDER_LOGN,
    /**< Global tree search (i!=j) */
    iBase_ALL_ORDER_N,
    /**< Global exhaustive search (i!=j) */
    iBase_SOME_ORDER_1,
    /**< Only some adjacency info, local (i!=j) */
    iBase_SOME_ORDER_LOGN,
    /**< Only some adjacency info, tree (i!=j) */
    iBase_SOME_ORDER_N,
    /**< Only some adjacency info, exhaustive (i!=j) */
    iBase_AVAILABLE,
    /**< ALL (intermediate) entities available. (i==j) */
    iBase_AdjacencyCost_MAX = iBase_AVAILABLE
      /**< facilitates iteration over all values */
  };

  enum iBase_CreationStatus {
    iBase_CreationStatus_MIN = 0,
    /**< facilitates iteration over all values */
    iBase_NEW = iBase_CreationStatus_MIN,
    /**< The entity was newly created */
    iBase_ALREADY_EXISTED,
    /**< The entity already existed and the handle for that
      already existing handle was returned */
    iBase_CREATED_DUPLICATE,
    /**< The entity already existed but a new, duplicate entity was
      nevertheless created */
    iBase_CREATION_FAILED,
    /**< Creation of the entity did not succeed */
    iBase_CreationStatus_MAX = iBase_CREATION_FAILED
      /**< facilitates iteration over all values */
  };

  enum iBase_ErrorType {
    iBase_ErrorType_MIN = 0,
    /**< facilitates iteration over all values */
    iBase_SUCCESS = iBase_ErrorType_MIN,
    iBase_MESH_ALREADY_LOADED,
    iBase_FILE_NOT_FOUND,
    iBase_FILE_WRITE_ERROR,
    iBase_NIL_ARRAY,
    iBase_BAD_ARRAY_SIZE,
    iBase_BAD_ARRAY_DIMENSION,
    iBase_INVALID_ENTITY_HANDLE,
    iBase_INVALID_ENTITY_COUNT,
    iBase_INVALID_ENTITY_TYPE,
    iBase_INVALID_ENTITY_TOPOLOGY,
    iBase_BAD_TYPE_AND_TOPO,
    iBase_ENTITY_CREATION_ERROR,
    iBase_INVALID_TAG_HANDLE,
    iBase_TAG_NOT_FOUND,
    iBase_TAG_ALREADY_EXISTS,
    iBase_TAG_IN_USE,
    iBase_INVALID_ENTITYSET_HANDLE,
    iBase_INVALID_ITERATOR_HANDLE,
    iBase_INVALID_ARGUMENT,
    iBase_MEMORY_ALLOCATION_FAILED,
    iBase_NOT_SUPPORTED,
    iBase_FAILURE,
    iBase_ErrorType_MAX = iBase_FAILURE
      /**< facilitates iteration over all values */
  };

  /***************************************************************************//**
   * \details
   * Many of the functions in iMesh can return arrays of tuples; that is, arrays
   * of multi-valued type. For example, the function iMesh_getVtxArrCoords,
   * returns an array of xyz coordinate 3-tuples (or, perhaps for geometrically
   * 2D meshes, xy 2-tuples). In these situations, there are multiple ways the
   * data can be organized in memory. For example, it could be stored xyz,xyz,xyz
   * or xxx...,yyy...,zzz.... These two different storage orders are referred
   * to as INTERLEAVED and BLOCKED, respectively. For some functions in iMesh,
   * the storage order is explicitly specified as an argument to the function.
   * For other functions, the storage order is not explicitly specified. And,
   * in these cases, it shall always be implicitly assumed to be INTERLEAVED.
   * This fact will be mentioned in the documentation for each specific function
   * where it applies. For example, in case of iMesh_getEntArrAdj, the returned
   * array of adjacent entities is multi-valued in that it stores for each
   * entity queried, all its adjacent entities. Such an array will be stored
   * INTERLEAVED with all adjacent entities for the first entity in the query
   * followed by all adjacent entities for the second entity in the query and
   * so forth.
   ******************************************************************************/
  enum iBase_StorageOrder {
    iBase_StorageOrder_MIN = 0,
    /**< facilitates iteration over all values */
    iBase_BLOCKED = iBase_StorageOrder_MIN,
    /**< xxx...yyy...zzz... */
    iBase_INTERLEAVED,
    /**< xyzxyzxyz... */
    iBase_StorageOrder_MAX = iBase_INTERLEAVED
      /**< facilitates iteration over all values */
  };

  enum iBase_TagValueType {
    iBase_TagValueType_MIN = 0,
    /**< facilitates iteration over all values */
    iBase_BYTES = iBase_TagValueType_MIN,
    /**< An opaque sequence of bytes, size always measured in bytes */
    iBase_INTEGER,
    /**< A value of type \c int */
    iBase_DOUBLE,
    /**< A value of type \c double */
    iBase_ENTITY_HANDLE,
    /**< A value of type \c iBase_EntityHandle */
    iBase_ENTITY_SET_HANDLE,
    /**< A value of type \c iBase_EntitySetHandle */
    iBase_TagValueType_MAX = iBase_ENTITY_SET_HANDLE
      /**< facilitates iteration over all values */
  };

  /***************************************************************************//**
   * \page ibase iBase: ITAPS Base Interface
   ******************************************************************************/

  /***************************************************************************//**
   * \defgroup iBase iBase
   ******************************************************************************/

  /***************************************************************************//**
   * \defgroup VersionNumbers Version Numbers
   * \ingroup iBase
   ******************************************************************************/

  /***************************************************************************//**
   * \defgroup ErrorHandling Error Handling
   * \ingroup iBase
   ******************************************************************************/

  /***************************************************************************//**
   * \defgroup Datatypes Datatypes
   * \ingroup iBase
   ******************************************************************************/

  /***************************************************************************//**
   * \mainpage The ITAPS Interfaces
   *
   * \subpage ibase
   *
   * \subpage imesh
   *
   * \subpage imeshp
   * 
   * \subpage igeom
   *
   * \subpage error
   *
   * \subpage trio
   *
   * \subpage strlen
   *
   * \subpage options
   *
   * \subpage numhops
   *
   * \page error Error Handling
   *
   * With few exceptions, every iMesh function includes an output argument,
   * 'int *err', which returns an error code indicating if the function call
   * may have failed. If the value returned for the 'err' argument is NOT
   * iBase_SUCCESS, the caller should NOT attempt to interpret (read the
   * values in) any of the other return arguments of the call. While some
   * implementations may actually return valid/useful results in other
   * return arguments of a call that has failed, there is no guarentee that
   * ALL implementations will do similarly and so depending on such behavior
   * is neither portable nor safe. This is true even if the returned values
   * are different from the values of the arguments before the call was
   * made.
   *
   * \page trio Array pointer, allocated and occupied sizes argument trio
   *
   * Many of the functions in iMesh have arguments corresponding to lists of 
   * objects.  In-type arguments for lists consist of a pointer to an array and
   * a list size.  Lists returned from functions are passed in three arguments,
   * a pointer to the array representing the list, and pointers to the
   * allocated and occupied lengths of the array.  These three arguments are 
   * inout-type arguments, because they can be allocated by the application and
   * passed into the interface to hold the results of the function.  Lists
   * which are pre-allocated must be large enough to hold the results of the
   * function; if this is not the case, an error is generated.  Otherwise, the
   * occupied size is changed to the size output from the function.  If a list
   * argument is unallocated (the list pointer points to a NULL value) or if
   * the incoming value of the allocated size is zero, the list storage will be
   * allocated by the implementation.
   *
   * IN ALL CASES, MEMORY ALLOCATED BY ITAPS INTERFACE IMPLEMENTATIONS IS DONE
   * USING THE C MALLOC FUNCTION, AND MUST BE DE-ALLOCATED USING THE C FREE
   * FUNCTION.
   *
   * \page strlen String Length Arguments
   *
   * Many of the functions in iMesh involve passing a string and also the length
   * of that string. How is the null character is handled?
   * For users of the iMesh interface calling iMesh functions, it is optional
   * as to whether or not to include the null character in computing the length
   * of the string. So, for example, calling iMesh from a C program, users could
   * pass strlen(my_string) or strlen(my_string)+1 as the length of the string.
   *
   * <em>Note to implementors</em>: However, it should be noted that the situation
   * is different for implementers of the iMesh interface. In implementing an
   * iMesh interface function, there can be no assumption that the string is
   * indeed null terminated. The length argument the caller passes in may or may
   * NOT include the null character and implementations must be coded to
   * accommodate this. This requirement is primarily due to differences in how
   * Fortran and C/C++ handle passing of strings as function arguments.
  *
    * \page numhops Indirection in Set-Inclusion and Parent-Child structures
    *
    * Various functions to query entities, entity sets and parent or child sets 
    * as well as the numbers of these involve a num_hops argument. If the set
    * upon which the query is originated is the root set, the num_hops argument
    * is irrelevant and is ignored by the implementation. Otherwise, the num_hops
    * argument represents the maximum number of levels of indirection employed in
    * satisfying the query not including the originating set. For example, using
    * value for num_hops of 0 (zero) in iMesh_getEntSets will return all the 
    * entity sets that are immediately contained in a given set. Likewise, a
    * value for num_hops of 1 (one) will return all entity sets that are 
    * immediately contained in the given set plus all entity sets that
    * are contained in those immediately contained sets (e.g. one level of
        * indirection). Using a value of -1 for num_hops will return results for
    * all possible levels of indirection. In other words, using a value of
    * -1 for num_hops is equivalent to setting the maximum number of levels
    * of indirection to infinity.
    *
    * \page options Option Strings
    *
    * A few of the functions in iMesh support arbitrary options passed as a
    * character string, called an 'Option String'. The format of and handling
    * of an Option String is as follows...
    *
    * 1. Option Strings are INsensitive to case.
    *
    * 2. Each option in an Option String is pre-pended with the implementation
    * name followed by a special character called the separator character.
    *
    * 3. The separator is a colon, ':'.
    *
    * 4. Multiple options existing in a single Option String are separated by a
    * special character called the delimiter character.
    *
    * 5. The delimiter character is a space, ' '.
    *
    * 6. The effect of multiple options in a single Option String is 
    * INsensitive to order of occurrence in the string.
    *
    * 7. By default, implementations silently ignore any options that
    * do not match on the implementation name part (everything before
        * the separator character). This way, a caller may included options
    * in a single string intended for multiple different implementations.
    *
    * 8. Implementations may (or may not) warn or error for option strings
    * that match on implementation name part but are found to be in error
    * for other reasons the implementation decides.
    *
    * 9. Whenever either the separator character, ':', or delimiter character,
    * ' ', need to appear in an option, they must be escaped with the
      * backslash character, '\'.
      *
      * For example, consider the Options String
      *
      *     "grummp:silant FMDB:TwoPhaseIO moab:mpiio_hints\ foo\:bar"
      *
      * In the above example, the space serves as the delimiter character
      * between multiple options in the string. The colon serves as the
      * implementation-name/option separator character. Because options are
      * required to be insensitive to case, the caller is free to use case as a
      * word separator as in 'TwoPhaseIO' and even in the implementation name,
    * as in 'FMDB:', although 'fmdb:twophaseio' and 'fmdb:TWOPHASEIO' would
      * all have the same effect. In the moab option, both the separator
      * character and delimiter character appear in the option and so are
      * pre-pended (e.g. escaped) with the backslash character.

      * GRUMMP will silently ignore the FMDB: and moab: options because they do
      * NOT match on the implementation name part. However, GRUMMP may
      * optionally error out, or warn or silently ignore 'grummp:silant' (it was
          * supposed to be spelled 'silent') as an invalid option.
      *
      * Note that iMesh itself currently does not define any options. In order
      * to discover options a given implementation defines, users are directed
      * to the developers of the respective implementations.
      ******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif /* #ifndef _ITAPS_iBase */

#ifndef IGEOM_FUNCS_HPP
#define IGEOM_FUNCS_HPP
#include "RefGroup.hpp"

/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

/**\file CATag.hpp
 *
 *\author Tim Tautges
 *\author Jason Kraftcheck
 *
 * Original file from SNL TSTT repository was named CATSTT.
 *
 * Renamed CATag and added to ANL ITAPS repository by J.Kraftcheck,
 * 2007-6-15
 */

class CATag;

class CGMTagManager
{
public:
#ifdef ITAPS_SHIM
//  iGeom_vtable *vtable;
#endif

  friend class CATag;

//  If the destructor is there, the plugin isn't recognized.
//  ~CGMTagManager();

  struct TagInfo
  {
    int tagLength;
    std::string tagName;
    int tagType;
    char *defaultValue;
    bool isActive;
  };

  static CubitAttrib* CATag_creator(RefEntity* entity, CubitSimpleAttrib &p_csa);

  iBase_ErrorType createTag (/*in*/ const char *tag_name,
                             /*in*/ const int tag_size,
                             /*in*/ const int tag_type,
                             /*in*/ char* default_value,
                             /*out*/ long *tag_handle);

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
  CGMTagManager();

  int CATag_att_type;

  std::vector<TagInfo> tagInfo;

  static TagInfo* const presetTagInfo;
  static const int numPresetTag;

  std::map<std::string, long> tagNameMap;
  static const char *CATag_NAME;
  static const char *CATag_NAME_INTERNAL;

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

  CubitStatus actuate() {return CUBIT_SUCCESS;}

//  the destructor will make the plugin not be recognized.
//  virtual ~CATag();

  CubitStatus update();

  CubitStatus reset();

  CubitSimpleAttrib cubit_simple_attrib();

  int int_attrib_type() {return myManager->CATag_att_type;}

  void add_csa_data(CubitSimpleAttrib *csa_ptr);

//  void print();

  iBase_ErrorType set_tag_data(long tag_num, const void *tag_data,
                               const bool can_shallow_copy = false);
};

extern "C" {
  /** \mainpage The ITAPS Geometry Interface iGeom
   *
   * The ITAPS Geometry Interface iGeom provides a common interface for
   * accessing geometry and data associated with a mesh.  Applications written
   * to use this interface can use a variety of implementations, choosing
   * the one that best meets its needs.  They can also use tools written
   * to this interface.
   *
   * \section ITAPS Data Model
   *
   * The ITAPS interfaces use a data model composed of four basic data types:\n
   * \em Entity: basic topological entities in a model, e.g. vertices,
   * edges, faces, regions. \n
   * \em Entity \em Set: arbitrary grouping of other entities and sets.
   * Entity sets also support parent/child relations with other sets which
   * are distinct from entities contained in those sets.  Parent/child links
   * can be used to embed graph relationships between sets, e.g. to
   * represent topological relationships between the sets. \n
   * \em Interface: the object with which model is associated and on which
   * functions in iGeom are called. \n
   * \em Tag: application data associated with objects of any of the other
   * data types.  Each tag has a designated name, size, and data type.
   *
   * \section JTAPS Entity Type
   * Each entity has a specific Entity Type.  The Entity
   * Type is one of VERTEX, EDGE, FACE, and REGION, and is synonymous with
   * the topological dimension of the entity.  Entity Type is an enumerated
   * type in the iBase_EntityType enumeration.
   *
   * \section KTAPS Entity-, Array-, and Iterator-Based Access
   *
   * The iGeom interface provides functions for accessing entities
   * individually, as arrays of entities, or using iterators.  These access
   * methods have different memory versus execution time tradeoffs,
   * depending on the implementation.
   *
   * \section LTAPS Lists Passed Through Interface
   *
   * Many of the functions in iGeom have arguments corresponding to lists of
   * objects.  In-type arguments for lists consist of a pointer to an array and
   * a list size.  Lists returned from functions are passed in three arguments,
   * a pointer to the array representing the list, and pointers to the
   * allocated and occupied lengths of the array.  These three arguments are
   * inout-type arguments, because they can be allocated by the application and
   * passed into the interface to hold the results of the function.  Lists
   * which are pre-allocated must be large enough to hold the results of the
   * function; if this is not the case, an error is generated.  Otherwise, the
   * occupied size is changed to the size output from the function.  If a list
   * argument is unallocated (the list pointer points to a NULL value) or if
   * the incoming value of the allocated size is zero, the list storage will be
   * allocated by the implementation.  IN ALL CASES, MEMORY ALLOCATED BY ITAPS
   * INTERFACE IMPLEMENTATIONS IS DONE USING THE C MALLOC FUNCTION, AND CAN BE
   * DE-ALLOCATED USING THE C FREE FUNCTION.
   *
   */

  /**\brief  Type used to store iGeom interface handle
   *
   * Type used to store iGeom interface handle
   */

typedef struct iGeom_Instance_Private* iGeom_Instance;

    /**\brief  Add an entity to a set
     *
     * Add an entity to a set
     * \param entity_to_add The entity being added
     * \param entity_set_handle the set being added to
     * \param *err Pointer to error type returned from function
     */

void iGeom_addEntToSet(iGeom_Instance instance,
                        iBase_EntityHandle entity_to_add,
                        iBase_EntitySetHandle entity_set_handle,
                        int* err);

    /**\brief  Create an entity set
     *
     * Create an entity set, either ordered (isList=1) or unordered
     * (isList=0).  Unordered entity sets can contain a given entity or
     * set only once.
     * \param entity_set Pointer to entity set created by function
     * \param *err Pointer to error type returned from function
     */

void iGeom_createEntSet(iGeom_Instance instance,
                        int isList,
                        iBase_EntitySetHandle *entity_set,
                        int* err);

    /**\brief Create a sphere
     *
     * Create a sphere of the specified radius centered on the origin.
     * \param radius radius of the sphere
     * \param geom_entity Pointer to new entity handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_createSphere(iGeom_Instance instance,
                        double radius,
                        iBase_EntityHandle *geom_entity,
                        int* err);

    /**\brief  Create an axis-oriented box
     *
     * Create an axis-oriented box of the given dimensions, centered at the
     * origin.
     * \param x x dimension of new box
     * \param y y dimension of new box
     * \param z z dimension of new box
     * \param geom_entity Pointer to new entity handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_createBrick(iGeom_Instance instance,
                       double x,
                       double y,
                       double z,
                       iBase_EntityHandle *geom_entity,
                       int* err);

    /**\brief  Create a cylinder
     *
     * Create a cylinder parallel to the z-axis and centered at the origin (so
     * that its z-coordinate extents are +height/2 and -height/2).
     * \param height The height of the cylinder.
     * \param major_rad The x-axis radius
     * \param minor_rad The y-axis radius. If minor_rad is 0, the cylinder will
     *        be circular (as if minor_rad == major_rad).
     * \param geom_entity Pointer to new entity handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_createCylinder(iGeom_Instance instance,
                          double height,
                          double major_rad,
                          double minor_rad,
                          iBase_EntityHandle *geom_entity,
                          int* err);

  /**\brief  Create a cone or tapered cylinder
   *
   * Create a cone parallel to the z-axis and centered at the origin (so that
   * its z-coordinate extents are +height/2 and -height/2). The 'base' of the
   * cylinder is at z = -height/2, and the top is at +height/2.
   * \param height The height of the cone.
   * \param major_rad_base The x-axis radius at the base of the cylinder
   * \param minor_rad_base The y-axis radius at the base.  If minor_rad_base
   *        is 0, the cylinder will be circular (as if minor_rad_base ==
   *        major_rad_base)
   * \param rad_top The x-axis radius at the top of the cone.  The y-axis
   *        radius at the top of the cone will be inferred to keep the aspect
   *        ratio of the top of the cone the same as the bottom. If rad_top is
   *        0, the cone terminates at a point.
   * \param geom_entity Pointer to new entity handle returned from function
   * \param *err Pointer to error type returned from function
   */

void iGeom_createCone(iGeom_Instance instance,
                      double height,
                      double major_rad_base,
                      double minor_rad_base,
                      double rad_top,
                      iBase_EntityHandle *geom_entity,
                      int* err);

    /**\brief  Create a torus
     *
     * Create a torus centered on the origin and encircling the z-axis.
     * \param major_rad The distance from the origin to the center of the
     *        torus's circular cross-section.
     * \param minor_rad The radius of the cross-section.
     * \param geom_entity Pointer to new entity handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_createTorus(iGeom_Instance instance,
                       double major_rad,
                       double minor_rad,
                       iBase_EntityHandle *geom_entity,
                       int* err);

    /**\brief  Get the bounding box of the specified entity
     *
     * Get the bounding box of the specified entity
     * \param entity_handle Entity being queried
     * \param min_x Minimum coordinate of bounding box
     * \param min_y Minimum coordinate of bounding box
     * \param min_z Minimum coordinate of bounding box
     * \param max_x Maximum coordinate of bounding box
     * \param max_y Maximum coordinate of bounding box
     * \param max_z Maximum coordinate of bounding box
     * \param *err Pointer to error type returned from function
     */

void iGeom_getEntBoundBox(iGeom_Instance instance,
                          iBase_EntityHandle entity_handle,
                          double* min_x,
                          double* min_y,
                          double* min_z,
                          double* max_x,
                          double* max_y,
                          double* max_z,
                          int* err);

    /**\brief  Move an entity by the given vector
     *
     * Move an entity by translating it along the given vector.
     * \param geom_entity the entity to move
     * \param x x coordinate of the vector
     * \param y y coordinate of the vector
     * \param z z coordinate of the vector
     * \param *err Pointer to error type returned from function
     */

void iGeom_moveEnt(iGeom_Instance instance,
                   iBase_EntityHandle geom_entity,
                   double x, double y, double z,
                   int* err);

    /**\brief  Rotate an entity about an axis
     *
     * Rotate an entity by the given angle about the given axis.
     * \param geom_entity the entity to rotate
     * \param angle the rotational angle, in degrees
     * \param axis_normal_x x coordinate of the axis
     * \param axis_normal_y y coordinate of the axis
     * \param axis_normal_z z coordinate of the axis
     * \param *err Pointer to error type returned from function
     */

void iGeom_rotateEnt(iGeom_Instance instance,
                     iBase_EntityHandle geom_entity,
                     double angle,
                     double axis_normal_x,
                     double axis_normal_y,
                     double axis_normal_z,
                     int* err);

    /**\brief  Make a copy of the specified entity
     *
     * Make a copy of the specified entity
     * \param geom_entity entity to be copied
     * \param geom_entity2 Pointer to the newly-created entity
     * \param *err Pointer to error type returned from function
     */

void iGeom_copyEnt(iGeom_Instance instance,
                   iBase_EntityHandle geom_entity,
                   iBase_EntityHandle *geom_entity2,
                   int* err);

   /**\brief  Scale an entity in the x, y, and z directions
     *
     * Scale an entity in the x, y, and z directions.
     * \param geom_entity the entity to scale
     * \param point_x  x coordinate of the scaling center
     * \param point_y  y coordinate of the scaling center
     * \param point_z  z coordinate of the scaling center
     * \param scale_x factor to scale by in the x direction
     * \param scale_y factor to scale by in the y direction
     * \param scale_z factor to scale by in the z direction
     * \param *err Pointer to error type returned from function
     */

void iGeom_scaleEnt(iGeom_Instance instance,
                    iBase_EntityHandle geom_entity,
                    double point_x,
                    double point_y,
                    double point_z,
                    double scale_x,
                    double scale_y,
                    double scale_z,
                    int* err);

    /**\brief  Geometrically unite entities
     *
     * Geometrically unite the specified entities.
     * \param geom_entities pointer to Array of entity handles being united
     * \param geom_entities_size Number of entities in geom_entities array
     * \param geom_entity Pointer to new entity handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_uniteEnts(iGeom_Instance instance,
                     iBase_EntityHandle const* geom_entities,
                     int geom_entities_size,
                     iBase_EntityHandle *geom_enttiy,
                     int* err);

    /**\brief  Geometrically subtract one entity from another
     *
     * Geometrically subtract the entity tool from the entity blank.
     * \param blank The entity to subtract from
     * \param tool The entity to subtract
     * \param geom_entity Pointer to new entity handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_subtractEnts(iGeom_Instance instance,
                        iBase_EntityHandle blank,
                        iBase_EntityHandle tool,
                        iBase_EntityHandle *geom_entity,
                        int* err);

    /**\brief  Delete specified entity
     *
     * Delete specified entity
     * \param geom_entity Entity to be deleted
     * \param *err Pointer to error type returned from function
     */

void iGeom_deleteEnt(iGeom_Instance instance,
                     iBase_EntityHandle geom_entity,
                     int* err);

    /**\brief  Geometrically intersect a pair of entities
     *
     * Geometrically intersect a pair of entities.
     * \param ent1 The entity to intersect
     * \param ent2 The entity to intersect
     * \param geom_entity Pointer to new entity handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_intersectEnts(iGeom_Instance instance,
                         iBase_EntityHandle ent1,
                         iBase_EntityHandle ent2,
                         iBase_EntityHandle *geom_entity,
                         int* err);

   /**\brief  Reflect an entity across a plane
     *
     * Reflect an entity across the given plane
     * \param geom_entity the entity to reflect
     * \param point_x  x coordinate of the point that the reflecting plane goes though
     * \param point_y  y coordinate of the point that the reflecting plane goes
though
     * \param point_z  z coordinate of the point that the reflecting plane goes
though
     * \param plane_normal_x x coordinate of the plane's normal
     * \param plane_normal_y y coordinate of the plane's normal
     * \param plane_normal_z z coordinate of the plane's normal
     * \param *err Pointer to error type returned from function
     */

void iGeom_reflectEnt(iGeom_Instance instance,
                      iBase_EntityHandle geom_entity,
                      double point_x,
                      double point_y,
                      double point_z,
                      double plane_normal_x,
                      double plane_normal_y,
                      double plane_normal_z,
                      int* err);

    /**\brief  Section (cut) a region with a plane
     *
     * Section (cut) a region with a plane, retaining one of the pieces and
     * discarding the other.
     * \param geom_entity The entity to section
     * \param plane_normal_x x coordinate of the plane's normal
     * \param plane_normal_y y coordinate of the plane's normal
     * \param plane_normal_z z coordinate of the plane's normal
     * \param offset Distance of the plane from the origin
     * \param reverse Keep the piece on the normal's side (=0) or not (=1)
     * \param geom_entity2 Pointer to new entity handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_sectionEnt(iGeom_Instance instance,
                      iBase_EntityHandle geom_entity,
                      double plane_normal_x,
                      double plane_normal_y,
                      double plane_normal_z,
                      double offset,
                      int reverse,
                      iBase_EntityHandle *geom_entity2,
                      int* err);

    /**\brief  Imprint entities
     *
     * Imprint entities by merging coincident surfaces.
     * \param gentity_handles Pointer to Array of entity handles being imprinted
     * \param gentity_handles_size Number of entities in geom_entities array
     * \param *err Pointer to error type returned from function
     */

void iGeom_imprintEnts(iGeom_Instance instance,
                       iBase_EntityHandle const* gentity_handles,
                       int gentity_handles_size,
                       int* err);

    /**\brief  Get a description of the error returned from the last iGeom call
     *
     * Get a description of the error returned from the last iGeom function
     * \param description_buffer Pointer to a character string to be filled with
     *        a description of the error from the last iGeom function
     * \param description_buffer_length Length of the character string
     *        pointed to by descr
     */

void iGeom_getDescription(iGeom_Instance instance,
                          char* description_buffer,
                          int description_buffer_length);

    /**\brief  Get entities of specific type in set or instance
     *
     * Get entities of specific type in set or instance.  All entities are
     * requested by specifying iBase_ALL_TYPES.  Specified type must be a value
     * in the iBase_EntityType enumeration.
     * \param set_handle Entity set being queried
     * \param gentity_type Type of entities being requested
     * \param entity_topology Topology of entities being requested
     * \param **gentity_handles Pointer to Pointer to array of entity handles
     *        returned from function
     * \param *gentity_handles_allocated Pointer to allocated size of
     *        gentity_handles array
     * \param *gentity_handles_size Pointer to occupied size of gentity_handles
     *        array
     * \param *err Pointer to error type returned from function
     */

void iGeom_getEntities(iGeom_Instance instance,
                       iBase_EntitySetHandle set_handle,
                       int gentity_type,
                       iBase_EntityHandle **gentity_handles,
                       int *gentity_handles_allocated,
                       int *gentity_handles_size,
                       int* err);

    /**\brief  Get the number of entities with the specified type in the
     *         instance or set
     *
     * Get the number of entities with the specified type in the instance
     * or set.  If entity set handle is zero, return information for instance,
     * otherwise for set.  Value of entity type must be from the
     * iBase_EntityType enumeration.  If iBase_ALL_TYPES is specified, total
     * number of entities (excluding entity sets) is returned.
     * \param set_handle Entity set being queried
     * \param gentity_type Type of entity requested
     * \param count Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_getNumOfType(iGeom_Instance instance,
                        iBase_EntitySetHandle set_handle,
                        int gentity_type,
                        int* count,
                        int* err);

    /**\brief  Get handle of the root set for this instance
     *
     * Get handle of the root set for this instance.  All geom in
     * this instance can be accessed from this set.
     * \param instance iGeom instance handle
     * \param root Pointer to set handle returned from function
     * \param *err Pointer to error type returned from function
     */

void iGeom_getRootSet(iGeom_Instance instance,
                      iBase_EntitySetHandle* root,
                      int* err);

    /**\brief  Get a the handle of an existing tag with the specified name
     *
     * Get a the handle of an existing tag with the specified name
     * \param instance iGeom instance handle
     * \param tag_name Pointer to Name of tag being queried
     * \param tag_handle Pointer to tag handle returned from function
     * \param *err Pointer to error type returned from function
     * \param tag_name_len Length of tag name string
     */

void iGeom_getTagHandle(iGeom_Instance instance,
                        const char *tag_name,
                        iBase_TagHandle* tag_handle,
                        int* err,
                        int tag_name_len);

    /**\brief  Get a the handle of an existing tag with the specified name
     *
     * Get a the handle of an existing tag with the specified name
     * \param instance iGeom instance handle
     * \param tag_handle Pointer to tag handle returned from function
     * \param tag_size Size of tag name
     * \param *err Pointer to error type returned from function
     */

void iGeom_getTagSizeBytes(iGeom_Instance instance,
                           iBase_TagHandle tag_handle,
                           int* tag_size,
                           int* err);

    /**\brief  Construct a new iGeom instance
     *
     * Construct a new iGeom instance, using implementation-specific options
     * \param options Pointer to implementation-specific options string
     * \param instance_out Pointer to iGeom instance handle returned from function
     * \param *err Pointer to error type returned from function
     * \param options_size Length of the character string pointed to by options
     */

void iGeom_newGeom(const  char* options,
                   iGeom_Instance* instance_out,
                   int* err,
                   const int options_size);

    /**\brief  Set a tag value of arbitrary type on an entity
     *
     * Set a tag value of arbitrary type on an entity.  Tag data
     * is passed as void*. tag_value_size specifies the size of the
     * memory pointed to by tag_value in terms of bytes. Applications may
     * use this function to set data of any type, not just iBase_BYTES.
     * However, because this function supports data of arbitrary type, in
     * all cases the size specified by tag_value_size is always in terms
     * of bytes.
     *
     * \param instance iGeom instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value_tmp Pointer to tag data being set on entity
     * \param *err Pointer to error type returned from function
     */

void iGeom_setData(iGeom_Instance instance,
                   iBase_EntityHandle entity_handle,
                   iBase_TagHandle tag_handle,
                   const void *tag_value_tmp,
                   int tag_value_size,
                   int* err);

    /**\brief  Set a tag value of arbitrary type on an entity set
     *
     * Set a tag value of arbitrary type on an entity set. The tag data
     * is passed as void*. tag_value_size specifies the size of the memory
     * pointed to by tag_value in terms of bytes. Applications are free to
     * use this function to set data of any type, not just iBase_BYTES.
     * However, in all cases, the size specified by tag_value_size is
     * always in terms of bytes.
     *
     * \param instance iGeom instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value_tmp Pointer to tag data being set on entity set
     * \param *err Pointer to error type returned from function
     */

void iGeom_setEntSetData(iGeom_Instance instance,
                         iBase_EntitySetHandle entity_set,
                         iBase_TagHandle tag_handle,
                         const void *tag_value_tmp,
                         const int tag_value_size,
                         int* err);

    /**\brief  Merge ents
     *
     * Merge entities of corresponding topology/geometry within the specified
     * tolerance.
     * \param gentit_handles Pointer to Array of entity handles being imprinted
     * \param gentity_handle_size Number of entities in geom_entities array
     * \param tolerance Tolerance within which entities are considered the same
     * \param *err Pointer to error type returned from function
     */

void iGeom_mergeEnts(iGeom_Instance instance,
                     iBase_EntityHandle const* gentity_handles,
                     int gentity_handles_size,
                     double tolerance,
                     int* err);

    /**\brief  Save a geom to a file
     *
     * Save a geom to a file.  If entity set is specified, save only the
     * geom contained in that set.
     * \param instance iGeom instance handle
     * \param entity_set_handle Entity set being saved
     * \param name File name to which geom is to be saved
     * \param options Pointer to implementation-specific options string
     * \param *err Pointer to error type returned from function
     * \param name_len Length of the file name character string
     * \param options_len Length of the options character string
     */
  void iGeom_save(iGeom_Instance instance,
                  char const* name,
                  char const* options,
                  int* err,
                  int name_len,
                  int options_len);

//Helper Functions

void iGeom_getErrorType(iGeom_Instance& instance,
                        int *error_type);

} //extern C

static void tokenize(const std::string& str,
                     std::vector<std::string>& tokens);

// Expect option of the form "NAME=VALUE".
// If NAME portion matches, pass back VALUE and return true.
// Otherwise, leave 'value' unchanged and return false.
static bool match_option(const std::string& opt,
                         const char* name,
                         std::string& value);

CubitStatus CubitCompat_export_solid_model(DLIList<RefEntity*>& ref_entity_list,
                                           const char* filename,
                                           const char * filetype,
                                           int &num_ents_exported,
                                           const CubitString &cubit_version,
                                           const char* logfile_name = NULL);

static CubitStatus iGeom_bounding_box(RefEntity* entity,
                                      CubitVector& minc,
                                      CubitVector& maxc);

static void copy_ibase_type(int ibase_type,
                            const DLIList<CubitEntity*>& list,
                            iBase_EntityHandle** entity_handles,
                            int* entith_handles_alloc,
                            int* entity_handles_size,
                            int* err);

static void append_all_ibase_type(int ibase_type,
                                  DLIList<RefEntity*>& target_list,
                                  int* err);

static int count_ibase_type(int ibase_type,
                            const DLIList<CubitEntity*>& list,
                            int* err);

template<typename SKIP_TYPE> static
int append_not_type(const DLIList<CubitEntity*>& source_list,
                    iBase_EntityHandle* array,
                    int array_size);

template <typename TARGET_TYPE> static
int append_type(const DLIList<CubitEntity*>& source_list,
                iBase_EntityHandle* array,
                int array_size);

template <typename T> static
int count_type(const DLIList<CubitEntity*>& list);

#endif
