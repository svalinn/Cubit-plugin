//#include "iGeom.h"
#include "iGeomError.h"
#include "../iGeom_funcs.hpp"
#include <string>
#include <assert.h>

#define iBase_SUCCESS_DESC "No Error"

std::string lastErrorDesc = iBase_SUCCESS_DESC;
iBase_ErrorType lastErrorType = iBase_SUCCESS;

#ifdef __cplusplus
extern "C" {
#endif

void CGM_iGeom_clearLastError()
{
  if (lastErrorType != iBase_SUCCESS) { // don't copy string not needed
    lastErrorType = iBase_SUCCESS;
    lastErrorDesc = iBase_SUCCESS_DESC;
  }
}

void CGM_iGeom_setLastError( int error_type, const char* description )
{
    // don't do string copies for non-errors
  if (error_type == iBase_SUCCESS && lastErrorType == iBase_SUCCESS)
    return;

  lastErrorType = static_cast<iBase_ErrorType>(error_type);
  if (description) {
    lastErrorDesc = description;
    return;
  }
  
#define ERROR_DESC( A, B ) \
  case iBase_ ## A : lastErrorDesc = B ; break;

  switch (error_type) {
    ERROR_DESC( SUCCESS                 , iBase_SUCCESS_DESC );
    ERROR_DESC( MESH_ALREADY_LOADED     , "Mesh already loaded" );
    ERROR_DESC( FILE_NOT_FOUND          , "Could not read file" );
    ERROR_DESC( FILE_WRITE_ERROR        , "File write failed" );
    ERROR_DESC( NIL_ARRAY               , "NULL or empty array" );
    ERROR_DESC( BAD_ARRAY_SIZE          , "Invalid array size" );
    ERROR_DESC( BAD_ARRAY_DIMENSION     , "Invalid array dimension" );
    ERROR_DESC( INVALID_ENTITY_HANDLE   , "Invalid entity handle" );
    ERROR_DESC( INVALID_ENTITY_COUNT    , "Invalid entity count" );
    ERROR_DESC( INVALID_ENTITY_TYPE     , "Invalid entity type" );
    ERROR_DESC( INVALID_ENTITY_TOPOLOGY , "Invalid entity topology" );
    ERROR_DESC( BAD_TYPE_AND_TOPO       , "Bad type and/or topology" );
    ERROR_DESC( ENTITY_CREATION_ERROR   , "Entity creation failed" );
    ERROR_DESC( INVALID_TAG_HANDLE      , "Invalid Tag" );
    ERROR_DESC( TAG_NOT_FOUND           , "Tag does not exist" );
    ERROR_DESC( TAG_ALREADY_EXISTS      , "Tag name conflict" );
    ERROR_DESC( TAG_IN_USE              , "Tag name conflict" );
    ERROR_DESC( INVALID_ENTITYSET_HANDLE, "Invalid entity set handle" );
    ERROR_DESC( INVALID_ITERATOR_HANDLE , "Invalid iterator handle" ); 
    ERROR_DESC( INVALID_ARGUMENT        , "Invalid argument" );
    ERROR_DESC( MEMORY_ALLOCATION_FAILED, "Out of memory" );
    ERROR_DESC( NOT_SUPPORTED           , "Feature not supported" );
    ERROR_DESC( FAILURE                 , "Unknown failure or internal error" );
    default:
      assert(false);
      lastErrorDesc = "IVALID OR UNKNOWN ERROR CODE";
  }
}

int CGM_iGeom_getLastErrorType()
{
  return lastErrorType;
}

void CGM_iGeom_getLastErrorDesc(char* description_buffer,
                                int description_buffer_length )
{
  if (description_buffer && description_buffer_length > 0) {
    lastErrorDesc.copy( description_buffer, description_buffer_length );
    if (lastErrorDesc.length() < (unsigned)description_buffer_length)
      description_buffer[lastErrorDesc.length()] = '\0';
  }
}
#ifdef __cplusplus
 } // extern "C"
#endif
