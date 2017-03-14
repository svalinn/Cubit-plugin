#ifndef I_GEOM_ERROR_H
#define I_GEOM_ERROR_H

#ifdef __cplusplus
extern "C" {
#endif

void CGM_iGeom_clearLastError();

void CGM_iGeom_setLastError( int error_type, const char* description = 0 );

int CGM_iGeom_getLastErrorType();

void CGM_iGeom_getLastErrorDesc( char* description_buffer,
                                 int description_buffer_length );



#ifdef __cplusplus
 } // extern "C"
#endif

#endif
