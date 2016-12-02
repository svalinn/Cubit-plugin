#ifndef IGEOM_FUNCS_HPP
#define IGEOM_FUNCS_HPP
#include "iBase.h"

void iGeom_createSphere( double radius,
                         iBase_EntityHandle *geom_entity );
void iGeom_createBrick( double x, 
                        double y, 
                        double z, 
                        iBase_EntityHandle *geom_entity );
void iGeom_createCylinder( double height, 
                           double major_rad, 
                           double minor_rad,
                           iBase_EntityHandle *geom_entity ); 
void iGeom_createCone( double height,
                       double major_rad_base,
                       double minor_rad_base,
                       double rad_top,
                       iBase_EntityHandle *geom_entity );
                       
void iGeom_createTorus( double major_rad, 
                        double minor_rad,
                        iBase_EntityHandle *geom_entity );
void iGeom_moveEnt( iBase_EntityHandle geom_entity, 
                    double x, double y, double z );

#endif
