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

void iGeom_rotateEnt( iBase_EntityHandle geom_entity,
                      double angle,
                      double axis_normal_x,
                      double axis_normal_y,
                      double axis_normal_z );

void iGeom_copyEnt( iBase_EntityHandle geom_entity,
                    iBase_EntityHandle *geom_entity2 );

void iGeom_scaleEnt( iBase_EntityHandle geom_entity,
                     double point_x,
                     double point_y,
                     double point_z,
                     double scale_x,
                     double scale_y,
                     double scale_z );

void iGeom_uniteEnts( iBase_EntityHandle const* geom_entities,
                      int geom_entities_size,
                      iBase_EntityHandle *geom_enttiy );

void iGeom_subtractEnts( iBase_EntityHandle blank,
                         iBase_EntityHandle tool,
                         iBase_EntityHandle *geom_entity );

void iGeom_deleteEnt( iBase_EntityHandle geom_entity );

void iGeom_intersectEnts( iBase_EntityHandle ent1,
                          iBase_EntityHandle ent2,
                          iBase_EntityHandle *geom_entity );

void iGeom_reflectEnt( iBase_EntityHandle geom_entity,
                        double point_x,
                        double point_y,
                        double point_z,
                        double plane_normal_x,
                        double plane_normal_y,
                        double plane_normal_z );
#endif
