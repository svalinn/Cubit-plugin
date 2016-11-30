#ifndef IGEOM_FUNCS_HPP
#define IGEOM_FUNCS_HPP

void iGeom_createSphere( double radius );
void iGeom_createBrick( double x, double y, double z );
void iGeom_createCylinder( double height, double major_rad, double minor_rad); 
void iGeom_createCone( double height,
                       double major_rad_base,
                       double minor_rad_base,
                       double rad_top );
                       
void iGeom_createTorus( double major_rad, double minor_rad );

#endif
