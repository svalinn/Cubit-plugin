#ifndef MCNP2CAD_HPP
#define MCNP2CAD_HPP

#include <fstream>

extern std::ofstream record;

bool import_mcnp(std::string filename);
class GeometryContext;


#endif // MCNP2CAD_HPP
