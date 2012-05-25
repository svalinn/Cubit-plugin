#ifndef READDAMSEL_HPP
#define READDAMSEL_HPP

#include "moab/Forward.hpp"
#include "DebugOutput.hpp"

#include "damsel.h"

namespace moab {

class DamselUtil 
{
public:
  //! needs to be a constructor to initialize dtom_data_type
  DamselUtil();

  static damsel_data_type mtod_data_type[MB_MAX_DATA_TYPE];

  std::map<damsel_data_type,DataType> dtom_data_type;

  static enum damsel_entity_type mtod_entity_type[MBMAXTYPE];

  static enum EntityType dtom_entity_type[DAMSEL_ENTITY_ALL_TYPES];
};

damsel_entity_type DamselUtil::mtod_entity_type[] = {
    DAMSEL_ENTITY_TYPE_VERTEX,      //  MBVERTEX
    DAMSEL_ENTITY_TYPE_EDGE,  // MBEDGE
    DAMSEL_ENTITY_TYPE_TRI, // MBTRI
    DAMSEL_ENTITY_TYPE_QUAD, //   MBQUAD
    DAMSEL_ENTITY_TYPE_POLYGON, //   MBPOLYGON
    DAMSEL_ENTITY_TYPE_TET,//  MBTET
    DAMSEL_ENTITY_TYPE_PYRAMID,  //   MBPYRAMID
    DAMSEL_ENTITY_TYPE_PRISM,  //   MBPRISM
    DAMSEL_ENTITY_TYPE_UNDEFINED, // MBKNIFE
    DAMSEL_ENTITY_TYPE_HEX,  //   MBHEX,
    DAMSEL_ENTITY_TYPE_POLYHEDRON, //   MBPOLYHEDRON
    DAMSEL_ENTITY_TYPE_UNDEFINED,   //   MBENTITYSET
    DAMSEL_ENTITY_TYPE_ALL_TYPES // MBMAXTYPE  /**<
};

EntityType DamselUtil::dtom_entity_type[] = {
    MBVERTEX,      //  MBVERTEX
    MBEDGE,  // MBEDGE
    MBTRI, // MBTRI
    MBQUAD, //   MBQUAD
    MBPOLYGON, //   MBPOLYGON
    MBTET,//  MBTET
    MBPRISM,  //   MBPRISM
    MBPYRAMID,  //   MBPYRAMID
    MBHEX,  //   MBHEX,
    MBPOLYHEDRON, //   MBPOLYHEDRON
    MBMAXTYPE,   //   MBENTITYSET
    MBMAXTYPE // MBMAXTYPE
};

inline DamselUtil::DamselUtil() {
  mtod_data_type[] = {
    DAMSEL_DATA_TYPE_BYTES, // MB_TYPE_OPAQUE
    DAMSEL_DATA_TYPE_INTEGER, // MB_TYPE_INTEGER
    DAMSEL_DATA_TYPE_DOUBLE, // MB_TYPE_DOUBLE
    DAMSEL_DATA_TYPE_INVALID, // MB_TYPE_BIT
    DAMSEL_DATA_TYPE_HANDLE // MB_TYPE_HANDLE
  };

  mtod_data_type[MB_TYPE_OPAQUE] = DAMSEL_DATA_TYPE_BYTES;
  mtod_data_type[MB_TYPE_INTEGER] = DAMSEL_DATA_TYPE_INTEGER;
  mtod_data_type[MB_TYPE_DOUBLE] = DAMSEL_DATA_TYPE_DOUBLE;
  mtod_data_type[MB_TYPE_BIT] = DAMSEL_DATA_TYPE_INVALID;
  mtod_data_type[MB_TYPE_HANDLE] = DAMSEL_DATA_TYPE_HANDLE;
}

}

#endif
