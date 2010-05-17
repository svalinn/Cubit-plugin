#include <stdlib.h>	// For exit()
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <algorithm>

#include "moab/CN.hpp"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"
#include "MBTagConventions.hpp"
#include "Internals.hpp"
#include "moab/ReadUtilIface.hpp"
#include "FileOptions.hpp"
#include "ReadCCMIO.hpp"
#include "moab/MeshTopoUtil.hpp"

#include "ccmio.h"

/*
 * CCMIO file structure
 *
 * Root
 *   State(kCCMIOState)
 *     Processor*
 *       Vertices
 *         ->ReadVerticesx, ReadMap
 *       Topology
 *         Boundary faces*(kCCMIOBoundaryFaces)
 *            ->ReadFaces, ReadFaceCells, ReadMap
 *         Internal faces(kCCMIOInternalFaces)
 *         Cells (kCCMIOCells)
 *            ->ReadCells (mapID), ReadMap, ReadCells (cellTypes)
 *       Solution
 *         Phase
 *           Field
 *             FieldData
 *   Problem(kCCMIOProblemDescription)
 *     CellType* (kCCMIOCellType)
 *       Index (GetEntityIndex), MaterialId(ReadOpti), MaterialType(ReadOptstr),
 *         PorosityId(ReadOpti), SpinId(ReadOpti), GroupId(ReadOpti)
 *
 * MaterialType (CCMIOReadOptstr in readexample)
 * constants (see readexample)
 * lagrangian data (CCMIOReadLagrangianData)
 * vertices label (CCMIOEntityDescription)
 * restart info: char solver[], iteratoins, time, char timeUnits[], angle
 *      (CCMIOReadRestartInfo, kCCMIORestartData), reference data?
 * phase:
 *   field: char name[], dims, CCMIODataType datatype, char units[]
 *       dims = kCCMIOScalar (CCMIOReadFieldDataf), 
 *              kCCMIOVector (CCMIOReadMultiDimensionalFieldData),
 *              kCCMIOTensor
 * MonitoringSets: num, name (CellSet, VertexSet, BoundarySet, BlockSet, SplineSet, CoupleSet)
 *      CCMIOGetProstarSet, CCMIOReadOpt1i,
 */

enum DataType { kScalar, kVector, kVertex, kCell, kInternalFace, kBoundaryFace,
                kBoundaryData, kBoundaryFaceData, kCellType };

namespace moab 
{
    
static int const kNValues = 10;	// Number of values of each element to print
static char const kDefaultState[] = "default";
static char const kUnitsName[] = "Units";
static int const kVertOffset = 2;
static int const kCellInc = 4;

#define CHKERR(a, b)                                 \
    {if (MB_SUCCESS != a) {if (b) readMeshIface->report_error(b); return a;}}

#define CHKCCMERR(a, b)                                 \
    {if (kCCMIONoErr != a) {if (b) readMeshIface->report_error(b); return MB_FAILURE;}}
  

ReaderIface* ReadCCMIO::factory( Interface* iface )
{ return new ReadCCMIO( iface ); }

ReadCCMIO::ReadCCMIO(Interface* impl)
    : mbImpl(impl)
{
  assert(impl != NULL);
  
  void* ptr = 0;
  impl->query_interface( "ReadUtilIface", &ptr );
  readMeshIface = reinterpret_cast<ReadUtilIface*>(ptr);

  // initialize in case tag_get_handle fails below
  mMaterialSetTag  = 0;
  mDirichletSetTag = 0;
  mNeumannSetTag   = 0;
  mHasMidNodesTag  = 0;
  mGlobalIdTag     = 0;

  //! get and cache predefined tag handles
  int dum_val = 0;
  ErrorCode result = impl->tag_get_handle(MATERIAL_SET_TAG_NAME,  mMaterialSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(MATERIAL_SET_TAG_NAME, 
                              sizeof(int), 
                              MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER,
                              mMaterialSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(DIRICHLET_SET_TAG_NAME, mDirichletSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(DIRICHLET_SET_TAG_NAME, 
                              sizeof(int), 
                              MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER,
                              mDirichletSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(NEUMANN_SET_TAG_NAME,   mNeumannSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(NEUMANN_SET_TAG_NAME, 
                              sizeof(int), 
                              MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER,
                              mNeumannSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(HAS_MID_NODES_TAG_NAME, mHasMidNodesTag);
  if (MB_TAG_NOT_FOUND == result) {
    int dum_val_array[] = {0, 0, 0, 0};
    result = impl->tag_create(HAS_MID_NODES_TAG_NAME, 
                              4*sizeof(int), 
                              MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER,
                              mHasMidNodesTag,
                              dum_val_array);
  }
  
  result = impl->tag_get_handle(GLOBAL_ID_TAG_NAME, mGlobalIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER, mGlobalIdTag, &dum_val);
  
}

ReadCCMIO::~ReadCCMIO() 
{}

ErrorCode ReadCCMIO::load_file(const char *file_name,
                                 const EntityHandle* file_set,
                                 const FileOptions& opts,
                                 const ReaderIface::IDTag* subset_list,
                                 int subset_list_length,
                                 const Tag* file_id_tag)
{
  CCMIOID rootID, problemID, stateID, processorID;
  CCMIOError error = kCCMIONoErr;

  CCMIOOpenFile(&error, file_name, kCCMIORead, &rootID);
  CHKCCMERR(error, "Problem opening file.");

    // get the file state
  ErrorCode rval = get_state(rootID, problemID, stateID);
  CHKERR(rval,NULL);

    // get processors
  std::set<CCMIOSize_t> procs;
  rval = get_processors(stateID, processorID, procs);
  CHKERR(rval,NULL);

  std::set<CCMIOSize_t>::iterator sit;
  Range new_ents, *new_ents_ptr = NULL;
  if (file_set) new_ents_ptr = &new_ents;
  
  for (sit = procs.begin(); sit != procs.end(); sit++) {
    rval = read_processor(stateID, problemID, processorID, *sit, new_ents_ptr);
    CHKERR(rval,NULL);
  }

    // load some meta-data
  rval = load_metadata(rootID, problemID, file_set);
  CHKERR(rval,NULL);

    // now, put all this into the file set, if there is one
  if (file_set) {
    rval = mbImpl->add_entities(*file_set, new_ents);
    CHKERR(rval, "Failed to add new entities to file set.");
  }
  
  return rval;
}

ErrorCode ReadCCMIO::load_metadata(CCMIOID rootID, CCMIOID problemID,
                                     const EntityHandle *file_set) 
{
    // Read the simulation title.
  CCMIOError error = kCCMIONoErr;
  ErrorCode rval = MB_SUCCESS;
  CCMIONode rootNode;

  if (kCCMIONoErr == CCMIOGetEntityNode(&error, rootID, &rootNode)) {
    char *name = NULL;
    CCMIOGetTitle(&error, rootNode, &name);

    if (NULL != name && strlen(name) != 0) {
        // make a tag for it and tag the read set
      Tag simname;
      rval = mbImpl->tag_get_handle("SimulationName", simname);
      if (MB_TAG_NOT_FOUND == rval) {
        rval = mbImpl->tag_create("SimulationName", strlen(name), MB_TAG_SPARSE, 
                                  MB_TYPE_OPAQUE, simname, NULL);
        CHKERR(rval, "Simulation name tag not found or created.");
      }
      EntityHandle tag_set = (NULL != file_set ? *file_set : 0);
      rval = mbImpl->tag_set_data(simname, &tag_set, 1, name);
      CHKERR(rval, "Problem setting simulation name tag.");

    }
    if (name) free(name);
  }

  rval = load_matset_data(problemID);
  CHKERR(rval, NULL);
  
  return rval;
}

ErrorCode ReadCCMIO::load_matset_data(CCMIOID problemID) 
{
    // make sure there are matsets
  if (newMatsets.empty()) return MB_SUCCESS;
  
    // ... walk through each cell type
  CCMIOSize_t i = CCMIOSIZEC(0);
  CCMIOID next;
  std::vector<char> mtype_name;
  CCMIOError error = kCCMIONoErr;
  Tag matNameTag = 0, matPorosityTag = 0, matSpinTag = 0, matGroupTag = 0;
  ErrorCode rval = create_matset_tags(matNameTag, matPorosityTag, 
                                        matSpinTag, matGroupTag);
  CHKERR(rval, NULL);
  
  while (CCMIONextEntity(NULL, problemID, kCCMIOCellType, &i, &next)
         == kCCMIONoErr) {
    int csize, mindex, idum;

    CCMIOGetEntityIndex(&error, next, &mindex);
    assert(mindex > 0 && mindex <= (int)newMatsets.size()+1);
    mindex--;

      // do Id first, since it should be set to something to identify set as a matset
    if (kCCMIONoErr != CCMIOReadOpti(NULL, next, "MaterialId", &idum))
      idum = mindex + 1;
    EntityHandle dum_ent = newMatsets[mindex];
    rval = mbImpl->tag_set_data(mMaterialSetTag, &dum_ent, 1, &idum);
    CHKERR(rval, "Failed to set material set id tag.");

      // MaterialType
    if (kCCMIONoErr == CCMIOReadOptstr(NULL, next, "MaterialType", &csize, NULL)) {
      mtype_name.resize(csize+1, '\0');
      CCMIOReadOptstr(&error, next, "MaterialType", &csize, &mtype_name[0]);

      rval = mbImpl->tag_set_data(matNameTag, &dum_ent, 1, &mtype_name[0]);
      CHKERR(rval, "Failed to set matNameTag.");
      
    }

      // porosity
    if (kCCMIONoErr == CCMIOReadOpti(NULL, next, "PorosityId", &idum)) {
      rval = mbImpl->tag_set_data(matPorosityTag, &dum_ent, 1, &idum);
      CHKERR(rval, "Failed to set matPorosityTag.");
      
    }

      // Spin
    if (kCCMIONoErr == CCMIOReadOpti(NULL, next, "SpinId", &idum)) {
      rval = mbImpl->tag_set_data(matSpinTag, &dum_ent, 1, &idum);
      CHKERR(rval, "Failed to set matSpinTag.");
      
    }

      // group
    if (kCCMIONoErr == CCMIOReadOpti(NULL, next, "GroupId", &idum)) {
      rval = mbImpl->tag_set_data(matGroupTag, &dum_ent, 1, &idum);
      CHKERR(rval, "Failed to set matGroupTag.");
    }
  }

  return MB_SUCCESS;
}

ErrorCode ReadCCMIO::create_matset_tags(Tag &matNameTag, Tag &matPorosityTag, 
                                          Tag &matSpinTag, Tag &matGroupTag)
{
  ErrorCode rval = mbImpl->tag_create("MaterialName", NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE,
                                        matNameTag, NULL, true);
  CHKERR(rval, "Failed to create matNameTag.");

  rval = mbImpl->tag_create("MaterialPorosity", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER,
                            matPorosityTag, NULL, true);
  CHKERR(rval, "Failed to create matPorosityTag.");

  rval = mbImpl->tag_create("MaterialSpin", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER,
                            matSpinTag, NULL, true);
  CHKERR(rval, "Failed to create matSpinTag.");

  rval = mbImpl->tag_create("MaterialGroup", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER,
                            matGroupTag, NULL, true);
  CHKERR(rval, "Failed to create matGroupTag.");

  return MB_SUCCESS;
}

ErrorCode ReadCCMIO::read_processor(CCMIOID stateID, CCMIOID problemID,
                                      CCMIOID processorID, CCMIOSize_t proc,
                                      Range *new_ents) 
{
  CCMIOError error = kCCMIONoErr;
  ErrorCode rval;
  bool has_solution = true;
  CCMIOID verticesID, topologyID, solutionID;
  
    // read the vertices, topology, and solution ids
  proc = CCMIOSIZEC(0);
  CCMIONextEntity(&error, stateID, kCCMIOProcessor, &proc, &processorID);
  CCMIOReadProcessor(&error, processorID, &verticesID, &topologyID, NULL, &solutionID);
  if(kCCMIONoErr != error) {
      // might not be a solution; try reading just verts & topology
    error = kCCMIONoErr;
    CCMIOReadProcessor(&error, processorID, &verticesID, &topologyID, NULL, NULL);
    if(kCCMIONoErr == error)
        hasSolution = false;
    else CHKERR(MB_FAILURE, "Couldn't get vertices and topology.");
  }

    // vert_map fields: s: none, i: gid, ul: vert handle, r: none
    //TupleList vert_map(0, 1, 1, 0, 0);
  TupleList vert_map;
  rval = read_vertices(proc, processorID, verticesID, topologyID, 
                       solutionID, has_solution, new_ents, vert_map);
  CHKERR(rval, NULL);
  
  rval = read_cells(proc, problemID, verticesID, topologyID, 
                    solutionID, has_solution, vert_map, new_ents);
  CHKERR(rval, NULL);

  return rval;
}

ErrorCode ReadCCMIO::read_cells(CCMIOSize_t proc, CCMIOID problemID,
                                  CCMIOID verticesID, CCMIOID topologyID,
                                  CCMIOID solutionID, bool has_solution,
                                  TupleList &vert_map, Range *new_ents) 
{

    // read the faces.
    // face_map fields: s:forward/reverse, i: cell id, ul: face handle, r: none
  ErrorCode rval;
#ifdef TUPLE_LIST
  TupleList face_map(1, 1, 1, 0, 0); 
#else
  TupleList face_map;
  SenseList sense_map;
#endif
  rval = read_all_faces(topologyID, vert_map, face_map
#ifndef TUPLE_LIST
                        , sense_map
#endif
                        , new_ents);
  CHKERR(rval, NULL);
  
    // now construct the cells; sort the face map by cell ids first
#ifdef TUPLE_LIST  
  rval = face_map.sort(1);
  CHKERR(rval, "Couldn't sort face map by cell id.");
#endif
  std::vector<EntityHandle> new_cells;
  rval = construct_cells(face_map, 
#ifndef TUPLE_LIST
                         sense_map,
#endif
                         vert_map, new_cells);
  CHKERR(rval, NULL);
  if (new_ents) {
    Range::iterator rit = new_ents->end();
    std::vector<EntityHandle>::reverse_iterator vit;
    for (vit = new_cells.rbegin(); vit != new_cells.rend(); vit++)
      rit = new_ents->insert(rit, *vit);
  }
  
  rval = read_gids_and_types(problemID, topologyID, new_cells);
  CHKERR(rval, NULL);
  
  return MB_SUCCESS;
}

ErrorCode ReadCCMIO::read_gids_and_types(CCMIOID problemID,
                                           CCMIOID topologyID,
                                           std::vector<EntityHandle> &cells) 
{
    // get the cells entity and number of cells
  CCMIOSize_t num_cells;
  CCMIOError error = kCCMIONoErr;
  CCMIOID cellsID, mapID;
  CCMIOGetEntity(&error, topologyID, kCCMIOCells, 0, &cellsID);
  CCMIOEntitySize(&error, cellsID, &num_cells, NULL);

    // check the number of cells against how many are in the cell array
  if (num_cells != (int)cells.size())
    CHKERR(MB_FAILURE, "Number of cells doesn't agree.");

    // read the gid map and set global ids
  std::vector<int> cell_gids(GETINT32(num_cells));
  CCMIOReadCells(&error, cellsID, &mapID, NULL,
                 CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CCMIOReadMap(&error, mapID, &cell_gids[0], 
               CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CHKCCMERR(error, "Couldn't read cells or cell id map.");

  ErrorCode rval = mbImpl->tag_set_data(mGlobalIdTag, &cells[0], 
                                          cells.size(), &cell_gids[0]);
  CHKERR(rval, "Couldn't set gids tag.");

    // now read cell types; first get the number of different types, so we can make material
    // sets
  int num_matsets = 0;
  CCMIOSize_t dum = CCMIOSIZEC(0);
  // if we don't have a problem description, return
  if (!CCMIOIsValidEntity(problemID)) return MB_SUCCESS;

    // get the number of cells
  CCMIOID next;
  while (kCCMIONoErr == CCMIONextEntity(NULL, problemID, kCCMIOCellType, &dum, &next))
    num_matsets++;

    // create the matsets
  for (int i = 0; i < num_matsets; i++) {
    EntityHandle matset;
    rval = mbImpl->create_meshset(MESHSET_SET, matset);
    CHKERR(rval, "Couldn't create material set.");
    newMatsets.insert(matset);
  }
  
    // read in cell types and use to assign cells to matsets; reuse cell_gids
  CCMIOReadCells(&error, cellsID, NULL, &cell_gids[0],
                 CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CHKCCMERR(error, "Trouble reading cell types.");

    // check min/max
  int min = *(std::min_element(cell_gids.begin(), cell_gids.end())),
      max = *(std::max_element(cell_gids.begin(), cell_gids.end()));
  if (1 != min || num_matsets != max) CHKERR(MB_FAILURE, "Wrong range for cell types.");
  
  if (min == max && num_matsets == 1) {
      // no streaming necessary, all elements are in this set
    rval = mbImpl->add_entities(newMatsets[0], &cells[0], cells.size());
    CHKERR(rval, "Trouble adding cells to matset.");
  }
  else {
      // now stream through cell types, gather cells for each type, add to matsets
    Range mcells;
    for (int i = 1; i <= num_matsets; i++) {
      mcells.clear();
        // go backwards, since it's more efficient when adding to the end of a range
      for (int j = num_cells-1; j >= 0; j--)
        if (cell_gids[j] == i) mcells.insert(cells[j]);
    
      rval = mbImpl->add_entities(newMatsets[i-1], mcells);
      CHKERR(rval, "Trouble adding cells to matset.");

    }
  }
  
  return MB_SUCCESS;
}


ErrorCode ReadCCMIO::construct_cells(TupleList &face_map, 
#ifndef TUPLE_LIST
                                       SenseList &sense_map, 
#endif
                                       TupleList &vert_map,
                                       std::vector<EntityHandle> &new_cells) 
{
  std::vector<EntityHandle> facehs;
  std::vector<int> senses;
  EntityHandle cell;
  ErrorCode tmp_rval, rval = MB_SUCCESS;
#ifdef TUPLE_LIST
  unsigned int i = 0;
  while (i < face_map.n) {
      // pull out face handles bounding the same cell
    facehs.clear();
    int this_id = face_map.get_int(i);
    unsigned int inext = i;
    while (face_map.get_int(inext) == this_id && inext <= face_map.n) {
      inext++;
      EntityHandle face = face_map.get_ulong(inext);
      facehs.push_back(face);
      senses.push_back(face_map.get_short(inext));
    }
#else
      
  std::map<int,std::vector<EntityHandle> >::iterator fmit;
  std::map<int,std::vector<int> >::iterator smit;
  for (fmit = face_map.begin(), smit = sense_map.begin();
       fmit != face_map.end(); fmit++, smit++) {
      // pull out face handles bounding the same cell
    facehs.clear();
    int this_id = (*fmit).first;
    facehs = (*fmit).second;
    senses.clear();
    senses = (*smit).second;
#endif
    tmp_rval = create_cell_from_faces(facehs, senses, cell);
    if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
    else {
      new_cells.push_back(cell);
        // tag cell with global id
      tmp_rval = mbImpl->tag_set_data(mGlobalIdTag, &cell, 1, &this_id);
      if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
    }
  }
    
  return MB_SUCCESS;
}

ErrorCode ReadCCMIO::create_cell_from_faces(std::vector<EntityHandle> &facehs,
                                              std::vector<int> &senses,
                                              EntityHandle &cell) 
{
    // test to see if they're one type
  EntityType this_type = mbImpl->type_from_handle(facehs[0]);
  bool same_type = true;
  for (std::vector<EntityHandle>::iterator vit = facehs.begin(); vit != facehs.end(); vit++) {
    if (this_type != mbImpl->type_from_handle(*vit)) {
      same_type = false;
      break;
    }
  }
  
    // if different, we can quit here, we'll consider this a polyhedron
  ErrorCode rval = MB_SUCCESS;
  if (!same_type || 
      (MBTRI == this_type && facehs.size() != 4) ||
      (MBQUAD == this_type && facehs.size() != 6) ||
      (MBQUAD != this_type && MBTRI != this_type)) {
    rval = mbImpl->create_element(MBPOLYHEDRON, &facehs[0], facehs.size(), cell);
    CHKERR(rval, "Couldn't make polyhedron.");
    return rval;
  }
  
    // try tet and hex elements; get connectivity of first face
  std::vector<EntityHandle> verts;
  rval = mbImpl->get_connectivity(&facehs[0], 1, verts);
  CHKERR(rval, "Couldn't get connectivity.");
  bool match = false;

    // reverse connectivity if sense is forward, since base face always points
    // into entity
  if (senses[0] > 0) std::reverse(verts.begin(), verts.end());

  std::vector<EntityHandle> storage;
  MeshTopoUtil mtu(mbImpl);
  if (MBTRI == this_type) {
      // get the 4th vertex through the next tri
    const EntityHandle *conn; int conn_size;
    rval = mbImpl->get_connectivity(facehs[1], conn, conn_size, true, &storage);
    CHKERR(rval, "Couldn't get connectivity.");
    int i = 0;
    while (std::find(verts.begin(), verts.end(), conn[i]) != verts.end() && i < conn_size) i++;
    if (conn_size == i) 
      CHKERR(MB_FAILURE, "Didn't find apex vertex.");

    match = true;
    this_type = MBTET;
  }
  else if (MBQUAD == this_type) {
      // build hex from quads
      // algorithm:
      // - verts = vertices from 1st quad
      // - find quad q1 sharing verts[0] and verts[1]
      // - find quad q2 sharing other 2 verts in q1
      // - find v1 = opposite vert from verts[1] in q1 , v2 = opposite from verts[0]
      // - get i = offset of v1 in verts2 of q2, rotate verts2 by i
      // - if verts2[i+1%4] != v2, flip verts2 by switching verts2[1] and verts2[3]
      // - append verts2 to verts


      // get the other vertices for this hex; need to find the quad with no common vertices
    Range tmp_faces, tmp_verts;

      // get q1, which shares 2 vertices with q0
    std::copy(facehs.begin(), facehs.end(), range_inserter(tmp_faces));
    rval = mbImpl->get_adjacencies(&verts[0], 2, 2, false, tmp_faces);
    if (MB_SUCCESS != rval || tmp_faces.size() != 2)
      CHKERR(MB_FAILURE, "Couldn't get adj face.");
    tmp_faces.erase(facehs[0]);
    EntityHandle q1 = *tmp_faces.begin();
      // get other 2 verts of q1
    rval = mbImpl->get_connectivity(&q1, 1, tmp_verts);
    CHKERR(rval, "Couldn't get adj verts.");
    tmp_verts.erase(verts[0]); tmp_verts.erase(verts[1]);
      // get q2
    std::copy(facehs.begin(), facehs.end(), range_inserter(tmp_faces));
    rval = mbImpl->get_adjacencies(tmp_verts, 2, false, tmp_faces);
    if (MB_SUCCESS != rval || tmp_faces.size() != 2)
      CHKERR(MB_FAILURE, "Couldn't get adj face.");
    tmp_faces.erase(q1);
    EntityHandle q2 = *tmp_faces.begin();
      // get verts in q2
    rval = mbImpl->get_connectivity(&q2, 1, storage);
    CHKERR(rval, "Couldn't get adj vertices.");

      // get verts in q1 opposite from v[1] and v[0] in q0
    EntityHandle v0 = 0, v1 = 0;
    rval = mtu.opposite_entity(q1, verts[1], v0);
    rval = mtu.opposite_entity(q1, verts[0], v1);
    if (!v0 || !v1)
      CHKERR(MB_FAILURE, "Trouble finding opposite vertices.");

      // offset of v0 in q2, then rotate and flip
    unsigned int ioff = std::find(storage.begin(), storage.end(), v0) - storage.begin();
    if (4 == ioff)
      CHKERR(MB_FAILURE, "Trouble finding offset.");

    if (storage[(ioff+1)%4] != v1) {
      std::reverse(storage.begin(), storage.end());
      ioff = std::find(storage.begin(), storage.end(), v0) - storage.begin();
    }
    if (0 != ioff)
      std::rotate(storage.begin(), storage.begin()+ioff, storage.end());

      // copy into verts, and make hex
    std::copy(storage.begin(), storage.end(), std::back_inserter(verts));
    match = true;
    this_type = MBHEX;
  }
  if (!match) 
    CHKERR(MB_FAILURE, "Couldn't find vertices for hex.");
  
    // now make the element
  rval = mbImpl->create_element(this_type, &verts[0], verts.size(), cell);
  CHKERR(rval, "create_element failed.");
  
  return MB_SUCCESS;
}
  
ErrorCode ReadCCMIO::read_all_faces(CCMIOID topologyID, TupleList &vert_map, 
                                      TupleList &face_map
#ifndef TUPLE_LIST
                                      ,SenseList &sense_map
#endif
                                      , Range *new_faces) 
{
  CCMIOSize_t index = CCMIOSIZEC(0);
  CCMIOID faceID;
  ErrorCode rval;

    // get total # internal/bdy faces, size the face map accordingly
  int nint_faces = 0, nbdy_faces = 0;
  CCMIOSize_t nf;
  CCMIOError error = kCCMIONoErr;
  while (kCCMIONoErr == CCMIONextEntity(NULL, topologyID, kCCMIOBoundaryFaces, &index, 
                                        &faceID))
  {
    CCMIOEntitySize(&error, faceID, &nf, NULL);
    nbdy_faces = nbdy_faces + nf;
  }
  CCMIOGetEntity(&error, topologyID, kCCMIOInternalFaces, 0, &faceID);
  CCMIOEntitySize(&error, faceID, &nf, NULL);
  nint_faces = nint_faces + nf;
#ifdef TUPLE_LIST
  face_map.resize(2*nint_faces + nbdy_faces);
#endif
  
    // get multiple blocks of bdy faces
  index = CCMIOSIZEC(0);
  while (kCCMIONoErr == CCMIONextEntity(NULL, topologyID, kCCMIOBoundaryFaces, &index, 
                                        &faceID))
  {
    rval = read_faces(faceID, kCCMIOBoundaryFaces, vert_map, face_map
#ifndef TUPLE_LIST
                      , sense_map
#endif
                      , new_faces);
    CHKERR(rval, "Trouble reading boundary faces.");
  }
  
    // now get internal faces
  CCMIOGetEntity(&error, topologyID, kCCMIOInternalFaces, 0, &faceID);

  rval = read_faces(faceID, kCCMIOInternalFaces, vert_map,face_map
#ifndef TUPLE_LIST
                    , sense_map
#endif
                    , new_faces);
  CHKERR(rval, "Trouble reading internal faces.");

  return rval;
}

ErrorCode ReadCCMIO::read_faces(CCMIOID faceID, CCMIOEntity bdy_or_int,
                                  TupleList &vert_map,
                                  TupleList &face_map
#ifndef TUPLE_LIST
                                  ,SenseList &sense_map
#endif
                                  , Range *new_faces)
{
  if (kCCMIOInternalFaces != bdy_or_int && kCCMIOBoundaryFaces != bdy_or_int)
    CHKERR(MB_FAILURE, "Face type isn't boundary or internal.");

  CCMIOSize_t num_faces;
  CCMIOError error = kCCMIONoErr;
  CCMIOEntitySize(&error, faceID, &num_faces, NULL);

    // get the size of the face connectivity array (not really a straight connect
    // array, has n, connect(n), ...)
  CCMIOSize_t farray_size = CCMIOSIZEC(0);
  CCMIOReadFaces(&error, faceID, bdy_or_int, NULL, &farray_size, NULL,
                 CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CHKCCMERR(error, "Trouble reading face connectivity length.");
    

    // allocate vectors for holding farray and cells for each face; use new for finer
    // control of de-allocation
  int num_sides = (kCCMIOInternalFaces == bdy_or_int ? 2 : 1);
  int *farray = new int[GETINT32(farray_size)];

    // read farray and make the faces
  CCMIOID mapID;
  CCMIOReadFaces(&error, faceID, bdy_or_int, &mapID, NULL,
                 farray, CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CHKCCMERR(error, "Trouble reading face connectivity.");

  std::vector<EntityHandle> face_handles(GETINT32(num_faces), 0);
  ErrorCode rval = make_faces(farray, vert_map, face_handles);
  CHKERR(rval, NULL);

    // read face cells and make tuples
  int *face_cells;
  if (num_sides*num_faces < farray_size) face_cells = new int[num_sides*GETINT32(num_faces)];
  else face_cells = farray;
  CCMIOReadFaceCells(&error, faceID, bdy_or_int, face_cells,
                     CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CHKCCMERR(error, "Trouble reading face cells.");

  int *tmp_ptr = face_cells;
  for (int i = 0; i < num_faces; i++) {
#ifdef TUPLE_LIST
    short forward = 1, reverse = -1;
    face_map.push_back(&forward, tmp_ptr++, &face_handles[i], NULL);
    if (2 == num_sides)
      face_map.push_back(&reverse, tmp_ptr++, &face_handles[i], NULL);
#else
    face_map[*tmp_ptr].push_back(face_handles[i]);
    sense_map[*tmp_ptr++].push_back(1);
    if (2 == num_sides) {
      face_map[*tmp_ptr].push_back(face_handles[i]);
      sense_map[*tmp_ptr++].push_back(-1);
    }
#endif
  }

    // now read & set face gids, reuse face_cells 'cuz we know it's big enough
  CCMIOReadMap(&error, mapID, face_cells, CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CHKCCMERR(error, "Trouble reading face gids.");

  rval = mbImpl->tag_set_data(mGlobalIdTag, &face_handles[0], face_handles.size(), face_cells);
  CHKERR(rval, "Couldn't set face global ids.");

    // ok, now sort face handles, and add to range
  std::sort(face_handles.begin(), face_handles.end());
  if (new_faces) {
    Range::iterator rit = new_faces->end();
    for (std::vector<EntityHandle>::reverse_iterator vit = face_handles.rbegin(); 
         vit != face_handles.rend(); vit++)
      rit = new_faces->insert(rit, *vit);
  }
  
  return MB_SUCCESS;
}
  

ErrorCode ReadCCMIO::make_faces(int *farray, 
                                  TupleList &vert_map,
                                  std::vector<EntityHandle> &new_faces) 
{
  unsigned int num_faces = new_faces.size();
  std::vector<EntityHandle> verts;
  ErrorCode tmp_rval = MB_SUCCESS, rval = MB_SUCCESS;
  
  for (unsigned int i = 0; i < num_faces; i++) {
    int num_verts = *farray++;
    verts.resize(num_verts);

      // fill in connectivity by looking up by gid in vert tuple_list
    for (int j = 0; j < num_verts; j++) {
#ifdef TUPLE_LIST
      int tindex = vert_map.find(1, farray[j]);
      if (-1 == tindex) {
        tmp_rval = MB_FAILURE;
        break;
      }
      verts[j] = vert_map.get_ulong(tindex, 0);
#else
      verts[j] = (vert_map[farray[j]])[0];
#endif      
    }
    farray += num_verts;

    if (MB_SUCCESS == tmp_rval) {
    
        // make face
      EntityType ftype = (3 == num_verts ? MBTRI :
                            (4 == num_verts ? MBQUAD : MBPOLYGON));
      EntityHandle faceh;
      tmp_rval = mbImpl->create_element(ftype, &verts[0], num_verts, faceh);
      if (faceh) new_faces[i] = faceh;
    }
    
    if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
  }
  
  return rval;
}

ErrorCode ReadCCMIO::read_vertices(CCMIOSize_t proc, CCMIOID processorID, CCMIOID verticesID,
                                     CCMIOID topologyID, CCMIOID solutionID, bool has_solution,
                                     Range *verts, TupleList &vert_map) 
{
  CCMIOError error = kCCMIONoErr;
  
    // pre-read the number of vertices, so we can pre-allocate & read directly in
  CCMIOSize_t nverts = CCMIOSIZEC(0);
  CCMIOEntitySize(&error, verticesID, &nverts, NULL);
  CHKCCMERR(error, "Couldn't get number of vertices.");

    // get # dimensions
  CCMIOSize_t dims;
  float scale;
  CCMIOReadVerticesf(&error, verticesID, &dims, NULL, NULL, NULL, CCMIOINDEXC(0), CCMIOINDEXC(1));
  CHKCCMERR(error, "Couldn't get number of dimensions.");

    // allocate vertex space
  EntityHandle node_handle = 0;
  std::vector<double*> arrays;
  readMeshIface->get_node_coords(3, GETINT32(nverts), MB_START_ID, node_handle, arrays);

    // read vertex coords
  CCMIOID mapID;
  std::vector<double> tmp_coords(GETINT32(dims)*GETINT32(nverts));
  CCMIOReadVerticesd(&error, verticesID, &dims, &scale, &mapID, &tmp_coords[0], 
                     CCMIOINDEXC(0), CCMIOINDEXC(0+nverts));
  CHKCCMERR(error, "Trouble reading vertex coordinates.");

    // copy interleaved coords into moab blocked coordinate space
  int i = 0, threei = 0;
  for (; i < nverts; i++) {
    arrays[0][i] = tmp_coords[threei++];
    arrays[1][i] = tmp_coords[threei++];
    if (3 == GETINT32(dims)) arrays[2][i] = tmp_coords[threei++];
    else arrays[2][i] = 0.0;
  }

    // scale, if necessary
  if (1.0 != scale) {
    for(i = 0; i < nverts; i++) {
      arrays[0][i] *= scale;
      arrays[1][i] *= scale;
      if (3 == GETINT32(dims)) arrays[2][i] *= scale;
    }
  }

    // put new vertex handles into range
  if (verts) verts->insert(node_handle, node_handle+nverts);

    // pack vert_map with global ids and handles for these vertices
  std::vector<int> gids(GETINT32(nverts));
  CCMIOReadMap(&error, mapID, &gids[0], CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CHKCCMERR(error, "Trouble reading vertex global ids.");
#ifdef TUPLE_LIST
  vert_map.resize(GETINT32(nverts));
  for (i = 0; i < GETINT32(nverts); i++) {
    vert_map.push_back(NULL, &gids[i], &node_handle, NULL);
#else
  for (i = 0; i < GETINT32(nverts); i++) {
    (vert_map[gids[i]]).push_back(node_handle);
#endif
    node_handle += 1;
  }
  
  return MB_SUCCESS;
}
  
ErrorCode ReadCCMIO::get_processors(CCMIOID stateID, CCMIOID &processorID,
                                      std::set<CCMIOSize_t> &procs) 
{
  CCMIOSize_t proc = CCMIOSIZEC(0);
  while (kCCMIONoErr == 
         CCMIONextEntity(NULL, stateID, kCCMIOProcessor, &proc, &processorID))
    procs.insert(proc);

  return MB_SUCCESS;
}

ErrorCode ReadCCMIO::get_state(CCMIOID rootID, CCMIOID &problemID, CCMIOID &stateID) 
{
  CCMIOError error = kCCMIONoErr;
  
    // first try default
  CCMIOGetState(&error, rootID, "default", &problemID, &stateID);
  if (kCCMIONoErr != error) {
    CCMIOSize_t i = CCMIOSIZEC(0);
    CCMIOError tmp_error = kCCMIONoErr;
    CCMIONextEntity(&tmp_error, rootID, kCCMIOState, &i, &stateID);
    if (kCCMIONoErr ==  tmp_error)
      CCMIONextEntity(&error, rootID, kCCMIOProblemDescription, 
                      &i, &problemID);
  }
  CHKCCMERR(error, "Couldn't find state.");

  return MB_SUCCESS;
}

ErrorCode ReadCCMIO::read_tag_values( const char* file_name,
                                        const char* tag_name,
                                        const FileOptions& opts,
                                        std::vector<int>& tag_values_out,
                                        const IDTag* subset_list,
                                        int subset_list_length) 
{
  return MB_FAILURE;
}

} // namespace moab
