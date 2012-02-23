#include "WriteDamsel.hpp"

#include "damsel.h"
#include "assert.h"
#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/Error.hpp"
#include "moab/WriteUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "EntitySequence.hpp"
#include "Internals.hpp"
#include "DenseTag.hpp"

namespace moab {

  // Some macros to handle error checking (cribbed from WriteHDF5).  The
  // CHK_MB_ERR_* check the value of an ErrorCode.
  // The *_0 macros accept no other arguments. The *_1
  // macros accept a single damsel handle to close on error.
  // All macros contain a "return" statement.  These macros are coded with a do if while
  // to allow statements calling them to be terminated with a ;
#define CHK_MB_ERR( A, B )                                    \
do if (MB_SUCCESS != (A)) { \
mError->set_last_error(B);\
return error(A);} while(false)

#define CHK_MB_ERR_2( A, B, C )                   \
do if (MB_SUCCESS != (A)) { \
mError->set_last_error(B, C);                 \
return error(A);} while(false)

#define CHK_MB_ERR_FINALIZE( A, B )       \
do if (MB_SUCCESS != (A)) {             \
  DMSLlib_finalize(dmslLib); \
  dmslLib = 0;          \
  mError->set_last_error(B);\
  return error(A);     \
} while(false)

#define CHK_DMSL_ERR( A, B )                 \
do if (DMSL_OK.id != (A.id)) {             \
mError->set_last_error(B);\
return error(MB_FAILURE);             \
} while(false)

#define CHK_DMSL_ERR_2( A, B, C )            \
do if (DMSL_OK.id != (A.id)) {             \
mError->set_last_error(B, C);            \
return error(MB_FAILURE);             \
} while(false)

#define CHK_DMSL_ERR_FINALIZE( A, B )        \
do if (DMSL_OK.id != (A.id)) {             \
  DMSLlib_finalize(dmslLib); \
  dmslLib = 0;          \
  mError->set_last_error(B);\
  return error(MB_FAILURE);                    \
} while(false)

// This function doesn't do anything useful.  It's just a nice
// place to set a break point to determine why the reader fails.
    static inline ErrorCode error( ErrorCode rval )
  { return rval; }

static damsel_entity_type moab_to_damsel_entity_type[] = {
    DAMSEL_ENTITY_TYPE_VERTEX,
    DAMSEL_ENTITY_TYPE_EDGE,
    DAMSEL_ENTITY_TYPE_TRI,
    DAMSEL_ENTITY_TYPE_QUAD,
    DAMSEL_ENTITY_TYPE_POLYGON,
    DAMSEL_ENTITY_TYPE_TET,
    DAMSEL_ENTITY_TYPE_PYRAMID,
    DAMSEL_ENTITY_TYPE_PRISM,
    DAMSEL_ENTITY_TYPE_UNDEFINED,
    DAMSEL_ENTITY_TYPE_HEX,
    DAMSEL_ENTITY_TYPE_POLYHEDRON,
    DAMSEL_ENTITY_TYPE_UNDEFINED,
    DAMSEL_ENTITY_TYPE_UNDEFINED
};

static damsel_data_type moab_to_damsel_data_type[] = {
    DAMSEL_DATA_TYPE_BYTES,
    DAMSEL_DATA_TYPE_INTEGER,
    DAMSEL_DATA_TYPE_DOUBLE,
    DAMSEL_DATA_TYPE_INVALID,
    DAMSEL_DATA_TYPE_HANDLE,
    DAMSEL_DATA_TYPE_INVALID
};

WriterIface* WriteDamsel::factory( Interface* iface )
  { return new WriteDamsel( iface ); }

WriteDamsel::WriteDamsel(Interface *impl) 
    : mbImpl(impl)
{
  assert(impl != NULL);

  impl->query_interface( mWriteIface );
  assert(mWriteIface);
  
  sequenceManager = dynamic_cast<Core*>(impl)->sequence_manager();
  assert(sequenceManager);

  impl->query_interface(mError);
  assert(mError);
  
  impl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
                       mGlobalIdTag, MB_TAG_SPARSE|MB_TAG_CREAT);

  moabHandleType = (sizeof(EntityHandle) == 64 ? DAMSEL_HANDLE_TYPE_HANDLE64 :
                    DAMSEL_HANDLE_TYPE_HANDLE32);
  
  dmslLib = DMSLlib_init();
}

WriteDamsel::~WriteDamsel() 
{
  mbImpl->release_interface(mWriteIface);
  DMSLlib_finalize(dmslLib);
}

ErrorCode WriteDamsel::write_file(const char *file_name, 
                                  const bool /* overwrite */,
                                  const FileOptions& opts,
                                  const EntityHandle *meshset_list,
                                  const int num_sets,
                                  const std::vector<std::string>& /* qa_records */,
                                  const Tag* /* tag_list */,
                                  int /* num_tags */,
                                  int /* requested_output_dimension */)
{
    // gather all entities into one big range
  Range all_ents;
  ErrorCode rval;
  damsel_err_t err;
  
  rval = mWriteIface->gather_entities(all_ents, meshset_list, num_sets);
  CHK_MB_ERR(rval, "Gather entities failed in WriteDamsel.");

  if (all_ents.empty()) return MB_SUCCESS;

  rval = init_dense_tag_info();
  CHK_MB_ERR(rval, NULL);
  
    // create a damsel model
  dmslModel = DMSLmodel_create(sizeof(EntityHandle) == 64 ? DAMSEL_HANDLE_TYPE_HANDLE64 : 
                               DAMSEL_HANDLE_TYPE_HANDLE32);
  
    // iterate through the groups of contiguous sequences of handles
  RangeSeqIntersectIter rsi(sequenceManager);
  rval = rsi.init(all_ents.begin(), all_ents.end());
  
  while (rval == MB_SUCCESS) {
    EntityHandle start = rsi.get_start_handle();
    
    if (MBVERTEX == mbImpl->type_from_handle(start)) 
      rval = write_vertices(rsi);
    
    else if (MBENTITYSET > mbImpl->type_from_handle(start)) 
      rval = write_entities(rsi);

    else
      rval = write_sets(rsi);
    
    rval = rsi.step();
    while (MB_ENTITY_NOT_FOUND == rval)
      rval = rsi.step();
  }
  
    // now tell Damsel to actually write it
  MPI_Comm comm = 0;
  err = DMSLmodel_attach(dmslModel, file_name, comm, NULL);
  CHK_DMSL_ERR(err, "DMSLmodel_attach failed.");
  
  damsel_request_t request;
  err = DMSLmodel_transfer_async(dmslModel, DAMSEL_TRANSFER_TYPE_WRITE, &request);
  CHK_DMSL_ERR(err, "DMSLmodel_transfer_asynch failed.");
  
  damsel_status_t status;
  err = DMSLmodel_wait(request, &status);
  CHK_DMSL_ERR(err, "DMSLmodel_wait failed.");
  
    // we should be done
  return MB_SUCCESS;
}

ErrorCode WriteDamsel::write_sets(RangeSeqIntersectIter &rsi) 
{
    // write the sets

  ErrorCode rval = MB_SUCCESS;
  std::vector<EntityHandle> ents;
  damsel_id_t dseth;
  damsel_err_t err;
  for (EntityHandle seth = rsi.get_start_handle(); seth <= rsi.get_end_handle(); seth++) {
    // get the set type (range or set)
    unsigned int opts;
    rval = mbImpl->get_meshset_options(seth, opts);
    CHK_MB_ERR_2(rval, "Failed to get options for meshset %lu.", seth);
    damsel_collection_type coll_type = (opts&MESHSET_SET ? DAMSEL_HANDLE_COLLECTION_TYPE_SET :
                                        DAMSEL_HANDLE_COLLECTION_TYPE_VECTOR);

      // make a damsel collection
    dseth = DMSLcoll_create(dmslModel, (damsel_handle_ptr)&seth, coll_type);
    if (!DMSLlib_ok(dmslLib, dseth)) CHK_MB_ERR_2(MB_FAILURE, "Bad handle returned by Damsel for meshset %lu.", seth);

      // get all the entities & add
    ents.clear();
    rval = mbImpl->get_entities_by_handle(seth, ents);
    CHK_MB_ERR_2(rval, "get_entities_by_handle failed for set %lu.", seth);
    if (!ents.empty()) {
      err = DMSLcoll_add_fast(dmslModel, (damsel_handle_ptr)&ents[0], ents.size());
      CHK_DMSL_ERR_2(err, "DMSLcoll_add_fast failed for meshset %lu.", seth);
    }

    // parents/children...
  }
  
    // dense tags
  rval = write_dense_tags(rsi);

  return rval;
}

ErrorCode WriteDamsel::write_entities(RangeSeqIntersectIter &rsi) 
{
    // write the entities; these entities will be in the same sequence and will be contiguous, guaranteed
  EntityHandle start_ent = rsi.get_start_handle(), end_ent = rsi.get_end_handle();

    // create a damsel container for these entity handles
  damsel_id_t ent_cont = DMSLhandle_create_sequence(dmslModel, moabHandleType, (damsel_handle_ptr)&start_ent,
                                                       end_ent-start_ent+1, 1);
  if (!DMSLlib_ok(dmslLib, ent_cont))
    CHK_MB_ERR(MB_FAILURE, "Bad sequence returned by Damsel.");

    // get # verts per entity and entity type
  EntityType etype = mbImpl->type_from_handle(start_ent);
  assert(MBMAXTYPE != etype);
  int num_connect = rsi.get_sequence()->values_per_entity();
  assert(0 < num_connect);
  
    // define the entities to damsel
  damsel_err_t err = DMSLentity_define(ent_cont, moab_to_damsel_entity_type[etype], num_connect);
  CHK_DMSL_ERR_2(err, "DMSLentity_define failed for entities starting with handle %lu.", rsi.get_start_handle());
  
    // get the connectivity storage location and pass to damsel
  Range ent_range(start_ent, end_ent);
  int count;
  EntityHandle *connect;
  int verts_per_ent;
  ErrorCode rval = mbImpl->connect_iterate(ent_range.begin(), ent_range.end(), connect, verts_per_ent, count);
  CHK_MB_ERR_2(rval, "Failed to get connect iterator for entities starting with handle %lu.", rsi.get_start_handle());
  if (count != (int)ent_range.size())
    CHK_MB_ERR_2(MB_FAILURE, "Entity subrange not in the same sequence for entities starting with handle %lu.", 
               rsi.get_start_handle());
  
  err = DMSLentity_set_connectivity_fast(ent_cont, (damsel_handle_ptr)connect);
  CHK_DMSL_ERR_2(err, "Failed in DMSLentity_set_connectivity_fast for entities starting with handle %lu.", 
               rsi.get_start_handle());

    // write dense tags
  rval = write_dense_tags(rsi, ent_cont);
  CHK_MB_ERR(rval, NULL);
  
  return MB_SUCCESS;
}

ErrorCode WriteDamsel::write_vertices(RangeSeqIntersectIter &rsi) 
{
    // write the vertices; these vertices will be in the same sequence and will be contiguous, guaranteed
  EntityHandle start_vert = rsi.get_start_handle(), end_vert = rsi.get_end_handle();

    // create a damsel container for these vertex handles
  damsel_id_t vertex_cont = DMSLhandle_create_sequence(dmslModel, moabHandleType, (damsel_handle_ptr)&start_vert,
                                                       end_vert-start_vert+1, 1);
  if (!DMSLlib_ok(dmslLib, vertex_cont)) 
    CHK_MB_ERR_2(MB_FAILURE, "Failed to create vertex sequence for vertices starting with handle %lu.", 
               rsi.get_start_handle());
  
    // define the entities to damsel
  damsel_err_t err = DMSLentity_define(vertex_cont, DAMSEL_ENTITY_TYPE_VERTEX, 1);
  CHK_DMSL_ERR_2(err, "Failure in DMSLentity_define for vertices starting with handle %lu.", 
               rsi.get_start_handle());
  
    // get the vertex coordinates storage locations and pass to damsel
  Range vert_range(start_vert, end_vert);
  double *xcoords = NULL, *ycoords = NULL, *zcoords = NULL;
  int count;
  ErrorCode rval = mbImpl->coords_iterate(vert_range.begin(), vert_range.end(),
                                          xcoords, ycoords, zcoords, count);
  CHK_MB_ERR_2(rval, "Failed to get coordinate iterator for vertices starting with handle %lu.", 
             rsi.get_start_handle());
  if (count != (int)vert_range.size()) {
    CHK_MB_ERR_2(MB_FAILURE, "Vertex subrange not in the same sequence for vertices starting with handle %lu.", 
               rsi.get_start_handle());
  }
  
  damsel_id_t xcoords_dtag, ycoords_dtag, zcoords_dtag;
  
  if (xcoords && !ycoords && !zcoords) {
      // interleaved
      // get/define the damsel tag for coordinates
    rval = damsel_coords_tags(xcoords_dtag, true);
    if (MB_SUCCESS != rval || !DMSLlib_ok(dmslLib, xcoords_dtag))
      CHK_MB_ERR(MB_FAILURE, "Bad Damsel tag returned.");

      // map the data to damsel
    err = DMSLtag_assign(xcoords_dtag, vertex_cont, xcoords);
    CHK_DMSL_ERR_2(err, "Failed to assign vertex coordinates tag for vertices starting with handle %lu.", 
                 rsi.get_start_handle());
  }
  else {
    rval = damsel_coords_tags(xcoords_dtag, ycoords_dtag, zcoords_dtag, true);
    if (MB_SUCCESS != rval || !DMSLlib_ok(dmslLib, xcoords_dtag) || !DMSLlib_ok(dmslLib, ycoords_dtag) ||
        !DMSLlib_ok(dmslLib, zcoords_dtag))
      CHK_DMSL_ERR_2(err, "Failed to assign vertex coordinates tag for vertices starting with handle %lu.", 
                   rsi.get_start_handle());
      // map the data to damsel
    err = DMSLtag_assign(xcoords_dtag, vertex_cont, xcoords);
    CHK_DMSL_ERR_2(err, "Failed to assign vertex coordinates tag for vertices starting with handle %lu.", 
                 rsi.get_start_handle());
    err = DMSLtag_assign(ycoords_dtag, vertex_cont, ycoords);
    CHK_DMSL_ERR_2(err, "Failed to assign vertex coordinates tag for vertices starting with handle %lu.", 
                 rsi.get_start_handle());
    err = DMSLtag_assign(zcoords_dtag, vertex_cont, zcoords);
    CHK_DMSL_ERR_2(err, "Failed to assign vertex coordinates tag for vertices starting with handle %lu.", 
                 rsi.get_start_handle());
  }

    // write dense tags
  rval = write_dense_tags(rsi, vertex_cont);
  CHK_MB_ERR(rval, NULL);
  
  return MB_SUCCESS;
}

ErrorCode WriteDamsel::damsel_coords_tags(damsel_id_t &xcoords_dtag, damsel_id_t &ycoords_dtag, 
                                          damsel_id_t &zcoords_dtag, bool create_if_missing) 
{
  xcoords_dtag = dmslXcoord; ycoords_dtag = dmslYcoord; zcoords_dtag = dmslZcoord;
  if ((!DMSLlib_ok(dmslLib, xcoords_dtag) || !DMSLlib_ok(dmslLib, xcoords_dtag) || !DMSLlib_ok(dmslLib, xcoords_dtag)) &&
      !create_if_missing) return MB_SUCCESS;

  EntityHandle tmp_handle;
  ErrorCode rval = MB_SUCCESS;
  if (!DMSLlib_ok(dmslLib, xcoords_dtag)) {
    tmp_handle = CREATE_HANDLE(MBMAXTYPE, 1);
    xcoords_dtag = dmslXcoord = DMSLtag_define(dmslModel, (damsel_handle_ptr)&tmp_handle,
                                               DAMSEL_DATA_TYPE_DOUBLE, "XCOORDS");
    if (!DMSLlib_ok(dmslLib, xcoords_dtag)) rval = MB_FAILURE;
  }
  
  if (!DMSLlib_ok(dmslLib, ycoords_dtag)) {
    tmp_handle = CREATE_HANDLE(MBMAXTYPE, 2);
    ycoords_dtag = dmslYcoord = DMSLtag_define(dmslModel, (damsel_handle_ptr)&tmp_handle,
                                               DAMSEL_DATA_TYPE_DOUBLE,
                                               "YCOORDS");
    if (!DMSLlib_ok(dmslLib, ycoords_dtag)) rval = MB_FAILURE;
  }
  
  if (!DMSLlib_ok(dmslLib, zcoords_dtag)) {
    tmp_handle = CREATE_HANDLE(MBMAXTYPE, 3);
    zcoords_dtag = dmslZcoord = DMSLtag_define(dmslModel, (damsel_handle_ptr)&tmp_handle,
                                               DAMSEL_DATA_TYPE_DOUBLE,
                                               "ZCOORDS");
    if (!DMSLlib_ok(dmslLib, zcoords_dtag)) rval = MB_FAILURE;
  }

  CHK_MB_ERR(rval, "Failed to get Damsel tag for blocked coordinates.");

  return rval;
}

ErrorCode WriteDamsel::damsel_coords_tags(damsel_id_t &coords_dtag, bool create_if_missing) 
{
  assert(false && "NEED TO DEFINE COMPOUND DATA TYPE TO STORE COORDS INTERLEAVED, AND DAMSEL DOESN'T HAVE THAT YET.");
  coords_dtag = dmslXcoord;
  if (!DMSLlib_ok(dmslLib, coords_dtag) && !create_if_missing) return MB_SUCCESS;

  EntityHandle tmp_handle;
  ErrorCode rval = MB_SUCCESS;
  if (!DMSLlib_ok(dmslLib, coords_dtag)) {
    tmp_handle = CREATE_HANDLE(MBMAXTYPE, 1);
    coords_dtag = dmslXcoord = DMSLtag_define(dmslModel, (damsel_handle_ptr)&tmp_handle,
                                               DAMSEL_DATA_TYPE_DOUBLE,
                                               "XYZCOORDS");
    if (!DMSLlib_ok(dmslLib, coords_dtag)) rval = MB_FAILURE;
  }
  
  CHK_MB_ERR(rval, "Failed to create interleaved coordinates tag.");
  return rval;
}

ErrorCode WriteDamsel::write_dense_tags(RangeSeqIntersectIter &rsi, damsel_id_t &ent_cont) 
{
    // all dense tags should have been initialized before this, so here we just go through
    // them and register data if there is any
  const unsigned char *val_ptr;
  ErrorCode rval = MB_SUCCESS;
  std::vector<Tag>::iterator tagit;
  std::vector<damsel_id_t>::iterator did_it;
  for (tagit = denseTags.begin(), did_it = dmslDenseTags.begin(); 
       tagit != denseTags.end(); tagit++, did_it++) {
      // get a ptr to memory for this tag/sequence
    DenseTag *dtag = dynamic_cast<DenseTag*>(*tagit);
    assert(dtag);
    rval = dtag->get_array(rsi.get_sequence(), val_ptr);
    CHK_MB_ERR_2(rval, "Failed to get tag coordinates pointer for vertices starting with handle %lu.",
               rsi.get_start_handle());

      // if ptr is NULL, no data for this tag in this sequence
    if (!val_ptr) continue;
    
    // else, register with damsel
    damsel_err_t err = DMSLtag_assign(*did_it, ent_cont, (void*)val_ptr);
    CHK_DMSL_ERR_2(err, "Failed to write coordinates tag for vertices starting with handle %lu.",
                 rsi.get_start_handle());
  }
  
  return rval;
}

ErrorCode WriteDamsel::write_dense_tags(RangeSeqIntersectIter &rsi) 
{
    // all dense tags should have been initialized before this, so here we just go through
    // them and register data if there is any
  const unsigned char *val_ptr;
  ErrorCode rval = MB_SUCCESS;
  std::vector<Tag>::iterator tagit;
  std::vector<damsel_id_t>::iterator did_it;
  damsel_err_t err;
  damsel_id_t ent_cont = DAMSEL_ID_INVALID;
  
  for (tagit = denseTags.begin(), did_it = dmslDenseTags.begin(); 
       tagit != denseTags.end(); tagit++, did_it++) {
      // get a ptr to memory for this tag/sequence
    DenseTag *dtag = dynamic_cast<DenseTag*>(*tagit);
    assert(dtag);
    rval = dtag->get_array(rsi.get_sequence(), val_ptr);
    CHK_MB_ERR_2(rval, "Failed to write coordinates tag for entities starting with handle %lu.",
               rsi.get_start_handle());

      // if ptr is NULL, no data for this tag in this sequence
    if (!val_ptr) continue;
    
    if (!DMSLlib_ok(dmslModel, ent_cont)) {
      EntityHandle starth = rsi.get_start_handle();
      ent_cont = DMSLhandle_create_sequence(dmslModel, moabHandleType, (damsel_handle_ptr)&starth, 
          rsi.get_end_handle()-starth+1, 1);
      if (!DMSLlib_ok(dmslLib, ent_cont)) 
        CHK_MB_ERR_2(MB_FAILURE, "Failure creating entity sequence when writing tag for entities starting with handle %lu.",
                   rsi.get_start_handle());
    }
    
    // else, register with damsel
    err = DMSLtag_assign(*did_it, ent_cont, val_ptr);
    CHK_DMSL_ERR_2(err, "Failed to write tag for entities starting with handle %lu.", rsi.get_start_handle());
  }
  
  return rval;
}

ErrorCode WriteDamsel::init_dense_tag_info() 
{
    // initialize allTags and tagIndices
  std::vector<Tag> tmp_tags;
  ErrorCode rval = dynamic_cast<Core*>(mbImpl)->tag_get_tags(tmp_tags);
  CHK_MB_ERR(rval, "Failed initializing all dense tags.");
  
  std::string tag_name;
  for (std::vector<Tag>::iterator vit = tmp_tags.begin(); vit != tmp_tags.end(); vit++) {
    if ((*vit)->get_storage_type() != MB_TAG_DENSE) continue;

    denseTags.push_back(*vit);

      // get a damsel id for this tag
    Tag thandle = *vit;
    damsel_id_t dtag = DMSLtag_define(dmslModel, (damsel_handle_ptr)&thandle, 
                                      moab_to_damsel_data_type[thandle->get_data_type()],
                                      tag_name.c_str());
    if (!DMSLlib_ok(dmslModel, dtag)) CHK_MB_ERR_2(MB_FAILURE, "Failure to get Damsel tag for MOAB tag %s.", 
                                                   tag_name.c_str());
    dmslDenseTags.push_back(dtag);
  }

  return MB_SUCCESS;
}

} // namespace moab

  
