/* WriteDamsel.cpp
 * The Damsel library provides mesh-aware parallel I/O; see http://trac.mcs.anl.gov/projects/damsel/wiki for details,
 * though for now that site is restricted to project participants.  Damsel uses a data model that's very similar to
 * that used in MOAB and ITAPS.  It uses the same basic data model concepts of entities, sets, tags, and interface.
 * In theory, we should be able to completely save/restore to/from Damsel any data that can be saved/restored to/from
 * our native HDF5-based reader/writer.
 *
 * Mapping between MOAB-Damsel data models
 * =======================================
 * Basic data model entities, MOAB <--> Damsel:
 *        Entity <--> Entity
 *     EntitySet <--> Collection
 *           Tag <--> Tag
 * API/data strutures:
 *          Range (n1) <--> Sequence container
 *         std::vector <--> Vector container
 *          Range (n2) <--> Tree container
 *
 * n1: single contiguous subrange
 * n2: multiple subranges
 *
 * Conventions
 * ===========
 * There are parts of MOAB data structures that need to be stored to Damsel that aren't represented in the Damsel data model,
 * e.g. dense vs. sparse storage type, set tracking flags.  
 * - We need to store these as tags in Damsel.  
 * - Since Damsel tags need to have a MOAB counterpart, we have to create those as tag data in MOAB too 
 *   (duplicating the data in the data structures, bummer).
 * - Because we may want to use these tags for multiple Damsel writes/models, we create the MOAB-side tags in the WriteDamsel
 *   constructor, not in the init_tags function that's called for every write
 * - Conventional tags have names prefixed with mbdmsl_ to avoid name conflicts with other MOAB tags.
 * Here we list the conventional tags used by MOAB's Damsel reader/writer.
 * 
 * Tag name                   Tag char's (storage type, data type, length, def val)       Values, used for what
 * --------                   -----------------------------------------------------       --------------------
 * mbdmsl_XCOORDS             dense; double[1]; 0.0                                       MOAB vertex x coordinate
 * mbdmsl_YCOORDS             dense; double[1]; 0.0                                       MOAB vertex y coordinate
 * mbdmsl_ZCOORDS             dense; double[1]; 0.0                                       MOAB vertex z coordinate
 * mbdmsl_COLL_FLAGS          sparse; char; 1; 0x0                                        bit 0: 0=set-type, 1=vector-type
 *                                                                                            1: 1=tracking, 0=not tracking
 * mbdmsl_PARENTS      sparse; handle; var;                                        (list of parent sets)
 * mbdmsl_CHILDS       sparse; handle; var;                                        (list of child sets)
 *
 * 
 * 
 */

#include "WriteDamsel.hpp"

#include "DamselUtil.hpp"
#include "damsel.h"
#include "damsel-internal.h"
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
#include "SparseTag.hpp"

namespace moab {

WriterIface* WriteDamsel::factory( Interface* iface )
  { return new WriteDamsel( iface ); }

WriteDamsel::WriteDamsel(Interface *impl) 
        : mbImpl(impl), mWriteIface(NULL), mError(NULL), sequenceManager(NULL),
          dmslLib(DAMSEL_LIBRARY_INVALID), dmslModel(DAMSEL_MODEL_INVALID),
          moabHandleType(DAMSEL_HANDLE_TYPE_INVALID)
{
  assert(impl != NULL);

  impl->query_interface( mWriteIface );
  assert(mWriteIface);
  
  sequenceManager = dynamic_cast<Core*>(impl)->sequence_manager();
  assert(sequenceManager);

  impl->query_interface(mError);
  assert(mError);
  
  ErrorCode rval = mbImpl->tag_get_handle("mbdmsl_XCOORDS", 1, MB_TYPE_DOUBLE,
                                         xcoordsTagPair.first, MB_TAG_DENSE|MB_TAG_CREAT);
  CHK_MB_ERR_NR(rval, "Failed to create_tag mbdmsl_XCOORDS.");
  tagMaps[2][xcoordsTagPair.first] = xcoordsTagPair.second = DAMSEL_TAG_INVALID;
  
  
  rval = mbImpl->tag_get_handle("mbdmsl_YCOORDS", 1, MB_TYPE_DOUBLE,
                               ycoordsTagPair.first, MB_TAG_DENSE|MB_TAG_CREAT);
  CHK_MB_ERR_NR(rval, "Failed to create_tag mbdmsl_YCOORDS.");
  tagMaps[2][ycoordsTagPair.first] = ycoordsTagPair.second = DAMSEL_TAG_INVALID;
  
  rval = mbImpl->tag_get_handle("mbdmsl_ZCOORDS", 1, MB_TYPE_DOUBLE,
                               zcoordsTagPair.first, MB_TAG_DENSE|MB_TAG_CREAT);
  CHK_MB_ERR_NR(rval, "Failed to create_tag mbdmsl_ZCOORDS.");
  tagMaps[2][zcoordsTagPair.first] = zcoordsTagPair.second = DAMSEL_TAG_INVALID;
  
  rval = mbImpl->tag_get_handle("mbdmsl_COLL_FLAGS", 1, MB_TYPE_OPAQUE,
                               collFlagsTagPair.first, MB_TAG_DENSE|MB_TAG_CREAT);
  CHK_MB_ERR_NR(rval, "Failed to create_tag mbdmsl_COLL_FLAGS.");
  tagMaps[2][collFlagsTagPair.first] = collFlagsTagPair.second = DAMSEL_TAG_INVALID;
  
  rval = mbImpl->tag_get_handle("mbdmsl_PARENTS", 1, MB_TYPE_HANDLE,
                               parentsTagPair.first, MB_TAG_DENSE|MB_TAG_CREAT|MB_TAG_VARLEN);
  CHK_MB_ERR_NR(rval, "Failed to create_tag mbdmsl_PARENTS.");
  tagMaps[2][parentsTagPair.first] = parentsTagPair.second = DAMSEL_TAG_INVALID;
  
  rval = mbImpl->tag_get_handle("mbdmsl_CHILDREN", 1, MB_TYPE_HANDLE,
                               childrenTagPair.first, MB_TAG_DENSE|MB_TAG_CREAT|MB_TAG_VARLEN);
  CHK_MB_ERR_NR(rval, "Failed to create_tag mbdmsl_CHILDREN.");
  tagMaps[2][childrenTagPair.first] = childrenTagPair.second = DAMSEL_TAG_INVALID;
  
  moabHandleType = (sizeof(EntityHandle) == 64 ? DAMSEL_HANDLE_TYPE_HANDLE64 :
                    DAMSEL_HANDLE_TYPE_HANDLE32);
  
}

WriteDamsel::~WriteDamsel() 
{
  if (mWriteIface) mbImpl->release_interface(mWriteIface);
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

  dmslLib = DMSLlib_init();

    // create a damsel model
  dmslModel = DMSLmodel_create(sizeof(EntityHandle) == 8 ? DAMSEL_HANDLE_TYPE_HANDLE64 : 
                               DAMSEL_HANDLE_TYPE_HANDLE32);
  
    // attach to a file, since we need it for creating containers
  MPI_Comm comm = MPI_COMM_WORLD;
  unlink(file_name);
  err = DMSLmodel_attach(dmslModel, file_name, comm, NULL);
  CHK_DMSL_ERR(err, "DMSLmodel_attach failed.");
  
  rval = mWriteIface->gather_entities(all_ents, meshset_list, num_sets);
  CHK_MB_ERR(rval, "Gather entities failed in WriteDamsel.");

  if (all_ents.empty()) return MB_SUCCESS;

    // create damsel tags for MOAB dense, sparse, and conventional tags
  rval = init_tag_info();
  CHK_MB_ERR(rval, NULL);

    // iterate through the groups of contiguous sequences of handles
  RangeSeqIntersectIter rsi(sequenceManager);
  rval = rsi.init(all_ents.begin(), all_ents.end());
  
  while (rval == MB_SUCCESS) {
      // write subrange of things to damsel: map handles, map entity definition data (connectivity/coords/set contents),
      // map dense tags
    rval = write_subrange(rsi);
    CHK_MB_ERR(rval, "Failed to write subrange.");
    
    rval = rsi.step();
    while (MB_ENTITY_NOT_FOUND == rval)
      rval = rsi.step();
  }

    // write sparse tags
  rval = map_sparse_tags();
  CHK_MB_ERR(rval, "Failed to write sparse tags.");

    //damsel_request_t request;
    //err = DMSLmodel_transfer_async(dmslModel, DAMSEL_TRANSFER_TYPE_WRITE, &request);
  err = DMSLmodel_transfer_sync(dmslModel, DAMSEL_TRANSFER_TYPE_WRITE);
  CHK_DMSL_ERR(err, "DMSLmodel_transfer_asynch failed.");
  
    //damsel_status_t status;
    //err = DMSLmodel_wait(request, &status);
  CHK_DMSL_ERR(err, "DMSLmodel_wait failed.");

  DMSLmodel_close(dmslModel);

  DMSLlib_finalize(dmslLib);
  
    // we should be done
  return MB_SUCCESS;
}

ErrorCode WriteDamsel::init_tag_info() 
{
    // initialize allTags and tagIndices
  std::vector<Tag> tmp_mtags;
  ErrorCode rval = mbImpl->tag_get_tags(tmp_mtags);
  CHK_MB_ERR(rval, "Failed to get all tag handles.");
  damsel_tag dtag;
  int dum_size;
  
    // define damsel tag handles for all dense/sparse tags
  for (std::vector<Tag>::iterator vit = tmp_mtags.begin(); vit != tmp_mtags.end(); vit++) {
    if (((*vit)->get_storage_type() != MB_TAG_DENSE && (*vit)->get_storage_type() != MB_TAG_SPARSE) ||
        mbImpl->tag_get_length(*vit, dum_size) == MB_VARIABLE_DATA_LENGTH ||
        dum_size != 1) {
      std::cerr << "Warning: tag " << (*vit)->get_name() 
                << "is not of type dense or sparse, and is not currently supported by the damsel writer." 
                << std::endl;
      continue;
    }
    else if (tagMaps[2].find(*vit) != tagMaps[2].end()) 
        // conventional tag - skip
      continue;
    
      // create a damsel counterpart for this tag
    Tag thandle = *vit;
    dtag = DMSLtag_define(dmslModel, (damsel_handle_ptr)&thandle, 
                          DamselUtil::mtod_data_type[(*vit)->get_data_type()],
                          (*vit)->get_name().c_str());
    if (DAMSEL_TAG_INVALID == dtag) CHK_MB_ERR_2(MB_FAILURE, "Failure to get Damsel tag for MOAB tag %s.", 
                                                 (*vit)->get_name().c_str());
    if ((*vit)->get_storage_type() == MB_TAG_DENSE) tagMaps[0][*vit] = dtag;
    else if ((*vit)->get_storage_type() == MB_TAG_SPARSE) tagMaps[1][*vit] = dtag;
  }

    // do the same for conventional tags: 
    // XCOORDS
  tagMaps[2][xcoordsTagPair.first] = xcoordsTagPair.second = 
      DMSLtag_define(dmslModel, (damsel_handle_ptr)&(xcoordsTagPair.first), 
                     DamselUtil::mtod_data_type[(xcoordsTagPair.first)->get_data_type()],
                     (xcoordsTagPair.first)->get_name().c_str());
  if (DAMSEL_TAG_INVALID == dtag) CHK_MB_ERR_2(MB_FAILURE, "Failure to get Damsel tag for MOAB tag %s.", 
                                               (xcoordsTagPair.first)->get_name().c_str());

    // YCOORDS
  tagMaps[2][ycoordsTagPair.first] = ycoordsTagPair.second = 
      DMSLtag_define(dmslModel, (damsel_handle_ptr)&(ycoordsTagPair.first), 
                     DamselUtil::mtod_data_type[(ycoordsTagPair.first)->get_data_type()],
                     (ycoordsTagPair.first)->get_name().c_str());
  if (DAMSEL_TAG_INVALID == dtag) CHK_MB_ERR_2(MB_FAILURE, "Failure to get Damsel tag for MOAB tag %s.", 
                                               (ycoordsTagPair.first)->get_name().c_str());

    // ZCOORDS
  tagMaps[2][zcoordsTagPair.first] = zcoordsTagPair.second = 
      DMSLtag_define(dmslModel, (damsel_handle_ptr)&(zcoordsTagPair.first), 
                     DamselUtil::mtod_data_type[(zcoordsTagPair.first)->get_data_type()],
                     (zcoordsTagPair.first)->get_name().c_str());
  if (DAMSEL_TAG_INVALID == dtag) CHK_MB_ERR_2(MB_FAILURE, "Failure to get Damsel tag for MOAB tag %s.", 
                                               (zcoordsTagPair.first)->get_name().c_str());

    // COLL_FLAGS
  tagMaps[2][collFlagsTagPair.first] = collFlagsTagPair.second = 
      DMSLtag_define(dmslModel, (damsel_handle_ptr)&(collFlagsTagPair.first), 
                     DamselUtil::mtod_data_type[(collFlagsTagPair.first)->get_data_type()],
                     (collFlagsTagPair.first)->get_name().c_str());
  if (DAMSEL_TAG_INVALID == dtag) CHK_MB_ERR_2(MB_FAILURE, "Failure to get Damsel tag for MOAB tag %s.", 
                                               (collFlagsTagPair.first)->get_name().c_str());

    /*
      SKIP PARENTS/CHILDREN FOR NOW, UNTIL WE HAVE VAR LENGTH TAGS IN DAMSEL

    // PARENTS
  parentsTagPair.second = DMSLtag_define(dmslModel, (damsel_handle_ptr)&(collFlagsTagPair.first), 
                                           DamselUtil::mtod_data_type[(collFlagsTagPair.first)->get_data_type()],
                                           (parentsTagPair.first)->get_name().c_str());
  if (DAMSEL_TAG_INVALID == dtag) CHK_MB_ERR_2(MB_FAILURE, "Failure to get Damsel tag for MOAB tag %s.", 
                                               (parentsTagPair.first)->get_name());

    // CHILDREN
  childrenTagPair.second = DMSLtag_define(dmslModel, (damsel_handle_ptr)&(collFlagsTagPair.first), 
                                           DamselUtil::mtod_data_type[(collFlagsTagPair.first)->get_data_type()],
                                           (childrenTagPair.first)->get_name().c_str());
  if (DAMSEL_TAG_INVALID == dtag) CHK_MB_ERR_2(MB_FAILURE, "Failure to get Damsel tag for MOAB tag %s.", 
                                               (childrenTagPair.first)->get_name());

    */

    // map the tag handles in one big call
  int num_tags = tagMaps[0].size() + tagMaps[1].size() + tagMaps[2].size();
  std::vector<Tag> moab_taghs;
  moab_taghs.reserve(num_tags);
  for (int i = 0; i < 3; i++) {
    for (std::map<Tag, damsel_tag>::iterator mit = tagMaps[i].begin(); mit != tagMaps[i].end(); mit++) {
      if (!(*mit).first || (*mit).second == DAMSEL_TAG_INVALID) continue;
      moab_taghs.push_back((*mit).first);
    }
  }
    
  damsel_container mtags = DMSLcontainer_create_vector(dmslModel, (damsel_handle_ptr)&moab_taghs[0], moab_taghs.size());
  std::cerr << "MOAB: created model container: mtags = " << mtags <<std::endl;
  
  damsel_err_t err = DMSLmodel_map_handles_inventing_file_handles(mtags);
  CHK_DMSL_ERR(err, "Failed to map tag handles.");

  err = DMSLcontainer_release(mtags);
  CHK_DMSL_ERR(err, "Problem releasing tag handle container.");
  
  return MB_SUCCESS;
}

ErrorCode WriteDamsel::write_vertices(RangeSeqIntersectIter &rsi) 
{
    // write the vertices; these vertices will be in the same sequence and will be contiguous, guaranteed
  EntityHandle start_vert = rsi.get_start_handle(), end_vert = rsi.get_end_handle();

    // create a damsel container for these vertex handles
  damsel_container vertex_cont = DMSLcontainer_create_sequence(dmslModel, start_vert, (int)(end_vert-start_vert+1), 1);
  std::cerr << "MOAB: created model container: vertex_cont = " << vertex_cont <<std::endl;
  if (DAMSEL_CONTAINER_INVALID == vertex_cont) 
    CHK_MB_ERR_2(MB_FAILURE, "Failed to create vertex sequence for vertices starting with handle %lu.", 
               rsi.get_start_handle());
  
  damsel_err_t err = DMSLmodel_map_handles_inventing_file_handles(vertex_cont);
  CHK_DMSL_ERR(err, "Failed to map handles.");

    // define the entities to damsel
  err = DMSLentity_define(vertex_cont, DAMSEL_ENTITY_TYPE_VERTEX, 1, vertex_cont);
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
  
  if (xcoords && !ycoords && !zcoords) {
      // interleaved

      // map the data to damsel
    err = DMSLmodel_map_tag(xcoords, vertex_cont, (damsel_handle_ptr)&xcoordsTagPair.first);
    CHK_DMSL_ERR_2(err, "Failed to assign vertex coordinates tag for vertices starting with handle %lu.", 
                 rsi.get_start_handle());
  }
  else {
      // map the data to damsel
    err = DMSLmodel_map_tag(xcoords, vertex_cont, (damsel_handle_ptr)&xcoordsTagPair.first);
    CHK_DMSL_ERR_2(err, "Failed to assign vertex x coordinates tag for vertices starting with handle %lu.", 
                   rsi.get_start_handle());
    err = DMSLmodel_map_tag(ycoords, vertex_cont, (damsel_handle_ptr)&ycoordsTagPair.first);
    CHK_DMSL_ERR_2(err, "Failed to assign vertex y coordinates tag for vertices starting with handle %lu.", 
                   rsi.get_start_handle());
    err = DMSLmodel_map_tag(zcoords, vertex_cont, (damsel_handle_ptr)&zcoordsTagPair.first);
    CHK_DMSL_ERR_2(err, "Failed to assign vertex z coordinates tag for vertices starting with handle %lu.", 
                   rsi.get_start_handle());
  }

    // write/map dense tags
  rval = map_dense_tags(rsi, vertex_cont);
  CHK_MB_ERR(rval, NULL);
  
  err = DMSLcontainer_release(vertex_cont);
  CHK_DMSL_ERR(err, "Problem releasing vertex handle container.");

  return MB_SUCCESS;
}

ErrorCode WriteDamsel::write_entities(RangeSeqIntersectIter &rsi) 
{
    // write the entities; these entities will be in the same sequence and will be contiguous, guaranteed
  EntityHandle start_ent = rsi.get_start_handle(), end_ent = rsi.get_end_handle();

    // create a damsel container for these entity handles
  damsel_container ent_cont;
  ent_cont = DMSLcontainer_create_sequence(dmslModel, start_ent, (int)(end_ent-start_ent+1), 1);
  std::cerr << "MOAB: created model container: ent_cont = " << ent_cont <<std::endl;
  if (DAMSEL_CONTAINER_INVALID == ent_cont)
    CHK_MB_ERR(MB_FAILURE, "Bad sequence returned by Damsel.");

  damsel_err_t err = DMSLmodel_map_handles_inventing_file_handles(ent_cont);
  CHK_DMSL_ERR(err, "Failed to map handles.");
  
    // get # verts per entity and entity type
  EntityType etype = mbImpl->type_from_handle(start_ent);
  assert(MBMAXTYPE != etype);
  int num_connect = rsi.get_sequence()->values_per_entity();
  assert(0 < num_connect);
  
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
  
    // define the entities to damsel
  err = DMSLentity_define_fast(ent_cont, DamselUtil::mtod_entity_type[etype], num_connect, (damsel_handle*)connect);
  CHK_DMSL_ERR_2(err, "DMSLentity_define failed for entities starting with handle %lu.", rsi.get_start_handle());
  
    // write dense tags
  rval = map_dense_tags(rsi, ent_cont);
  CHK_MB_ERR(rval, NULL);
  
  err = DMSLcontainer_release(ent_cont);
  CHK_DMSL_ERR(err, "Problem releasing entity handle container.");

  return MB_SUCCESS;
}

ErrorCode WriteDamsel::map_dense_tags(RangeSeqIntersectIter &rsi, damsel_container &ent_cont) 
{
    // all dense_tags have been initialized before this, so here we just go through
    // them and map data if there is any
  const unsigned char *val_ptr;
  ErrorCode rval = MB_SUCCESS;
  std::map<Tag, damsel_tag>::iterator tagit;
  damsel_err_t err;
  for (tagit = tagMaps[0].begin(); tagit != tagMaps[0].end(); tagit++) {
      // get a ptr to memory for this tag/sequence
    DenseTag *dtag = dynamic_cast<DenseTag*>((*tagit).first);
    assert(dtag);
    rval = dtag->get_array(rsi.get_sequence(), val_ptr);
    CHK_MB_ERR_2(rval, "Failed to get tag coordinates pointer for vertices starting with handle %lu.",
               rsi.get_start_handle());

      // if ptr is NULL, no data for this tag in this sequence
    if (!val_ptr) continue;
    
      // else, register with damsel
    err = DMSLmodel_map_tag((void*)val_ptr, ent_cont, (damsel_handle_ptr)&dtag);
    CHK_DMSL_ERR_2(err, "Failed to write coordinates tag for vertices starting with handle %lu.",
                   rsi.get_start_handle());
  }
  
  return rval;
}

ErrorCode WriteDamsel::map_sparse_tags() 
{
    // all sparse_tags have been initialized before this, so here we just go through
    // them and map data if there is any
  ErrorCode rval = MB_SUCCESS;
  std::map<Tag, damsel_tag>::iterator tagit;
  damsel_err_t err;
  std::vector<unsigned char> tag_values;
  std::vector<EntityHandle> tagged_ents;
  damsel_container ent_cont;
  for (tagit = tagMaps[1].begin(); tagit != tagMaps[1].end(); tagit++) {
      // get a ptr to memory for this tag/sequence
    SparseTag *stag = dynamic_cast<SparseTag*>((*tagit).first);
    assert(stag);
    Range output_ents;
    rval = stag->get_tagged_entities(sequenceManager, output_ents);
    CHK_MB_ERR_2(rval, "Trouble getting tagged entities for tag %s.", stag->get_name().c_str());

      // if no entities have this tag set, don't map it
    if (output_ents.empty()) continue;
    
      // else, register with damsel
      // allocate space for and get values
    tag_values.resize(stag->get_size() * output_ents.size());
    rval = mbImpl->tag_get_data(stag, output_ents, &tag_values[0]);
    CHK_MB_ERR_2(rval, "Trouble getting tag values for tag %s.", stag->get_name().c_str());

      // build a vector of entity handles from the range, and a container from that
    tagged_ents.resize(output_ents.size());
    std::copy(output_ents.begin(), output_ents.end(), tagged_ents.begin());
    ent_cont = DMSLcontainer_create_vector(dmslModel, (damsel_handle_ptr)&tagged_ents[0], tagged_ents.size());
    std::cerr << "MOAB: created model container: sparse_tag_ent_cont = " << ent_cont <<std::endl;
    if (ent_cont == DAMSEL_CONTAINER_INVALID) 
      CHK_MB_ERR_2(MB_FAILURE, "Trouble creating entity handle container for tag %s.", stag->get_name().c_str());

      // now map it
    err = DMSLmodel_map_tag((void*)&tag_values[0], ent_cont, (damsel_handle_ptr)&stag);
    CHK_DMSL_ERR_2(err, "Failed to write tag %s.", stag->get_name().c_str());

    err = DMSLcontainer_release(ent_cont);
    CHK_DMSL_ERR(err, "Problem releasing entity handle container.");
  }
  
  return rval;
}

ErrorCode WriteDamsel::write_sets(RangeSeqIntersectIter &rsi) 
{
    // write the sets
  ErrorCode rval = MB_SUCCESS;
  std::vector<EntityHandle> ents;
  damsel_container mcont;
  damsel_err_t err;
  unsigned int i, num_sets = rsi.get_end_handle() - rsi.get_start_handle() + 1;
  std::vector<unsigned char> set_flags(num_sets);
  EntityHandle seth;
  for (seth = rsi.get_start_handle(), i = 0; seth <= rsi.get_end_handle(); seth++, i++) {
      // get all the entities in the set
    ents.clear();
    rval = mbImpl->get_entities_by_handle(seth, ents);
    CHK_MB_ERR_2(rval, "get_entities_by_handle failed for set %lu.", seth);
    if (!ents.empty()) {
      mcont = DMSLcontainer_create_vector(dmslModel, (damsel_handle*)&ents[0], ents.size());
    }
    else {
      mcont = DMSLcontainer_create_vector(dmslModel, (damsel_handle*)NULL, 0);
    }
    std::cerr << "MOAB: created model container: sets_cont = " << mcont <<std::endl;

    // get the set type (range or set)
    unsigned int opts;
    rval = mbImpl->get_meshset_options(seth, opts);
    CHK_MB_ERR_2(rval, "Failed to get options for meshset %lu.", seth);
    damsel_collection_type coll_type = (opts&MESHSET_SET ? DAMSEL_HANDLE_COLLECTION_TYPE_SET :
                     DAMSEL_HANDLE_COLLECTION_TYPE_VECTOR);

      // parents/children...

      // set flags
    if (opts & MESHSET_TRACK_OWNER)
      set_flags[i] |= MESHSET_TRACK_OWNER;
    else
      set_flags[i] |= 0;

    err = DMSLcontainer_release(mcont);
    CHK_DMSL_ERR(err, "Problem releasing set entity handle container.");
  }

    // set the COLL_FLAGS tag, using assign (direct)
    // make a container of set handles...
  mcont = DMSLcontainer_create_sequence(dmslModel, rsi.get_start_handle(), num_sets, 1);
  std::cerr << "MOAB: created model container: sets_cont = " << mcont <<std::endl;
    // assign the tags on them
  err = DMSLmodel_map_tag(&set_flags[0], mcont, (damsel_handle_ptr)&(collFlagsTagPair.first));
  CHK_DMSL_ERR(err, "Failed to assign COLL_FLAGS tag for sets.");

  err = DMSLmodel_map_handles_inventing_file_handles(mcont);
  CHK_DMSL_ERR(err, "Failed to map set handles.");
  
    // map other dense tags
  rval = map_dense_tags(rsi, mcont);
  CHK_MB_ERR(rval, "Failed to map dense tags for sets.");

  return rval;
}


} // namespace moab

  
