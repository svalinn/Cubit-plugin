#include "ReadDamsel.hpp"

#include "assert.h"
#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/Error.hpp"
#include "moab/ReadUtilIface.hpp"
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
do if (MB_SUCCESS != (A   )) { \
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
do if (DMSL_OK.id != A.id) {             \
mError->set_last_error(B);\
return error(MB_FAILURE);             \
} while(false)

#define CHK_DMSL_ERR_2( A, B, C )            \
do if (DMSL_OK.id != A.id) {             \
mError->set_last_error(B, C);            \
return error(MB_FAILURE);             \
} while(false)

#define CHK_DMSL_ERR_FINALIZE( A, B )        \
do if (DMSL_OK.id != A.id) {             \
  DMSLlib_finalize(dmslLib); \
  dmslLib = 0;          \
  mError->set_last_error(B);\
  return error(MB_FAILURE);                    \
} while(false)


ReadDamsel::ReadDamsel() 
{
  init();
}

ErrorCode ReadDamsel::init()
{
  mdbImpl->query_interface(readMeshIface);
}

ErrorCode ReadDamsel::parse_options(FileOptions &opts,
				    bool &parallel) 
{
    // Handle parallel options
  std::string junk;
  bool use_mpio = (MB_SUCCESS == opts.get_null_option("USE_MPIO"));
  rval = opts.match_option("PARALLEL", "READ_PART");
  bool parallel = (rval != MB_ENTITY_NOT_FOUND);
  nativeParallel = (rval == MB_SUCCESS);
  if (use_mpio && !parallel) {
    readUtil->report_error( "'USE_MPIO' option specified w/out 'PARALLEL' option" );
    return MB_NOT_IMPLEMENTED;
  }

  return MB_SUCCESS;
}

// ASSUMPTIONS:
// Partition collection is a *flat* collection of handles for entities and other collections that
// will be represented on a part

ErrorCode ReadDamsel::load_file( const char* filename, 
                                 const EntityHandle* file_set, 
                                 const FileOptions& opts,
                                 const ReaderIface::SubsetList* subset_list,
                                 const Tag* file_id_tag )
{
  ErrorCode rval;
 
  rval = process_options( filename, opts );
  if (MB_SUCCESS != rval)
    return rval;

    // initialize damsel
  dmslLib = DMSLlib_init();
  
    // figure out handle size
  moab_to_damsel_data_type[MB_TYPE_OPAQUE] = DAMSEL_DATA_TYPE_BYTES;
  moab_to_damsel_data_type[MB_TYPE_INTEGER] = DAMSEL_DATA_TYPE_INTEGER;
  moab_to_damsel_data_type[MB_TYPE_DOUBLE] = DAMSEL_DATA_TYPE_DOUBLE;
  moab_to_damsel_data_type[MB_TYPE_BIT] = DAMSEL_DATA_TYPE_INVALID;
  moab_to_damsel_data_type[MB_TYPE_HANDLE] = DAMSEL_DATA_TYPE_HANDLE;
  
    // create a damsel model
  dmslModel = DMSLmodel_create(sizeof(EntityHandle) == 8 ? DAMSEL_HANDLE_TYPE_HANDLE64 : 
                               DAMSEL_HANDLE_TYPE_HANDLE32);
  
    // model attach - need model id from make model, filename
  int proc_rank = 0, proc_size = 1;
#ifdef USE_MPI
  MPI_Comm comm = 0;
  if (nativeParallel) {
    proc_rank = myPcomm->proc_config().proc_rank();
    proc_size = myPcomm->proc_config().proc_size();
    comm = myPcomm->proc_config().proc_comm();
  }
#endif

  err = DMSLmodel_attach(dmslModel, fileName, comm, NULL);
  CHK_DMSL_ERR(err, "DMSLmodel_attach failed.");
  
    // STEP 0: GET COLLECTION, TAG, ENTITY INFOS FOR GLOBAL MODEL
  int num_containers = 0, num_tag_infos = 0, num_ent_infos = 0;
  damsel_err_t err;
  err = DMSLget_tuple_count(dmslModel, &num_containers, &num_tag_infos);
  CHK_DMSL_ERR(err, "DMSLget_tuple_count failed.");
  err = DMSLentity_get_count(dmslModel, &num_ent_infos);
  CHK_DMSL_ERR(err, "DMSLentity_get_count failed.");
  num_coll_infos = DMSLmodel_get_collection_count(dmslModel);
  CHK_DMSL_ERR(err, "DMSLmodel_get_collection_count failed.");

  std::vector<damsel_entity_buf_type> ent_infos(num_ent_infos);
  std::vector<damsel_collection_buf_type> coll_infos(num_coll_infos);
  std::vector<damsel_tag_buf_type> tag_infos(num_tag_infos);
  std::vector<damsel_container_buf_type> cont_infos(num_containers);
  err = DMSLmodel_get_entity_infos(dmslModel, &ent_infos[0]);
  CHK_DMSL_ERR(err, "Failure getting entity infos.");
  err = DMSLmodel_get_collection_infos(dmslModel, &coll_infos[0]);
  CHK_DMSL_ERR(err, "Failure getting collection infos.");
  err = DMSLmodel_get_tag_infos(dmslModel, &tag_infos[0]);
  CHK_DMSL_ERR(err, "Failure getting tag infos.");
  err = DMSLmodel_get_container_infos(dmslModel, &cont_infos[0]);
  CHK_DMSL_ERR(err, "Failure getting container infos.");

    // create MOAB-side tags for all damsel tags except pre-defined ones
  rval = process_tags(tag_infos); RR;
  
  if (parallel) {
      // STEP 1: GET COLLECTION(S) REPRESENTING PARTITION: 
      // input: tag name, optionally value; 
      // output: container with file-side handles of collections satisfying those criteria
      // - get all handles/values for tag
      // - select handles matching criteria for tag value (will be collection handles)
    std::string partn_tag_name("PARALLEL_PARTITION");
    int *partn_vals = NULL;
    damsel_handle partn_tag = DMSLselect_tag_by_name(dmslModel, partn_tag_name.c_str());
      // get all the parts with that tag regardless of value
    damsel_container part_handles = DLSLselect_handles_with_values(dmslModel, partn_tag);

      // STEP 2: GET HANDLES FOR TAGS WE NEED FOR THIS READER:
      // - "SET_CHARACTERISTIC"
    damsel_handle setchar_tag = DMSLselect_tag_by_name(dmslModel, "SET_CHARACTERISTIC");
      // - "PARENT_LIST"
    damsel_handle plist_tag = DMSLselect_tag_by_name(dmslModel, "PARENT_LIST");
      // - "CHILD_LIST"
    damsel_handle clist_tag = DMSLselect_tag_by_name(dmslModel, "CHILD_LIST");

      // STEP 3: GET VALUES FOR "SET_CHARACTERISTIC" TAG ON PARTITION COLLECTIONS,
      //         GET VECTOR- OR SET-TYPE FLAGS FOR PARTITION COLLECTIONS
      // (gives tracking flag for moab)
    int num_handles = DMSLcont_get_count(dmslModel, part_handles);
    std::vector<unsigned> char_tagvals(num_handles);
      // map the set chars
    err = DMSLmodel_map_tag(&char_tagvals[0], part_handles, setchar_tag); OK;
      // execute the transfer
    err = DMSLmodel_transfer_sync(dmslModel, DAMSEL_TRANSFER_TYPE_READ); OK;

      // STEP 4: READ/PROCESS PARTITION COLLECTION(S)
      // decide the handles I am responsible using round-robin for now
    // - GET TYPE, CONTENTS OF COLLECTION CONTENTS CONTAINER
    // - allocate moab-side container (using count from container)
    // - MAP storage TO CONTAINER 
    // - EXECUTE
    // ==> have list of all handles (entities + collections) represented on this proc

    int tmp_num = num_handles / proc_size, extra = num_handles % proc_size;
    if (extra) tmp_num++;
    int my_num_handles = tmp_num;
    if (proc_rank >= extra) my_num_handles--;
    int first_ind = std::min(proc_rank,extra) * tmp_num + 
        std::max(proc_rank-extra,0) * (tmp_num-1);
    int end_ind = first_ind + my_num_handles;

      // - create moab entity sets for partition collection(s)
    EntityHandle start_handle;
    rval = readMeshIface->create_entity_sets(my_num_handles, &char_tagvals[first_ind], 0, start_handle); RR;
  }
  else {
      // initialize just by entity; each call to process_ent_info will:
      // a. create moab-side representation to read into
      // b. map those handles to damsel handles
      // c. map coords / connectivity storage to damsel equivalent
      // d. for each tag, map moab storage to damsel storage
    std::vector<damsel_entity_buf_type>::iterator eiit =  ent_infos.begin();
    for (; eiit != ent_enfos.end(); eiit++) {
      rval = process_ent_info(*eiit); RR;
    }
  }

    // process collections
  rval = process_coll_infos(coll_infos); RR;

    // STEP 5: process into list of local info structs, each represents file-side struct and
    // portion of that struct
    // ASSUMPTION: each local info struct represents single entity type & # vertices or collection
  
    // STEP 6: For each local info struct:

    // STEP 6b: READ CONTAINER INTO LOCAL BUFFER
    // STEP 6c: create app representation of entities/vertices/collection, and damsel container for them,
    //    and MAP APP HANDLE CONTAINER TO DAMSEL CONTAINER
    // STEP 6d: process vertices/entities/collection
    //    6d1: if vertices, continue
    //    6d2: if entities:
    //    - MAP LOCAL CONNECTIVITY REP'N TO DAMSEL (might be tag, don't know yet)
    //    6d3: if collection:
    //    - (everything in STEP 4 for this collection except for EXECUTE)
    //    - might need to filter out things not represented on this rank
    //    6d4: if sparse tag:
    //    - create app-side representation of sparse tag
    //    - READ CONTAINER OF MODEL HANDLES INTO LOCAL BUFFER
    //    - allocate app-side storage for tag values
    //    - MAP APP STORAGE TO MODEL TAG + (implicit?) CONTAINER

    // STEP 6e: process dense tags for the local info struct; for each dense tag:
    //   - get app tag handle for model tag handle
    //   - get app storage for app tag handle + container
    //   - MAP APP STORAGE TO MODEL TAG + CONTAINER

    // STEP 7: EXECUTE
    //   - assign all mapped data
    //   - translate all damsel handles to app handles
    // uninit damsel

}

ErrorCode ReadDamsel::insert_in_id_map(EntityHandle starth, int num, damsel_container cont, unsigned incr_size) 
{
  int ind = 0, count = DMSLcontainer_count(cont);
  while (ind < count) {
    long start = DMSLcontainer_handle_at_position(cont, ind), end = start+incr_size;
    while (ind < count-1 && DMSLcontainer_handle_at_position(cont, ind+1) == end+incr_size) 
      ind++, end+=incr_size;
    idMap.insert(start, end-start, starth);
    starth += end-start;
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadDamsel::process_ent_info(const damsel_entity_buf_type &einfo) 
{
    // create this chunk of entities
  EntityHandle *connect, start_handle;
  ErrorCode rval;
  damsel_err_t err;
  damsel_container app_cont;
  Range these_ents;
  
  if (einfo.entity_type != DAMSEL_VERTEX) {
      // create the moab entities
    rval = readMeshIface->get_element_connect(einfo.count, einfo.vertices_per_entity,
                                              moab_to_damsel_entity_type[einfo.entity_type],
                                              0, start_handle, connect);
    RR;
    these_ents.insert(start_handle, start_handle+einfo.count-1);

      // create an app-side sequence and map to file-side container
    app_cont = DMSLhandle_create_sequence(dmslModel, einfo.count, start_handle, 1);
    err = DMSLmodel_map_handles(app_cont, einfo.entity_container);
    if (DMSL_OK != err) return MB_FAILURE;

      // map connectivity
    err = DMSLmodel_map_connectivity(connect, app_cont); OK;
  }
  
  else {
      // get the number of coordinate arrays
    int num_ctags = 0;
    damsel_handle xcoord_dtag = DMSL_select_tag_by_name(dmslModel, "XCOORDS");
    if (xcoord_dtag) num_ctags++;
    damsel_handle ycoord_dtag = DMSL_select_tag_by_name(dmslModel, "YCOORDS");
    if (ycoord_dtag) num_ctags++;
    damsel_handle zcoord_dtag = DMSL_select_tag_by_name(dmslModel, "ZCOORDS");
    if (zcoord_dtag) num_ctags++;
    
      // should have one vertex per entity
    assert(einfo.vertices_per_entity == 1);
    std::vector<double*> coord_arrays;
    rval = readMeshIface->get_node_coords(num_ctags, einfo.count, 0, start_handle, coord_arrays);
    RR;

    these_ents.insert(start_handle, start_handle+einfo.count-1);

      // create an app-side sequence and map to file-side container
    app_cont = DMSLhandle_create_sequence(dmslModel, einfo.count, start_handle, 1);
    err = DMSLmodel_map_handles(app_cont, einfo.entity_container);
    if (DMSL_OK != err) return MB_FAILURE;

      // map the coords storage
    if (xcoord_dtag != 0) {err = DMSLmodel_map_tag(coord_arrays[0], app_cont, xcoord_dtag); OK;}
    if (ycoord_dtag != 0) {err = DMSLmodel_map_tag(coord_arrays[1], app_cont, ycoord_dtag); OK;}
    if (zcoord_dtag != 0) {err = DMSLmodel_map_tag(coord_arrays[2], app_cont, zcoord_dtag); OK;}
  }

    // save mapping from moab entity to einfo
  entityMap[start_handle] = einfo;

  rval = process_entity_tags(einfo.tag_count, einfo.tag_handle_container, app_cont, these_ents);
  
  return rval;
}

ErrorCode ReadDamsel::process_entity_tags(int count, damsel_container tag_container, 
                                          damsel_container app_cont, Range &these_ents) 
{
    // process tags on these entities
  for (int i = 0; i < tag_count; i++) {
    damsel_handle dtagh = DMSLcontainer_handle_at_position(tag_container, i);

      // don't do coords tags here, was done earlier
    if (xcoordDtag == dtagh || ycoordDtag == dtagh || zcoordDtag == dtagh) continue;
        
    Tag tagh = tagMap[dtagh];
    assert(tagh);
    void *tag_data;
    int ecount = these_ents.size();
    rval = mbImpl->tag_iterate(tagh, these_ents.begin(), these_ents.end(), ecount, tag_data); RR;
    assert(ecount == these_ents.size());
    err = DMSLmodel_map_tag(tag_data, app_cont, &dtagh); OK;
  }

  return MB_SUCCESS;
}
  
ErrorCode ReadDamsel::process_tags(std::vector<damsel_tag_buf_type> &tag_infos) 
{
  Tag tagh;
  ErrorCode rval = MB_SUCCESS, tmp_rval;
  for (std::vector<damsel_tag_buf_type>::iterator tit = tag_infos.begin(); tit != tag_infos.end(); tit++) {
    if (!strncmp(*tit.name, "##", 2)) {
        // predefined tag name, store the handle
      if (!strcmp(*tit.name, "##child_list")) childListTag = *tit.tag_handle;
      else if (!strcmp(*tit.name, "##parent_list")) parentListTag = *tit.tag_handle;
      else {
        rval = MB_FAILURE;
        continue;
      }
    }
    
    else if (damsel_to_moab_data_type(*tit.tag_datatype) == MB_TYPE_OPAQUE) {
      rval = MB_FAILURE;
      tagMap[*tit.tag_handle] = 0;
      continue;
    }
      
    tmp_rval = mbImpl->tag_get_handle(*tit.name, 1, damsel_to_moab_data_type[*tit.tag_datatype],
                                      tagh, MB_TAG_CREAT || MB_TAG_DENSE);
    if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
    else tagMap[*tit.tag_handle] = tagh;
  }
  
  return rval;
}

ErrorCode ReadDamsel::process_coll_infos(std::vector<collection_info_buf_type> &coll_infos) 
{
  ErrorCode rval = MB_SUCCESS;
  EntityHandle seth;
  for (std::vector<collection_info_buf_type>::iterator cit = coll_infos.begin(); cit != coll_infos.end(); cit++) {
      // make the set
    tmp_rval = mbImpl->create_meshset((*cit.type ==  DAMSEL_HANDLE_COLLECTION_TYPE_SET ? MESHSET_SET : MESHSET_ORDERED),
                                  seth);
    if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
      // make datastructures to pass things to process_entity_tags
    Range tmp_range(seth);
    damsel_container ch = DMSLhandle_create_sequence(dmslModel, 1, seth, 1);

      // get the tags on this set
    tmp_rval = process_entity_tags(*cit.tag_count, *cit.tag_container, ch, tmp_range);
    if (MB_SUCCESS != tmp_rval) rval = tmp_rval;

      // process the set contents
    if (*cit.type == DAMSEL_HANDLE_COLLECTION_TYPE_SET) {
      Range ents;
      tmp_rval = get_container_range(*cit.contents, ents);
      if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
      else if (!ents.empty()) {
	tmp_rval = mbImpl->add_entities(seth, ents);
	if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
      }
    }
    else {
      std::vector<EntityHandle> ents;
      tmp_rval = get_container_range(*cit.contents, ents);
      if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
      else if (!ents.empty()) {
	tmp_rval = mbImpl->add_entities(seth, &ents[0], ents.size());
	if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
      }
    }
  }
  
  return rval;
}

ErrorCode ReadDamsel::get_contents(damsel_model m, damsel_container c, Range &ents) {
  EntityHandle eh;
  if (DMSLcontainer_get_type(c) == DAMSEL_CONTAINER_TYPE_SEQUENCE) {
    damsel_handle start;
    size_t count, stride;
    damsel_err_t err = DMSLcontainer_sequence_get_contents(m, c, &start, &count, &stride); OK;
    if (stride == 1) {
      while (count) {
	// get start in rangemap
	RangeMap::iterator beg = entityMap.lower_bound(start);
	if (beg == entityMap.end()) return MB_SUCCESS;
	unsigned int diff = std::max(*beg.begin-start, 0);
	unsigned int num = std::min(count-diff, *beg.count);
	ents.insert(*beg.begin+diff, *beg.begin+diff+num-1);
	count -= (diff + num);
	beg++;
      }
    }
    else {
      for (int i = count-1; i >= 0; i--) {
	if (entityMap.find(start+i, eh)) ents.insert(eh);
      }
    }
  }
  else if {DMSLcontainer_get_type(c) == DAMSEL_CONTAINER_TYPE_VECTOR) {
    EntityHandle *handle_ptr;
    size_t count;
    damsel_err_t err = DMSLcontainer_vector_get_contents(m, c, &handle_ptr, &count); OK;
    for (int i = count-1; i >= 0; i--) {
      if (entityMap.find(handle_ptr[i], eh)) ents.insert(eh);
    }
  }
  else if {DMSLcontainer_get_type(c) == DAMSEL_CONTAINER_TYPE_TREE) {
    damsel_handle_ptr node_ptr = NULL;
    damsel_container cont = NULL;
    while (DMSLcontainer_tree_get_contents(m, c, &node_ptr, &cont) == DMSL_OK &&
	   cont) {
      rval = get_contents(m, c, ents);
      if (MB_SUCCESS != rval) return rval;
    }
  }

  return MB_SUCCESS;
}
