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
  int num_seq_infos = 0, num_vec_infos = 0, num_tag_infos = 0, num_ent_infos = 0;
  damsel_err_t err;
  err = DMSLget_tuple_count(dmslModel, &num_containers, &num_tags);
  CHK_DMSL_ERR(err, "DMSLget_tuple_count failed.");
  err = DMSLentity_get_count(dmslModel, &num_ent_infos);
  CHK_DMSL_ERR(err, "DMSLentity_get_count failed.");
  num_coll_infos = DMSLmodel_get_collection_count(dmslModel);
  CHK_DMSL_ERR(err, "DMSLmodel_get_collection_count failed.");

  std::vector<damsel_entity_buf_type> ent_infos(num_ent_infos);
  std::vector<damsel_collection_buf_type> coll_infos(num_coll_infos);
  std::vector<damsel_tag_buf_type> tag_infos(num_seq_infos+num_vec_infos);
  err = DMSLmodel_get_entity_infos(dmslModel, &ent_infos[0]);
  CHK_DMSL_ERR(err, "Failure getting entity infos.");
  err = DMSLmodel_get_collection_infos(dmslModel, &coll_infos[0]);
  CHK_DMSL_ERR(err, "Failure getting entity infos.");
  err = DMSLmodel_get_tag_infos(dmslModel, &tag_infos[0]);
  CHK_DMSL_ERR(err, "Failure getting entity infos.");

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
  int tmp_num = num_handles / proc_size, extra = num_handles % proc_size;
  if (extra) tmp_num++;
  int my_num_handles = tmp_num;
  if (proc_rank >= extra) my_num_handles--;
  int first_ind = std::min(proc_rank,extra) * tmp_num + 
      std::max(proc_rank-extra,0) * (tmp_num-1);
  int end_ind = first_ind + my_num_handles;

    // - create moab entity sets for partition collection(s)
  EntityHandle start_handle;
  rval = mbImpl->create_entity_sets(my_num_handles, &char_tagvals[first_ind], 0, start_handle); RR;
  
    // - GET TYPE, CONTENTS OF COLLECTION CONTENTS CONTAINER
    // - allocate moab-side container (using count from container)
    // - MAP storage TO CONTAINER 
    // - EXECUTE
    // ==> have list of all handles (entities + collections) represented on this proc

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

ErrorCode parse_options(FileOptions &opts,
                        bool &parallel, 
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

                        
