#ifndef READCCMIO_HPP
#define READCCMIO_HPP

#ifndef IS_BUILDING_MB
#error "ReadCCMIO.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <set>
#include <map>
#include <string>

#include "MBForward.hpp"
#include "MBRange.hpp"
#include "ExoIIInterface.hpp"
#include "MBReaderIface.hpp"
#include "TupleList.hpp"
#include "ccmio.h"

#define TupleList std::map<int, std::vector<MBEntityHandle> > 
#define SenseList std::map<int, std::vector<int> > 
#undef TUPLE_LIST

class MBReadUtilIface;

class MB_DLL_EXPORT ReadCCMIO : public MBReaderIface
{
 
public:

    //! Constructor
  ReadCCMIO(MBInterface *impl);

   //! Destructor
  virtual ~ReadCCMIO();
  
  static MBReaderIface* factory( MBInterface* );

  MBErrorCode load_file(const char *file_name,
                        const MBEntityHandle* file_set,
                        const FileOptions& opts,
                        const MBReaderIface::IDTag* subset_list,
                        int subset_list_length,
                        const MBTag* file_id_tag);

private:
  
  MBErrorCode read_processor(CCMIOID rootID, CCMIOID problemID,
                             CCMIOID processorID, CCMIOSize_t proc,
                             MBRange *new_ents);


  MBErrorCode read_cells(CCMIOSize_t proc, CCMIOID processorID,
                         CCMIOID verticesID, CCMIOID topologyID,
                         CCMIOID solutionID, bool has_solution,
                         TupleList &vert_map, MBRange *new_cells);


  MBErrorCode construct_cells(TupleList &face_map, 
#ifndef TUPLE_LIST
                              SenseList &sense_map,
#endif
                              TupleList &vert_map, 
                              std::vector<MBEntityHandle> &new_cells);


  MBErrorCode create_cell_from_faces(std::vector<MBEntityHandle> &facehs,
                                     std::vector<int> &senses,
                                     MBEntityHandle &cell);

  MBErrorCode read_gids_and_types(CCMIOID problemID,
                                  CCMIOID topologyID,
                                  std::vector<MBEntityHandle> &cells);

  MBErrorCode read_all_faces(CCMIOID topologyID, TupleList &vert_map, 
                             TupleList &face_map
#ifndef TUPLE_LIST
                             ,SenseList &sense_map
#endif
                             , MBRange *new_faces);


  MBErrorCode read_faces(CCMIOID faceID, CCMIOEntity bdy_or_int,
                         TupleList &vert_map,
                         TupleList &face_map
#ifndef TUPLE_LIST
                         ,SenseList &sense_map
#endif
                         , MBRange *new_faces);

  MBErrorCode make_faces(int *farray, 
                         TupleList &vert_map,
                         std::vector<MBEntityHandle> &new_faces);

  MBErrorCode read_vertices(CCMIOSize_t proc, CCMIOID processorID, CCMIOID verticesID,
                            CCMIOID topologyID, CCMIOID solutionID, bool has_solution,
                            MBRange *verts, TupleList &vert_map);


  MBErrorCode get_processors(CCMIOID stateID, CCMIOID &processorID,
                             std::set<CCMIOSize_t> &procs);


  MBErrorCode get_state(CCMIOID rootID, CCMIOID &problemID, CCMIOID &stateID);


  virtual MBErrorCode read_tag_values( const char* file_name,
                                       const char* tag_name,
                                       const FileOptions& opts,
                                       std::vector<int>& tag_values_out,
                                       const IDTag* subset_list = 0,
                                       int subset_list_length = 0 );
  
  MBErrorCode load_matset_data(CCMIOID problemID);
  
  MBErrorCode load_metadata(CCMIOID rootID, CCMIOID problemID, 
                            const MBEntityHandle *file_set);
  
  MBErrorCode create_matset_tags(MBTag &matNameTag, MBTag &matPorosityTag, 
                                 MBTag &matSpinTag, MBTag &matGroupTag);

    //! Cached tags for reading.  Note that all these tags are defined when the
    //! core is initialized.
  MBTag mMaterialSetTag;
  MBTag mDirichletSetTag;
  MBTag mNeumannSetTag;
  MBTag mHasMidNodesTag;
  MBTag mGlobalIdTag;
  
  MBInterface *mbImpl;

  MBReadUtilIface* readMeshIface;

  MBRange newMatsets;
  
  bool hasSolution;
};

#endif
