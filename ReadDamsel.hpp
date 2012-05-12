//-------------------------------------------------------------------------
// Filename      : ReadDamsel.hpp
//
// Purpose       : Damsel file reader
//
// Creator       : Tim Tautges
//-------------------------------------------------------------------------

#ifndef READDAMSEL_HPP
#define READDAMSEL_HPP

#ifndef IS_BUILDING_MB
//#error "ReadDamsel.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <map>
#include <string>

#include "moab/Forward.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/Range.hpp"
#include "DebugOutput.hpp"

#ifdef USE_MPI
#  include "moab_mpi.h"
#  include "moab/ParallelComm.hpp"
#endif 

#include "damsel.h"

namespace moab {

class ReadUtilIface;
class ParallelComm;

class ReadDamsel : public ReaderIface
{
   
public:
  
  static ReaderIface* factory( Interface* );
  
    //! load an NC file
  ErrorCode load_file( const char* file_name,
                       const EntityHandle* file_set,
                       const FileOptions& opts,
                       const SubsetList* subset_list = 0,
                       const Tag* file_id_tag = 0 );

   //! Constructor
   ReadDamsel(Interface* impl = NULL);

   //! Destructor
  virtual ~ReadDamsel();

  virtual ErrorCode read_tag_values( const char* file_name,
                                     const char* tag_name,
                                     const FileOptions& opts,
                                     std::vector<int>& tag_values_out,
                                     const SubsetList* subset_list = 0 );
private:

//------------member variables ------------//

    //! interface instance
  Interface* mbImpl;

    //! utilIface
  ReadUtilIface *readMeshIface;
  
    //! file name
  std::string fileName;

    //! whether this reader natively supports parallel semantics
  bool nativeParallel;

    //! parallel info
  ParallelComm *myPcomm;
  
    //! Used to track file handles
  Tag mGlobalIdTag;
  
    //! file name
  std::string fileName;

    //! damsel library id
  damsel_library dmslLib;
  
    //! damsel model id
  damsel_model dmslModel;
  
    //! damsel coordinates tag ids (only first is used if interleaved)
  damsel_tag dmslXcoord, dmslYcoord, dmslZcoord;

    //! all dense tag handles in model
  std::vector<Tag> denseTags;
  
    //! damsel ids for dense tags
  std::vector<damsel_tag> dmslDenseTags;

    //! Damsel handle type used in (this build of) MOAB
  damsel_handle_type moabHandleType;

  damsel_data_type moab_to_damsel_data_type[MB_MAX_DATA_TYPE];
};

} // namespace moab

#endif
