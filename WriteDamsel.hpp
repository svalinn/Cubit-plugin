//-------------------------------------------------------------------------
// Filename      : WriteDamsel.hpp
//
// Purpose       : ExodusII writer
//
// Special Notes : Lots of code taken from verde implementation
//
// Creator       : Corey Ernst 
//
// Date          : 8/02
//
// Owner         : Corey Ernst 
//-------------------------------------------------------------------------

#ifndef WRITEDAMSEL_HPP
#define WRITEDAMSEL_HPP

#ifndef IS_BUILDING_MB
#error "WriteDamsel.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <string>

#include "moab/Forward.hpp"
#include "moab/Range.hpp"
#include "moab/WriterIface.hpp"
#include "RangeSeqIntersectIter.hpp"
#include "FileOptions.hpp"

#include "damsel.h"

namespace moab {

class WriteUtilIface;
class SequenceManager;
class Error;

class WriteDamsel : public WriterIface
{
 
public:
    //! factory function, for ReaderWriter
  static WriterIface* factory( Interface* iface );

    //! Constructor
  WriteDamsel(Interface *impl);

    //! Destructor
  virtual ~WriteDamsel();

    //! Primary interface function
    //! \param file_name Filename being written
    //! \param overwrite If true and the file exists, an error is returned
    //! \param opts File options, e.g. whether and how to write in parallel
    //! \param meshset_list If non-NULL, a vector of sets defining what is to be written
    //! \param num_sets The number of sets to be written, only used if meshset_list is non-NULL
    //! \param qa_records Strings defining provenance information
    //! \param tag_list If non-NULL, only these tags should be written
    //! \param num_tags The number of tag handles in tag_list, used only if tag_list is non-NULL
    //! \param requested_output_dimension Dimension used for coordinates
  ErrorCode write_file(const char *file_name, 
                       const bool /* overwrite */,
                       const FileOptions& opts,
                       const EntityHandle *meshset_list,
                       const int num_sets,
                       const std::vector<std::string>& /* qa_records */,
                       const Tag* /* tag_list */,
                       int /* num_tags */,
                       int /* requested_output_dimension */);

  enum {DAMSEL_IS_TRACKING = 0x1
  } DAMSEL_FLAGS;
  
private:

    //! Write the sets in the model, for the handles in the specified RangeSeqIntersectIter
    //! \param rsi Range sequence iterator defining range of entities/sets to be written
  ErrorCode write_sets(RangeSeqIntersectIter &rsi);
  

    //! Write the entities in the model, for the handles in the specified RangeSeqIntersectIter
    //! \param rsi Range sequence iterator defining range of entities/sets to be written
  ErrorCode write_entities(RangeSeqIntersectIter &rsi);
  

    //! Write the vertices in the model, for the handles in the specified RangeSeqIntersectIter
    //! \param rsi Range sequence iterator defining range of entities/sets to be written
  ErrorCode write_vertices(RangeSeqIntersectIter &rsi);
  
    //! Get the damsel tag ids for the corresponding MOAB coordinates
    //! This version assumes coordinates are stored in three separate Damsel tags
    //! \param xcoords_dtag Damsel id for tag used to store x coordinates for blocked storage
    //! \param ycoords_dtag Damsel id for tag used to store y coordinates for blocked storage
    //! \param zcoords_dtag Damsel id for tag used to store z coordinates for blocked storage
    //! \param create_if_missing If true and damsel id for a tag is not yet initialized, make one
  ErrorCode damsel_coords_tags(damsel_tag &xcoords_dtag, damsel_tag &ycoords_dtag, 
                               damsel_tag &zcoords_dtag, bool create_if_missing);
  
    //! Get the damsel tag id for the corresponding MOAB coordinates
    //! This version assumes coordinates are stored in the same Damsel tag
    //! \param xcoords_dtag Damsel id for tag used to store x coordinates for interleaved storage
    //! \param ycoords_dtag Damsel id for tag used to store y coordinates for interleaved storage
    //! \param zcoords_dtag Damsel id for tag used to store z coordinates for interleaved storage
    //! \param create_if_missing If true and damsel id for a tag is not yet initialized, make one
  ErrorCode damsel_coords_tags(damsel_tag &coords_dtag, bool create_if_missing);
  
    //! Write dense tags for the specified entities, using the specified damsel entity container
  ErrorCode write_dense_tags(RangeSeqIntersectIter &rsi, damsel_container &ent_cont);

    //! Write dense tags for the specified entities; a new entity container is created if necessary
  ErrorCode write_dense_tags(RangeSeqIntersectIter &rsi);
  
    //! Initialize global information about dense tags, once for entire write_file call
  ErrorCode init_dense_tag_info();

    //! interface instance
  Interface *mbImpl;

    //! WriteUtil object used in this writer
  WriteUtilIface* mWriteIface;

    //! Error object used to register errors
  Error *mError;
  
    //! Used to initialize the RangeSeqIntersectIter
  SequenceManager *sequenceManager;

    //! Used to track file handles
  Tag mGlobalIdTag;
  
    //! Used for set flags
  Tag mSetFlagsTag;
  
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

  damsel_data_type moab_to_damsel_data_type[MB_MAX_DATA_TYPE];}
;

} // namespace moab

#endif
