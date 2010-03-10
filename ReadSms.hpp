#ifndef READ_SMS_HPP
#define READ_SMS_HPP

#include "MBForward.hpp"
#include "MBReaderIface.hpp"
#include "MBRange.hpp"
#include <vector>

class MBReadUtilIface;

// Base class for binary and ASCII readers
class ReadSms : public MBReaderIface
{
   
public:

    //! factory method 
  static MBReaderIface* factory( MBInterface* );

  MBErrorCode load_file( const char *file_name,
                         const MBEntityHandle* file_set,
                         const FileOptions& opts,
                         const MBReaderIface::IDTag* subset_list = 0,
                         int subset_list_length = 0,
                         const MBTag* file_id_tag = 0 );

  MBErrorCode read_tag_values( const char* file_name,
                               const char* tag_name,
                               const FileOptions& opts,
                               std::vector<int>& tag_values_out,
                               const IDTag* subset_list = 0,
                               int subset_list_length = 0 );
  
    //! Constructor
  ReadSms(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadSms();

private:

  MBErrorCode add_entities( MBEntityHandle start, 
                            MBEntityHandle count,
                            const MBTag* file_id_tag );

  MBErrorCode load_file_impl( FILE* file, const MBTag* file_id_tag );
  
  MBErrorCode get_set(std::vector<MBEntityHandle> *sets,
                      int set_type, int set_id,
                      MBTag set_tag,
                      MBEntityHandle &this_set,
                      const MBTag* file_id_tag );

  MBErrorCode read_parallel_info(FILE *file_ptr);

  MBReadUtilIface* readMeshIface;

    //! interface instance
  MBInterface* mdbImpl;
  
  MBTag globalId, paramCoords, geomDimension;
  
  int setId;
};


#endif
