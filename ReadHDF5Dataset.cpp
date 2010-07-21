/** \file   ReadHDF5Dataset.cpp
 *  \author Jason Kraftcheck 
 *  \date   2010-07-09
 */

#include <assert.h>
#include <H5Dpublic.h>

#include "ReadHDF5Dataset.hpp"

namespace moab {

ReadHDF5Dataset::ReadHDF5Dataset( hid_t data_set_handle,
                                  const Range& file_ids,
                                  EntityHandle start_id,
                                  hid_t data_type,
                                  bool close_data_set )
  : closeDataSet(close_data_set),
    dataSet( data_set_handle ),
    dataSpace(-1),
    dataType(data_type)
{ init( &file_ids, start_id ); }

ReadHDF5Dataset::ReadHDF5Dataset( hid_t data_set_handle,
                                  const Range& file_ids,
                                  EntityHandle start_id,
                                  hid_t data_type,
                                  int column,
                                  bool close_data_set )
  : closeDataSet(close_data_set),
    dataSet( data_set_handle ),
    dataSpace(-1),
    dataType(data_type)
{ init( &file_ids, start_id, column ); }

ReadHDF5Dataset::ReadHDF5Dataset( hid_t data_set_handle,
                                  hid_t data_type,
                                  bool close_data_set )
  : closeDataSet(close_data_set),
    dataSet( data_set_handle ),
    dataSpace(-1),
    dataType(data_type)
{ init( 0, 0 ); }

ReadHDF5Dataset::ReadHDF5Dataset( hid_t data_set_handle,
                                  hid_t data_type,
                                  int column,
                                  bool close_data_set )
  : closeDataSet(close_data_set),
    dataSet( data_set_handle ),
    dataSpace(-1),
    dataType(data_type)
{ init( 0, 0, column ); }

void ReadHDF5Dataset::init( const Range* file_ids, EntityHandle start_id, int column )
{
  dataSpace = H5Dget_space( dataSet );
  if (dataSpace < 0)
    throw Exception(__LINE__);
  
  dataSpaceRank = H5Sget_simple_extent_dims( dataSpace, dataSetCount, dataSetOffset );
  if (dataSpaceRank < 0) 
    throw Exception(__LINE__);
  
  for (int i = 0; i < dataSpaceRank; ++i)
    dataSetOffset[i] = 0;
  if (column == 0 && dataSpaceRank == 1) {
    // do nothing
  }
  if (column >= 0) {
    set_column( column );
  }
  
  if (!file_ids) {
    start_id = 1;
    internalRange.insert( (EntityHandle)1, (EntityHandle)(dataSetCount[0]) );
    file_ids = &internalRange;
  }
  
  set_file_ids( *file_ids, start_id );
}

unsigned ReadHDF5Dataset::columns() const
{
  if (dataSpaceRank == 1)
    return 1;
  else if (dataSpaceRank == 2)
    return dataSetCount[1];
  
  throw Exception(__LINE__);
}

void ReadHDF5Dataset::set_column( unsigned column )
{
  if (dataSpaceRank != 2 || column < 0 || column >= dataSetCount[1])
    throw Exception(__LINE__);
  dataSetCount[1] = 1;
  dataSetOffset[1] = column;
}

void ReadHDF5Dataset::set_file_ids( const Range& file_ids, EntityHandle start_id )
{
  startID = start_id;
  currOffset = file_ids.begin();
  rangeEnd = file_ids.end();
}

ReadHDF5Dataset::~ReadHDF5Dataset() 
{
  if (dataSpace >= 0)
    H5Sclose( dataSpace );
  if (closeDataSet && dataSet >= 0)
    H5Dclose( dataSet );
  dataSpace = dataSet = -1;
}

void ReadHDF5Dataset::read( void* buffer,
                            size_t max_rows,
                            size_t& rows_read,
                            hid_t io_prop )
{
  herr_t err;
  rows_read = 0;
  if (done()) {
    rows_read = 0;
    return;
  }
  
  assert(startID <= *currOffset);

    // Build H5S hyperslab selection describing the portions of the
    // data set to read
  H5S_seloper_t sop = H5S_SELECT_SET;
  size_t avail = max_rows;
  while (avail && currOffset != rangeEnd) {
    size_t count = *(currOffset.end_of_block()) - *currOffset + 1;
    if (count > avail)
      count = avail;
    avail -= count;
    rows_read += count;
    
    dataSetOffset[0] = *currOffset - startID;
    dataSetCount[0] = count;
    err = H5Sselect_hyperslab( dataSpace, sop, dataSetOffset, NULL, dataSetCount, 0 );
    if (err < 0)
      throw Exception(__LINE__);
    sop = H5S_SELECT_OR; // subsequent calls to select_hyperslab append
  
    currOffset += count;
  }
  
    // Create a data space describing the memory in which to read the data
  dataSetCount[0] = rows_read;
  hid_t mem_id = H5Screate_simple( dataSpaceRank, dataSetCount, NULL );
  if (mem_id < 0)
    throw Exception(__LINE__);
  
    // Do the actual read
  err = H5Dread( dataSet, dataType, mem_id, dataSpace, io_prop, buffer );
  H5Sclose( mem_id );
  if (err < 0)
    throw Exception(__LINE__);
}

void ReadHDF5Dataset::null_read( hid_t io_prop )
{
  herr_t err;
  err = H5Sselect_none( dataSpace );
  if (err < 0)
    throw Exception(__LINE__);
  
#if H5_VERS_MAJOR > 1 || H5_VERS_MINOR >= 8
  hid_t mem_id = H5Screate(H5S_NULL);
  if (mem_id < 0)
    throw Exception(__LINE__);
#else
  hid_t mem_id = H5Screate(H5S_SIMPLE);
  if (mem_id < 0)
    throw Exception(__LINE__);
  err = H5Sselect_none( mem_id );
  if (err < 0)
    throw Exception(__LINE__);
#endif

  err = H5Dread( dataSet, dataType, mem_id, dataSpace, io_prop, 0 );
  if (err < 0)
    throw Exception(__LINE__);
}

} // namespace moab
