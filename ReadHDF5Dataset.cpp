/** \file   ReadHDF5Dataset.cpp
 *  \author Jason Kraftcheck 
 *  \date   2010-07-09
 */

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

#include "ReadHDF5Dataset.hpp"

#include <H5Dpublic.h>
#include <H5Tpublic.h>
#include <H5Ppublic.h>
#ifdef HDF5_PARALLEL
#  include <H5FDmpi.h>
#  include <H5FDmpio.h>
#endif

#define HDF5_16API (H5_VERS_MAJOR < 2 && H5_VERS_MINOR < 8)

namespace moab {

//const size_t MAX_POINT_BUFFER_SIZE = 8*1024*1024/sizeof(hsize_t);
const size_t READALL_BUFFER_SIZE = 1u<<12; // 4M

ReadHDF5Dataset::ReadHDF5Dataset( hid_t data_set_handle,
                                  hid_t data_type,
                                  bool close_data_set )
  : ioMode(HYPERSLAB),
    closeDataSet(close_data_set),
    dataSet( data_set_handle ),
    dataSpace(-1),
    dataType(data_type),
    ioProp(H5P_DEFAULT)
{ 
  dataSpace = H5Dget_space( dataSet );
  if (dataSpace < 0)
    throw Exception(__LINE__);
  
  dataSpaceRank = H5Sget_simple_extent_dims( dataSpace, dataSetCount, dataSetOffset );
  if (dataSpaceRank < 0) 
    throw Exception(__LINE__);
  rowsInTable = dataSetCount[0];
  
  for (int i = 0; i < dataSpaceRank; ++i)
    dataSetOffset[i] = 0;

  currOffset = rangeEnd = internalRange.end();
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

const char MODE_ENV_VAR[] = "H5M_SPARSE_READ_MODE";

static const char* mode_name( ReadHDF5Dataset::Mode mode ) 
{
  static const char hyperslab[] = "HYPERSLAB";
  static const char point[] = "POINT";
  static const char contig[] = "CONTIGUOUS";
  static const char unknown[] = "<unknown>";
  switch (mode) {
    case ReadHDF5Dataset::HYPERSLAB:  return hyperslab;
    case ReadHDF5Dataset::POINT:      return point;
    case ReadHDF5Dataset::CONTIGUOUS: return contig;
  }
  return unknown;
}

const char* ReadHDF5Dataset::get_mode_str() const
  { return mode_name( ioMode ); }

void ReadHDF5Dataset::set_file_ids( const Range& file_ids, 
                                    EntityHandle start_id,
                                    hid_t io_prop,
                                    const Comm* communicator )
{
  startID = start_id;
  currOffset = file_ids.begin();
  rangeEnd = file_ids.end();
  ioProp = io_prop;

  if (file_ids.empty()) {
    ioMode = HYPERSLAB;
    return;
  }
  
  // Calculate mode to use.
  // We should default to HYPERSLAB if reading a single contiguous
  // block.  READALL and HYPERSLAB are techincally the same, but we
  // cannot honor the clients request for collective IO when doing
  // READALL.
  
  // Begin by calculating some statistics 
  
  //const size_t maxsize = std::numerical_limits<size_t>::max();
  size_t /*min_size = maxsize, max_size = 0,*/ sum_size = 0;
  size_t /*min_gap = maxsize, max_gap = 0,*/ sum_gap = 0;
  // NOTE: sums of squares will overflow for large sizes.  check
  //       max_size/max_gap before relying on these values.
  //size_t sqr_size = 0, sqr_gap = 0;
  size_t num_blocks = 0;
  
  Range::const_iterator j, i = file_ids.begin();
  for (;;) {
    j = i.end_of_block();
    size_t n = *j - *i + 1;
    //if (n < min_size)
    //  min_size = n;
    //if (n > max_size)
    //  max_size = n;
    sum_size += n;
    //sqr_size += n*n;
    ++num_blocks;
    
    i = j;
    ++i;
    if (i == file_ids.end())
      break;
      
    n = *i - *j - 1;
    //if (n < min_gap)
    //  min_gap = n;
    //if (n > max_gap)
    //  max_gap = n;
    sum_gap += n;
    //sqr_gap += n*n;
  }
  
  size_t avg_size = (sum_size + num_blocks/2) / num_blocks;
  //size_t var_size = (sqr_size - sum_size*sum_size + num_blocks/2) / num_blocks;
  size_t avg_gap/*, var_gap*/;
  if (num_blocks > 1) {
    avg_gap = (sum_gap + (num_blocks-1)/2) / (num_blocks-1);
    //var_gap = (sqr_gap - sum_gap*sum_gap + (num_blocks-1)/2) / (num_blocks-1);
  }
  else {
    avg_gap = 0;
    // var_gap = 0;
  }
  
    // If the blocks are relatively large, then use hyperslabs.
  if (num_blocks < 1000 || (num_blocks < 10000 && avg_size > 3)) {
    ioMode = HYPERSLAB;
  }
  
    // If the gap size is reasonably close to the data size, then
    // just read everything and pick out the values we want
  else if (avg_gap < 3*avg_size) {
    ioMode = CONTIGUOUS;
  }
  
    // If the average data size is less than 6, do point-wise
    // selection.  We cannot do point-wise selection for data 
    // sets with rank greater than 2 (because it isn't implemented)
  else if (dataSpaceRank <= 2) {
    ioMode = POINT;
  }
  
    // Otherwise just use hyperslab
  else {
    ioMode = HYPERSLAB;
  }
  
  const char* e = getenv( MODE_ENV_VAR );
  if (e && *e) {
    Mode old_mode = ioMode;
    ioMode = (Mode)atoi(e);
    switch (ioMode) {
      case HYPERSLAB: break;
      case POINT:     break;
      case CONTIGUOUS:   break;
      default: std::cerr << "INVALID MODE sppecified in \"" << MODE_ENV_VAR
                         << "\" env var: " << (int)ioMode << std::endl;
             ioMode = old_mode;
    }
  }

  // If doing collective IO, make sure all procs have compatible schemes.
  // Everything should work OK if some procs are doing HYPERSLAB and some
  // are doing POINT, but if any proc is doing CONTIGUOUS, then all procs
  // need to do that scheme because it forces independent IO.
#ifdef HDF5_PARALLEL
  H5FD_mpio_xfer_t mode = H5FD_MPIO_INDEPENDENT;
  if (io_prop != H5P_DEFAULT) {
    herr_t err = H5Pget_dxpl_mpio( io_prop, &mode );
    if (err < 0)
      mode = H5FD_MPIO_INDEPENDENT;
  }
  if (mode == H5FD_MPIO_COLLECTIVE) {
    if (!communicator)
      throw Exception(__LINE__);
    int int_mode = ioMode, all_mode = ioMode;
    int rval = MPI_Allreduce( &int_mode, &all_mode, 1, MPI_INT, MPI_MAX, *communicator );
    if (rval != MPI_SUCCESS)
      throw Exception(__LINE__);
    if (all_mode == CONTIGUOUS)
      ioMode = (Mode)all_mode;
  }
#endif
}

void ReadHDF5Dataset::set_all_file_ids( hid_t io_prop, const Comm* communicator )
{
  internalRange.clear();
  internalRange.insert( (EntityHandle)1, (EntityHandle)(rowsInTable) );
  set_file_ids( internalRange, 1, io_prop, communicator );
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
                            size_t& rows_read  )
{
  assert(startID <= *currOffset);

  switch (ioMode) {
    case HYPERSLAB:
      read_hyperslab( buffer, max_rows, rows_read );
      break;
    case POINT:
      read_point( buffer, max_rows, rows_read );
      break;
    case CONTIGUOUS:
      read_contig( buffer, max_rows, rows_read );
      break;
  }
}

void ReadHDF5Dataset::read_hyperslab( void* buffer,
                                      size_t max_rows,
                                      size_t& rows_read )
{
  herr_t err;
  rows_read = 0;
  if (done()) {
    return;
  }

#ifndef NDEBUG
  size_t N = rangeEnd -  currOffset;
  Range::const_iterator beg = currOffset;
#endif

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
  
#ifndef NDEBUG
    assert( N >= count );
    N -= count;
#endif
    currOffset += count;
  }

#ifndef NDEBUG
  size_t dist = currOffset - beg;
  assert(dist == rows_read);
  assert((size_t)(rangeEnd - currOffset) == N);
#endif
  
    // Create a data space describing the memory in which to read the data
  dataSetCount[0] = rows_read;
  hid_t mem_id = H5Screate_simple( dataSpaceRank, dataSetCount, NULL );
  if (mem_id < 0)
    throw Exception(__LINE__);
  
    // Do the actual read
  err = H5Dread( dataSet, dataType, mem_id, dataSpace, ioProp, buffer );
  H5Sclose( mem_id );
  if (err < 0)
    throw Exception(__LINE__);
}

void ReadHDF5Dataset::read_point( void* buffer,
                                  size_t max_rows,
                                  size_t& rows_read )
{
  herr_t err;
  rows_read = 0;
  if (done()) {
    return;
  }
  
  if (dataSpaceRank == 1) {
    // limit buffer size to keep memory use nominal
//    if (max_rows > MAX_POINT_BUFFER_SIZE)
//      max_rows = MAX_POINT_BUFFER_SIZE;
    selectData.clear();
    while (currOffset != rangeEnd && rows_read < max_rows) {
      selectData.push_back( *currOffset - startID );
      ++rows_read;
      ++currOffset;
    }
  }
  else if (dataSpaceRank == 2) {
    // limit buffer size to keep memory use nominal
//    if (max_rows*dataSetCount[1]*2 > MAX_POINT_BUFFER_SIZE)
//      max_rows = MAX_POINT_BUFFER_SIZE/(dataSetCount[1]*2);
    selectData.clear();
    while (currOffset != rangeEnd && rows_read < max_rows) {
      for (size_t j = 0; j < dataSetCount[1]; ++j) {
        selectData.push_back( *currOffset - startID );
        selectData.push_back( dataSetOffset[1] + j );
      }
      ++rows_read;
      ++currOffset;
    }
  }
  else {
    // not implemented for rank > 2;
    assert(false);
    throw Exception( __LINE__ );
  }

    // because we may have limited max_rows to prevent the coordiante
    // buffer from growing w/out bound on any processor, we cannot do 
    // collective reads because different procs may do differnt numbers
    // of reads.  Force independent IO.
//  io_prop = H5P_DEFAULT;

  const hsize_t count = selectData.size() / dataSpaceRank;
#if HDF5_16API
//  selectVector.clear();
//  for (size_t i = 0; i < selectData.size(); i += dataSpaceRank)
//    selectVector.push_back( &selectData[i] );
  err = H5Sselect_elements( dataSpace, H5S_SELECT_SET, count, (const hsize_t**)&selectData[0] );
#else
  err = H5Sselect_elements( dataSpace, H5S_SELECT_SET, count, &selectData[0] );
#endif
  if (err < 0)
    throw Exception(__LINE__);
  
    // Create a data space describing the memory in which to read the data
  hid_t mem_id = H5Screate_simple( 1, &count, NULL );
  if (mem_id < 0)
    throw Exception(__LINE__);
  
    // Do the actual read
  err = H5Dread( dataSet, dataType, mem_id, dataSpace, ioProp, buffer );
  H5Sclose( mem_id );
  if (err < 0)
    throw Exception(__LINE__);
}

void ReadHDF5Dataset::read_contig( void* buffer,
                                   size_t max_rows,
                                   size_t& rows_read )
{
  herr_t err;
  rows_read = 0;
  if (done()) {
    return;
  }

    // advance I until it is at the end or the start of the next pair
    // beyond what we want to read
  Range::const_iterator i = currOffset;
  while (i != rangeEnd && *(i.end_of_block()) - *currOffset  < max_rows) {
    i = i.end_of_block();
    ++i;
  }
  if (i != rangeEnd && *i - *currOffset < max_rows) {
    assert( *(i.end_of_block()) - *currOffset  >= max_rows );
    i += max_rows - (*i - *currOffset);
  }
  --i;  // point to last included value, rather than one past it
  
    // read a contiguous block of data  
  dataSetOffset[0] = *currOffset - startID;
  dataSetCount[0] = *i - *currOffset + 1;
  assert( dataSetCount[0] <= max_rows );
  //assert( i+1 == rangeEnd || dataSetCount[0] == max_rows );
  err = H5Sselect_hyperslab( dataSpace, H5S_SELECT_SET, dataSetOffset, NULL, dataSetCount, 0 );
  if (err < 0)
    throw Exception(__LINE__);
  
    // Create a data space describing the memory in which to read the data
  dataSetCount[0] = dataSetCount[0];
  hid_t mem_id = H5Screate_simple( dataSpaceRank, dataSetCount, NULL );
  if (mem_id < 0)
    throw Exception(__LINE__);
  
    // Do the actual read
  err = H5Dread( dataSet, dataType, mem_id, dataSpace, H5P_DEFAULT, buffer );
  H5Sclose( mem_id );
  if (err < 0)
    throw Exception(__LINE__);
  
    // figure out the size of one row
  size_t bytes = H5Tget_size( dataType );
  for (int i = 1; i < dataSpaceRank; ++i) 
    bytes *= dataSetCount[i];
  
    // remove unwanted values
  ++i;
  rows_read = 0;
  EntityHandle start = *currOffset;
  char* dst = (char*)buffer;
  while (currOffset != i) {
    size_t n = *currOffset.end_of_block() - *currOffset + 1;
    if (i != rangeEnd && *currOffset.end_of_block() >= *i)
      n = *i - *currOffset;
    void* src = (char*)buffer + (*currOffset - start)*bytes;
    memmove( dst, src, n*bytes );
    dst += n*bytes;
    rows_read += n;
    currOffset += n;
  }
  assert( dst - (char*)buffer == (ptrdiff_t)(rows_read * bytes) );
}


void ReadHDF5Dataset::null_read()
{
  herr_t err;
  err = H5Sselect_none( dataSpace );
  if (err < 0)
    throw Exception(__LINE__);
  
#if HDF5_16API
  hid_t mem_id = H5Screate(H5S_SIMPLE);
  if (mem_id < 0)
    throw Exception(__LINE__);
  err = H5Sselect_none( mem_id );
  if (err < 0)
    throw Exception(__LINE__);
#else
  hid_t mem_id = H5Screate(H5S_NULL);
  if (mem_id < 0)
    throw Exception(__LINE__);
#endif

  err = H5Dread( dataSet, dataType, mem_id, dataSpace, ioProp, 0 );
  if (err < 0)
    throw Exception(__LINE__);
}

} // namespace moab
