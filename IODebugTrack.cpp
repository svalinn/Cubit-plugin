#include "IODebugTrack.hpp"
#include <iostream>
#include <vector>
#include <assert.h>

#ifdef USE_MPI
#  include "moab_mpi.h"
#endif

const char PFX[] = ">>> ";

namespace moab {

IODebugTrack::IODebugTrack( bool enabled,
                            const std::string name,
                            std::ostream output_stream,
                            unsigned long table_size )
          : enableOutput(enabled),
            tableName(name),
            ostr(output_stream),
            maxSize(table_size) 
{
#ifdef USE_MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank );
#else
  mpiRank = 0;
#endif
}


IODebugTrack::IODebugTrack( bool enabled,
                            const std::string name,
                            unsigned long table_size )
          : enableOutput(enabled),
            tableName(name),
            ostr(std::cerr),
            maxSize(table_size) 
{
#ifdef USE_MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank );
#else
  mpiRank = 0;
#endif
}

IODebugTrack::~IODebugTrack()
{
  if (!enableOutput || mpiRank) // only root prints gap summary
    return;
  
  if (dataSet.empty()) {
    ostr << PFX << tableName << " : No Data Written!!!!" << std::endl;
    return;
  }
  
  std::list<DRange>::const_iterator j, i = dataSet.begin();
  if (i->begin > 0) {
    ostr << PFX << tableName << " : [0," << i->begin-1 << "] : No Data Written!" << std::endl;
    ostr.flush();
  }
  j = i;
  for (++j; j != dataSet.end(); ++i, ++j) 
    if (i->end + 1 != j->begin) {
      ostr << PFX << tableName << " : [" << i->end + 1 << "," << j->begin-1 << "] : No Data Written!" << std::endl;
      ostr.flush();
    }
  if (maxSize && i->end + 1 < maxSize) {
    ostr << PFX << tableName << " : [" << i->end + 1 << "," << maxSize - 1 << "] : No Data Written!" << std::endl;
    ostr.flush();
  }
}

void IODebugTrack::record_io( unsigned long begin, unsigned long count )
{
  if (enableOutput && count) {
    DRange ins = { begin, begin+count-1, mpiRank };
    record_io( ins );
  }
}

void IODebugTrack::record_io( DRange ins )
{
  if (!enableOutput)
    return;

    // only root should get non-local data
  assert(!mpiRank || ins.rank == (unsigned)mpiRank);

    // test for overlap with all existing ranges
  std::list<DRange>::iterator i;
  for (i = dataSet.begin(); i != dataSet.end(); ++i) {
    if (i->end >= ins.begin && i->begin <= ins.end) { // if overlap
      ostr << PFX << tableName;
      if (i->rank == ins.rank) {
        if (mpiRank == (int)ins.rank) 
          ostr << ": Local overwrite on rank " << mpiRank;
        
        // otherwise should have been logged on remote proc, do nothing here
      }
      else 
        ostr << ": Conflicting write for ranks " << i->rank << " and " << ins.rank;
      
      ostr << ": [" << i->begin << "," << i->end << "] and [" << ins.begin
           << "," << ins.end << "]" << std::endl;
      ostr.flush();
    }
  }

  dataSet.push_back( ins );
}

void IODebugTrack::all_reduce()
{
#ifdef USE_MPI
  if (!enableOutput)
    return;

  int commsize;
  MPI_Comm_size( MPI_COMM_WORLD, &commsize);
  int count = dataSet.size();
  std::vector<int> displs(commsize);
  MPI_Gather( &count, 1, MPI_INT, 
              &displs[0], 1, MPI_INT,
              0, MPI_COMM_WORLD );
  int total = 0;
  for (int i = 0; i < commsize; ++i) {
    int tmp = displs[i];
    displs[i] = 3*total;
    total += tmp;
  }
              
  std::vector<int> counts(commsize);
  std::vector<DRange> send(dataSet.size()), recv(total);
  std::copy( dataSet.begin(), dataSet.end(), send.begin() );
  MPI_Gatherv( &send[0], 3*send.size(), MPI_UNSIGNED_LONG,
               &recv[0], &counts[0], &displs[0], MPI_UNSIGNED_LONG,
               0, MPI_COMM_WORLD );
  
  if (0 == mpiRank) {
    for (int i = count; i < total; ++i)
      record_io( recv[i] );
  }
  dataSet.clear();  
#endif
}
    


} // namespace moab
