/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#include <stdlib.h>
#include <string.h>
#include <H5Ipublic.h>
#include "file-handle.h"
#include "status.h"
#include "util.h"

#define FILE_HANDLE_MAGIC 0xFEEDFEED

int mhdf_check_valid_file( FileHandle* handle, mhdf_Status* status )
{
  if (!handle)
  {
    mhdf_setFail( status, "NULL file handle." );
    return 0;
  }
  
  if (handle->magic !=  FILE_HANDLE_MAGIC)
  {
    mhdf_setFail( status, "Invalid file handle." );
    return 0;
  }
  
  return 1;
}

FileHandle* mhdf_alloc_FileHandle( hid_t hdf_table, mhdf_Status* status )
{
  FileHandle* rval = (FileHandle*)mhdf_malloc( sizeof(FileHandle), status );
  if (!rval) return NULL;
  
  rval->magic = FILE_HANDLE_MAGIC;
  rval->hdf_handle = hdf_table;
  rval->open_handle_count = 0;
  rval->max_id = 0L;
  return rval;
}
