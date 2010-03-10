/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**
 * \class ReadSmf
 * \brief SMF reader from QSLIM
 * \author Michael Garland 
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "ReadSmf.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"
#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "FileOptions.hpp"
#include "MBAffineXform.hpp"

inline int streq(const char *a,const char *b) { return strcmp(a,b)==0; }

ReadSmf::cmd_entry ReadSmf::read_cmds[] = {
    { "v", &ReadSmf::vertex },
    { ":vn", &ReadSmf::v_normal },
    { ":vc", &ReadSmf::v_color },
    { ":fc", &ReadSmf::f_color },
    { "t", &ReadSmf::face },
    { "f", &ReadSmf::face },

    { "begin", &ReadSmf::begin },
    { "end", &ReadSmf::end },
    { "set", &ReadSmf::set },
    { "inc", &ReadSmf::inc },
    { "dec", &ReadSmf::dec },

    { "mmult", &ReadSmf::mload },
    { "mload", &ReadSmf::mmult },
    { "trans", &ReadSmf::trans },
    { "scale", &ReadSmf::scale },
    { "rot", &ReadSmf::rot },

    { NULL, NULL }
};

MBAffineXform mat_from_args(std::vector<std::string> & argv)
{
    double m3[9], offset[3];
    for (int i=0; i<9; i++)
	m3[i] = atof(argv[i].c_str());
    for (int j=0; j<3; j++)
	offset[j] = atof(argv[j+9].c_str());// only the first 12 are used, the last row (or column?) is 0001?
    MBAffineXform M(m3, offset);
    return M;
}

void bad_annotation(char *cmd)
{
    std::cerr << "SMF: Malformed annotation ["<< cmd << "]" << std::endl;
}

MBReaderIface* ReadSmf::factory( MBInterface* iface )
  { return new ReadSmf( iface ); }

ReadSmf::ReadSmf(MBInterface* impl)
    : mdbImpl(impl)
{
  void* ptr = 0;
  mdbImpl->query_interface("MBReadUtilIface", &ptr);
  readMeshIface = reinterpret_cast<MBReadUtilIface*>(ptr);
  _numNodes= _numFaces = 0;
  _numNodesInFile = _numElementsInFile = 0;
}

ReadSmf::~ReadSmf()
{
  if (readMeshIface) {
    mdbImpl->release_interface("MBReadUtilIface", readMeshIface);
    readMeshIface = 0;
  }
}


MBErrorCode ReadSmf::read_tag_values( const char* /* file_name */,
                                      const char* /* tag_name */,
                                      const FileOptions& /* opts */,
                                      std::vector<int>& /* tag_values_out */,
                                      const IDTag* /* subset_list */,
                                      int /* subset_list_length */ )
{
  return MB_NOT_IMPLEMENTED;
}


MBErrorCode ReadSmf::load_file( const char *filename,
                                const MBEntityHandle* ,
                                const FileOptions& opts,
                                const MBReaderIface::IDTag* subset_list,
                                int subset_list_length,
                                const MBTag* file_id_tag) 
{
  MBErrorCode result;
  
  if (subset_list && subset_list_length) {
    readMeshIface->report_error( "Reading subset of files not supported for VTK." );
    return MB_UNSUPPORTED_OPERATION;
  }

  // Does the caller want a field to be used for partitioning the entities?
  // If not, we'll assume any scalar integer field named MATERIAL_SET specifies partitions.
  std::string partition_tag_name;
  result = opts.get_option( "PARTITION", partition_tag_name );
  if ( result == MB_SUCCESS )
    mPartitionTagName = partition_tag_name;

  std::ifstream smfFile ( filename );
  if (!smfFile)
  {
    return MB_FILE_DOES_NOT_EXIST;
  }
  init();
  

  while( !smfFile.eof() )
    {
	if( smfFile.getline(line, SMF_MAXLINE, '\n').good() )
        {
	    result = parse_line(line);
            if (MB_SUCCESS != result)
    		return result;
	}
    }

  // at this point we have _numNodesInFile vertices and _numElementsInFile triangles
  // the coordinates are in _coords, and connectivities in _connec
  // std::vector<double> _coords; // 3*numNodes; we might not know the number of nodes
  // std::vector<int> _connec; // 3*num of elements; we might not know them;
  
    // Create vertices
  std::vector<double*> arrays;
  MBEntityHandle start_handle_out;
  start_handle_out = 0;
  result = readMeshIface->get_node_arrays( 3, _numNodesInFile, MB_START_ID,
                                           start_handle_out, arrays );

  if (MB_SUCCESS != result)
    return result;

  
  // fill the arrays with data from _coords
  for (int i=0; i<_numNodesInFile; i++)
    {
       int i3 = 3*i;
       arrays[0][i] = _coords[i3];
       arrays[1][i] = _coords[i3+1];
       arrays[2][i] = _coords[i3+2];
    }
  // elements
  
  MBEntityHandle start_handle_elem_out;
  start_handle_elem_out = 0;
  MBEntityHandle* conn_array_out;
  result = readMeshIface->get_element_array( _numElementsInFile,
                                             3,
                                             MBTRI , // MBEntityType
                                             MB_START_ID,
                                             start_handle_elem_out,
                                             conn_array_out );
  if (MB_SUCCESS != result)
    return result;
  for (int j=0; j<_numElementsInFile*3; j++)
  {
     conn_array_out[j] = _connec[j];
  }

  // notify MOAB of the new elements
  result = readMeshIface->update_adjacencies(start_handle_elem_out, _numElementsInFile,
                                               3, conn_array_out);

  if (MB_SUCCESS != result)
    return result;

  MBRange range(start_handle_out, start_handle_out+_numElementsInFile-1);
 
 
  return MB_SUCCESS;
}

void ReadSmf::init( )
{
    ivar.next_face = 1;
    ivar.next_vertex = 1;
    state = new SMF_State(ivar);
    line = new char [4096];
    return;
}



void ReadSmf::annotation(char *cmd,  std::vector<std::string> & argv)
{
   // Skip over the '#$' prefix
    cmd+=2;

    if( streq(cmd,"SMF") ) {
	if( atof(argv[0].c_str() ) != 1.0 )
	    std::cerr << "SMF: Version mismatch ("
		 << argv[0] << " instead of "
		 << "1.0" << ")" << std::endl;
    }
    else if( streq(cmd,"vertices") )
    {
	if( argv.size() == 1 )
	    _numNodes = atoi(argv[0].c_str() );
	else
	    bad_annotation(cmd);
    }
    else if( streq(cmd, "faces") )
    {
	if( argv.size() == 1 )
	    _numFaces = atoi(argv[0].c_str() );
	else
	    bad_annotation(cmd);

    }
    else if( streq(cmd, "BBox") )
    {
    }
    else if( streq(cmd, "BSphere") )
    {
    }
    else if( streq(cmd, "PXform") )
    {
	if( argv.size() == 16 )
	    mat_from_args(argv);
	else
	    bad_annotation(cmd);
    }
    else if( streq(cmd, "MXform") )
    {
	if( argv.size() == 16 )
	    mat_from_args(argv);
	else
	    bad_annotation(cmd);
    }
    
}

MBErrorCode ReadSmf::parse_line(char *line)
{
    char *cmd,*s;
    std::vector<std::string>  argv;

    while( *line==' ' || *line=='\t' ) line++;  // skip initial white space

    // Ignore empty lines
    if( line[0]=='\n' || line[0]=='\0' ) return MB_SUCCESS;

    // Ignore comments
    if( line[0]=='#' && line[1]!='$' ) return MB_SUCCESS;

    //
    // First, split the line into tokens
    cmd = strtok(line, " \t\n");

    while( (s=strtok(NULL, " \t\n")) )
    {
	std::string stg(s);
	argv.push_back(stg);
    }

    //
    // Figure out what command it is and execute it
    if( cmd[0]=='#' && cmd[1]=='$' )
	annotation(cmd,argv);
    else
    {
	cmd_entry *entry = &read_cmds[0];
	bool handled = 0;

	while( entry->name && !handled )
	    if( streq(entry->name, cmd) )
	    {
		(this->*(entry->cmd))(argv);
		handled = 1;
	    }
	    else
		entry++;

	if( !handled  )
	{
	    // Invalid command:
	    std::cerr << "SMF: Illegal command [" << cmd << "]" << std::endl;
	    return MB_UNSUPPORTED_OPERATION;
	}
    }
    return MB_SUCCESS;
}


void ReadSmf::vertex(std::vector<std::string> & argv)
{
    double  v[3];

    for (int i=0; i<3; i++)
	v[i] = atof(argv[i].c_str() );

    state->vertex(v);
    ivar.next_vertex++;
    _numNodesInFile++;
    for (int j=0; j<3; j++)
	_coords.push_back(v[j]);
    //model->in_Vertex(v);
}
void ReadSmf::v_normal( std::vector<std::string> & argv )
{
}
void ReadSmf::v_color(std::vector<std::string> & argv)
{
}
void ReadSmf::f_color(std::vector<std::string> & argv)
{
}
void ReadSmf::face(std::vector<std::string> & argv )
{
    int vert[3];

    for(int i=0; i<argv.size(); i++)
	vert[i] = atoi(argv[i].c_str());

    state->face(vert, ivar);
    ivar.next_face++;
    for (int j=0; j<3; j++)
	_connec.push_back(vert[j]);
    _numElementsInFile++;

}

void ReadSmf::begin(std::vector<std::string> & argv)
{
   state = new SMF_State(ivar,state);
}
void ReadSmf::end(std::vector<std::string> & argv)
{
    SMF_State *old = state;
    state=state->pop();
    delete old;
}
void ReadSmf::set(std::vector<std::string> & argv)
{
    state->set(argv);
}
void ReadSmf::inc(std::vector<std::string> & argv)
{
    std::cerr << "SMF: INC not yet implemented." << std::endl;
}
void ReadSmf::dec(std::vector<std::string> & argv)
{
    std::cerr << "SMF: DEC not yet implemented." << std::endl;
}

void ReadSmf::trans(std::vector<std::string> & argv)
{
    double v3[3];
    for (int i=0; i<3; i++)
	v3[i] = atof(argv[i].c_str());
    MBAffineXform M = MBAffineXform::translation(v3);
    //Mat4 M = Mat4::trans(atof(argv(0)), atof(argv(1)), atof(argv(2)));
    state->mmult(M);
}
void ReadSmf::scale(std::vector<std::string> & argv)
{
    double v3[3];
    for (int i=0; i<3; i++)
	v3[i] = atof(argv[i].c_str());
    MBAffineXform M = MBAffineXform::scale(v3);
    //Mat4 M = Mat4::scale(atof(argv(0)), atof(argv(1)), atof(argv(2)));
    state->mmult(M);
}
void ReadSmf::rot(std::vector<std::string> & argv)
{
    double angle =  atof(argv[1].c_str() ) * M_PI /180.0;
    double axis[3] = {0., 0., 0.};
    switch( argv[0][0] )
    {
	case 'x':
	 axis[0] = 1.;
	break;
	case 'y':
	 axis[1] = 1.;
	break;
	case 'z':
	 axis[2] = 1.; 
	break;
	default:
	  std::cerr << "SMF: Malformed rotation command" << std::endl;
	break;
    }
    MBAffineXform M = MBAffineXform::rotation( angle, axis );
    state->mmult(M);
}
void ReadSmf::mmult(std::vector<std::string> & argv)
{
    state->mmult(mat_from_args(argv));
}
void ReadSmf::mload(std::vector<std::string> & argv)
{
    state->mload(mat_from_args(argv));
}
