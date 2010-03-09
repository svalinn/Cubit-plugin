

#include "SMF_State.hpp"
#include <cstring>
#include <cstdlib>

inline int streq(const char *a,const char *b) { return std::strcmp(a,b)==0; }


SMF_State::SMF_State(const SMF_ivars& ivar, SMF_State *link)
{
    next = link;
    first_vertex = ivar.next_vertex;
    if( next )
    {
	vertex_correction = next->vertex_correction;
	xform = next->xform;
    }
    else
    {
	vertex_correction = 0;
        MBAffineXform identity;
	xform = identity;
    }

}

void SMF_State::vertex(double v[3])
{
    xform.xform_point(v);
}

void SMF_State::normal(double normal[3])
{
    xform.xform_vector(normal);
}

void SMF_State::face( int * verts, const SMF_ivars& ivar)
{
    for(int i=0; i<3; i++)
    {
	if( verts[i] < 0 )
	    verts[i] += ivar.next_vertex;
	else
	    verts[i] += vertex_correction + (first_vertex - 1);
    }
}

void SMF_State::set( std::vector<std::string> & argv)
{
    const char *cmd = argv[0].c_str();

    if( streq(cmd, "vertex_correction") )
	vertex_correction = atoi(argv[1].c_str());
}

void SMF_State::mmult(const MBAffineXform &M)
{
    // initially, we tried this:
    // xform.accumulate(M);
    // maybe we should do M.accumulate(xform) 
    MBAffineXform tmp=M;
    tmp.accumulate(xform);
    xform = tmp;
}

void SMF_State::mload(const MBAffineXform& M)
{
    xform = M;
}
