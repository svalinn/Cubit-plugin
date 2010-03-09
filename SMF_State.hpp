#ifndef GFXSMF_STATE_INCLUDED // -*- C++ -*-
#define GFXSMF_STATE_INCLUDED

#include "MBAffineXform.hpp"
#include <string>
#include <vector>

struct SMF_ivars
{
    int next_vertex;
    int next_face;
};

class SMF_State
{
private:
    SMF_State *next;

    //
    // Standard state variables
    int first_vertex;
    int vertex_correction;
    MBAffineXform xform;

public:
    SMF_State(const SMF_ivars& ivar,SMF_State *link=NULL);
    SMF_State *pop() { return next; }

    void set(std::vector<std::string> & argv);
    void inc(const char *var, int delta=1);
    void dec(const char *var, int delta=1);

    void mmult(const MBAffineXform&);
    void mload(const MBAffineXform&);

    void vertex(double v[3]);
    void normal(double n[3]);
    void face(int *, const SMF_ivars& ivar);
};


// GFXSMF_STATE_INCLUDED
#endif
