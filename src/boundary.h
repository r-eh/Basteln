/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is boundary of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#define EMPTY  -1
#define PASSIVE 1
#define ACTIVE 10
#define MOVING 20


#include"increment.h"

class lexer;
class ghostcell;
class fdm_nhf;

class boundary : public increment
{
public:
    boundary(lexer*, ghostcell *);
    ~boundary();
    
// functions
    // ini
    void ini_storage(lexer*,ghostcell*);
    void fill(lexer*,ghostcell*);
    
    // add
    void add(lexer*,ghostcell*,double,double,double,double,double);
    
    void resize(lexer*,int);

    // remove
    void remove(int);
    void erase_all();
    
    // update
    void update_nhflow(lexer*,fdm_nhf*,ghostcell*);

    
// data arrays
    int *iloc,*jloc,*kloc;
    int *cellside;
    int *bc_type;
    double *ks;
    
    int grid_level;
    
// iterators
    int index; // replace loopindex
    int index_empty,index_empty0; //
    int numactive,numempty;   // number of active boundaryicles
    int capacity; // length of allocated array
    int capacity_para,maxnum;
    
    int n,q;
    
private:

    
};

#endif
