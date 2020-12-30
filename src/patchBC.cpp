/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"patchBC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

patchBC::patchBC(lexer *p) 
{
    
}

patchBC::~patchBC()
{
}

void patchBC::patchBC_ini(lexer *p, ghostcell *pgc)
{
    patchBC_gcb_convert(p,pgc);
    patchBC_IDcount(p,pgc);

    
    // creat patch objects
    patch = new patch_obj*[obj_count];
    
    for(qn=0; qn<obj_count;++qn)
    patch[qn] = new patch_obj(p);
    
    // fill patch objects
    

} 
