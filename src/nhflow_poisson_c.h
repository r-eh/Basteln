
/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

This file is part of REEF3D.

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"increment.h"

class lexer;
class fdm_nhf;
class ghostcell;
class ioflow;
class poisson;
class solver;

#ifndef NHFLOW_POISSON_C_H_
#define NHFLOW_POISSON_C_H_

using namespace std;


class nhflow_poisson_c : public increment
{

public:

	nhflow_poisson_c (lexer *);
	virtual ~nhflow_poisson_c();

	virtual void start(lexer *,fdm_nhf*,double*);

private:

	int count,n,q;
    double teta;
    
};

#endif


