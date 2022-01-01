/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"geotopo.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

int geotopo::conv(double a)
{

int b,c;
double d,diff;

c= int( a);
d=double(c);
diff=a-d;

b=c;

if(diff>0.5)
b=c+1;

if(diff<=-0.5)
b=c-1;


return b;

}
