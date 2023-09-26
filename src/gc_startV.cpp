/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"lexer.h"
#include"fdm_nhf.h"
#include"sliceint.h"
#include"field.h"
#include"ghostcell.h"

void ghostcell::start1V(lexer *p, double *f, int gcv)
{
    //  MPI Boundary Swap
    starttime=timer();
    gcparaxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    p->xtime+=timer()-starttime;
    
    // 10 U
    // 11 V
    // 12 W
    // 14 E
    starttime=timer();
    ULOOP
    {  
    // s
        // Fx
        if(p->flag1[Im1JK]<0 && gcv==10 && p->B98<3)
        {
        f[Im1JK] = 0.5*fabs(p->W22)*d->eta(i-1,j)*d->eta(i-1,j) + fabs(p->W22)*d->eta(i-1,j)*d->dfx(i,j);
        }
        
        if(p->flag1[Im1JK]<0 && gcv==10 && p->B98>=3)
        {
        f[Im1JK] = d->UH[Im1JK]*d->U[Im1JK] + 0.5*fabs(p->W22)*d->eta(i-1,j)*d->eta(i-1,j) + fabs(p->W22)*d->eta(i-1,j)*d->dfx(i,j);
        }
        
        // Gx
        if(p->flag1[Im1JK]<0 && gcv==11 && p->B98<3)
        {
        f[Im1JK] = 0.0;
        }
        
        if(p->flag1[Im1JK]<0 && gcv==11 && p->B98>=3)
        {
        f[Im1JK] = d->VH[Im1JK]*d->U[Im1JK];
        }
        
        // Hx
        if(p->flag1[Im1JK]<0 && gcv==12 && p->B98<3)
        {
        f[Im1JK] = 0.0;
        }
        
        if(p->flag1[Im1JK]<0 && gcv==12 && p->B98>=3)
        {
        f[Im1JK] = d->WH[Im1JK]*d->U[Im1JK];
        }
         
        // Ex
        if(p->flag1[Im1JK]<0 && gcv==14 && p->B98<3)
        {
        f[Im1JK] = 0.0;
        }
        
        if(p->flag1[Im1JK]<0 && gcv==14 && p->B98>=3)
        {
        f[Im1JK] = d->UH[Im1JK];
        }
        
        
    // n
        // Fx
        if(p->flag1[Ip1JK]<0 && gcv==10 && p->B99<3)
        {
        f[Ip1JK] = 0.5*fabs(p->W22)*d->eta(i+1,j)*d->eta(i+1,j) + fabs(p->W22)*d->eta(i+1,j)*d->dfx(i,j);
        }
        
        if(p->flag1[Ip1JK]<0 && gcv==10 && p->B99>=3)
        {
        f[Ip1JK] = d->UH[Ip2JK]*d->U[Ip2JK] + 0.5*fabs(p->W22)*d->eta(i+1,j)*d->eta(i+1,j) + fabs(p->W22)*d->eta(i+1,j)*d->dfx(i,j);
        }
        
        // Gx
        if(p->flag1[Ip1JK]<0 && gcv==11 && p->B99<3)
        {
        f[Ip1JK] = 0.0;
        }
        
        if(p->flag1[Ip1JK]<0 && gcv==11 && p->B99>=3)
        {
        f[Ip1JK] = d->VH[Ip1JK]*d->U[Ip1JK];
        }
        
        // Hx
        if(p->flag1[Ip1JK]<0 && gcv==14 && p->B99<3)
        {
        f[Ip1JK] = 0.0;
        }
        
        if(p->flag1[Ip1JK]<0 && gcv==12 && p->B99>=3)
        {
        f[Ip1JK] = d->WH[Ip1JK]*d->U[Ip1JK];
        }
        
        // Ex
        if(p->flag1[Ip1JK]<0 && gcv==14 && p->B99<3)
        {
        f[Ip1JK] = 0.0;
        }
        
        if(p->flag1[Ip1JK]<0 && gcv==14 && p->B99>=3)
        {
        f[Ip1JK] = d->UH[Ip1JK];
        }
        
    // e
        if(p->flag1[IJm1K]<0 && p->j_dir==1)
        {
        f[IJm1K] = 0.0;
        }
    
    // w
        if(p->flag1[IJp1K]<0 && p->j_dir==1)
        {
        f[IJp1K] = 0.0;
        }
        
    // b
        if(p->flag1[IJKm1]<0)
        {
        f[IJKm1] = 0.0;
        }
        
    // t
        if(p->flag1[IJKp1]<0)
        {
        f[IJKp1] = 0.0;
        }
    }
    p->gctime+=timer()-starttime;
}

void ghostcell::start2V(lexer *p, double *f, int gcv)
{
    //  MPI Boundary Swap
    starttime=timer();
    gcparaxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    p->xtime+=timer()-starttime;
    
    starttime=timer();
    VLOOP
    {  
    // s
        if(p->flag2[Im1JK]<0)
        {
        f[Im1JK] = 0.0;
        }
          
    // n
        if(p->flag2[Ip1JK]<0)
        {
        f[Ip1JK] = 0.0;
        }
        
    // e
        if(p->flag2[IJm1K]<0 &&  gcv==11 && p->j_dir==1)
        {
        f[IJm1K] = 0.5*fabs(p->W22)*d->eta(i,j-1)*d->eta(i,j-1) + fabs(p->W22)*d->eta(i,j-1)*d->dfy(i,j);
        }
        
        if(p->flag2[IJm1K]<0 &&  gcv==14 && p->j_dir==1)
        {
        f[IJm1K] = 0.0;
        }
        
        if(p->flag2[IJm1K]<0 && gcv!=11 && gcv!=14 && p->j_dir==1)
        {
        f[IJm1K] = 0.0;
        }
        
    // w
        if(p->flag2[IJp1K]<0 &&  gcv==11 && p->j_dir==1)
        {
        f[IJp1K] = 0.5*fabs(p->W22)*d->eta(i,j+1)*d->eta(i,j+1) + fabs(p->W22)*d->eta(i,j+1)*d->dfy(i,j);
        }
        
        if(p->flag2[IJp1K]<0 &&  gcv==14 && p->j_dir==1)
        {
        f[IJp1K] = 0.0;
        }
        
        if(p->flag2[IJp1K]<0 && gcv!=11 && gcv!=14 && p->j_dir==1)
        {
        f[IJp1K] = 0.0;
        }
        
    // b
        if(p->flag2[IJKm1]<0 && p->j_dir==1)
        {
        f[IJKm1] = 0.0;
        }
        
    // t
        if(p->flag2[IJKp1]<0 && p->j_dir==1)
        {
        f[IJKp1] = 0.0;
        }
    }
    p->gctime+=timer()-starttime;
}

void ghostcell::start3V(lexer *p, double *f, int gcv)
{
    starttime=timer();
    gcparaxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    gcparacoxV1(p, f, gcv);
    p->xtime+=timer()-starttime;

    starttime=timer();
    WLOOP
    {  
        if(p->flag3[Im1JK]<0)
        {
        f[Im1JK] = 0.0;
        }
          
        if(p->flag3[Ip1JK]<0)
        {
        f[Ip1JK] = 0.0;
        }
        
        if(p->flag3[IJm1K]<0)
        {
        f[IJm1K] = 0.0;
        }
        
        if(p->flag3[IJp1K]<0)
        {
        f[IJp1K] = 0.0;
        }
        
        if(p->flag3[IJKm1]<0)
        {
        f[IJKm1] = 0.0;
        }
        
        if(p->flag3[IJKp1]<0)
        {
        f[IJKp1] = 0.0;
        }
    }
    p->gctime+=timer()-starttime;
}

void ghostcell::start4V(lexer *p, double *f, int gcv)
{
    starttime=timer();
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    p->xtime+=timer()-starttime;
    
    
    starttime=timer();
    LOOP
    {  

    // xxxxxxx
        if(p->flag4[Im1JK]<0 && gcv==10 && p->B98<3)
        {
        f[Im1JK] = 0.0;
        f[Im2JK] = 0.0;
        f[Im3JK] = 0.0;
        }
        
        if(p->flag4[Im1JK]<0 && gcv!=10 && p->B98<3)
        {
        f[Im1JK] = 0.0;
        f[Im2JK] = 0.0;
        f[Im3JK] = 0.0;
        }
          
        if(p->flag4[Ip1JK]<0 && gcv==10 && p->B99<3)
        {
        f[Ip1JK] = 0.0;
        f[Ip2JK] = 0.0;
        f[Ip3JK] = 0.0;
        }
        
        if(p->flag4[Ip1JK]<0 && gcv!=10 && p->B99<3)
        {
        f[Ip1JK] = 0.0;
        f[Ip2JK] = 0.0;
        f[Ip3JK] = 0.0;
        }
        
        
    // yyyyy
        if(p->flag4[IJm1K]<0 && p->j_dir==1 && gcv==11)
        {
        f[IJm1K] = 0.0;
        f[IJm2K] = 0.0;
        f[IJm3K] = 0.0;
        }
        
        if(p->flag4[IJm1K]<0 && p->j_dir==1 && gcv!=11)
        {
        f[IJm1K] = 0.0;
        f[IJm2K] = 0.0;
        f[IJm3K] = 0.0;
        }
        
        if(p->flag4[IJp1K]<0 && p->j_dir==1 && gcv==11)
        {
        f[IJp1K] = 0.0;
        f[IJp2K] = 0.0;
        f[IJp3K] = 0.0;
        }
        
        if(p->flag4[IJp1K]<0 && p->j_dir==1 && gcv!=11)
        {
        f[IJp1K] = 0.0;
        f[IJp2K] = 0.0;
        f[IJp3K] = 0.0;
        }
        
    // zzzzz
        if(p->flag4[IJKp1]<0 && gcv!=12)
        {
        f[IJKp1] = 0.0;
        f[IJKp2] = 0.0;
        f[IJKp3] = 0.0;
        }
        
        if(p->flag4[IJKm1]<0 && gcv!=12)
        {
        f[IJKm1] = 0.0;
        f[IJKm2] = 0.0;
        f[IJKm3] = 0.0;
        }
    }
    p->gctime+=timer()-starttime;
}

void ghostcell::start5V(lexer *p, double *f, int gcv)
{    
    starttime=timer();
	gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
	p->xtime+=timer()-starttime;
}

void ghostcell::start7V(lexer *p, double *f, sliceint &bc, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
	p->xtime+=timer()-starttime;
    }
    
    starttime=timer();
    if(gcv==250)
    fivec(p,f,bc);
    
    if(gcv==150)
    fivec2D(p,f,bc);
    
    if(gcv==210)
    fivec_vel(p,f,bc);
    
    if(gcv==110)
    fivec2D_vel(p,f,bc);
    p->gctime+=timer()-starttime;
}

void ghostcell::start7P(lexer *p, double *f, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
	p->xtime+=timer()-starttime;
    }
    
    /*
    starttime=timer();
    FLOOP
    {  
        //if(p->B98!=3||bc(i-1,j)==0)
    // xxxxxxx
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK];
        f[FIm2JK] = f[FIJK];
        f[FIm3JK] = f[FIJK];
        }

        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK];
        f[FIp2JK] = f[FIJK];
        f[FIp3JK] = f[FIJK];
        }
        
        
    // yyyyy
        if(p->flag7[FIJm1K]<0 && p->j_dir==1)
        {
        f[FIJm1K] = f[FIJK];
        f[FIJm2K] = f[FIJK];
        f[FIJm3K] = f[FIJK];
        }
        
        if(p->flag7[FIJp1K]<0 && p->j_dir==1)
        {
        f[FIJp1K] = f[FIJK];
        f[FIJp2K] = f[FIJK];
        f[FIJp3K] = f[FIJK];
        }
        
    // zzzzz
        if(p->flag7[FIJKp1]<0)
        {
        f[FIJK] = 0.0;
        f[FIJKp1] = 0.0;
        f[FIJKp2] = 0.0;
        f[FIJKp3] = 0.0;
        }
        
        if(p->flag7[FIJKm1]<0)
        {
        f[FIJKm1] = f[FIJK];
        f[FIJKm2] = f[FIJK];
        f[FIJKm3] = f[FIJK];
        }
    }
    p->gctime+=timer()-starttime;*/
}

void ghostcell::start7S(lexer *p, double *f, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
	p->xtime+=timer()-starttime;
    }
}

