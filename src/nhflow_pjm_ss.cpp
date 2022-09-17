
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
Author: Hans Bihs
--------------------------------------------------------------------*/

#define HX (fabs(a->hx(i,j))>1.0e-20?a->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(a->WL(i,j)+a->WL(i+1,j)))>1.0e-20?0.5*(a->WL(i,j)+a->WL(i+1,j)):1.0e20)
#define HY (fabs(a->hy(i,j))>1.0e-20?a->hy(i,j):1.0e20)

#include"pjm_sigss.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"hypre_struct.h"
#include"hypre_sstruct_fnpf.h"
 
pjm_sigss::pjm_sigss(lexer* p, fdm_nhf *d, ghostcell *pgc, heat *&pheat, concentration *&pconc)
{
    if((p->F80==0||p->A10==5||p->A10==55) && p->H10==0 && p->W30==0)
	pd = new density_f(p);
	
	if(p->F80==0 && p->H10==0 && p->W30==1)
	pd = new density_comp(p);
	
	if(p->F80==0 && p->H10>0)
	pd = new density_heat(p,pheat);
	
	if(p->F80==0 && p->C10>0)
	pd = new density_conc(p,pconc);

    gcval_press=540;  

	gcval_u=10;
	gcval_v=11;
	gcval_w=12;

    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Darray(M,vecsize*15);
    p->Darray(x,vecsize);
    p->Darray(rhs,vecsize);
    
    teta=0.5;
}

pjm_sigss::~pjm_sigss()
{
}

void pjm_sigss::start(fdm_nhf *d,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, ioflow *pflow, double *U, double *V, double *W, double alpha)
{    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
				
    rhscalc(p,a,pgc,uvel,vvel,wvel,alpha);
    
    if(p->j_dir==0)
    poisson2D(p,a,a->press);
    
    if(p->j_dir==1)
    poisson3D(p,a,a->press);
    
    fillvec(p,a,pgc);
	
        starttime=pgc->timer();

    psolv->startM(p,a,pgc,x,rhs,M,5);

        endtime=pgc->timer();
        
    fillvec_back(p,a,pgc);
    
	pgc->start4(p,a->press,gcval_press);
    
	ucorr(p,a,uvel,alpha);
	vcorr(p,a,vvel,alpha);
	wcorr(p,a,wvel,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
    
    pgc->start4(p,a->test,1);
}

void pjm_sigss::ucorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{	
	if(p->D37==1)
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,a,1,0,0))*((a->press(i+1,j,k)-a->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1))-0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(p->D37==2)
    ULOOP
    {       
     check=0;
    
        if(p->flag1[IJKp1]<0)
        check=1;        
    
    if(check==1)
    uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,a,1,0,0))*((a->press(i+1,j,k)-a->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*((1.0 - 1.0/teta)*(a->press(i,j,k)+a->press(i+1,j,k)))-0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,a,1,0,0))*((a->press(i+1,j,k)-a->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1))-0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void pjm_sigss::vcorr(lexer* p, fdm_nhf *d, double *V, double alpha)
{	 
    if(p->D37==1)
    VLOOP
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((a->press(i,j+1,k)-a->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(a->press(i,j,k+1)+a->press(i,j+1,k+1))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
                
                
    if(p->D37==2)
    VLOOP
    {       
     check=0;
    
        if(p->flag2[IJKp1]<0)
        check=1;        
    
    if(check==1)
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((a->press(i,j+1,k)-a->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(1.0 - 1.0/teta)*(a->press(i,j,k)+a->press(i,j+1,k))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((a->press(i,j+1,k)-a->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(a->press(i,j,k+1)+a->press(i,j+1,k+1))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void pjm_sigss::wcorr(lexer* p, fdm_nhf *d, double *W, double alpha)
{
    if(p->D37==1)
    WLOOP 	
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
    
    
    if(p->D37==2)
	WLOOP
    {
    check=0;
    
        if(p->flag3[IJKp1]<0)
        check=1;

    if(check==1)    
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*(((1.0 - 1.0/teta)*a->press(i,j,k)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
    
    if(check==0)
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
    }
}
 
void pjm_sigss::rhscalc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    NLOOP4
	rhs[n]=0.0;
	
    pip=p->Y50;

    n=0;
    KJILOOP
    {
    PCHECK
    rhs[n] =        -  ((u(i,j,k)-u(i-1,j,k))/p->DXN[IP]
                            + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(0.5*(u(i,j,k+1)+u(i-1,j,k+1))-0.5*(u(i,j,k-1)+u(i-1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                            
                            + (v(i,j,k)-v(i,j-1,k))/p->DYN[JP] 
                            + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(0.5*(v(i,j,k+1)+v(i,j-1,k+1))-0.5*(v(i,j,k-1)+v(i,j-1,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                           
                            + p->sigz[IJ]*(w(i,j,k)-w(i,j,k-1))/p->DZN[KP])/(alpha*p->dt);
                           
    SCHECK
    rhs[n] = 0.0;
                                        
    ++n;
    }
    pip=0;
}
 
void pjm_sigss::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, ouble *U, double *V, double *W, double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
}

void pjm_sigss::upgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->D38==1 && p->A540==1)
    ULOOP
	a->F(i,j,k) -= PORVAL1*fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DXP[IP];
    
    if(p->D38==1 && p->A540==2)
    ULOOP
	a->F(i,j,k) -= PORVAL1*fabs(p->W22)*(a->eta(i+1,j) - a->eta(i,j))/p->DXP[IP];
    
   
    if(p->D38==2 && p->A540==1)
    ULOOP
	a->F(i,j,k) -= PORVAL1*fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(eta(i+1,j),2.0) - pow(eta(i,j),2.0))/p->DXP[IP]
                    
                    + ((p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j))*a->depth(i+1,j) - (p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j))*a->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*((p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j)) + (p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j)))*(a->depth(i+1,j)-a->depth(i,j))/p->DXP[IP]);
    
    if(p->D38==2 && p->A540==2)
    ULOOP
	a->F(i,j,k) -= PORVAL1*fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(a->eta(i+1,j),2.0) - pow(a->eta(i,j),2.0))/p->DXP[IP]
                    
                    + (a->eta(i+1,j)*a->depth(i+1,j) - a->eta(i,j)*a->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*(a->eta(i,j) + a->eta(i+1,j))*(a->depth(i+1,j)-a->depth(i,j))/p->DXP[IP]);
    
    // fx = 1/2 g (eta^2 - 2* eta *z_b)
    // Sx = -g * eta * eta * Bx
}

void pjm_sigss::vpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
   if(p->D38==1 && p->A540==1)
    VLOOP
	a->G(i,j,k) -= PORVAL2*fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DYP[JP];
    
    if(p->D38==1 && p->A540==2)
    VLOOP
	a->G(i,j,k) -= PORVAL2*fabs(p->W22)*(a->eta(i,j+1) - a->eta(i,j))/p->DYP[JP];
}

void pjm_sigss::wpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
}

void pjm_sigss::fillvec(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    n=0;
    KJILOOP
    {
    PCHECK
    x[n] = a->press(i,j,k);
    
    SCHECK
    x[n] = 0.0;
    
    ++n;
    }
}

void pjm_sigss::fillvec_back(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    n=0;
    KJILOOP
    {
    PCHECK
    a->press(i,j,k) = x[n];
    
    ++n;
    }
}





