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

#include"nhflow_momentum_RK2.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"nhflow_reconstruct.h"
#include"nhflow_fsf_reconstruct.h"
#include"nhflow_convection.h"
#include"nhflow_signal_speed.h"
#include"diffusion.h"
#include"nhflow_pressure.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"nhflow.h"
#include"nhflow.h"
#include"nhflow_fsf.h"
#include"nhflow_turbulence.h"
#include"vrans.h"
#include"nhflow_weno_flux.h"

nhflow_momentum_RK2::nhflow_momentum_RK2(lexer *p, fdm_nhf *d, ghostcell *pgc)
                                                    : bcmom(p), nhflow_sigma(p), etark1(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    
    p->Darray(UHRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHRK1,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(UDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WDIFF,p->imax*p->jmax*(p->kmax+2));
    
    pweno = new nhflow_weno_flux(p);
}

nhflow_momentum_RK2::~nhflow_momentum_RK2()
{
}

void nhflow_momentum_RK2::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, nhflow_signal_speed *pss, nhflow_fsf_reconstruct *pfsfrecon, 
                                     nhflow_reconstruct *precon, nhflow_convection *pconvec, diffusion *pdiff, 
                                     nhflow_pressure *ppress, solver *psolv, nhflow *pnhf, nhflow_fsf *pfsf, nhflow_turbulence *pnhfturb, vrans *pvrans)
{	
    pflow->discharge_nhflow(p,d,pgc);
    pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1);
		
//Step 1
//--------------------------------------------------------
    
    reconstruct(p,d,pgc,pss,pfsfrecon,precon,d->eta,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    pconvec->start(p,d,d->U,4,d->U,d->V,d->W,d->eta);
    pfsf->rk2_step1(p, d, pgc, pflow, d->U, d->V, d->W, etark1, etark1, 1.0);
    
    sigma_update(p,d,pgc,etark1,d->eta,1.0);
    omega_update(p,d,pgc,d->U,d->V,d->W);
    
    
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	pflow->isource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,d,etark1,d->eta);
	irhs(p,d,pgc);
	pconvec->start(p,d,d->UH,1,d->U,d->V,d->W,etark1);
    //pweno->start(p,d,d->U,1,d->U,d->V,d->W,etark1);
	//pdiff->diff_u(p,a,pgc,psolv,udiff,a->u,a->v,a->w,1.0);

	LOOP
	UHRK1[IJK] = d->UH[IJK]
				+ p->dt*CPORNH*d->F[IJK];

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	pflow->jsource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,d,etark1,d->eta);
	jrhs(p,d,pgc);
    pconvec->start(p,d,d->VH,2,d->U,d->V,d->W,etark1);
	//pdiff->diff_v(p,a,pgc,psolv,vdiff,a->u,a->v,a->w,1.0);

	LOOP
	VHRK1[IJK] = d->VH[IJK]
				+ p->dt*CPORNH*d->G[IJK];

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,d,etark1,d->eta);
	krhs(p,d,pgc);
	pconvec->start(p,d,d->WH,3,d->U,d->V,d->W,etark1);
    //pweno->start(p,d,d->W,3,d->U,d->V,d->W,etark1);
	//pdiff->diff_w(p,a,pgc,psolv,wdiff,a->u,a->v,a->w,1.0);

	LOOP
	WHRK1[IJK] = d->WH[IJK]
				+ p->dt*CPORNH*d->H[IJK];
	
    p->wtime=pgc->timer()-starttime;
    
    velcalc(p,d,pgc,UHRK1,VHRK1,WHRK1);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    //pflow->pressure_io(p,a,pgc);
	ppress->start(p,d,psolv,pgc,pflow,UHRK1,VHRK1,WHRK1,1.0);
	
    velcalc(p,d,pgc,UHRK1,VHRK1,WHRK1);

    pflow->U_relax(p,pgc,d->U,UHRK1);
    pflow->V_relax(p,pgc,d->V,VHRK1);
    pflow->W_relax(p,pgc,d->W,WHRK1);

	pflow->P_relax(p,pgc,d->P);

	pgc->start1V(p,UHRK1,gcval_u);
    pgc->start2V(p,VHRK1,gcval_v);
    pgc->start3V(p,WHRK1,gcval_w);
    
    clearrhs(p,d,pgc);
    
//Step 2
//--------------------------------------------------------
 
    reconstruct(p,d,pgc,pss,pfsfrecon,precon,etark1,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,etark1);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    pconvec->start(p,d,UHRK1,4,d->U,d->V,d->W,etark1);
    pfsf->rk2_step2(p, d, pgc, pflow, d->U,d->V,d->W, etark1, etark1, 0.5);
    
    sigma_update(p,d,pgc,d->eta,etark1,0.5);
    omega_update(p,d,pgc,d->U,d->V,d->W);
    
    
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	//pflow->isource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,d,d->eta,etark1);
	irhs(p,d,pgc);
    pconvec->start(p,d,UHRK1,1,d->U,d->V,d->W,d->eta);
    //pweno->start(p,d,d->U,1,d->U,d->V,d->W,etark1);
	//pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->UH[IJK] = 0.5*d->UH[IJK] + 0.5*UHRK1[IJK]
				+ 0.5*p->dt*CPORNH*d->F[IJK];
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	//pflow->jsource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,d,d->eta,etark1);
	jrhs(p,d,pgc);
	pconvec->start(p,d,VHRK1,2,d->U,d->V,d->W,d->eta);
	//pdiff->diff_v(p,a,pgc,psolv,vdiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->VH[IJK] = 0.5*d->VH[IJK] + 0.5*VHRK1[IJK]
				+ 0.5*p->dt*CPORNH*d->G[IJK];
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,d,d->eta,etark1);
	krhs(p,d,pgc);
	pconvec->start(p,d,WHRK1,3,d->U,d->V,d->W,d->eta);
    //pweno->start(p,d,d->W,3,d->U,d->V,d->W,d->eta);
	//pdiff->diff_w(p,a,pgc,psolv,wdiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->WH[IJK] = 0.5*d->WH[IJK] + 0.5*WHRK1[IJK]
				+ 0.5*p->dt*CPORNH*d->H[IJK];
	
    p->wtime+=pgc->timer()-starttime;
    
    velcalc(p,d,pgc,d->UH,d->VH,d->WH);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
	//pflow->pressure_io(p,a,pgc);
    ppress->start(p,d,psolv,pgc,pflow,d->UH,d->VH,d->WH,0.5);
    
    velcalc(p,d,pgc,d->UH,d->VH,d->WH);
	
	pflow->U_relax(p,pgc,d->U,d->UH);
    pflow->V_relax(p,pgc,d->V,d->VH);
    pflow->W_relax(p,pgc,d->W,d->WH);

	pflow->P_relax(p,pgc,d->P);

	pgc->start1V(p,d->UH,gcval_u);
    pgc->start2V(p,d->VH,gcval_v);
    pgc->start3V(p,d->WH,gcval_w);
    
    clearrhs(p,d,pgc);
    
    velcalc(p,d,pgc,d->UH,d->VH,d->WH);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
}

void nhflow_momentum_RK2::reconstruct(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_signal_speed *pss, nhflow_fsf_reconstruct *pfsfrecon, 
                                     nhflow_reconstruct *precon, slice &eta, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    // reconstruct eta
    pfsfrecon->reconstruct_2D(p, pgc, d, eta, d->ETAs, d->ETAn, d->ETAe, d->ETAw);
    pfsfrecon->reconstruct_2D_WL(p, pgc, d);
    
    // reconstruct U 
    precon->reconstruct_3D_x(p, pgc, d, U, d->Us, d->Un);
    precon->reconstruct_3D_y(p, pgc, d, U, d->Ue, d->Uw);
    precon->reconstruct_3D_z(p, pgc, d, U, d->Ub, d->Ut);
    
    // reconstruct  V
    precon->reconstruct_3D_x(p, pgc, d, V, d->Vs, d->Vn);
    precon->reconstruct_3D_y(p, pgc, d, V, d->Ve, d->Vw);
    precon->reconstruct_3D_z(p, pgc, d, V, d->Vb, d->Vt);
    
    // reconstruct  W
    precon->reconstruct_3D_x(p, pgc, d, W, d->Ws, d->Wn);
    precon->reconstruct_3D_y(p, pgc, d, W, d->We, d->Ww);
    precon->reconstruct_3D_z(p, pgc, d, W, d->Wb, d->Wt);
    
    // reconstruct UH
    precon->reconstruct_3D_x(p, pgc, d, UH, d->UHs, d->UHn);
    precon->reconstruct_3D_y(p, pgc, d, UH, d->UHe, d->UHw);
    
    // reconstruct  VH
    precon->reconstruct_3D_x(p, pgc, d, VH, d->VHs, d->VHn);
    precon->reconstruct_3D_y(p, pgc, d, VH, d->VHe, d->VHw);
    
    // reconstruct  WH
    precon->reconstruct_3D_x(p, pgc, d, WH, d->WHs, d->WHn);
    precon->reconstruct_3D_y(p, pgc, d, WH, d->WHe, d->WHw);
    
    pss->signal_speed_update(p, pgc, d, d->Us, d->Un, d->Ve, d->Vw, d->Ds, d->Dn, d->De, d->Dw);
}

void nhflow_momentum_RK2::velcalc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *UH, double *VH, double *WH)
{
    LOOP
    {
    d->U[IJK] = UH[IJK]/d->WL(i,j);
    d->V[IJK] = VH[IJK]/d->WL(i,j);
    d->W[IJK] = WH[IJK]/d->WL(i,j);       
    }
    
    pgc->start1V(p,d->U,gcval_u);
    pgc->start2V(p,d->V,gcval_v);
    pgc->start3V(p,d->W,gcval_w);
}

void nhflow_momentum_RK2::irhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
	n=0;
	LOOP
	{
    d->maxF=MAX(fabs(d->rhsvec.V[n]),d->maxF);
	d->F[IJK] += (d->rhsvec.V[n])*PORVALNH;
	d->rhsvec.V[n]=0.0;
	++n;
	}*/
}

void nhflow_momentum_RK2::jrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]),a->maxG);
	a->G[IJK] += (a->rhsvec.V[n])*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}*/
}

void nhflow_momentum_RK2::krhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
	n=0;
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n]),a->maxH);
    
    if(p->D38==0)
    a->H[IJK] += (a->rhsvec.V[n])*PORVAL3;
    
    if(p->D38>0)
	a->H[IJK] += (a->rhsvec.V[n])*PORVAL3;
    
    
	a->rhsvec.V[n]=0.0;
	++n;
	}
    */
}

void nhflow_momentum_RK2::clearrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
	n=0;
	LOOP
	{
	d->rhsvec.V[n]=0.0;
	++n;
	}
    
}

void nhflow_momentum_RK2::inidisc(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_fsf *pfsf)
{
    sigma_ini(p,d,pgc,d->eta);
    sigma_update(p,d,pgc,d->eta,d->eta,1.0);
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
}
     



