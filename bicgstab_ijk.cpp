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

#include"bicgstab_ijk.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver_void.h"
#include"jacobi_scaling.h"
#include"jacobi_block.h"
#include"sip.h"

bicgstab_ijk::bicgstab_ijk(lexer* p, fdm *a, ghostcell *pgc):epsi(1e-19)
{

    p->Darray(sj,p->imax*p->jmax*p->kmax);
    p->Darray(rj,p->imax*p->jmax*p->kmax);
    p->Darray(r0,p->imax*p->jmax*p->kmax);
    p->Darray(vj,p->imax*p->jmax*p->kmax);
    p->Darray(tj,p->imax*p->jmax*p->kmax);
    p->Darray(pj,p->imax*p->jmax*p->kmax); 
    p->Darray(ph,p->imax*p->jmax*p->kmax);
    p->Darray(sh,p->imax*p->jmax*p->kmax);
    p->Darray(aii,p->imax*p->jmax*p->kmax);
    p->Darray(x,p->imax*p->jmax*p->kmax);
    p->Darray(rhs,p->imax*p->jmax*p->kmax);

}

bicgstab_ijk::~bicgstab_ijk()
{
}

void bicgstab_ijk::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
}

void bicgstab_ijk::start(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{
	p->preconiter=0;
    
	
	if(var==1)
    {
    flag = p->flag1;
    ulast=p->ulast;
    vlast=0;
    wlast=0;
    }
	
	if(var==2)
    {
    flag = p->flag2;
    ulast=0;
    vlast=p->vlast;
    wlast=0;
    }
	
	if(var==3)
    {
    flag = p->flag3;
    ulast=0;
    vlast=0;
    wlast=p->wlast;
    }
	
	if(var==4||var==5)
    {
    flag = p->flag4;
    ulast=0;
    vlast=0;
    wlast=0;
    }
    
    fillxvec(p,a,f,rhsvec);
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,C);
	
	finalize(p,a,f);
}

void bicgstab_ijk::startF(lexer* p, fdm_fnpf* c, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var, int gcv, double stop_crit)
{
}
	
void bicgstab_ijk::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{
	solveriter=0;
	residual = 1.0e9;

	// -----------------
	precon_setup(p,a,pgc);
	// -----------------

 restart:
    r_j=norm_r0=0.0;	
	pgc->gcparaxijk_single(p,x,var);
	
	matvec_axb(p,a,x,rj);
	
	FLEXLOOP
	{
		r0[IJK]=pj[IJK]=rj[IJK];
		r_j += rj[IJK]*r0[IJK];
    }

    r_j=pgc->globalsum(r_j);
    norm_r0=sqrt(r_j);

    if((residual>=stop_crit) && (solveriter<maxiter))
	{

	do{
	    sigma=0.0;
	    norm_vj=0.0;
	    norm_rj=0.0;
		
		// -------------------------
		precon_solve(p,a,pgc,ph,pj);
		pgc->gcparaxijk_single(p,ph,var);				
		// -------------------------
		
		matvec_std(p,a,ph,vj);
		
		FLEXLOOP
		{
			sigma   += vj[IJK]*r0[IJK];
			norm_vj += vj[IJK]*vj[IJK];
			norm_rj += rj[IJK]*rj[IJK];
	    }
		
        sigma = pgc->globalsum(sigma);
		norm_vj = sqrt(pgc->globalsum(norm_vj));
		norm_rj = sqrt(pgc->globalsum(norm_rj));

	    alpha=r_j/sigma;

	if(fabs(sigma) <= (1.0e-12*(norm_vj*norm_r0)))
	{	
		residual=res_calc(p,a,pgc,x);
		++solveriter;

		goto restart;
	}

    if((fabs(alpha)*norm_vj/(norm_rj==0?1.0e-15:norm_rj))<=0.08)
	{
		residual=res_calc(p,a,pgc,x);
		++solveriter;

		goto restart;
	}

		norm_sj=0.0;
		
		FLEXLOOP
		{
		sj[IJK] = rj[IJK] - alpha*vj[IJK];
		norm_sj += sj[IJK]*sj[IJK];
		}

	    norm_sj=sqrt(pgc->globalsum(norm_sj));

    if(norm_sj>stop_crit)
	{
		// -------------------------
		precon_solve(p,a,pgc,sh,sj);
        pgc->gcparaxijk_single(p,sh,var);		
		// -------------------------

		matvec_std(p,a,sh,tj);
		
		w1=w2=0.0;
		
		FLEXLOOP
		{
		    w1 += tj[IJK]*sj[IJK];
		    w2 += tj[IJK]*tj[IJK];
		}

		w1=pgc->globalsum(w1);
		w2=pgc->globalsum(w2);

		w=w1/(w2==0?1.0e-15:w2);

		r_j1=0.0;
		
		FLEXLOOP
		{
		x[IJK] += alpha*ph[IJK] + w*sh[IJK];
		rj[IJK]  = sj[IJK]-w*tj[IJK];
		r_j1 += rj[IJK]*r0[IJK];
		}

		r_j1=pgc->globalsum(r_j1);

		beta=alpha*r_j1/(w*r_j==0?1.0e-15:(w*r_j));
		
		FLEXLOOP
		pj[IJK] = rj[IJK] + beta*(pj[IJK]-w*vj[IJK]);
	}


	if(norm_sj<=stop_crit)
	{
	r_j1=0.0;
		
		FLEXLOOP
		{
		x[IJK] += alpha*ph[IJK];
		rj[IJK]=sj[IJK];
		r_j1 += rj[IJK]*r0[IJK];
		}

    r_j1=pgc->globalsum(r_j1);
	}

	    r_j = r_j1 ;

	    residual=0.0;
		
		FLEXLOOP
		residual += rj[IJK]*rj[IJK];

	    residual = sqrt(pgc->globalsum(residual))/double(p->cellnumtot);
		
	    ++solveriter;

	}while((residual>=stop_crit) && (solveriter<maxiter));

    } 
		
	
	LOOP
	{
	ph[IJK]=0.0;
	sh[IJK]=0.0;
	}

}

void bicgstab_ijk::matvec_axb(lexer *p, fdm* a, double *x, double *y)
{
    n=0;
	FLEXLOOP
	{
	y[IJK]  = rhs[IJK]

			-(a->M.p[n]*x[IJK]
			+ a->M.n[n]*x[Ip1JK] 
			+ a->M.s[n]*x[Im1JK]
			+ a->M.w[n]*x[IJp1K]
			+ a->M.e[n]*x[IJm1K]
			+ a->M.t[n]*x[IJKp1]
			+ a->M.b[n]*x[IJKm1]);
    ++n;
	}
}

void bicgstab_ijk::matvec_std(lexer *p, fdm* a, double *x, double *y)
{
    n=0;
	FLEXLOOP
	{
	y[IJK]      = a->M.p[n]*x[IJK]
				+ a->M.n[n]*x[Ip1JK] 
				+ a->M.s[n]*x[Im1JK]
				+ a->M.w[n]*x[IJp1K]
				+ a->M.e[n]*x[IJm1K]
				+ a->M.t[n]*x[IJKp1]
				+ a->M.b[n]*x[IJKm1];
    ++n;
	}
}

double bicgstab_ijk::res_calc(lexer *p, fdm *a, ghostcell *pgc, double *x)
{
	double y;
	double resi=0.0;
    
    n=0;
	FLEXLOOP
	{	
	y  = rhs[IJK]

		-(a->M.p[n]*x[IJK]
		+ a->M.n[n]*x[Ip1JK] 
		+ a->M.s[n]*x[Im1JK]
		+ a->M.w[n]*x[IJp1K]
		+ a->M.e[n]*x[IJm1K]
		+ a->M.t[n]*x[IJKp1]
		+ a->M.b[n]*x[IJKm1]);

	resi+=y*y;
    
    ++n;
	}

	resi=sqrt(pgc->globalsum(resi));

	return resi/double(p->cellnumtot);	
}

void bicgstab_ijk::precon_setup(lexer* p,fdm* a, ghostcell* pgc)
{
    n=0;
	FLEXLOOP
    {
	aii[IJK]=-1.0/(a->M.p[n]+epsi);
    ++n;
    }
}

void bicgstab_ijk::precon_solve(lexer* p,fdm* a, ghostcell* pgc, double *f, double *b)
{
	FLEXLOOP
	f[IJK]=b[IJK]*aii[IJK];
}

void bicgstab_ijk::fillxvec(lexer* p, fdm* a, field& f, vec &rhsvec)
{
    n=0;
	FLEXLOOP
	{
	x[IJK] = f(i,j,k);
    
    /*
        if(flag[Im1JK]<0)
		x[Im1JK]=0.0;
		
		if(flag[Ip1JK]<0)
		x[Ip1JK]=0.0;
		
		if(flag[IJm1K]<0)
		x[IJm1K]=0.0;
		
		if(flag[IJp1K]<0)
		x[IJp1K]=0.0;
		
		if(flag[IJKm1]<0)
		x[IJKm1]=0.0;
		
		if(flag[IJKp1]<0)
		x[IJKp1]=0.0;*/
        
    rhs[IJK] = rhsvec.V[n];

    ++n;
    }
}


void bicgstab_ijk::finalize(lexer *p, fdm *a, field &f)
{  
        FLEXLOOP
        f(i,j,k)=x[IJK];
}

