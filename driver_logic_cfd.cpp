/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"driver.h"
#include"lexer.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"6DOF_header.h"
#include"waves_header.h"

void driver::logic()
{
	makegrid_cds();
	pini = new initialize(p);

    if(p->mpirank==0)
	cout<<"starting ini"<<endl;
	pini->start(a,p,pgc);

	if(p->mpirank==0)
    cout<<"creating objects"<<endl;

// time stepping
    if(p->N48==0)
	ptstep=new fixtimestep(p);

    if(p->N40<=10 && (p->N48==1 || p->N48==3 || p->N48==4)  && (p->D20!=0&&p->D20!=2))
	ptstep=new etimestep(p);
	
	if(p->N40<=10 && (p->N48==1 || p->N48==3 || p->N48==4) && (p->D20==0||p->D20>=2))
	ptstep=new ietimestep(p);

	if(p->N40>10  && (p->N48==1 || p->N48==3 || p->N48==4))
	ptstep=new itimestep(p);
 

	if(p->N48==2)
	ptstep=new lsmtimestep(p);
	
// Multiphase
	if(p->F300==0)
	pmp = new multiphase_v();
	
	if(p->F300>0)
	pmp = new multiphase_f(p,a,pgc);


//discretization scheme

    //Convection	
	if(p->D10==0)
	pconvec=new convection_void(p);

	if(p->D10==1 && p->N40<=10)
	pconvec=new fou(p);

	if(p->D10==2 && p->N40<=10)
	pconvec=new cds2(p);

	if(p->D10==3 && p->N40<=10)
	pconvec=new quick(p);

	if(p->D10==4 && p->N40<=10)
	pconvec=new weno_flux_nug(p);
	
	if(p->D10==5 && p->N40<=10)
	pconvec=new weno_hj_nug(p);
	
	if(p->D10==6 && p->N40<=10)
	pconvec=new cds4(p);
    
    if(p->D10==7 && p->N40<=10)
	pconvec=new weno3_flux(p);
    
    if(p->D10==8 && p->N40<=10)
	pconvec=new weno3_hj(p);
    
    if(p->D10==9 && p->N40<=10)
	pconvec=new weno_flux(p);
	
	if(p->D10>=10 && p->D10<30 && p->N40<=10)
	pconvec=new hires(p,p->D10);


	if(p->D10==1 && p->N40>10)
	pconvec=new ifou(p);

	if(p->D10==2 && p->N40>10)
	pconvec=new icds2(p);

	if(p->D10==3 && p->N40>10)
	pconvec=new iquick(p);

	if(p->D10==4 && p->N40>10)
	pconvec=new iweno_flux(p);

	if(p->D10==5 && p->N40>10)
	pconvec=new iweno_hj(p);
	
	if(p->D10==6 && p->N40>10)
	pconvec=new icds4(p);
	
	if(p->D10>=10 && p->D10<30 && p->N40>10)
	pconvec=new ihires(p,p->D10);
	

	// Convection Turbulence
	if(p->T12==0)
	pturbdisc=new convection_void(p);

	if(p->T12==1 && p->T11<=10)
	pturbdisc=new fou(p);

	if(p->T12==2 && p->T11<=10)
	pturbdisc=new cds2(p);

	if(p->T12==3 && p->T11<=10)
	pturbdisc=new quick(p);

	if(p->T12==4 && p->T11<=10)
	pturbdisc=new weno_flux(p);

	if(p->T12==5 && p->T11<=10)
	pturbdisc=new weno_hj(p);
	
	if(p->T12==6 && p->T11<=10)
	pturbdisc=new cds4(p);
	
	if(p->T12>=10 && p->T12<30 && p->T11<=10)
	pturbdisc=new hires(p,p->T12);


	if(p->T12==1 && p->T11>10)
	pturbdisc=new ifou(p);

	if(p->T12==2 && p->T11>10)
	pturbdisc=new icds2(p);

	if(p->T12==3 && p->T11>10)
	pturbdisc=new iquick(p);

	if(p->T12==4 && p->T11>10)
	pturbdisc=new iweno_flux(p);

	if(p->T12==5 && p->T11>10)
	pturbdisc=new iweno_hj(p);
	
	if(p->T12==6 && p->T11>10)
	pturbdisc=new icds4(p);

	if(p->T12>=10 && p->T12<30 && p->T11>10)
	pturbdisc=new ihires(p,p->T12);
	
	//  Convection FSF
	if(p->F35==0&&p->F85==0)
	pfsfdisc=new convection_void(p);
	
	if(p->F35==1 && p->F30<=10)
	pfsfdisc=new fou(p);

	if(p->F35==2 && p->F30<=10)
	pfsfdisc=new cds2_alt(p);

	if(p->F35==3 && p->F30<=10)
	pfsfdisc=new quick(p);

	if(p->F35==4 && p->F30<=10)
	pfsfdisc=new weno_flux(p);

	if(p->F35==5 && p->F30<=10)
	pfsfdisc=new weno_hj_nug(p);
	
	if(p->F35==6 && p->F30<=10)
	pfsfdisc=new cds4(p);
    
    if(p->F35==7 && p->F30<=10)
	pfsfdisc=new weno3_flux(p);
    
    if(p->F35==8 && p->F30<=10)
	pfsfdisc=new weno3_hj(p);
	
	if(p->F35>=10 && p->F35<30 && p->F30<=10)
	pfsfdisc=new hires(p,p->F35);
	
	if(p->F35>=40 && p->F35<50 && p->F30<=10)
	pfsfdisc=new hires(p,p->F35);


	if(p->F35==1 && p->F30>10)
	pfsfdisc=new ifou(p);

	if(p->F35==2 && p->F30>10)
	pfsfdisc=new icds2(p);

	if(p->F35==3 && p->F30>10)
	pfsfdisc=new iquick(p);

	if(p->F35==4 && p->F30>10)
	pfsfdisc=new iweno_flux(p);

	if(p->F35==5 && p->F30>10)
	pfsfdisc=new iweno_hj(p);
	
	if(p->F35==6 && p->F30>10)
	pfsfdisc=new icds4(p);
	
	if(p->F35>=10 && p->F35<30 && p->F30>10)
	pfsfdisc=new ihires(p,p->F35);
	
	
//  Convection Multiphase LSM
	if(p->F305==0)
	pmpconvec=new convection_void(p);
	
	if(p->F305==1 && p->F300<=10)
	pmpconvec=new fou(p);

	if(p->F305==2 && p->F300<=10)
	pmpconvec=new cds2(p);

	if(p->F305==3 && p->F300<=10)
	pmpconvec=new quick(p);

	if(p->F305==4 && p->F300<=10)
	pmpconvec=new weno_flux(p);

	if(p->F305==5 && p->F300<=10)
	pmpconvec=new weno_hj(p);
	
	if(p->F305==6 && p->F300<=10)
	pmpconvec=new cds4(p);
	
	if(p->F305>=10 && p->F305<30 && p->F300<=10)
	pmpconvec=new hires(p,p->F305);
	
	if(p->F305>=40 && p->F305<50 && p->F300<=10)
	pmpconvec=new hires(p,p->F305);


	if(p->F305==1 && p->F300>10)
	pmpconvec=new ifou(p);

	if(p->F305==2 && p->F300>10)
	pmpconvec=new icds2(p);

	if(p->F305==3 && p->F300>10)
	pmpconvec=new iquick(p);

	if(p->F305==4 && p->F300>10)
	pmpconvec=new iweno_flux(p);

	if(p->F305==5 && p->F300>10)
	pmpconvec=new iweno_hj(p);
	
	if(p->F305==6 && p->F300>10)
	pmpconvec=new icds4(p);
	
	if(p->F305>=10 && p->F305<30 && p->F300>10)
	pmpconvec=new ihires(p,p->F305);
	
//  Convection Concentration
	if(p->C15==0)
	pconcdisc=new convection_void(p);
	
	if(p->C15==1 && p->C10<=10)
	pconcdisc=new fou(p);

	if(p->C15==2 && p->C10<=10)
	pconcdisc=new cds2_alt(p);

	if(p->C15==3 && p->C10<=10)
	pconcdisc=new quick(p);

	if(p->C15==4 && p->C10<=10)
	pconcdisc=new weno_flux(p);

	if(p->C15==5 && p->C10<=10)
	pconcdisc=new weno_hj(p);
	
	if(p->C15==6 && p->C10<=10)
	pconcdisc=new cds4(p);
	
	if(p->C15>=10 && p->C15<30 && p->C10<=10)
	pconcdisc=new hires(p,p->C15);
	
	if(p->C15>=40 && p->C15<50 && p->C10<=10)
	pconcdisc=new hires(p,p->C15);


	if(p->C15==1 && p->C10>10)
	pconcdisc=new ifou(p);

	if(p->C15==2 && p->C10>10)
	pconcdisc=new icds2(p);

	if(p->C15==3 && p->C10>10)
	pconcdisc=new iquick(p);

	if(p->C15==4 && p->C10>10)
	pconcdisc=new iweno_flux(p);

	if(p->C15==5 && p->C10>10)
	pconcdisc=new iweno_hj(p);
	
	if(p->C15==6 && p->C10>10)
	pconcdisc=new icds4(p);
	
	if(p->C15>=10 && p->C15<30 && p->C10>10)
	pconcdisc=new ihires(p,p->C15);
	
	if(p->S60==11 || p->S60==12)
	pconcdisc=new iweno_hj(p);
	
	if(p->S60>0&&p->S60<10)
	pconcdisc=new weno_hj(p);

//Diffusion
	// momentum and scalars
	if(p->D20==0)
	pdiff=new diff_void;

	if(p->D20==1 && p->N40<=10)
	pdiff=new ediff2(p);
	
	if(p->D20>=2 && p->N40<=10)
	pdiff=new idiff2_FS(p);

	if(p->N40>10)
	pdiff=new idiff2(p);
	
	// turbulence
	if(p->D20==0 || p->T10==0)
	pturbdiff=new diff_void;

	if(p->D20==1 && p->T11<=10 && p->T10>0)
	pturbdiff=new ediff2(p);
	
	if(p->D20>=2 && p->T11<=10 && p->T10>0)
	pturbdiff=new idiff2_FS(p);
	
	if(p->T11>10)
	pturbdiff=new idiff2(p);	
	
	// concentration
	if(p->D20==0 || p->C10==0)
	pconcdiff=new diff_void;

	if(p->D20==1 && p->C10<=10 && p->C10>0)
	pconcdiff=new ediff2(p);
	
	if(p->D20>=2 && p->C10<=10 && p->C10>0)
	pconcdiff=new idiff2_FS(p);
	
	if(p->C10>10)
	pconcdiff=new idiff2(p);
	
	// susepdended 
	if(p->S60<11 && p->S60>0)
	psuspdiff=new ediff2(p);
	
	if(p->D20==1 && p->S60<=10 && p->S60>0)
	psuspdiff=new ediff2(p);
	
	if(p->D20>=2 && p->S60<=10 && p->S60>0)
	psuspdiff=new idiff2_FS(p);
	
	if(p->S60>10)
	psuspdiff=new idiff2(p);
	


//turbulence model
	if(p->T10==0)
	pturb = new kepsilon_void(p,a,pgc);

	//ke
	if((p->T10==1 || p->T10==21) && p->T11==1)
	pturb = new kepsilon_AB(p,a,pgc);

	if((p->T10==1 || p->T10==21) && p->T11==2)
	pturb = new kepsilon_RK3(p,a,pgc);

	if((p->T10==1 || p->T10==21) && p->T11==3)
	pturb = new kepsilon_RK3(p,a,pgc);

	if((p->T10==1 || p->T10==21) && p->T11==11)
	pturb = new kepsilon_IM1(p,a,pgc);

	if((p->T10==1 || p->T10==21) && p->T11==12)
	pturb = new kepsilon_IM2(p,a,pgc);

	if(p->T10==11 && p->T11==1)
	pturb = new wallin_AB(p,a,pgc);

	if(p->T10==11 && p->T11==3)
	pturb =new wallin_RK3(p,a,pgc);

	if(p->T10==11 && p->T11==11)
	pturb = new wallin_IM1(p,a,pgc);

	if(p->T10==11 && p->T11==12)
	pturb = new wallin_IM2(p,a,pgc);

    //kw
    if((p->T10==2 || p->T10==22) && p->T11==1)
	pturb = new komega_AB(p,a,pgc);

	if((p->T10==2 || p->T10==22) && p->T11==2)
	pturb = new komega_RK2(p,a,pgc);

	if((p->T10==2 || p->T10==22) && p->T11==3)
	pturb = new komega_RK3(p,a,pgc);

	if((p->T10==2 || p->T10==22) && p->T11==11)
	pturb = new komega_IM1(p,a,pgc,pmp);

	if((p->T10==2 || p->T10==22) && p->T11==12)
	pturb = new komega_IM2(p,a,pgc,pmp);

	if(p->T10==12 && p->T11==1)
	pturb =new wallin_ABkw(p,a,pgc);

	if(p->T10==12 && p->T11==3)
	pturb =new wallin_RK3kw(p,a,pgc);

	if(p->T10==12 && p->T11==11)
	pturb =new wallin_IM1kw(p,a,pgc,pmp);

	if(p->T10==12 && p->T11==12)
	pturb =new wallin_IM2kw(p,a,pgc,pmp);

    // LES
	if(p->T10==31)
	pturb =new LES_smagorinsky(p,a);

	if(p->T10==32)
	pturb =new LES_germano(p,a);

//Heat
    if(p->H10==0)
	pheat =  new heat_void(p,a,pgc);

	if(p->H10==1)
	pheat =  new heat_AB(p,a,pgc,pheat);

	if(p->H10==2)
	pheat =  new heat_RK2(p,a,pgc,pheat);

	if(p->H10==3)
	pheat =  new heat_RK3(p,a,pgc,pheat);

	if(p->H10==11)
	pheat =  new heat_IM1(p,a,pgc,pheat);

	if(p->H10==12)
	pheat =  new heat_IM2(p,a,pgc,pheat);
	
// Concentration
    if(p->C10==0 && p->F101==0)
	pconc =  new concentration_void(p,a,pgc);

	if(p->C10==1)
	pconc =  new concentration_AB(p,a,pgc);

	if(p->C10==2)
	pconc =  new concentration_RK2(p,a,pgc);

	if(p->C10==3)
	pconc =  new concentration_RK3(p,a,pgc);

	if(p->C10==11)
	pconc =  new concentration_IM1(p,a,pgc);

	if(p->C10==12)
	pconc =  new concentration_IM2(p,a,pgc);
	
// Air Entrainment
    if(p->F101==1)
	{
	pconc =  new air_entrainment_IM1(p,a,pgc);
	pconcdisc=new iweno_hj(p);
	pconcdiff=new idiff2(p);
	}
    
// Wave Models
    if(p->A10==5 || p->A10==0)
    pnse = new nsewave_v(p,a,pgc,pheat,pconc);
    
    if(p->A10==4)
    {
    if(p->A410==1)
    pnse = new nsewave_f(p,a,pgc,pheat,pconc);
    
    if(p->A410==2)
    pnse = new nsewave_geo(p,a,pgc,pheat,pconc);
    }
    

// Free Surface
    if(p->F10==1)
    poneph = new onephase_f(p,a,pgc);
    
    if(p->F10==2)
    poneph = new onephase_v(p,a,pgc);
    
    if((p->F30==0 && p->F80==0) || p->F11==1)
	pfsf = new levelset_void(p,a,pgc);

	if(p->F30==1)
	pfsf = new levelset_AB2(p,a,pgc,pheat,pconc);

	if(p->F30==2)
	pfsf = new levelset_RK2(p,a,pgc,pheat,pconc);

	if(p->F30==3 && p->F11==0)
	pfsf = new levelset_RK3(p,a,pgc,pheat,pconc);
	
	if(p->F30==4)
	pfsf = new levelset_RK4(p,a,pgc,pheat,pconc);	
	
	if(p->F30==5)
	pfsf = new levelset_AB3(p,a,pgc,pheat,pconc);

	if(p->F30==11)
	pfsf = new levelset_IM1(p,a,pgc,pheat,pconc,a->phi);

	if(p->F30==12)
	pfsf = new levelset_IM2(p,a,pgc,pheat,pconc,a->phi);
	
	if(p->F40==0)
	preini = new reini_void(p);
    
    if(p->F40==1)
	preini = new reinifluid_AB2(p,a);
	
	if(p->F40==2)
	preini = new reinifluid_AB3(p,a);
	
    if(p->F40==3||p->F40==33)
    preini = new reinifluid_RK3(p,1);
	
	if(p->F40==4)
	preini = new reinifluid_RK4(p,a);
    
	if(p->F40==21)
	preini = new reini_AB2(p,a);
	
	if(p->F40==22)
	preini = new reini_AB3(p,a);
	
	if(p->F40==23)
	preini = new reini_RK3(p,1);
	
	if(p->F40==14)
	preini = new reini_RK4(p,a);
	
	if(p->F40==5)
	preini = new reinivc_RK3(p);
	
	if(p->F40==7)
	preini = new reinigc_RK3(p,a);
	
	if(p->F40==8)
	preini = new reinigc_RK4(p,a);
	

	if(p->F40==11 || p->F40==13 || p->F40==14)
	preini = new directreini(p,a);



	if(p->F31==0)
	ppart = new particle_void();

	if(p->F31==1 || p->F31==2)
	ppart = new particle(p,a,pgc);
	
	
	if(p->F80==1)
	pfsf = new VOF_AB(p,a,pgc,pheat);

	if(p->F80==3)
	pfsf = new VOF_RK3(p,a,pgc,pheat);
	
	if(p->F80==4)
	pfsf = new VOF_PLIC(p,a,pgc,pheat);
	
	if(p->F80==11)
	pfsf = new VOF_IM1(p,a,pgc,pheat);
    
    
    //  Convection VOF
	if(p->F85==0 && p->F35==0)
	pfsfdisc=new convection_void(p);
	
	if(p->F85==1 && p->F80<=10)
	pfsfdisc=new fou(p);

	if(p->F85==2 && p->F80<=10)
	pfsfdisc=new cds2_alt(p);

	if(p->F85==3 && p->F80<=10)
	pfsfdisc=new quick(p);

	if(p->F85==4 && p->F80<=10)
	pfsfdisc=new weno_flux(p);

	if(p->F85==5 && p->F80<=10)
	pfsfdisc=new weno_hj(p);
	
	if(p->F85==6 && p->F80<=10)
	pfsfdisc=new cds4(p);
	
	if(p->F85>=10 && p->F85<30 && p->F80<=10)
	pfsfdisc=new hires(p,p->F85);
	
	if(p->F85>=40 && p->F85<50 && p->F80<=10)
	pfsfdisc=new hires(p,p->F85);
    
    if(p->F85==51 && p->F80<=10)
	pfsfdisc=new hric(p);
	
	if(p->F85==52 && p->F80<=10)
	pfsfdisc=new hric_mod(p);
	
	if(p->F85==53 && p->F80<=10)
	pfsfdisc=new cicsam(p);


	if(p->F85==1 && p->F80>10)
	pfsfdisc=new ifou(p);

	if(p->F85==2 && p->F80>10)
	pfsfdisc=new icds2(p);

	if(p->F85==3 && p->F80>10)
	pfsfdisc=new iquick(p);

	if(p->F85==4 && p->F80>10)
	pfsfdisc=new iweno_flux(p);

	if(p->F85==5 && p->F80>10)
	pfsfdisc=new iweno_hj(p);
	
	if(p->F85==6 && p->F80>10)
	pfsfdisc=new icds4(p);
	
	if(p->F85>=10 && p->F85<30 && p->F80>10)
	pfsfdisc=new ihires(p,p->F85);
	

//pressure scheme
	if(p->D30==0)
	ppress = new pressure_void(p);

	if(p->D30==1 && p->W30==0 && p->F10==2)
	ppress = new pjm(p,a);
    
    if(p->D30==1 && p->W30==1 && p->F10==2)
	ppress = new pjm_comp(p,a,pgc);
    
    if(p->D30==1 && p->F10==1)
	ppress = new pjm_nse(p,a);
    
    if(p->D30==2)
	ppress = new pjm_fsm(p,a);
    
    if(p->D30==3)
	ppress = new pjm_corr(p,a);
	
	if(p->D30==4)
	ppress = new pjm_fsi(p,a);
    
    if(p->D30==5)
	ppress = new pjm_4th(p,a);

	if(p->D30==11)
	ppress = new simple(p,a);

	if(p->D30==12)
	ppress = new piso(p,a);
	

//poisson scheme for pressure
	if(p->D30<5 && p->F10==2)
	ppois = new poisson_f(p);
    
    if(p->D30==5 && p->F10==2)
	ppois = new poisson_f(p);
    
    if(p->D30<9 && p->F10==1)
	ppois = new poisson_nse(p);

	if(p->D30>10 && p->D30<20)
	ppois = new presscorr(p);
	

//Solver
	if(p->N8==0)
	psolv = new solver_void(p,a,pgc);
	
	if(p->N8==1)
	psolv = new jacobi_block(p,a,pgc);
	
	if(p->N8==2)
	psolv = new sip(p,a,pgc);
	
	if(p->N8==3)
	psolv = new bicgstab(p,a,pgc,p->N9);
	
    
//Poison Solver	
	if(p->N10==0)
	ppoissonsolv = new solver_void(p,a,pgc);
	
	if(p->N10==1)
	ppoissonsolv = new jacobi_block(p,a,pgc);
	
	if(p->N10==2)
	ppoissonsolv = new sip(p,a,pgc);
	
	if(p->N10==3)
	ppoissonsolv = new bicgstab(p,a,pgc,p->N11);
	
	#ifdef HYPRE_COMPILATION
	if(p->N10>10 && p->N10<=20)
	ppoissonsolv = new hypre_struct(p,a,pgc);
	#endif
    
    #ifdef HYPRE_COMPILATION
	if(p->N10>20 && p->N10<=30)
	ppoissonsolv = new hypre_aij(p,a,pgc);
	#endif

    //psolv=ppoissonsolv;

//IOFlow
	if(p->B60==0 && p->B90==0 && p->B180==0 )
	pflow = new ioflow_v(p,pgc);

	if(p->B60>=1)
	pflow = new ioflow_f(p,pgc);

	if(p->B90>=1)
	pflow= new iowave(p,pgc);

	if(p->B180==1||p->B191==1||p->B192==1)
	pflow= new ioflow_gravity(p,pgc);

//Potential Flow Solver
    if(p->I11==0)
    potflow = new potential_v;

    if(p->I11==1 )
    potflow = new potential_f(p);

// Benchmark
    if(p->F150==0)
    pbench = new benchmark_void;

    if(p->F150==1)
    pbench = new benchmark_vortex(p,a);
	
	if(p->F150==2)
    pbench = new benchmark_disk(p,a);
	
	if(p->F150==3)
    pbench = new benchmark_vortex3D(p,a);
	
	if(p->F150==11)
    pbench = new benchmark_convection(p,a);

// Printer
	pprint = new vtu3D(p,a,pgc);
	
// Sediment
    if(p->S10==0)
    psed = new sediment_void();

    if(p->S10>0)
    psed = new sediment_f(p,pturb);

    if(p->S11==0)
    pbed = new bedload_void();

    if(p->S11==1)
    pbed = new bedload_VR(p,pturb);

    if(p->S11==2)
    pbed = new bedload_MPM(p,pturb);
	
	if(p->S11==3)
    pbed = new bedload_EF(p,pturb);
    
    if(p->S11==4)
    pbed = new bedload_einstein(p,pturb);

    if(p->S10==0)
    ptopo = new topo_void(p,a,pgc);
	
	if(p->S10==1)
    ptopo = new topo_direct(p,a,pgc,pturb);

    if(p->S10==0 && p->G1==0)
    preto = new reinitopo_void();

    if(p->S10==1 || p->G1==1)
    {
    if(p->G40==0)
    preto = new reinitopo_void();
    
    if(p->G40==1)
    preto = new reinitopo_AB2(p);
    
    if(p->G40==3)
    preto = new reinitopo_RK3(p);
    }

    if(p->S60==0)
    psusp = new suspended_void();

    if(p->S60==1)
    psusp = new suspended_AB(p,a,pturb);

    if(p->S60==2)
    psusp = new suspended_RK2(p,a,pturb);

    if(p->S60==3)
    psusp = new suspended_RK3(p,a,pturb);

    if(p->S60==11)
    psusp = new suspended_IM1(p,a,pturb);

    if(p->S60==12)
    psusp = new suspended_IM2(p,a,pturb);

// Velocities
	if(p->N40==1)
	pmom = new momentum_AB2(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);
	
	if(p->N40==2)
	pmom = new momentum_RK2(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);

	if(p->N40==3 && p->F11==0)
	pmom = new momentum_RK3(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);
    
    if(p->N40==3 && p->F11==1 && p->F30==3)
	pmom = new momentum_FSFC_RK3(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pfsfdisc,preini,pheat,pconc);

	if(p->N40==4)
	pmom = new momentum_RK4(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);

	if(p->N40==6 && p->F11==0)
	pmom = new momentum_FS3(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);
    
     if(p->N40==6 && p->F11==1 && p->F30==3)
	pmom = new momentum_FSFC_RK3(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pfsfdisc,preini,pheat,pconc);

	if(p->N40==7)
	pmom = new momentum_FS4(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);
    
    if(p->N40==8)
	pmom = new momentum_MK3(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);
    
    if(p->N40==9)
	pmom = new momentum_MF3(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);

	if(p->N40==11)
	pmom = new momentum_IM1(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);

	if(p->N40==12)
	pmom = new momentum_IM2(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);

	if(p->N40==0 && p->X10==1 && p->X13==1)
	pmom = new momentum_FSI(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);
	
	if(p->N40==0 && p->X13==0)
	pmom = new momentum_void();	
	
// 6DOF
	if(p->X10==0)
    p6dof = new sixdof_void;
	
	if(p->X10==1)
    p6dof = new sixdof_f(p,a,pgc,pmom,pflow,pfsf,pfsfdisc,psolv,preini,ppart);

	
// Start MAINLOOP
	if(p->A10==4)
    loop_nsewave(a);
    
    if(p->A10==5 && p->X10 == 1 && p->X13 >= 1 && p->N40==0) 
	{
		loop_cfd_fsi(a);
	}
    else if(p->A10==5)
	{
		loop_cfd(a);
	}
}
