/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::solve(lexer *p,fdm* a, ghostcell *pgc)
{
	update();
	
	dphi=dtheta=dpsi=0.0;
	dxg=dyg=dzg=0.0;
	
	if(p->X11_u==1)
	dUs = (Xs - p->X26_Cu*Us)/Mfb + Vs*Rs - Ws*Qs;
	if(p->X11_v==1)
	dVs = (Ys - p->X26_Cu*Vs)/Mfb + Ws*Ps - Us*Rs;
	if(p->X11_w==1)
	dWs = (Zs - p->X26_Cu*Ws)/Mfb + Us*Qs - Vs*Ps;
	
	if(p->X11_p==1)
	dPs = (Ks - p->X25_Cp*Ps - (Iz-Iy)*Qs*Rs)/Ix;
	if(p->X11_q==1)
	dQs = (Ms - p->X25_Cq*Qs - (Ix-Iz)*Rs*Ps)/Iy;
	if(p->X11_r==1)
	dRs = (Ns - p->X25_Cr*Rs - (Iy-Ix)*Ps*Qs)/Iz;
	
	if(p->X11_u==1)
	Us = Usn + 0.5*p->dt*(3.0*dUs - dUsn);
	if(p->X11_v==1)
	Vs = Vsn + 0.5*p->dt*(3.0*dVs - dVsn);
	if(p->X11_w==1)
	Ws = Wsn + 0.5*p->dt*(3.0*dWs - dWsn);
	
	if(p->X11_p==1)
	Ps = Psn + 0.5*p->dt*(3.0*dPs - dPsn);
	if(p->X11_q==1)
	Qs = Qsn + 0.5*p->dt*(3.0*dQs - dQsn);
	if(p->X11_r==1)
	Rs = Rsn + 0.5*p->dt*(3.0*dRs - dRsn);
	
	transform_vec_SE(Us,Vs,Ws,Ue,Ve,We);
	transform_angle_SE(Ps,Qs,Rs,Pe,Qe,Re);

	if(p->X11_u==1)
	dxg = 0.5*p->dt*(3.0*Ue - Uen);
	if(p->X11_v==1)
	dyg = 0.5*p->dt*(3.0*Ve - Ven);
	if(p->X11_w==1)
	dzg =  0.5*p->dt*(3.0*We - Wen);
	
	if(p->X11_p==1)
	dphi = 0.5*p->dt*(3.0*Pe - Pen);
	if(p->X11_q==1)
	dtheta = 0.5*p->dt*(3.0*Qe - Qen);
	if(p->X11_r==1)
	dpsi = 0.5*p->dt*(3.0*Re - Ren);
}

void sixdof_gc::update()
{
	dUsnnn=dUsnn;
	dUsnn=dUsn;
	dUsn=dUs;
	
	dVsnnn=dVsnn;
	dVsnn=dVsn;
	dVsn=dVs;
	
	dWsnnn=dWsnn;
	dWsnn=dWsn;
	dWsn=dWs;
	
	dPsnnn=dPsnn;
	dPsnn=dPsn;
	dPsn=dPs;
	
	dQsnnn=dQsnn;
	dQsnn=dQsn;
	dQsn=dQs;
	
	dRsnnn=dRsnn;
	dRsnn=dRsn;
	dRsn=dRs;
	
	Usnnn=Usnn;
	Usnn=Usn;
	Usn=Us;
	
	Vsnnn=Vsnn;
	Vsnn=Vsn;
	Vsn=Vs;
	
	Wsnnn=Wsnn;
	Wsnn=Wsn;
	Wsn=Ws;
	
	Psnnn=Psnn;
	Psnn=Psn;
	Psn=Ps;
	
	Qsnnn=Qsnn;
	Qsnn=Qsn;
	Qsn=Qs;
	
	Rsnnn=Rsnn;
	Rsnn=Rsn;
	Rsn=Rs;
	
	phi_s = phi_sn;
	theta_s = theta_sn;
	psi_s = psi_sn;
	
	
	dRsnnn=dRsnn;
	dRsnn=dRsn;
	dRsn=dRs;
	
	dUennn=dUenn;
	dUenn=dUen;
	dUen=dUe;	

	dVennn=dVenn;
	dVenn=dVen;
	dVen=dVe;

	dWennn=dWenn;
	dWenn=dWen;
	dWen=dWe;	
	
	Uennn=Uenn;
	Uenn=Uen;
	Uen=Ue;
	
	Vennn=Venn;
	Venn=Ven;
	Ven=Ve;
	
	Wennn=Wenn;
	Wenn=Wen;
	Wen=We;
	
	Pennn=Penn;
	Penn=Pen;
	Pen=Pe;
	
	Qennn=Qenn;
	Qenn=Qen;
	Qen=Qe;
	
	Rennn=Renn;
	Renn=Ren;
	Ren=Re;
	
	xgn=xg;
	ygn=yg;
	zgn=zg;
}


void sixdof_gc::solve_quaternion()
{
	// Symplectic integration using 2nd-order Strömer-Verlet scheme

	// Linear momentum
/*	
	e_[7] = e_[7] + p->dt/Mfb*e_[10] + p->dt*p->dt/(2.0*Mfb)*Xe;
	e_[8] = e_[8] + p->dt/Mfb*e_[11] + p->dt*p->dt/(2.0*Mfb)*Ye;
	e_[9] = e_[9] + p->dt/Mfb*e_[12] + p->dt*p->dt/(2.0*Mfb)*Ze;

	e_[10] = e_[10] + p->dt/Mfb*Xe;
	e_[11] = e_[11] + p->dt/Mfb*Ye;
	e_[12] = e_[12] + p->dt/Mfb*Ze;


	// Angular momentum
	
    std::vector<double> e1 = e_;
	std::vector<double> L1(13); 
    double norm, eold;
	
	for (int k = 0; k < 20; k++)
	{
		L_ = get_h(e_,e1);
		
		norm = 0.0;
		
		for (int i = 4; i < 7; i++)
		{
			eold = e1[i];

			e1[i] = e_[i] + 0.5*p->dt*L_[i];
			
			norm += fabs(eold - e1[i]); 
		}	
		
		if (norm < 1e-7) break;
	}

	L_ = get_e(e_,e1);
	
	for (int k = 0; k < 20; k++)
	{
		L1 = get_e(e1,e1);
		
		norm = 0.0;
		
		for (int i = 0; i < 4; i++)
		{
			eold = e1[i];

			e1[i] = e_[i] + 0.5*p->dt*(L_[i] + L1[i]);
			
			norm += fabs(eold - e1[i]); 
		}	
		
		if (norm < 1e-7) break;
	}
        
	L_ = get_h(e1,e1);

	for (int i = 0; i < 7; i++)
	{		
		e_[i] = i < 4 ? e1[i] : (e_[i] + p->dt*L_[i]);
	}

*/	

	
    // RK4 scheme

	std::vector<double> e1 = e_;
	std::vector<double> e2 = e_;
	std::vector<double> e3 = e_; 
	std::vector<double> ek = e_; 
		
	L_ = get_R(e_);
		
	for (int i=0; i<13; i++)
	{
		e1[i] = p->dt*L_[i];
			
		ek[i] = e_[i] + 0.5*e1[i];
	}   

	L_ = get_R(ek);

	for (int i=0; i<13; i++)
	{
		e2[i] = p->dt*L_[i]; 
		
		ek[i] = e_[i] + 0.5*e2[i];
	} 

	L_ = get_R(ek); 

	for (int i=0; i<13; i++)
	{
		e3[i] = p->dt*L_[i]; 

		ek[i] = e_[i] + e3[i]; 
	} 

	L_ = get_R(ek); 
	
	for (int i=0; i<13; i++)
	{
		e_[i] += 1.0/6.0*(e1[i] + 2.0*e2[i] + 2.0*e3[i] + p->dt*L_[i]); 
	}	

	update_quaternion(); 
}


std::vector<double> sixdof_gc::get_e
(
	const std::vector<double>& e1,
	const std::vector<double>& e2
)
{
	std::vector<double> L(13,0.0); 
	std::vector<double> e(7,0.0); 

	for (int i = 0; i < 7; i++)
	{		
		e[i] = i < 4 ? e1[i] : e2[i];
	}	
	
	L[0] = 
	 -(I_[0][0]*I_[1][1]*e[3]*e[6] - I_[0][0]*I_[1][2]*e[2]*e[6] - I_[0][1]*I_[1][0]*e[3]*e[6] 
	 + I_[0][1]*I_[1][2]*e[1]*e[6] + I_[0][2]*I_[1][0]*e[2]*e[6] - I_[0][2]*I_[1][1]*e[1]*e[6] 
	 - I_[0][0]*I_[2][1]*e[3]*e[5] + I_[0][0]*I_[2][2]*e[2]*e[5] + I_[0][1]*I_[2][0]*e[3]*e[5] 
	 - I_[0][1]*I_[2][2]*e[1]*e[5] - I_[0][2]*I_[2][0]*e[2]*e[5] + I_[0][2]*I_[2][1]*e[1]*e[5] 
	 + I_[1][0]*I_[2][1]*e[3]*e[4] - I_[1][0]*I_[2][2]*e[2]*e[4] - I_[1][1]*I_[2][0]*e[3]*e[4] 
	 + I_[1][1]*I_[2][2]*e[1]*e[4] + I_[1][2]*I_[2][0]*e[2]*e[4] - I_[1][2]*I_[2][1]*e[1]*e[4])
	 /(2*(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] 
	 + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]));

	L[1] = 
	 (I_[0][0]*I_[1][1]*e[2]*e[6] - I_[0][1]*I_[1][0]*e[2]*e[6] + I_[0][1]*I_[1][2]*e[0]*e[6] 
	 - I_[0][2]*I_[1][1]*e[0]*e[6] + I_[0][0]*I_[1][2]*e[3]*e[6] - I_[0][2]*I_[1][0]*e[3]*e[6] 
	 - I_[0][0]*I_[2][1]*e[2]*e[5] + I_[0][1]*I_[2][0]*e[2]*e[5] - I_[0][1]*I_[2][2]*e[0]*e[5] 
	 + I_[0][2]*I_[2][1]*e[0]*e[5] - I_[0][0]*I_[2][2]*e[3]*e[5] + I_[0][2]*I_[2][0]*e[3]*e[5] 
	 + I_[1][0]*I_[2][1]*e[2]*e[4] - I_[1][1]*I_[2][0]*e[2]*e[4] + I_[1][1]*I_[2][2]*e[0]*e[4] 
	 - I_[1][2]*I_[2][1]*e[0]*e[4] + I_[1][0]*I_[2][2]*e[3]*e[4] - I_[1][2]*I_[2][0]*e[3]*e[4])
	 /(2*(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] 
	 + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]));

	L[2] = 
	 -(I_[0][0]*I_[1][1]*e[1]*e[6] + I_[0][0]*I_[1][2]*e[0]*e[6] - I_[0][1]*I_[1][0]*e[1]*e[6] 
	 - I_[0][2]*I_[1][0]*e[0]*e[6] - I_[0][1]*I_[1][2]*e[3]*e[6] + I_[0][2]*I_[1][1]*e[3]*e[6] 
	 - I_[0][0]*I_[2][1]*e[1]*e[5] - I_[0][0]*I_[2][2]*e[0]*e[5] + I_[0][1]*I_[2][0]*e[1]*e[5] 
	 + I_[0][2]*I_[2][0]*e[0]*e[5] + I_[0][1]*I_[2][2]*e[3]*e[5] - I_[0][2]*I_[2][1]*e[3]*e[5] 
	 + I_[1][0]*I_[2][1]*e[1]*e[4] + I_[1][0]*I_[2][2]*e[0]*e[4] - I_[1][1]*I_[2][0]*e[1]*e[4] 
	 - I_[1][2]*I_[2][0]*e[0]*e[4] - I_[1][1]*I_[2][2]*e[3]*e[4] + I_[1][2]*I_[2][1]*e[3]*e[4])
	 /(2*(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] 
	 + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]));

	L[3] = 
	 (I_[0][0]*I_[1][1]*e[0]*e[6] - I_[0][1]*I_[1][0]*e[0]*e[6] - I_[0][0]*I_[1][2]*e[1]*e[6] 
	 + I_[0][2]*I_[1][0]*e[1]*e[6] - I_[0][1]*I_[1][2]*e[2]*e[6] + I_[0][2]*I_[1][1]*e[2]*e[6] 
	 - I_[0][0]*I_[2][1]*e[0]*e[5] + I_[0][1]*I_[2][0]*e[0]*e[5] + I_[0][0]*I_[2][2]*e[1]*e[5] 
	 - I_[0][2]*I_[2][0]*e[1]*e[5] + I_[0][1]*I_[2][2]*e[2]*e[5] - I_[0][2]*I_[2][1]*e[2]*e[5] 
	 + I_[1][0]*I_[2][1]*e[0]*e[4] - I_[1][1]*I_[2][0]*e[0]*e[4] - I_[1][0]*I_[2][2]*e[1]*e[4] 
	 + I_[1][2]*I_[2][0]*e[1]*e[4] - I_[1][1]*I_[2][2]*e[2]*e[4] + I_[1][2]*I_[2][1]*e[2]*e[4])
	 /(2*(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] 
	 + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]));  
	
	return L;
}


std::vector<double> sixdof_gc::get_h
(
	const std::vector<double>& e1,
	const std::vector<double>& e2
)
{
	std::vector<double> L(13,0.0); 
	std::vector<double> e(7,0.0); 
	
	for (int i = 0; i < 7; i++)
	{		
		e[i] = i < 4 ? e1[i] : e2[i];
	}	
	
	double Ke_, Me_, Ne_;

	// According to 	Shivarama and Schwab
	Ke_ = Me*(2*e[0]*e[3] + 2*e[1]*e[2]) - Ne*(2*e[0]*e[2] - 2*e[1]*e[3]) + Ke*(e[0]*e[0] + e[1]*e[1] - e[2]*e[2] - e[3]*e[3]);
	Me_ = Ne*(2*e[0]*e[1] + 2*e[2]*e[3]) - Ke*(2*e[0]*e[3] - 2*e[1]*e[2]) + Me*(e[0]*e[0] - e[1]*e[1] + e[2]*e[2] - e[3]*e[3]);
	Ne_ = Ke*(2*e[0]*e[2] + 2*e[1]*e[3]) - Me*(2*e[0]*e[1] - 2*e[2]*e[3]) + Ne*(e[0]*e[0] - e[1]*e[1] - e[2]*e[2] + e[3]*e[3]);
		
	
	L[4] = 
	 (I_[0][0]*I_[1][2]*e[0]*e[0]*e[6]*e[6] - I_[0][2]*I_[1][0]*e[0]*e[0]*e[6]*e[6] + I_[0][0]*I_[1][2]*e[1]*e[1]*e[6]*e[6] 
	 - I_[0][2]*I_[1][0]*e[1]*e[1]*e[6]*e[6] + I_[0][0]*I_[1][2]*e[2]*e[2]*e[6]*e[6] - I_[0][2]*I_[1][0]*e[2]*e[2]*e[6]*e[6] 
	 + I_[0][0]*I_[1][2]*e[3]*e[3]*e[6]*e[6] - I_[0][2]*I_[1][0]*e[3]*e[3]*e[6]*e[6] - I_[0][0]*I_[2][1]*e[0]*e[0]*e[5]*e[5] 
	 + I_[0][1]*I_[2][0]*e[0]*e[0]*e[5]*e[5] - I_[0][0]*I_[2][1]*e[1]*e[1]*e[5]*e[5] + I_[0][1]*I_[2][0]*e[1]*e[1]*e[5]*e[5] 
	 - I_[0][0]*I_[2][1]*e[2]*e[2]*e[5]*e[5] + I_[0][1]*I_[2][0]*e[2]*e[2]*e[5]*e[5] - I_[0][0]*I_[2][1]*e[3]*e[3]*e[5]*e[5] 
	 + I_[0][1]*I_[2][0]*e[3]*e[3]*e[5]*e[5] + I_[0][0]*I_[1][1]*I_[2][2]*Ke_ - I_[0][0]*I_[1][2]*I_[2][1]*Ke_ 
	 - I_[0][1]*I_[1][0]*I_[2][2]*Ke_ + I_[0][1]*I_[1][2]*I_[2][0]*Ke_ + I_[0][2]*I_[1][0]*I_[2][1]*Ke_ 
	 - I_[0][2]*I_[1][1]*I_[2][0]*Ke_ + I_[0][0]*I_[1][1]*e[0]*e[0]*e[5]*e[6] - I_[0][1]*I_[1][0]*e[0]*e[0]*e[5]*e[6] 
	 + I_[0][0]*I_[1][1]*e[1]*e[1]*e[5]*e[6] - I_[0][1]*I_[1][0]*e[1]*e[1]*e[5]*e[6] + I_[0][0]*I_[1][1]*e[2]*e[2]*e[5]*e[6] 
	 - I_[0][1]*I_[1][0]*e[2]*e[2]*e[5]*e[6] + I_[0][0]*I_[1][1]*e[3]*e[3]*e[5]*e[6] - I_[0][1]*I_[1][0]*e[3]*e[3]*e[5]*e[6] 
	 - I_[0][0]*I_[2][2]*e[0]*e[0]*e[5]*e[6] + I_[0][2]*I_[2][0]*e[0]*e[0]*e[5]*e[6] - I_[0][0]*I_[2][2]*e[1]*e[1]*e[5]*e[6] 
	 + I_[0][2]*I_[2][0]*e[1]*e[1]*e[5]*e[6] - I_[0][0]*I_[2][2]*e[2]*e[2]*e[5]*e[6] + I_[0][2]*I_[2][0]*e[2]*e[2]*e[5]*e[6] 
	 - I_[0][0]*I_[2][2]*e[3]*e[3]*e[5]*e[6] + I_[0][2]*I_[2][0]*e[3]*e[3]*e[5]*e[6] + I_[1][0]*I_[2][1]*e[0]*e[0]*e[4]*e[5] 
	 - I_[1][1]*I_[2][0]*e[0]*e[0]*e[4]*e[5] + I_[1][0]*I_[2][1]*e[1]*e[1]*e[4]*e[5] - I_[1][1]*I_[2][0]*e[1]*e[1]*e[4]*e[5] 
	 + I_[1][0]*I_[2][1]*e[2]*e[2]*e[4]*e[5] + I_[1][0]*I_[2][2]*e[0]*e[0]*e[4]*e[6] - I_[1][1]*I_[2][0]*e[2]*e[2]*e[4]*e[5] 
	 - I_[1][2]*I_[2][0]*e[0]*e[0]*e[4]*e[6] + I_[1][0]*I_[2][1]*e[3]*e[3]*e[4]*e[5] + I_[1][0]*I_[2][2]*e[1]*e[1]*e[4]*e[6] 
	 - I_[1][1]*I_[2][0]*e[3]*e[3]*e[4]*e[5] - I_[1][2]*I_[2][0]*e[1]*e[1]*e[4]*e[6] + I_[1][0]*I_[2][2]*e[2]*e[2]*e[4]*e[6] 
	 - I_[1][2]*I_[2][0]*e[2]*e[2]*e[4]*e[6] + I_[1][0]*I_[2][2]*e[3]*e[3]*e[4]*e[6] - I_[1][2]*I_[2][0]*e[3]*e[3]*e[4]*e[6])
	 /(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] 
	 + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);
	 
	L[5] = 
	 (I_[0][1]*I_[1][2]*e[0]*e[0]*e[6]*e[6] - I_[0][2]*I_[1][1]*e[0]*e[0]*e[6]*e[6] + I_[0][1]*I_[1][2]*e[1]*e[1]*e[6]*e[6] 
	 - I_[0][2]*I_[1][1]*e[1]*e[1]*e[6]*e[6] + I_[0][1]*I_[1][2]*e[2]*e[2]*e[6]*e[6] - I_[0][2]*I_[1][1]*e[2]*e[2]*e[6]*e[6] 
	 + I_[0][1]*I_[1][2]*e[3]*e[3]*e[6]*e[6] - I_[0][2]*I_[1][1]*e[3]*e[3]*e[6]*e[6] - I_[1][0]*I_[2][1]*e[0]*e[0]*e[4]*e[4] 
	 + I_[1][1]*I_[2][0]*e[0]*e[0]*e[4]*e[4] - I_[1][0]*I_[2][1]*e[1]*e[1]*e[4]*e[4] + I_[1][1]*I_[2][0]*e[1]*e[1]*e[4]*e[4] 
	 - I_[1][0]*I_[2][1]*e[2]*e[2]*e[4]*e[4] + I_[1][1]*I_[2][0]*e[2]*e[2]*e[4]*e[4] - I_[1][0]*I_[2][1]*e[3]*e[3]*e[4]*e[4] 
	 + I_[1][1]*I_[2][0]*e[3]*e[3]*e[4]*e[4] + I_[0][0]*I_[1][1]*I_[2][2]*Me_ - I_[0][0]*I_[1][2]*I_[2][1]*Me_ 
	 - I_[0][1]*I_[1][0]*I_[2][2]*Me_ + I_[0][1]*I_[1][2]*I_[2][0]*Me_ + I_[0][2]*I_[1][0]*I_[2][1]*Me_ 
	 - I_[0][2]*I_[1][1]*I_[2][0]*Me_ - I_[0][0]*I_[1][1]*e[0]*e[0]*e[4]*e[6] + I_[0][1]*I_[1][0]*e[0]*e[0]*e[4]*e[6] 
	 - I_[0][0]*I_[1][1]*e[1]*e[1]*e[4]*e[6] + I_[0][1]*I_[1][0]*e[1]*e[1]*e[4]*e[6] - I_[0][0]*I_[1][1]*e[2]*e[2]*e[4]*e[6] 
	 + I_[0][1]*I_[1][0]*e[2]*e[2]*e[4]*e[6] - I_[0][0]*I_[1][1]*e[3]*e[3]*e[4]*e[6] + I_[0][1]*I_[1][0]*e[3]*e[3]*e[4]*e[6] 
	 + I_[0][0]*I_[2][1]*e[0]*e[0]*e[4]*e[5] - I_[0][1]*I_[2][0]*e[0]*e[0]*e[4]*e[5] + I_[0][0]*I_[2][1]*e[1]*e[1]*e[4]*e[5] 
	 - I_[0][1]*I_[2][0]*e[1]*e[1]*e[4]*e[5] + I_[0][0]*I_[2][1]*e[2]*e[2]*e[4]*e[5] - I_[0][1]*I_[2][0]*e[2]*e[2]*e[4]*e[5] 
	 + I_[0][0]*I_[2][1]*e[3]*e[3]*e[4]*e[5] - I_[0][1]*I_[2][0]*e[3]*e[3]*e[4]*e[5] - I_[0][1]*I_[2][2]*e[0]*e[0]*e[5]*e[6] 
	 + I_[0][2]*I_[2][1]*e[0]*e[0]*e[5]*e[6] - I_[0][1]*I_[2][2]*e[1]*e[1]*e[5]*e[6] + I_[0][2]*I_[2][1]*e[1]*e[1]*e[5]*e[6] 
	 - I_[0][1]*I_[2][2]*e[2]*e[2]*e[5]*e[6] + I_[0][2]*I_[2][1]*e[2]*e[2]*e[5]*e[6] - I_[0][1]*I_[2][2]*e[3]*e[3]*e[5]*e[6] 
	 + I_[0][2]*I_[2][1]*e[3]*e[3]*e[5]*e[6] + I_[1][1]*I_[2][2]*e[0]*e[0]*e[4]*e[6] - I_[1][2]*I_[2][1]*e[0]*e[0]*e[4]*e[6] 
	 + I_[1][1]*I_[2][2]*e[1]*e[1]*e[4]*e[6] - I_[1][2]*I_[2][1]*e[1]*e[1]*e[4]*e[6] + I_[1][1]*I_[2][2]*e[2]*e[2]*e[4]*e[6] 
	 - I_[1][2]*I_[2][1]*e[2]*e[2]*e[4]*e[6] + I_[1][1]*I_[2][2]*e[3]*e[3]*e[4]*e[6] - I_[1][2]*I_[2][1]*e[3]*e[3]*e[4]*e[6])
	 /(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] 
	 + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);
	 
	 L[6] = 
	 (I_[0][1]*I_[2][2]*e[0]*e[0]*e[5]*e[5] - I_[0][2]*I_[2][1]*e[0]*e[0]*e[5]*e[5] + I_[0][1]*I_[2][2]*e[1]*e[1]*e[5]*e[5] 
	 - I_[0][2]*I_[2][1]*e[1]*e[1]*e[5]*e[5] + I_[0][1]*I_[2][2]*e[2]*e[2]*e[5]*e[5] - I_[0][2]*I_[2][1]*e[2]*e[2]*e[5]*e[5] 
	 + I_[0][1]*I_[2][2]*e[3]*e[3]*e[5]*e[5] - I_[0][2]*I_[2][1]*e[3]*e[3]*e[5]*e[5] - I_[1][0]*I_[2][2]*e[0]*e[0]*e[4]*e[4] 
	 + I_[1][2]*I_[2][0]*e[0]*e[0]*e[4]*e[4] - I_[1][0]*I_[2][2]*e[1]*e[1]*e[4]*e[4] + I_[1][2]*I_[2][0]*e[1]*e[1]*e[4]*e[4] 
	 - I_[1][0]*I_[2][2]*e[2]*e[2]*e[4]*e[4] + I_[1][2]*I_[2][0]*e[2]*e[2]*e[4]*e[4] - I_[1][0]*I_[2][2]*e[3]*e[3]*e[4]*e[4] 
	 + I_[1][2]*I_[2][0]*e[3]*e[3]*e[4]*e[4] + I_[0][0]*I_[1][1]*I_[2][2]*Ne_ - I_[0][0]*I_[1][2]*I_[2][1]*Ne_ 
	 - I_[0][1]*I_[1][0]*I_[2][2]*Ne_ + I_[0][1]*I_[1][2]*I_[2][0]*Ne_ + I_[0][2]*I_[1][0]*I_[2][1]*Ne_
	 - I_[0][2]*I_[1][1]*I_[2][0]*Ne_ - I_[0][0]*I_[1][2]*e[0]*e[0]*e[4]*e[6] + I_[0][2]*I_[1][0]*e[0]*e[0]*e[4]*e[6] 
	 - I_[0][0]*I_[1][2]*e[1]*e[1]*e[4]*e[6] + I_[0][2]*I_[1][0]*e[1]*e[1]*e[4]*e[6] - I_[0][0]*I_[1][2]*e[2]*e[2]*e[4]*e[6] 
	 - I_[0][1]*I_[1][2]*e[0]*e[0]*e[5]*e[6] + I_[0][2]*I_[1][0]*e[2]*e[2]*e[4]*e[6] + I_[0][2]*I_[1][1]*e[0]*e[0]*e[5]*e[6] 
	 - I_[0][0]*I_[1][2]*e[3]*e[3]*e[4]*e[6] - I_[0][1]*I_[1][2]*e[1]*e[1]*e[5]*e[6] + I_[0][2]*I_[1][0]*e[3]*e[3]*e[4]*e[6] 
	 + I_[0][2]*I_[1][1]*e[1]*e[1]*e[5]*e[6] - I_[0][1]*I_[1][2]*e[2]*e[2]*e[5]*e[6] + I_[0][2]*I_[1][1]*e[2]*e[2]*e[5]*e[6] 
	 - I_[0][1]*I_[1][2]*e[3]*e[3]*e[5]*e[6] + I_[0][2]*I_[1][1]*e[3]*e[3]*e[5]*e[6] + I_[0][0]*I_[2][2]*e[0]*e[0]*e[4]*e[5] 
	 - I_[0][2]*I_[2][0]*e[0]*e[0]*e[4]*e[5] + I_[0][0]*I_[2][2]*e[1]*e[1]*e[4]*e[5] - I_[0][2]*I_[2][0]*e[1]*e[1]*e[4]*e[5] 
	 + I_[0][0]*I_[2][2]*e[2]*e[2]*e[4]*e[5] - I_[0][2]*I_[2][0]*e[2]*e[2]*e[4]*e[5] + I_[0][0]*I_[2][2]*e[3]*e[3]*e[4]*e[5] 
	 - I_[0][2]*I_[2][0]*e[3]*e[3]*e[4]*e[5] - I_[1][1]*I_[2][2]*e[0]*e[0]*e[4]*e[5] + I_[1][2]*I_[2][1]*e[0]*e[0]*e[4]*e[5] 
	 - I_[1][1]*I_[2][2]*e[1]*e[1]*e[4]*e[5] + I_[1][2]*I_[2][1]*e[1]*e[1]*e[4]*e[5] - I_[1][1]*I_[2][2]*e[2]*e[2]*e[4]*e[5] 
	 + I_[1][2]*I_[2][1]*e[2]*e[2]*e[4]*e[5] - I_[1][1]*I_[2][2]*e[3]*e[3]*e[4]*e[5] + I_[1][2]*I_[2][1]*e[3]*e[3]*e[4]*e[5])
	 /(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] 
	 + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);	
	
	return L;	
}

std::vector<double> sixdof_gc::get_R
(
	const std::vector<double>& e
)
{
	std::vector<double> L(13);    
	
	double Ke_, Me_, Ne_;

	// According to 	Shivarama and Schwab
	Ke_ = Me*(2*e[0]*e[3] + 2*e[1]*e[2]) - Ne*(2*e[0]*e[2] - 2*e[1]*e[3]) + Ke*(e[0]*e[0] + e[1]*e[1] - e[2]*e[2] - e[3]*e[3]);
	Me_ = Ne*(2*e[0]*e[1] + 2*e[2]*e[3]) - Ke*(2*e[0]*e[3] - 2*e[1]*e[2]) + Me*(e[0]*e[0] - e[1]*e[1] + e[2]*e[2] - e[3]*e[3]);
	Ne_ = Ke*(2*e[0]*e[2] + 2*e[1]*e[3]) - Me*(2*e[0]*e[1] - 2*e[2]*e[3]) + Ne*(e[0]*e[0] - e[1]*e[1] - e[2]*e[2] + e[3]*e[3]);

	L[0] = 
	 -(I_[0][0]*I_[1][1]*e[3]*e[6] - I_[0][0]*I_[1][2]*e[2]*e[6] - I_[0][1]*I_[1][0]*e[3]*e[6] 
	 + I_[0][1]*I_[1][2]*e[1]*e[6] + I_[0][2]*I_[1][0]*e[2]*e[6] - I_[0][2]*I_[1][1]*e[1]*e[6] 
	 - I_[0][0]*I_[2][1]*e[3]*e[5] + I_[0][0]*I_[2][2]*e[2]*e[5] + I_[0][1]*I_[2][0]*e[3]*e[5] 
	 - I_[0][1]*I_[2][2]*e[1]*e[5] - I_[0][2]*I_[2][0]*e[2]*e[5] + I_[0][2]*I_[2][1]*e[1]*e[5] 
	 + I_[1][0]*I_[2][1]*e[3]*e[4] - I_[1][0]*I_[2][2]*e[2]*e[4] - I_[1][1]*I_[2][0]*e[3]*e[4] 
	 + I_[1][1]*I_[2][2]*e[1]*e[4] + I_[1][2]*I_[2][0]*e[2]*e[4] - I_[1][2]*I_[2][1]*e[1]*e[4])
	 /(2*(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] 
	 + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]));

	L[1] = 
	 (I_[0][0]*I_[1][1]*e[2]*e[6] - I_[0][1]*I_[1][0]*e[2]*e[6] + I_[0][1]*I_[1][2]*e[0]*e[6] 
	 - I_[0][2]*I_[1][1]*e[0]*e[6] + I_[0][0]*I_[1][2]*e[3]*e[6] - I_[0][2]*I_[1][0]*e[3]*e[6] 
	 - I_[0][0]*I_[2][1]*e[2]*e[5] + I_[0][1]*I_[2][0]*e[2]*e[5] - I_[0][1]*I_[2][2]*e[0]*e[5] 
	 + I_[0][2]*I_[2][1]*e[0]*e[5] - I_[0][0]*I_[2][2]*e[3]*e[5] + I_[0][2]*I_[2][0]*e[3]*e[5] 
	 + I_[1][0]*I_[2][1]*e[2]*e[4] - I_[1][1]*I_[2][0]*e[2]*e[4] + I_[1][1]*I_[2][2]*e[0]*e[4] 
	 - I_[1][2]*I_[2][1]*e[0]*e[4] + I_[1][0]*I_[2][2]*e[3]*e[4] - I_[1][2]*I_[2][0]*e[3]*e[4])
	 /(2*(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] 
	 + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]));

	L[2] = 
	 -(I_[0][0]*I_[1][1]*e[1]*e[6] + I_[0][0]*I_[1][2]*e[0]*e[6] - I_[0][1]*I_[1][0]*e[1]*e[6] 
	 - I_[0][2]*I_[1][0]*e[0]*e[6] - I_[0][1]*I_[1][2]*e[3]*e[6] + I_[0][2]*I_[1][1]*e[3]*e[6] 
	 - I_[0][0]*I_[2][1]*e[1]*e[5] - I_[0][0]*I_[2][2]*e[0]*e[5] + I_[0][1]*I_[2][0]*e[1]*e[5] 
	 + I_[0][2]*I_[2][0]*e[0]*e[5] + I_[0][1]*I_[2][2]*e[3]*e[5] - I_[0][2]*I_[2][1]*e[3]*e[5] 
	 + I_[1][0]*I_[2][1]*e[1]*e[4] + I_[1][0]*I_[2][2]*e[0]*e[4] - I_[1][1]*I_[2][0]*e[1]*e[4] 
	 - I_[1][2]*I_[2][0]*e[0]*e[4] - I_[1][1]*I_[2][2]*e[3]*e[4] + I_[1][2]*I_[2][1]*e[3]*e[4])
	 /(2*(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] 
	 + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]));

	L[3] = 
	 (I_[0][0]*I_[1][1]*e[0]*e[6] - I_[0][1]*I_[1][0]*e[0]*e[6] - I_[0][0]*I_[1][2]*e[1]*e[6] 
	 + I_[0][2]*I_[1][0]*e[1]*e[6] - I_[0][1]*I_[1][2]*e[2]*e[6] + I_[0][2]*I_[1][1]*e[2]*e[6] 
	 - I_[0][0]*I_[2][1]*e[0]*e[5] + I_[0][1]*I_[2][0]*e[0]*e[5] + I_[0][0]*I_[2][2]*e[1]*e[5] 
	 - I_[0][2]*I_[2][0]*e[1]*e[5] + I_[0][1]*I_[2][2]*e[2]*e[5] - I_[0][2]*I_[2][1]*e[2]*e[5] 
	 + I_[1][0]*I_[2][1]*e[0]*e[4] - I_[1][1]*I_[2][0]*e[0]*e[4] - I_[1][0]*I_[2][2]*e[1]*e[4] 
	 + I_[1][2]*I_[2][0]*e[1]*e[4] - I_[1][1]*I_[2][2]*e[2]*e[4] + I_[1][2]*I_[2][1]*e[2]*e[4])
	 /(2*(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] 
	 + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]));


	L[4] = 
	 (I_[0][0]*I_[1][2]*e[0]*e[0]*e[6]*e[6] - I_[0][2]*I_[1][0]*e[0]*e[0]*e[6]*e[6] + I_[0][0]*I_[1][2]*e[1]*e[1]*e[6]*e[6] 
	 - I_[0][2]*I_[1][0]*e[1]*e[1]*e[6]*e[6] + I_[0][0]*I_[1][2]*e[2]*e[2]*e[6]*e[6] - I_[0][2]*I_[1][0]*e[2]*e[2]*e[6]*e[6] 
	 + I_[0][0]*I_[1][2]*e[3]*e[3]*e[6]*e[6] - I_[0][2]*I_[1][0]*e[3]*e[3]*e[6]*e[6] - I_[0][0]*I_[2][1]*e[0]*e[0]*e[5]*e[5] 
	 + I_[0][1]*I_[2][0]*e[0]*e[0]*e[5]*e[5] - I_[0][0]*I_[2][1]*e[1]*e[1]*e[5]*e[5] + I_[0][1]*I_[2][0]*e[1]*e[1]*e[5]*e[5] 
	 - I_[0][0]*I_[2][1]*e[2]*e[2]*e[5]*e[5] + I_[0][1]*I_[2][0]*e[2]*e[2]*e[5]*e[5] - I_[0][0]*I_[2][1]*e[3]*e[3]*e[5]*e[5] 
	 + I_[0][1]*I_[2][0]*e[3]*e[3]*e[5]*e[5] + I_[0][0]*I_[1][1]*I_[2][2]*Ke_ - I_[0][0]*I_[1][2]*I_[2][1]*Ke_ 
	 - I_[0][1]*I_[1][0]*I_[2][2]*Ke_ + I_[0][1]*I_[1][2]*I_[2][0]*Ke_ + I_[0][2]*I_[1][0]*I_[2][1]*Ke_ 
	 - I_[0][2]*I_[1][1]*I_[2][0]*Ke_ + I_[0][0]*I_[1][1]*e[0]*e[0]*e[5]*e[6] - I_[0][1]*I_[1][0]*e[0]*e[0]*e[5]*e[6] 
	 + I_[0][0]*I_[1][1]*e[1]*e[1]*e[5]*e[6] - I_[0][1]*I_[1][0]*e[1]*e[1]*e[5]*e[6] + I_[0][0]*I_[1][1]*e[2]*e[2]*e[5]*e[6] 
	 - I_[0][1]*I_[1][0]*e[2]*e[2]*e[5]*e[6] + I_[0][0]*I_[1][1]*e[3]*e[3]*e[5]*e[6] - I_[0][1]*I_[1][0]*e[3]*e[3]*e[5]*e[6] 
	 - I_[0][0]*I_[2][2]*e[0]*e[0]*e[5]*e[6] + I_[0][2]*I_[2][0]*e[0]*e[0]*e[5]*e[6] - I_[0][0]*I_[2][2]*e[1]*e[1]*e[5]*e[6] 
	 + I_[0][2]*I_[2][0]*e[1]*e[1]*e[5]*e[6] - I_[0][0]*I_[2][2]*e[2]*e[2]*e[5]*e[6] + I_[0][2]*I_[2][0]*e[2]*e[2]*e[5]*e[6] 
	 - I_[0][0]*I_[2][2]*e[3]*e[3]*e[5]*e[6] + I_[0][2]*I_[2][0]*e[3]*e[3]*e[5]*e[6] + I_[1][0]*I_[2][1]*e[0]*e[0]*e[4]*e[5] 
	 - I_[1][1]*I_[2][0]*e[0]*e[0]*e[4]*e[5] + I_[1][0]*I_[2][1]*e[1]*e[1]*e[4]*e[5] - I_[1][1]*I_[2][0]*e[1]*e[1]*e[4]*e[5] 
	 + I_[1][0]*I_[2][1]*e[2]*e[2]*e[4]*e[5] + I_[1][0]*I_[2][2]*e[0]*e[0]*e[4]*e[6] - I_[1][1]*I_[2][0]*e[2]*e[2]*e[4]*e[5] 
	 - I_[1][2]*I_[2][0]*e[0]*e[0]*e[4]*e[6] + I_[1][0]*I_[2][1]*e[3]*e[3]*e[4]*e[5] + I_[1][0]*I_[2][2]*e[1]*e[1]*e[4]*e[6] 
	 - I_[1][1]*I_[2][0]*e[3]*e[3]*e[4]*e[5] - I_[1][2]*I_[2][0]*e[1]*e[1]*e[4]*e[6] + I_[1][0]*I_[2][2]*e[2]*e[2]*e[4]*e[6] 
	 - I_[1][2]*I_[2][0]*e[2]*e[2]*e[4]*e[6] + I_[1][0]*I_[2][2]*e[3]*e[3]*e[4]*e[6] - I_[1][2]*I_[2][0]*e[3]*e[3]*e[4]*e[6])
	 /(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] 
	 + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);
	 
	L[5] = 
	 (I_[0][1]*I_[1][2]*e[0]*e[0]*e[6]*e[6] - I_[0][2]*I_[1][1]*e[0]*e[0]*e[6]*e[6] + I_[0][1]*I_[1][2]*e[1]*e[1]*e[6]*e[6] 
	 - I_[0][2]*I_[1][1]*e[1]*e[1]*e[6]*e[6] + I_[0][1]*I_[1][2]*e[2]*e[2]*e[6]*e[6] - I_[0][2]*I_[1][1]*e[2]*e[2]*e[6]*e[6] 
	 + I_[0][1]*I_[1][2]*e[3]*e[3]*e[6]*e[6] - I_[0][2]*I_[1][1]*e[3]*e[3]*e[6]*e[6] - I_[1][0]*I_[2][1]*e[0]*e[0]*e[4]*e[4] 
	 + I_[1][1]*I_[2][0]*e[0]*e[0]*e[4]*e[4] - I_[1][0]*I_[2][1]*e[1]*e[1]*e[4]*e[4] + I_[1][1]*I_[2][0]*e[1]*e[1]*e[4]*e[4] 
	 - I_[1][0]*I_[2][1]*e[2]*e[2]*e[4]*e[4] + I_[1][1]*I_[2][0]*e[2]*e[2]*e[4]*e[4] - I_[1][0]*I_[2][1]*e[3]*e[3]*e[4]*e[4] 
	 + I_[1][1]*I_[2][0]*e[3]*e[3]*e[4]*e[4] + I_[0][0]*I_[1][1]*I_[2][2]*Me_ - I_[0][0]*I_[1][2]*I_[2][1]*Me_ 
	 - I_[0][1]*I_[1][0]*I_[2][2]*Me_ + I_[0][1]*I_[1][2]*I_[2][0]*Me_ + I_[0][2]*I_[1][0]*I_[2][1]*Me_ 
	 - I_[0][2]*I_[1][1]*I_[2][0]*Me_ - I_[0][0]*I_[1][1]*e[0]*e[0]*e[4]*e[6] + I_[0][1]*I_[1][0]*e[0]*e[0]*e[4]*e[6] 
	 - I_[0][0]*I_[1][1]*e[1]*e[1]*e[4]*e[6] + I_[0][1]*I_[1][0]*e[1]*e[1]*e[4]*e[6] - I_[0][0]*I_[1][1]*e[2]*e[2]*e[4]*e[6] 
	 + I_[0][1]*I_[1][0]*e[2]*e[2]*e[4]*e[6] - I_[0][0]*I_[1][1]*e[3]*e[3]*e[4]*e[6] + I_[0][1]*I_[1][0]*e[3]*e[3]*e[4]*e[6] 
	 + I_[0][0]*I_[2][1]*e[0]*e[0]*e[4]*e[5] - I_[0][1]*I_[2][0]*e[0]*e[0]*e[4]*e[5] + I_[0][0]*I_[2][1]*e[1]*e[1]*e[4]*e[5] 
	 - I_[0][1]*I_[2][0]*e[1]*e[1]*e[4]*e[5] + I_[0][0]*I_[2][1]*e[2]*e[2]*e[4]*e[5] - I_[0][1]*I_[2][0]*e[2]*e[2]*e[4]*e[5] 
	 + I_[0][0]*I_[2][1]*e[3]*e[3]*e[4]*e[5] - I_[0][1]*I_[2][0]*e[3]*e[3]*e[4]*e[5] - I_[0][1]*I_[2][2]*e[0]*e[0]*e[5]*e[6] 
	 + I_[0][2]*I_[2][1]*e[0]*e[0]*e[5]*e[6] - I_[0][1]*I_[2][2]*e[1]*e[1]*e[5]*e[6] + I_[0][2]*I_[2][1]*e[1]*e[1]*e[5]*e[6] 
	 - I_[0][1]*I_[2][2]*e[2]*e[2]*e[5]*e[6] + I_[0][2]*I_[2][1]*e[2]*e[2]*e[5]*e[6] - I_[0][1]*I_[2][2]*e[3]*e[3]*e[5]*e[6] 
	 + I_[0][2]*I_[2][1]*e[3]*e[3]*e[5]*e[6] + I_[1][1]*I_[2][2]*e[0]*e[0]*e[4]*e[6] - I_[1][2]*I_[2][1]*e[0]*e[0]*e[4]*e[6] 
	 + I_[1][1]*I_[2][2]*e[1]*e[1]*e[4]*e[6] - I_[1][2]*I_[2][1]*e[1]*e[1]*e[4]*e[6] + I_[1][1]*I_[2][2]*e[2]*e[2]*e[4]*e[6] 
	 - I_[1][2]*I_[2][1]*e[2]*e[2]*e[4]*e[6] + I_[1][1]*I_[2][2]*e[3]*e[3]*e[4]*e[6] - I_[1][2]*I_[2][1]*e[3]*e[3]*e[4]*e[6])
	 /(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] 
	 + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);
	 
	 L[6] = 
	 (I_[0][1]*I_[2][2]*e[0]*e[0]*e[5]*e[5] - I_[0][2]*I_[2][1]*e[0]*e[0]*e[5]*e[5] + I_[0][1]*I_[2][2]*e[1]*e[1]*e[5]*e[5] 
	 - I_[0][2]*I_[2][1]*e[1]*e[1]*e[5]*e[5] + I_[0][1]*I_[2][2]*e[2]*e[2]*e[5]*e[5] - I_[0][2]*I_[2][1]*e[2]*e[2]*e[5]*e[5] 
	 + I_[0][1]*I_[2][2]*e[3]*e[3]*e[5]*e[5] - I_[0][2]*I_[2][1]*e[3]*e[3]*e[5]*e[5] - I_[1][0]*I_[2][2]*e[0]*e[0]*e[4]*e[4] 
	 + I_[1][2]*I_[2][0]*e[0]*e[0]*e[4]*e[4] - I_[1][0]*I_[2][2]*e[1]*e[1]*e[4]*e[4] + I_[1][2]*I_[2][0]*e[1]*e[1]*e[4]*e[4] 
	 - I_[1][0]*I_[2][2]*e[2]*e[2]*e[4]*e[4] + I_[1][2]*I_[2][0]*e[2]*e[2]*e[4]*e[4] - I_[1][0]*I_[2][2]*e[3]*e[3]*e[4]*e[4] 
	 + I_[1][2]*I_[2][0]*e[3]*e[3]*e[4]*e[4] + I_[0][0]*I_[1][1]*I_[2][2]*Ne_ - I_[0][0]*I_[1][2]*I_[2][1]*Ne_ 
	 - I_[0][1]*I_[1][0]*I_[2][2]*Ne_ + I_[0][1]*I_[1][2]*I_[2][0]*Ne_ + I_[0][2]*I_[1][0]*I_[2][1]*Ne_
	 - I_[0][2]*I_[1][1]*I_[2][0]*Ne_ - I_[0][0]*I_[1][2]*e[0]*e[0]*e[4]*e[6] + I_[0][2]*I_[1][0]*e[0]*e[0]*e[4]*e[6] 
	 - I_[0][0]*I_[1][2]*e[1]*e[1]*e[4]*e[6] + I_[0][2]*I_[1][0]*e[1]*e[1]*e[4]*e[6] - I_[0][0]*I_[1][2]*e[2]*e[2]*e[4]*e[6] 
	 - I_[0][1]*I_[1][2]*e[0]*e[0]*e[5]*e[6] + I_[0][2]*I_[1][0]*e[2]*e[2]*e[4]*e[6] + I_[0][2]*I_[1][1]*e[0]*e[0]*e[5]*e[6] 
	 - I_[0][0]*I_[1][2]*e[3]*e[3]*e[4]*e[6] - I_[0][1]*I_[1][2]*e[1]*e[1]*e[5]*e[6] + I_[0][2]*I_[1][0]*e[3]*e[3]*e[4]*e[6] 
	 + I_[0][2]*I_[1][1]*e[1]*e[1]*e[5]*e[6] - I_[0][1]*I_[1][2]*e[2]*e[2]*e[5]*e[6] + I_[0][2]*I_[1][1]*e[2]*e[2]*e[5]*e[6] 
	 - I_[0][1]*I_[1][2]*e[3]*e[3]*e[5]*e[6] + I_[0][2]*I_[1][1]*e[3]*e[3]*e[5]*e[6] + I_[0][0]*I_[2][2]*e[0]*e[0]*e[4]*e[5] 
	 - I_[0][2]*I_[2][0]*e[0]*e[0]*e[4]*e[5] + I_[0][0]*I_[2][2]*e[1]*e[1]*e[4]*e[5] - I_[0][2]*I_[2][0]*e[1]*e[1]*e[4]*e[5] 
	 + I_[0][0]*I_[2][2]*e[2]*e[2]*e[4]*e[5] - I_[0][2]*I_[2][0]*e[2]*e[2]*e[4]*e[5] + I_[0][0]*I_[2][2]*e[3]*e[3]*e[4]*e[5] 
	 - I_[0][2]*I_[2][0]*e[3]*e[3]*e[4]*e[5] - I_[1][1]*I_[2][2]*e[0]*e[0]*e[4]*e[5] + I_[1][2]*I_[2][1]*e[0]*e[0]*e[4]*e[5] 
	 - I_[1][1]*I_[2][2]*e[1]*e[1]*e[4]*e[5] + I_[1][2]*I_[2][1]*e[1]*e[1]*e[4]*e[5] - I_[1][1]*I_[2][2]*e[2]*e[2]*e[4]*e[5] 
	 + I_[1][2]*I_[2][1]*e[2]*e[2]*e[4]*e[5] - I_[1][1]*I_[2][2]*e[3]*e[3]*e[4]*e[5] + I_[1][2]*I_[2][1]*e[3]*e[3]*e[4]*e[5])
	 /(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] 
	 + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);
     
     L[7] = e[10]/Mfb;
     L[8] = e[11]/Mfb;
     L[9] = e[12]/Mfb;
     L[10] = Xe;
     L[11] = Ye;
     L[12] = Ze;     
	
	return L;
}


void sixdof_gc::update_quaternion()
{
    Lnn_ = Ln_;
    Ln_ = L_;
    ennnn_ = ennn_;
    ennn_ = enn_;
    enn_ = en_;
    en_ = e_;
}
