/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF.h"
#include<vector>
#include<fstream>
#include<iostream>
#include <Eigen/Dense>
#include"increment.h"
#include"slice4.h"
#include"sliceint5.h"
#include"ddweno_f_nug.h"

class lexer;
class fdm2D;
class ghostcell;
class net;
class slice;

using namespace std;

#ifndef SIXDOF_SFLOW_H_
#define SIXDOF_SFLOW_H_

class sixdof_sflow : public sixdof, public increment, public ddweno_f_nug
{
public:
	
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    sixdof_sflow(lexer*, fdm2D*, ghostcell*);
	void ini(lexer*,fdm2D*,ghostcell*);
	
    virtual ~sixdof_sflow();
	virtual void start(lexer*,fdm*,ghostcell*,double,vrans*,vector<net*>&);
	void start(lexer*,fdm2D*,ghostcell*);
	virtual void initialize(lexer*,fdm*,ghostcell*,vector<net*>&);
    
    virtual void isource(lexer*,fdm*,ghostcell*);
    virtual void jsource(lexer*,fdm*,ghostcell*);
    virtual void ksource(lexer*,fdm*,ghostcell*);
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*);
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*);
    
private:
	
    void cylinder(lexer*,fdm2D*,ghostcell*);
    void box(lexer*,fdm2D*,ghostcell*);
    void ini_parameter(lexer*, fdm2D*, ghostcell*);
    void print_ini(lexer*, fdm2D*, ghostcell*);
    void print_parameter(lexer*,ghostcell*);
    void print_stl(lexer*,ghostcell*);
    
    void iniPosition_RBM(lexer*, fdm2D*, ghostcell*);
    void rotation_tri(lexer*,double,double,double,double&,double&,double&, const double&, const double&, const double&);
    void quat_matrices(const Eigen::Vector4d&);
   
    void ray_cast(lexer*, ghostcell*);
	void ray_cast_io_x(lexer*, ghostcell*,int,int);
	void ray_cast_io_ycorr(lexer*, ghostcell*,int,int);
    void ray_cast_x(lexer*, ghostcell*,int,int);
	void ray_cast_y(lexer*, ghostcell*,int,int);
    void reini(lexer*,ghostcell*,slice&);
    void disc(lexer*,ghostcell*,slice&);
    void time_preproc(lexer*);

    double Hsolidface(lexer*, int,int);
    void updateFSI(lexer*, fdm2D*, ghostcell*);
    void updatePosition(lexer*, fdm2D*, ghostcell*);
    void updateForcing_hemisphere(lexer*, fdm2D*, ghostcell*);
    void updateForcing_ship(lexer*, fdm2D*, ghostcell*);

    double phi, theta, psi;
    double Uext, Vext, Wext, Pext, Qext, Rext;
    Eigen::Matrix3d quatRotMat;
    int reiniter, tricount, n6DOF, printtime;

    slice4 press,frk1,frk2,L,dt,fb;
    
    Eigen::Vector4d e_;
    Eigen::Matrix<double, 3, 4> E_, G_;
    Eigen::Matrix3d R_, Rinv_;

    // Raycast
    sliceint5 cutl,cutr,fbio;
    double **tri_x,**tri_y,**tri_z,**tri_x0,**tri_y0,**tri_z0;
    double xs,xe,ys,ye,zs,ze;
    int entity_sum, count, rayiter;
    int *tstart,*tend;
    double epsifb;
    const double epsi; 

};

#endif