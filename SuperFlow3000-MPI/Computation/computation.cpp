/*
 * computation.cpp
 *
 *  Created on: 08.11.2013
 *      Author: felix
 */
#include<math.h>
#include"computation.h"
#include"../Stencil/stencil.h"
#include<iostream>
#include"../Grid/gridfunction.h"

Computation::Computation(Simparam param){
		this->param=param;
	}

RealType Computation::computeTimestep (RealType uMax, RealType vMax, const PointType& h){
    RealType minimum = param.RE/(2*(1/(h[0]*h[0])+1/(h[1]*h[1])));
    if (minimum > h[0]/abs(uMax)) {minimum = h[0]/abs(uMax);}
    if (minimum > h[1]/abs(vMax)) {minimum = h[1]/abs(vMax);}
    return param.tau*minimum;
}


void Computation::computeNewVelocities(GridFunction* u, GridFunction* v,
                                GridFunctionType& f, GridFunctionType& g,
                                GridFunctionType& p, const PointType& h,
                                RealType deltaT){
	//compute u
	MultiIndexType bb (1,1);
	MultiIndexType ee (u->griddimension[0]-3,u->griddimension[1]-2);
	u->SetGridFunction(bb,ee,1,f);
	RealType factor = -deltaT/h[0];
	MultiIndexType offset (1,0);
	u->AddToGridFunction (bb,ee, factor,p,offset);
	offset[0]=0;
	u->AddToGridFunction (bb,ee,-factor,p,offset);

	//compute v
	ee[0]= u->griddimension[0]-2; ee[1]= u->griddimension[1]-3;
	v->SetGridFunction(bb,ee,1,g);
	factor = -deltaT/h[1];
	offset[1]=1;
	v->AddToGridFunction (bb,ee, factor,p,offset);
	offset[1]=0;
	v->AddToGridFunction (bb,ee,-factor,p,offset);
}




void Computation::computeMomentumEquations(GridFunction* f, GridFunction* g,
                                GridFunctionType* u, GridFunctionType* v,
                                GridFunctionType& gx, GridFunctionType& gy,
                                const PointType& h, RealType& deltaT) {
	MultiIndexType dim = f->griddimension;

	Stencil sten(3,h);
	GridFunction derivative (dim);
	RealType factor;
	//  --  compute F  --
	MultiIndexType bread (0,0);
	MultiIndexType eread (dim[0]-2,dim[1]-1);
	MultiIndexType bwrite (1,1);
	MultiIndexType ewrite (dim[0]-3,dim[1]-2);
	derivative.SetGridFunction(bread,eread,0);  //set to zero
	//add u
	f->SetGridFunction(bwrite,ewrite,1,*u);
	//add derivatives:
	sten.ApplyFxxStencilOperator(bread,eread,bwrite,ewrite,*u,derivative);
	factor = deltaT/param.RE;
	GridFunctionType bla = derivative.GetGridFunction(); //ToDo -> fragen wieso?
	f->AddToGridFunction(bwrite,ewrite,factor,bla);

	sten.ApplyFyyStencilOperator(bread,eread,bwrite,ewrite,*u,derivative);
	bla = derivative.GetGridFunction();
	f->AddToGridFunction(bwrite,ewrite,factor,bla);

	factor=-deltaT;
	sten.ApplyUSqxStencilOperator(bread,eread,bwrite,ewrite,*u,derivative,param.alpha);
	bla = derivative.GetGridFunction();
	f->AddToGridFunction(bwrite,ewrite,factor,bla);

	sten.ApplyUVyStencilOperator(bread,eread,bwrite,ewrite,*u,*v,derivative,param.alpha);
	bla = derivative.GetGridFunction();
	f->AddToGridFunction(bwrite,ewrite,factor,bla);

	f->AddToGridFunction(bwrite,ewrite,-factor,gx);

	//  --  compute G  --
	derivative.SetGridFunction(bread,eread,0);  //set derivative to zero
	eread[0] =dim[0]-1; eread[1] =dim[1]-2;
	ewrite[0]=dim[0]-2; ewrite[1]=dim[1]-3;
	//add v
	g->SetGridFunction(bwrite,ewrite,1,*v);
	//add derivatives:
	sten.ApplyFxxStencilOperator(bread,eread,bwrite,ewrite,*v,derivative);
	factor = deltaT/param.RE;
    bla = derivative.GetGridFunction(); //ToDo -> fragen wieso?
	g->AddToGridFunction(bwrite,ewrite,factor,bla);

	sten.ApplyFyyStencilOperator(bread,eread,bwrite,ewrite,*v,derivative);
	bla = derivative.GetGridFunction();
	g->AddToGridFunction(bwrite,ewrite,factor,bla);

	factor=-deltaT;
	sten.ApplyVSqyStencilOperator(bread,eread,bwrite,ewrite,*v,derivative,param.alpha);
	bla = derivative.GetGridFunction();
	g->AddToGridFunction(bwrite,ewrite,factor,bla);

	sten.ApplyUVxStencilOperator(bread,eread,bwrite,ewrite,*u,*v,derivative,param.alpha);
	bla = derivative.GetGridFunction();
	g->AddToGridFunction(bwrite,ewrite,factor,bla);

	g->AddToGridFunction(bwrite,ewrite,-factor,gy);

}
void Computation::setBoundaryU(GridFunction& u){
    RealType value = 0;
    // left -> 0
    MultiIndexType bb(0,1);
    MultiIndexType ee(0,u.griddimension[1]-2);
    u.SetGridFunction(bb,ee,value);
    //right -> 0
    bb[0]= u.griddimension[0]-2; bb[1] = 1;
    ee[0]= u.griddimension[0]-2; ee[1] = u.griddimension[1]-2;
    u.SetGridFunction(bb,ee,value);

    //bottom
    bb[0]= 1; bb[1] = 0;
    ee[0]= u.griddimension[0]-2; ee[1] = 0;
    MultiIndexType offset(0,1);
    u.SetGridFunction(bb,ee,-1,offset);
    //ToDo: testen
    //top
    bb[0]= 1; bb[1] = u.griddimension[1]-1;
    ee[0]= u.griddimension[0]-2; ee[1] = u.griddimension[1]-1;
    offset[1]=-1;
    u.SetGridFunction(bb,ee,-1,offset);
}
void Computation::setBoundaryV(GridFunction& v){
    // left
    MultiIndexType bb (0,1);
    MultiIndexType ee (0,v.griddimension[1]-2);
    MultiIndexType offset(1,0);
    v.SetGridFunction(bb,ee,-1,offset);

    RealType value = 0;
    //bottom ->0
    bb[0] = 1; bb[1] = 0;
    ee[0] = v.griddimension[0]-2; ee[1] = 0;
    v.SetGridFunction(bb,ee,value);
    //top ->0
    bb[0] = 1; bb[1] = v.griddimension[1]-2;
    ee[0] = v.griddimension[0]-2; ee[1] = v.griddimension[1]-2;
    v.SetGridFunction(bb,ee,value);

    //right
    bb[0] = v.griddimension[0]-1; bb[1] = 1;
    ee[0] = v.griddimension[0]-1; ee[1] = v.griddimension[1]-2;
    offset[0]=-1;
    v.SetGridFunction(bb,ee,-1,offset);
}

void Computation::setBoundaryP(GridFunction& p){
	// left
    MultiIndexType bb (0,1);
    MultiIndexType ee (0,p.griddimension[1]-2);
    MultiIndexType offset(1,0);
    p.SetGridFunction(bb,ee,1,offset);
    //right
    bb[0] = p.griddimension[0]-1; bb[1] = 1;
    ee[0] = p.griddimension[0]-1; ee[1] = p.griddimension[1]-2;
    offset[0]=-1;
    p.SetGridFunction(bb,ee,1,offset);
    //bottom
    bb[0] = 1; bb[1] = 0;
    ee[0] = p.griddimension[0]-2; ee[1] = 0;
    offset[0]=0; offset[1]=1;
    p.SetGridFunction(bb,ee,1,offset);
    //top
    bb[0] = 1; bb[1] = p.griddimension[1]-1;
    ee[0] = p.griddimension[0]-2; ee[1] = p.griddimension[1]-1;
    offset[0]=0; offset[1]=-1;
    p.SetGridFunction(bb,ee,1,offset);
}
//ToDo: referenz reingeben?
void Computation::setBoundaryF(GridFunction& f, GridFunctionType& u){
	// left
	MultiIndexType bb (0,1);
	MultiIndexType ee (0,f.griddimension[1]-2);
	f.SetGridFunction(bb,ee,1,u);
    //right
    bb[0] = f.griddimension[0]-2; bb[1] = 1;
    ee[0] = f.griddimension[0]-2; ee[1] = f.griddimension[1]-2;
    f.SetGridFunction(bb,ee,1,u);
}

void Computation::setBoundaryG(GridFunction& g, GridFunctionType& v){
    //bottom
	MultiIndexType bb (1,0);
	MultiIndexType ee (g.griddimension[0]-2,0);
    g.SetGridFunction(bb,ee,1,v);
    //top
    bb[0] = 1; bb[1] = g.griddimension[1]-2;
    ee[0] = g.griddimension[0]-2; ee[1] = g.griddimension[1]-2;
    g.SetGridFunction(bb,ee,1,v);
}

void Computation::computeRighthandSide(GridFunction* rhs,
    		GridFunctionType& f,
    		GridFunctionType& g,
    		const PointType& delta,
    		RealType deltaT){

	MultiIndexType bb (0,0);
	MultiIndexType ee (rhs->griddimension[0]-1,rhs->griddimension[1]-1);
	rhs->SetGridFunction(bb,ee,0);
	bb[0]= 1; bb[1]= 1;
	ee[0]= rhs->griddimension[0]-2; ee[1]= rhs->griddimension[1]-2;
	RealType factor = 1/(deltaT*delta[0]);
	MultiIndexType offset (0,0);
	rhs->AddToGridFunction(bb,ee,factor,f,offset);
	offset[0]=-1;
	rhs->AddToGridFunction(bb,ee,-factor,f,offset);
	factor = (1/(deltaT*delta[1]));
	offset[0]=0;
	rhs->AddToGridFunction(bb,ee,factor,g,offset);
	offset[1]=-1;
	rhs->AddToGridFunction(bb,ee,-factor,g,offset);
}



