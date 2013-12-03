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
	MultiIndexType bb = u->bottomleft;
	MultiIndexType ee(u->upperright[0],u->upperright[1]);
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
                                GridFunction* u, GridFunction* v,
                                GridFunctionType& gx, GridFunctionType& gy,
                                const PointType& h, RealType& deltaT) {
	MultiIndexType dim = f->griddimension;

	Stencil sten(3,h);
	//ToDo welchen character uebergeben?? wahrscheinlich egal
	GridFunction derivative (dim,'p');
	RealType factor;
	//  --  compute F  --
	MultiIndexType bread;
	bread = f-> beginread;
	MultiIndexType eread;
	eread = f-> endread;
	MultiIndexType bwrite;
	bwrite = f-> beginwrite;
	MultiIndexType ewrite;
	ewrite = f-> endwrite;

	derivative.SetGridFunction(bread,eread,0);  //set to zero
	//add u

	GridFunctionType tmpu = u->GetGridFunction();
	f->SetGridFunction(bwrite,ewrite,1,tmpu);
	//add derivatives:
	sten.ApplyFxxStencilOperator(bread,eread,bwrite,ewrite,u->GetGridFunction(),derivative);
	factor = deltaT/param.RE;
	GridFunctionType bla = derivative.GetGridFunction(); //ToDo -> fragen wieso?
	f->AddToGridFunction(bwrite,ewrite,factor,bla);

	sten.ApplyFyyStencilOperator(bread,eread,bwrite,ewrite,u->GetGridFunction(),derivative);
	bla = derivative.GetGridFunction();
	f->AddToGridFunction(bwrite,ewrite,factor,bla);

	factor=-deltaT;
	sten.ApplyUSqxStencilOperator(bread,eread,bwrite,ewrite,(u->GetGridFunction()),derivative,param.alpha);
	bla = derivative.GetGridFunction();
	f->AddToGridFunction(bwrite,ewrite,factor,bla);

	sten.ApplyUVyStencilOperator(bread,eread,bwrite,ewrite,(u->GetGridFunction()),(v->GetGridFunction()),derivative,param.alpha);
	bla = derivative.GetGridFunction();
	f->AddToGridFunction(bwrite,ewrite,factor,bla);

	f->AddToGridFunction(bwrite,ewrite,-factor,gx);

	//  --  compute G  --
	bread = g->beginread;
	eread = g->endread;
	bwrite = g->beginwrite;
	ewrite = g->endwrite;

	derivative.SetGridFunction(bread,eread,0);  //set derivative to zero

	//add v
	GridFunctionType tmpv = v->GetGridFunction();
	g->SetGridFunction(bwrite,ewrite,1,tmpv);
	//add derivatives:
	sten.ApplyFxxStencilOperator(bread,eread,bwrite,ewrite,v->GetGridFunction(),derivative);
	factor = deltaT/param.RE;
    bla = derivative.GetGridFunction(); //ToDo -> fragen wieso?
	g->AddToGridFunction(bwrite,ewrite,factor,bla);

	sten.ApplyFyyStencilOperator(bread,eread,bwrite,ewrite,v->GetGridFunction(),derivative);
	bla = derivative.GetGridFunction();
	g->AddToGridFunction(bwrite,ewrite,factor,bla);

	factor=-deltaT;
	sten.ApplyVSqyStencilOperator(bread,eread,bwrite,ewrite,v->GetGridFunction(),derivative,param.alpha);
	bla = derivative.GetGridFunction();
	g->AddToGridFunction(bwrite,ewrite,factor,bla);

	sten.ApplyUVxStencilOperator(bread,eread,bwrite,ewrite,u->GetGridFunction(),v->GetGridFunction(),derivative,param.alpha);
	bla = derivative.GetGridFunction();
	g->AddToGridFunction(bwrite,ewrite,factor,bla);

	g->AddToGridFunction(bwrite,ewrite,-factor,gy);

}
void Computation::setBoundaryU(GridFunction& u){
    RealType value = 0;
    MultiIndexType bb;
    MultiIndexType ee;
    MultiIndexType offset;

    if(u.globalboundary[3]){
    // left -> 0
    	bb[0] = 1; bb[1] = 2;
    	ee[0] = 1; ee[1] = u.griddimension[1]-2;
    	u.SetGridFunction(bb,ee,value);
    }
    //right -> 0
    if(u.globalboundary[1]){
    	bb[0]= u.griddimension[0]-2; bb[1] = 2;
    	ee[0]= u.griddimension[0]-2; ee[1] = u.griddimension[1]-2;
    	u.SetGridFunction(bb,ee,value);
    }

    //bottom
    if(u.globalboundary[0]){
    	bb[0]= 1; bb[1] = 1;
    	ee[0]= u.griddimension[0]-2; ee[1] = 1;
    	offset[0] = 0;
    	offset[1] = 1;
    	u.SetGridFunction(bb,ee,-1,offset);
    }

    //top
    if(u.globalboundary[2]) {
    	bb[0]= 1; bb[1] = u.griddimension[1]-1;
    	ee[0]= u.griddimension[0]-2; ee[1] = u.griddimension[1]-1;
    	offset[0] = 0;
    	offset[1] = -1;
    	u.SetGridFunction(bb,ee,-1,offset);
    }
}


void Computation::setBoundaryV(GridFunction& v){
	MultiIndexType bb;
	MultiIndexType ee;
	MultiIndexType offset;
	RealType value = 0;

	// left
	if(v.globalboundary[3]) {
		bb[0] = 1; bb[1] = 1;
		ee[0] = 1; ee[1] = v.griddimension[1]-2;
		offset[0] = 1;
		offset[1] = 0;
		v.SetGridFunction(bb,ee,-1,offset);
	}

    //bottom ->0
    if(v.globalboundary[0]) {
    	bb[0] = 2; bb[1] = 1;
    	ee[0] = v.griddimension[0]-2; ee[1] = 1;
    	v.SetGridFunction(bb,ee,value);
    }
    //top ->0
    if(v.globalboundary[2]) {
    	bb[0] = 2; bb[1] = v.griddimension[1]-2;
    	ee[0] = v.griddimension[0]-2; ee[1] = v.griddimension[1]-2;
    	v.SetGridFunction(bb,ee,value);
    }

    //right
    if(v.globalboundary[1]) {
    	bb[0] = v.griddimension[0]-1; bb[1] = 1;
    	ee[0] = v.griddimension[0]-1; ee[1] = v.griddimension[1]-2;
    	offset[0]=-1;
    	offset[1] = 0;
    	v.SetGridFunction(bb,ee,-1,offset);
    }
}

void Computation::setBoundaryP(GridFunction& p){
	MultiIndexType bb;
	MultiIndexType ee;
	MultiIndexType offset;

	// left
	if(p.globalboundary[3]) {
		bb[0] = 1; bb[1] = 2;
		ee[0] = 1; ee[1] = p.griddimension[1]-2;
		offset[0] = 1;
		offset[1] = 0;
		p.SetGridFunction(bb,ee,1,offset);
	}
    //right
	if(p.globalboundary[1]) {
		bb[0] = p.griddimension[0]-1; bb[1] = 2;
		ee[0] = p.griddimension[0]-1; ee[1] = p.griddimension[1]-2;
		offset[0]=-1;
		offset[1] = 0;
		p.SetGridFunction(bb,ee,1,offset);
	}
    //bottom
	if(p.globalboundary[0]) {
		bb[0] = 2; bb[1] = 1;
		ee[0] = p.griddimension[0]-2; ee[1] = 1;
		offset[0]=0; offset[1]=1;
		p.SetGridFunction(bb,ee,1,offset);
	}
    //top
	if(p.globalboundary[2]) {
		bb[0] = 2; bb[1] = p.griddimension[1]-1;
		ee[0] = p.griddimension[0]-2; ee[1] = p.griddimension[1]-1;
		offset[0]=0; offset[1]=-1;
		p.SetGridFunction(bb,ee,1,offset);
	}
}

void Computation::setBoundaryF(GridFunction& f, GridFunctionType& u){
	MultiIndexType bb;
	MultiIndexType ee;

	// left
	if(f.globalboundary[2]) {
		bb[0] = 1; bb[1] = 2;
		ee[0] = 1; ee[1] = f.griddimension[1]-2;
		f.SetGridFunction(bb,ee,1,u);
	}
    //right
	if(f.globalboundary[1]) {
		bb[0] = f.griddimension[0]-2; bb[1] = 2;
		ee[0] = f.griddimension[0]-2; ee[1] = f.griddimension[1]-2;
		f.SetGridFunction(bb,ee,1,u);
	}
}

void Computation::setBoundaryG(GridFunction& g, GridFunctionType& v){
	MultiIndexType bb;
	MultiIndexType ee;

	//bottom
	if(g.globalboundary[0]) {
		bb[0] = 2; bb[1] = 1;
		ee[0] = g.griddimension[0]-2; ee[1] = 1;
		g.SetGridFunction(bb,ee,1,v);
	}
    //top
	if(g.globalboundary[2]) {
		bb[0] = 2; bb[1] = g.griddimension[1]-2;
		ee[0] = g.griddimension[0]-2; ee[1] = g.griddimension[1]-2;
		g.SetGridFunction(bb,ee,1,v);
	}
}

void Computation::computeRighthandSide(GridFunction* rhs,
    		GridFunctionType& f,
    		GridFunctionType& g,
    		const PointType& delta,
    		RealType deltaT){

	// alle Werte mit 0 initialisieren
	MultiIndexType bwrite (0,0);
	MultiIndexType ewrite (rhs->griddimension[0]-1,rhs->griddimension[1]-1);
	rhs->SetGridFunction(bwrite,ewrite,0);

	bwrite = rhs->beginwrite;
	ewrite = rhs->endwrite;

	RealType factor = 1/(deltaT*delta[0]);
	MultiIndexType offset (0,0);
	rhs->AddToGridFunction(bwrite,ewrite,factor,f,offset);
	offset[0]=-1;
	rhs->AddToGridFunction(bwrite,ewrite,-factor,f,offset);
	factor = (1/(deltaT*delta[1]));
	offset[0]=0;
	rhs->AddToGridFunction(bwrite,ewrite,factor,g,offset);
	offset[1]=-1;
	rhs->AddToGridFunction(bwrite,ewrite,-factor,g,offset);
}



