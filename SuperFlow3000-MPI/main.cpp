/*
 * main.cpp
 *
 *  Created on: Nov 5, 2013
 *      Author: ischinger
 */
#include <iostream>
#include "Misc/template.h"
#include "Misc/typedef.h"
#include "IO/IO.hpp"
#include "Stencil/stencil.h"
#include "Grid/gridfunction.h"
#include "Computation/computation.h"
#include "Solver/solver.h"
#include "mpi.h"

int main(int argc, char *argv[]){
	std::cout << "#### SuperFlow3000 ####\n";

	MPI_Init(&argc, &argv);

	int mpiRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	int mpiSize;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);




	char InputFileName[] = "inputvals.bin";
	char OutputFolderName[] = "output";  // output folder! -> be careful, if folder is not there, no data are saved
	// load simparam
	IO Reader(InputFileName,OutputFolderName);
	Simparam simparam = Reader.getSimparam();
	Computation pc (simparam);
	Solver sol (simparam);
	// initialize grids
	// different sizes for MPI

	// willkuerliche gewaehlt - vertikal immer nur zwei gebiete
	IndexType mpiSizeH = mpiSize/2;
	IndexType mpiSizeV = 2;


	//Grids sollen gleiche groesse haben!!
	IndexType il=2;
	IndexType ir=il+(simparam.iMax)/mpiSizeH-1;
	IndexType jb=2;
	IndexType jt=jb+(simparam.jMax)/mpiSizeV-1;


	MultiIndexType griddimension ((ir-il+4),(jt-jb+4));
    GridFunction p(griddimension,simparam.PI);
    p.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'p');

    GridFunction rhs(griddimension);
    rhs.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'r');

    GridFunction v(griddimension,simparam.VI);
    v.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'v');

    GridFunction g(griddimension);
    g.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'g');

    GridFunction u(griddimension,simparam.UI);
    u.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'u');

    GridFunction f(griddimension);
    f.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'f');

    // Durch ir, il, jb und jt abgedeckt
    /*
    MultiIndexType bb(1,1); //lower left
    MultiIndexType ee(simparam.iMax,simparam.jMax); //upper right
    */

	const PointType h(simparam.xLength/simparam.iMax , simparam.yLength/simparam.jMax);

	RealType deltaT = simparam.deltaT;
	RealType t = 0;
	int step = 0;

	// so gross wie u
	GridFunction gx(griddimension,simparam.GX);

	// so gross wie v
	GridFunction gy(griddimension,simparam.GY);

	//---- for Boundary condition ----
	//for driven cavity


	// steht jetzt in der gridfunction selber
	//MultiIndexType upperleft (1,                 griddimension[1]-1);
	//MultiIndexType upperright(griddimension[0]-2,griddimension[1]-1);
	MultiIndexType offset (0,-1);
	//evtl. zum testen noetig (einfach durchfliessen)
	//MultiIndexType linksunten (0,1);
	//MultiIndexType linksoben  (0,griddimension[1]-2);
	//MultiIndexType rechtsunten (griddimension[0]-2,1);
	//MultiIndexType rechtsoben  (griddimension[0]-2,griddimension[1]-2);
	// write first output data

	//wird hier ein vtk file geschrieben, ohne dass randwerte in matritzen geschrieben wurden?
	Reader.writeVTKFile(griddimension,u.GetGridFunction(),v.GetGridFunction(), p.GetGridFunction(), h, step);
	// start time loop
	while (t <= simparam.tEnd){

		// compute deltaT
		deltaT = 0.1;
		//deltaT = pc.computeTimestep(u.MaxValueGridFunction(bb,ee),v.MaxValueGridFunction(bb,ee),h);
		// set boundary
		pc.setBoundaryU(u); //First implementation: only no-flow boundaries-> everything is zero!
		pc.setBoundaryV(v);
		// driven cavity:
		if (u.globalboundary[2]){
			u.SetGridFunction(u.upperleft,u.upperright,-1.0,offset,2.0);
		}
		//einfach durchfliesen
		//u.SetGridFunction(linksunten,linksoben,1);
		//u.SetGridFunction(rechtsunten,rechtsoben,1);

		//u.PlotGrid();
		if (0 == (step % 5)) {
			Reader.writeVTKFile(griddimension,u.GetGridFunction(),v.GetGridFunction(), p.GetGridFunction(), h, step);
		}

	    // compute f / g
		GridFunctionType blgx = gx.GetGridFunction(); //ToDo: schoener machen!
		GridFunctionType blgy = gy.GetGridFunction();
		GridFunctionType blu  = u.GetGridFunction();
		GridFunctionType blv  = v.GetGridFunction();
		pc.computeMomentumEquations(&f,&g,&u,&v,blgx,blgy,h,deltaT);
		pc.setBoundaryF(f,blu);
		pc.setBoundaryG(g,blv);
		// set right side of pressure equation
		GridFunctionType blf = f.GetGridFunction();
		GridFunctionType blg = g.GetGridFunction();
		pc.computeRighthandSide(&rhs, blf, blg,h,deltaT);

		// solver
		//ToDo enventuell muss die iterationschleife hier rein!
		GridFunctionType blrhs = rhs.GetGridFunction();
		sol.SORCycle(&p,blrhs,h);

		//Update velocity
		GridFunctionType blp = p.GetGridFunction();
		pc.computeNewVelocities(&u, &v,blf,blg,blp,h,deltaT);
		// update time
		t += deltaT;
		step++;


		// write files
		std::cout<< step<<"  -  "<<t<<" / " <<simparam.tEnd<<std::endl;
	}

	GridFunction TestGridU(7,7);
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			TestGridU.SetGridFunction(i,j,i);
		}
	}

	GridFunction TestGridV(7,7);
		for (int i = 0; i < 7; i++) {
			for (int j = 0; j < 7; j++) {
				TestGridV.SetGridFunction(i,j,j);
			}
		}


	GridFunction DerivGrid(7,7);
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			DerivGrid.SetGridFunction(i,j,0);
		}
	}

	Stencil stenci(3,h);
	MultiIndexType beginread(0,0);
	MultiIndexType endread(6,6);
	MultiIndexType beginwrite(1,1);
	MultiIndexType endwrite(5,5);
	stenci.ApplyUVyStencilOperator(beginread,endread,beginwrite,endwrite, TestGridU.GetGridFunction(), TestGridV.GetGridFunction(),DerivGrid, simparam.alpha);
	std::cout << "bla " << std::endl;
	DerivGrid.PlotGrid();

	// ToDo Liste
	/*
	 * Stencil testen - kommen die richtigen Matrizen raus?
	 * Solver testen?
	 * Computation fertig machen
	 * Stencil andere Ableitungen
	 * grid wieder zum laufen bekommen
	 */

	/*// Teste Apply Grid Function
	GridFunction sourcegridfunction;

	GridFunction imagegridfunctinon;

	ApplyStencilOperator(gridreadbegin,
				gridreadend,
				gridwritebegin,
				gridwriteend,
				sourcegridfunction,
				imagegridfunction);
*/	std::cout<<"laeuft... immerhin ;)";


/*	GridFunction testgrid(griddimension,simparam.UI);
	for (int i = 0; i < griddimension[0]; i++) {
		for (int j=0; j < griddimension[1]; j++) {
			testgrid.SetGridFunction(i,j,i);
		}
	}

	testgrid.PlotGrid();

	GridFunction derivative(griddimension,0);
	Stencil sten(3,h);
	sten.ApplyPxStencilOperator(bb,ee,bb,ee,testgrid.GetGridFunction(),derivative);

	derivative.PlotGrid();
//bla
*/
	return 0;
}

