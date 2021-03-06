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
#include "Communication/communication.h"

int main(int argc, char *argv[]){
	std::cout << "#### SuperFlow3000 ####\n";

	char InputFileName[] = "inputvals.bin";
	char OutputFolderName[] = "output";  // output folder! -> be careful, if folder is not there, no data are saved
	// load simparam

	Simparam simparam;
	IO Reader(InputFileName,OutputFolderName);
	simparam = Reader.getSimparam();

	//std::cout << simparam.iMax << std::endl;
	MPI_Init(&argc, &argv);
	int mpiRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	int mpiSize;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

	// willkuerliche gewaehlt - vertikal immer nur zwei gebiete

	IndexType mpiSizeV = 2;
	IndexType mpiSizeH = 2;

	bool plot = false;

	//IndexType mpiSizeH = 1;
	//IndexType mpiSizeV = 1;




	Computation pc (simparam);
	Solver sol (simparam);


	//Grids sollen gleiche groesse haben!!
	IndexType il=2;

	IndexType ir=il+(simparam.iMax)/mpiSizeH-1;
	//std::cout << simparam.xLength << std::endl;
	IndexType jb=2;
	IndexType jt=jb+(simparam.jMax)/mpiSizeV-1;


	MultiIndexType griddimension ((ir-il+4),(jt-jb+4));

	GridFunction p(griddimension,simparam.PI,'p');
    p.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'p');

    GridFunction rhs(griddimension,'r');
    rhs.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'r');

    GridFunction v(griddimension,simparam.VI,'v');
    v.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'v');

    GridFunction g(griddimension,'g');
    g.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'g');

    GridFunction u(griddimension,simparam.UI,'u');
    u.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'u');

    GridFunction f(griddimension,'f');
    f.InitializeGlobalBoundaryPosition(mpiRank,mpiSizeH,mpiSizeV,'f');




    // Durch ir, il, jb und jt abgedeckt
    /*
    MultiIndexType bb(1,1); //lower left
    MultiIndexType ee(simparam.iMax,simparam.jMax); //upper right
    */

	const PointType h(simparam.xLength/simparam.iMax , simparam.yLength/simparam.jMax);


	RealType deltaT = simparam.deltaT;
	RealType local_deltaT = simparam.deltaT;
	RealType t = 0;
	int step = 0;


	int ibegin = p.beginwrite[0];
	int iend   = p.endwrite[0];
	int jbegin = p.beginwrite[1];
	int jend   = p.endwrite[1];
	int localgriddimensionX = iend-ibegin+1;
	int localgriddimensionY = jend-jbegin+1;

	std::cout << "Schritt 1" << std::endl;

	// so gross wie u
	GridFunction gx(griddimension,simparam.GX,'u');

	// so gross wie v
	GridFunction gy(griddimension,simparam.GY,'v');

	//---- for Boundary condition ----
	//for driven cavity


	// steht jetzt in der gridfunction selber
	//MultiIndexType upperleft (1,                 griddimension[1]-1);
	//MultiIndexType upperright(griddimension[0]-2,griddimension[1]-1);
	MultiIndexType offset (0,-1);
	//evtl. zum testen noetig (einfach durchfliessen)
	MultiIndexType linksunten (0,1);
	//MultiIndexType linksoben  (0,griddimension[1]-2);
	//MultiIndexType rechtsunten (griddimension[0]-2,1);
	MultiIndexType rechtsoben  (griddimension[0]-2,griddimension[1]-2);
	// write first output data

	//wird hier ein vtk file geschrieben, ohne dass randwerte in matritzen geschrieben wurden?

	//if(mpiRank == 2){
			//Reader.writeVTKFile(griddimension,u.GetGridFunction(),v.GetGridFunction(), p.GetGridFunction(), h, step, mpiRank);
		//}
	// parallele visualisierung

	if (mpiRank == 0) {
			 Reader.writeVTKMasterfile(mpiSizeH, mpiSizeV, step, localgriddimensionX, localgriddimensionY);
	}

	Reader.writeVTKSlavefile(u, v,  p, h, mpiSizeH, mpiSizeV, step,mpiRank);

	// start time loop
	Communication communicator(mpiRank, mpiSizeH, mpiSizeV, p.globalboundary);
	while (t <= simparam.tEnd){
		// compute deltaT
		//deltaT = 0.005;
		local_deltaT = pc.computeTimestep(u.MaxValueGridFunction(u.beginwrite,u.endwrite),v.MaxValueGridFunction(v.beginwrite,v.endwrite),h);
		MPI_Allreduce ( &local_deltaT, &deltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		// set boundary
		pc.setBoundaryU(u); //First implementation: only no-flow boundaries-> everything is zero!
		pc.setBoundaryV(v);

		// driven cavity:
		if (u.globalboundary[2]){
			//MultiIndexType UpperLeft(1,u.griddimension[1]-1);
			MultiIndexType UpperLeft(0,u.griddimension[1]-1);
			MultiIndexType UpperRight(u.griddimension[0]-2, u.griddimension[1]-1);
			u.SetGridFunction(UpperLeft,UpperRight,-1.0,offset,2.0);
		}

		//Testrandwerte
		/*if (mpiRank==2){
			MultiIndexType LowerLeft(u.beginwrite[0],u.beginwrite[1]);
			MultiIndexType UpperLeft(u.beginwrite[0], u.endwrite[1]);
			u.SetGridFunction(LowerLeft,UpperLeft,2.0);
		}*/
		//einfach durchfliesen
		//u.SetGridFunction(linksunten,linkusoben,1);
		//u.SetGridFunction(rechtsunten,rechtsoben,1);

		//u.PlotGrid();
		/*std::cout << "Rank: " << mpiRank << " 0: ";
		if (u.globalboundary[0]) std::cout << 1;
		else std::cout << 0;
		std::cout << " 1: ";
		if (u.globalboundary[1]) std::cout << 1;
		else std::cout << 0;
		std::cout << " 2: ";
		if (u.globalboundary[2]) std::cout << 1;
		else std::cout << 0;
		std::cout << " 3: ";
		if (u.globalboundary[3]) std::cout <<1;
		else std::cout << 0;
		std::cout << std::endl;
*/

		if((step%10) == 0){

			Reader.writeVTKFile(griddimension,u.GetGridFunction(),v.GetGridFunction(), p.GetGridFunction(), h, step, mpiRank);
		}

		// parallel visualisierung

		if (0 == (step % 10)) {
			//Reader.writeVTKFile(griddimension,u.GetGridFunction(),v.GetGridFunction(), p.GetGridFunction(), h, step);
			 if (mpiRank == 0) {
				 Reader.writeVTKMasterfile(mpiSizeH, mpiSizeV, step, localgriddimensionX, localgriddimensionY);
			 }
			 Reader.writeVTKSlavefile(u, v,  p, h, mpiSizeH, mpiSizeV, step,mpiRank);
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

		sol.SORCycle(&p,blrhs,h,&communicator);

		if (plot == true) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpiRank == 0) {
			std::cout << "Rank: " << mpiRank<< std::endl;
			u.PlotGrid();
			std::cout << "--------------------------" << std::endl;

		}

		MPI_Barrier(MPI_COMM_WORLD);
			if (mpiRank == 1) {
				std::cout << "Rank: " << mpiRank<< std::endl;
				u.PlotGrid();
				std::cout << "--------------------------" << std::endl;

			}
		MPI_Barrier(MPI_COMM_WORLD);
				if (mpiRank == 2) {
					std::cout << "Rank: " << mpiRank<< std::endl;
					u.PlotGrid();
					std::cout << "--------------------------" << std::endl;
				}
		MPI_Barrier(MPI_COMM_WORLD);
				if (mpiRank == 3) {
					std::cout << "Rank: " << mpiRank<< std::endl;
					u.PlotGrid();
					std::cout << "--------------------------" << std::endl;
				}
		}
		//Update velocity
		GridFunctionType blp = p.GetGridFunction();

		pc.computeNewVelocities(&u, &v,blf,blg,blp,h,deltaT);

		if (plot == true) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpiRank == 0) {
			std::cout << "Rank: " << mpiRank<< std::endl;
			u.PlotGrid();
			std::cout << "--------------------------" << std::endl;

		}

		MPI_Barrier(MPI_COMM_WORLD);
			if (mpiRank == 1) {
				std::cout << "Rank: " << mpiRank<< std::endl;
				u.PlotGrid();
				std::cout << "--------------------------" << std::endl;

			}
		MPI_Barrier(MPI_COMM_WORLD);
				if (mpiRank == 2) {
					std::cout << "Rank: " << mpiRank<< std::endl;
					u.PlotGrid();
					std::cout << "--------------------------" << std::endl;
				}
		MPI_Barrier(MPI_COMM_WORLD);
				if (mpiRank == 3) {
					std::cout << "Rank: " << mpiRank<< std::endl;
					u.PlotGrid();
					std::cout << "--------------------------" << std::endl;
				}

		}

		communicator.ExchangePValues(u);
		communicator.ExchangePValues(v);



		if (plot == true) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (mpiRank == 0) {
			std::cout << "Rank: " << mpiRank<< std::endl;
			u.PlotGrid();
			std::cout << "--------------------------" << std::endl;

		}

		MPI_Barrier(MPI_COMM_WORLD);
			if (mpiRank == 1) {
				std::cout << "Rank: " << mpiRank<< std::endl;
				u.PlotGrid();
				std::cout << "--------------------------" << std::endl;

			}
		MPI_Barrier(MPI_COMM_WORLD);
				if (mpiRank == 2) {
					std::cout << "Rank: " << mpiRank<< std::endl;
					u.PlotGrid();
					std::cout << "--------------------------" << std::endl;
				}
		MPI_Barrier(MPI_COMM_WORLD);
				if (mpiRank == 3) {
					std::cout << "Rank: " << mpiRank<< std::endl;
					u.PlotGrid();
					std::cout << "--------------------------" << std::endl;
				}
		}



		// update time
		t += deltaT;
		step++;

		// write files
		std::cout<< step<<"  -  "<<t<<" / " <<simparam.tEnd<<std::endl;
	}
	MPI_Finalize();
	return 0;
}

