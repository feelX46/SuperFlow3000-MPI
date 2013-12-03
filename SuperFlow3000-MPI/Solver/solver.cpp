/*----------------------------------------------------------------------------------------------
 *  AUTHOR: Aaron Kraemer, and others..?
 *  Version: 1
 *
 *  File: solver.cpp
 *----------------------------------------------------------------------------------------------
 * This is the class file to class Solver
 *
 * the Solver contains the function to compute the current residual and a function
 * in which the residual is variated, in the way of SOR
 *
 */
#include "solver.h"
#include <math.h>
#include "../IO/IO.hpp"
#include "../Stencil/stencil.h"
#include "../Computation/computation.h"

Solver::Solver(Simparam param){
		this->param=param;
	}

RealType Solver::computeResidual(GridFunction& sourcegridfunction,
    				     GridFunctionType& rhs,
    					 const PointType& h){
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // The pre-value to be returned (return sqrt(doubleSum)):
    RealType doubleSum = 0.0;

    /* We need to compute the derivatives p_xx and p_yy, therefore the stencil has to be applied.
     */

    MultiIndexType dim = sourcegridfunction.GetGridDimension();
    MultiIndexType bread;
    bread = sourcegridfunction.beginread;
    MultiIndexType eread;
    eread = sourcegridfunction.endread;
    MultiIndexType bwrite;
    bwrite = sourcegridfunction.beginwrite;
    MultiIndexType ewrite;
    ewrite = sourcegridfunction.endwrite;

    //Compute the needed derivations for the whole (inner?) area
    Stencil stencil(3,h); 					// bzw. Kann man einfach const weitergeben? /Wie?
    //Get the values for derivative in x-direction:
    // ToDo welchen character uebergeben?? Wahrscheinlich egal, da interne beginread.... nicht benutzt werden
    GridFunction Fxx(dim,'p');
    stencil.ApplyFxxStencilOperator(bread, eread, bwrite, ewrite, sourcegridfunction.GetGridFunction(), Fxx);
    //Get the values for derivative in y-direction:
    GridFunction Fyy(dim,'p');
    stencil.ApplyFyyStencilOperator(bread, eread, bwrite, ewrite, sourcegridfunction.GetGridFunction(), Fyy);

    // Compute the residual: res = sqrt(Sum_i^I(Sum_j^J((p_xx+p_yy-rightHandSide)�/(I*J))))
    RealType derivator;
    for (IndexType i = bwrite[0]; i <= ewrite[0]; i++)
    {
    	for (IndexType j = bwrite[1]; j <= ewrite[1]; j++)
    	{
    		derivator = Fxx.GetGridFunction()[i][j]+ Fyy.GetGridFunction()[i][j] - rhs[i][j];
            doubleSum +=  derivator*derivator / (ewrite[0]-bwrite[0]+1) / (ewrite[1]-bwrite[1]+1);
    	}
    }
    //std::cout<<doubleSum<<std::endl;
    return sqrt(doubleSum);
}

void Solver::SORCycle(GridFunction* gridfunction,
			  GridFunctionType& rhs,
			  const PointType& h){

	MultiIndexType dim = gridfunction->GetGridDimension();
	MultiIndexType bread;
	bread = gridfunction->beginread;
	MultiIndexType eread;
	eread = gridfunction->endread;
	MultiIndexType bwrite;
	bwrite = gridfunction->beginwrite;
	MultiIndexType ewrite;
	ewrite = gridfunction-> endwrite;

	Computation pc (param);

	//Initialization of the residual. Just choose a value, that should be a bad error.
	RealType res = 10e20;
	int iterationCounter = 0;
	// SOR-cycling until error is small enough, or the number of iterations gets to high:
	RealType neighbours_x, neighbours_y;
	while (iterationCounter < param.iterMax && res > param.eps )
	{

		pc.setBoundaryP(*gridfunction);
		 for (IndexType i = bwrite[0]; i <= ewrite[0]; i++)
		{
			for (IndexType j = bwrite[1]; j <= ewrite[1]; j++)
			{
				//help-values "neighbours_x" and "neighbours_y" for better overview
				neighbours_x = (gridfunction->GetGridFunction()[i+1][j] + gridfunction->GetGridFunction()[i-1][j])
									  / h[0] / h[0];
				neighbours_y = (gridfunction->GetGridFunction()[i][j+1] + gridfunction->GetGridFunction()[i][j-1])
											  / h[1] / h[1];
				//SOR-iteration
				gridfunction->SetGridFunction(i,j,(1.0 - param.omg)*gridfunction->GetGridFunction()[i][j]
							 + param.omg /(2.0*(1/h[0]/h[0]+1/h[1]/h[1]))
							 * (neighbours_x + neighbours_y - rhs[i][j]));
			}
		}
		iterationCounter++;
		// Hier muessen Druckwerte zwischen Prozessoren ausgetauscht werden!!!


		res = computeResidual(*gridfunction, rhs, h);
	}
	if (iterationCounter >= param.iterMax)
		std::cout<<"iteration abort with error res = "<<res<<std::endl;

}
