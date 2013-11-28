/*
 * solver.h
 *
 *  Created on: 08.11.2013
 *      Author: Aaron
 */

#ifndef SOLVER_H_
#define SOLVER_H_


#include"../Misc/typedef.h"
#include"../Grid/gridfunction.h"
#include "../IO/IO.hpp"
/*! @class Class Solver runs the Successive-Over-Relaxation-Cycle
 *
 */
class Solver {
public:
	Solver(Simparam param);
	/*! @brief Function to compute the global residual
		   * @param sourcegridfunction ?The discretized solution
		   * @param rhs The right hand side of the discretized local PDE
		   * @param h ?what are these two RealTypes for?
		   */
	//Solver::~Solver();


    RealType computeResidual(GridFunction& sourcegridfunction,
    						 GridFunctionType& rhs,
    						 const PointType& h);
    //-------------------------------------------------------------------------------

    /*! @brief Function to compute the global residual
    		   * @param gridfunction Pointer on the discretized solution
    		   * @param rhs The right hand side of the discretized local PDE
    		   * @param delta delta_i are the gridwidths in x- and y-direction
    		   * @param omega The relaxation-parameter of the SOR-cycle
    		   */
    void SORCycle(GridFunction* gridfunction,
    	    	  GridFunctionType& rhs,
    	    	  const PointType& h);


private:
    Simparam param;
};


#endif /* SOLVER_H_ */
