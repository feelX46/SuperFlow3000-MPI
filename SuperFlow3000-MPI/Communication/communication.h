/*
 * communication.h
 *
 *  Created on: 29.11.2013
 *      Author: felix
 */

#ifndef COMMUNICATION_H_
#define COMMUNICATION_H_

#include"../Misc/typedef.h"
#include"../Grid/gridfunction.h"
#include "mpi.h"

class Communication {
public:
	Communication(int rank, int mpiSizeH, int mpiSizeV, bool *globalboundary);

	//~Communication();

	void ExchangePValues(GridFunction& p);
	void ExchangeValues(GridFunction& p, int rank);

	IndexType neighbors[4];

	bool red;
	int mpiRank;


};




#endif /* COMMUNICATION_H_ */
