/*
 * communication.cpp
 *
 *  Created on: 29.11.2013
 *      Author: felix
 */

#include"communication.h"
#include <iostream>



Communication::Communication(int rank, int mpiSizeH, int mpiSizeV, bool *globalboundary) {
	neighbors[0] = rank - mpiSizeH;
	neighbors[1] = rank + 1;
	neighbors[2] = rank + mpiSizeH;
	neighbors[3] = rank - 1;

	if (globalboundary[0]) {
		neighbors[0] = -1;
	}

	if (globalboundary[1]) {
		neighbors[1] = -1;
	}

	if (globalboundary[2]) {
		neighbors[2] = -1;
	}

	if (globalboundary[3]) {
		neighbors[3] = -1;
	}


	//ToDo  hier muss noch allgemeiner formuliert werden - gerade nur gueltig fuer mpiSizeV = 2 und mpiSizeH ungerade
	if (mpiSizeH%2 == 1 && rank % 2 == 0) {
		red = false;
	}
	else {
		red = true;
	}

	// Fehlermeldung, falls red aufteilung schief laeuft!
	if (mpiSizeH % 2 != 0) {std::cout << "ERROR : mpiSizeH ist nicht gerade!"; exit(0);}



}
void Communication::ExchangePValues(GridFunction& p){
	IndexType sizeV = p.endwrite[1]-p.beginwrite[1]+1;

	double sliceV[sizeV];



	if (red == false) {
		// rechter Nachbar vorhanden
		if (neighbors[1] >= 0) {
			// kopiere rechter rand in array
			for (int i=0; i < sizeV; i++) {
				sliceV[i] = p.GetGridFunction()[p.endwrite[0]][p.beginwrite[1]+i];
			}
			// sende nach rechts
			MPI_Send(&sliceV, sizeV, MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD);

			// warte auf receive von rechts (zweig red == true)
			MPI_Recv(&sliceV, sizeV, MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// schreibe sliceV in gridfunction
			for (int i = 0; i < sizeV; i++) {
				p.SetGridFunction(p.endread[0],p.beginwrite[1]+i,sliceV[i]);
			}

		}
		// linker nachbar vorhanden
		if (neighbors[3] >= 0) {
			// schreibe array mit linkem rand
			for (int i=0; i < sizeV; i++) {
				sliceV[i] = p.GetGridFunction()[p.beginwrite[0]][p.beginwrite[1]+i];
			}
			// sende nach links
			MPI_Send(&sliceV, sizeV, MPI_DOUBLE, neighbors[3], 0, MPI_COMM_WORLD);

			// warte auf receive von links (zweig red == true)
			MPI_Recv(&sliceV, sizeV, MPI_DOUBLE, neighbors[3], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// schreibe sliceV in gridfunction
			for (int i = 0; i < sizeV; i++) {
				p.SetGridFunction(p.beginread[0],p.beginwrite[1]+i,sliceV[i]);
			}

		}
	}
	else {
		// linker nachbar vorhanden
		if(neighbors[3] >= 0) {
			// warte auf receive von links (zweig red == false)
			MPI_Recv(&sliceV, sizeV, MPI_DOUBLE, neighbors[3], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// schreibe empfangene Daten in aktuelle gridfunction
			// maybe ToDo gleiche for schleife benutzen
			// schreibe array in gridfunction
			for (int i = 0; i < sizeV; i++) {
				p.SetGridFunction(p.beginread[0],p.beginwrite[1]+i,sliceV[i]);
			}
			// schreibe array fuer kopieren nach links
			for (int i=0; i < sizeV; i++) {
				sliceV[i] = p.GetGridFunction()[p.beginwrite[0]][p.beginwrite[1]+i];
			}
			// sende array nach links
			MPI_Send(&sliceV, sizeV, MPI_DOUBLE, neighbors[3], 0, MPI_COMM_WORLD);
			}
		// rechter nachbar vorhanden
		if(neighbors[1] >= 0) {
			// warte auf receive von rechts (zweig red == false)
			MPI_Recv(&sliceV, sizeV, MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// schreibe array in gridfunction
			for (int i = 0; i < sizeV; i++) {
				p.SetGridFunction(p.endread[0],p.beginwrite[1]+i,sliceV[i]);
			}

			// schreibe array fuer kopieren nach rechts
			for (int i=0; i < sizeV; i++) {
				sliceV[i] = p.GetGridFunction()[p.endwrite[0]][p.beginwrite[1]+i];
			}
			// sende array nach rechts
			MPI_Send(&sliceV, sizeV, MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD);

		}
		}
	}

