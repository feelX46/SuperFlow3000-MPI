/*
 * communication.cpp
 *
 *  Created on: 29.11.2013
 *      Author: felix
 */

#include"communication.h"
#include <iostream>



Communication::Communication(int rank, int mpiSizeH, int mpiSizeV, bool *globalboundary) {
	mpiRank = rank;
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


	// wegen fehlermeldungen bei ausfuehrung mit 6 prozessoren jetzt erstmal hart auf 4 prozessoren geschrieben
	if(rank == 0) {red = true;}
	if(rank == 1) {red = false;}
	if(rank == 2) {red = false;}
	if(rank == 3) {red = true;}


	//std::cout << "Rank: " << rank << " n0 " << neighbors[0] << " n1 " << neighbors[1] << " n2 " << neighbors[2] << " n3 " << neighbors[3] << std::endl;

	/*
	//ToDo  hier muss noch allgemeiner formuliert werden - gerade nur gueltig fuer mpiSizeV = 2 und mpiSizeH ungerade
	if (mpiSizeH%2 == 1 ) {
		if (rank % 2 == 0) {
			red = false;
		}
		else {
			red = true;
		}
		std::cout << "rank: " << rank << " color red = " << red << std::endl;
	}
	else {
		// Fehlermeldung, falls red aufteilung schief laeuft!
		std::cout << "ERROR : mpiSizeH ist gerade!";
		exit(0);
	}
*/


}
void Communication::ExchangePValues(GridFunction& p){

	// !!! LINKS RECHTS TAUSCHEN !!!
	IndexType sizeV = p.endwrite[1]-p.beginwrite[1]+1;
	//IndexType sizeV = p.endwrite[1]-p.beginwrite[1]+3;

	double sliceV[sizeV];

	if (red == false) {
		// rechter Nachbar vorhanden
		if (neighbors[1] >= 0) {
			// kopiere rechter rand in array
			for (int i=0; i < sizeV; i++) {
				sliceV[i] = p.GetGridFunction()[p.endwrite[0]][p.beginwrite[1]+i];
				//sliceV[i] = p.GetGridFunction()[p.endwrite[0]][p.beginwrite[1]-1+i];
			}
			// sende nach rechts
			MPI_Send(&sliceV, sizeV, MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD);

			// warte auf receive von rechts (zweig red == true)
			MPI_Recv(&sliceV, sizeV, MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// schreibe sliceV in gridfunction
			for (int i = 0; i < sizeV; i++) {
				p.SetGridFunction(p.endread[0],p.beginwrite[1]+i,sliceV[i]);
				//p.SetGridFunction(p.endread[0],p.beginwrite[1]-1+i,sliceV[i]);
			}

		}
		// linker nachbar vorhanden
		if (neighbors[3] >= 0) {
			// schreibe array mit linkem rand
			for (int i=0; i < sizeV; i++) {
				sliceV[i] = p.GetGridFunction()[p.beginwrite[0]][p.beginwrite[1]+i];
				//sliceV[i] = p.GetGridFunction()[p.beginwrite[0]][p.beginwrite[1]-1+i];
			}
			// sende nach links
			MPI_Send(&sliceV, sizeV, MPI_DOUBLE, neighbors[3], 0, MPI_COMM_WORLD);

			// warte auf receive von links (zweig red == true)
			MPI_Recv(&sliceV, sizeV, MPI_DOUBLE, neighbors[3], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// schreibe sliceV in gridfunction
			for (int i = 0; i < sizeV; i++) {
				p.SetGridFunction(p.beginread[0],p.beginwrite[1]+i,sliceV[i]);
				//p.SetGridFunction(p.beginread[0],p.beginwrite[1]-1+i,sliceV[i]);
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
				//p.SetGridFunction(p.beginread[0],p.beginwrite[1]-1+i,sliceV[i]);
			}
			// schreibe array fuer kopieren nach links
			for (int i=0; i < sizeV; i++) {
				sliceV[i] = p.GetGridFunction()[p.beginwrite[0]][p.beginwrite[1]+i];
				//sliceV[i] = p.GetGridFunction()[p.beginwrite[0]][p.beginwrite[1]-1+i];
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
				p.SetGridFunction(p.endread[0],p.beginwrite[1]-1+i,sliceV[i]);
			}

			// schreibe array fuer kopieren nach rechts
			for (int i=0; i < sizeV; i++) {
				sliceV[i] = p.GetGridFunction()[p.endwrite[0]][p.beginwrite[1]-1+i];
			}
			// sende array nach rechts
			MPI_Send(&sliceV, sizeV, MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD);

		}
		}


	// !! OBEN UNTEN TAUSCHEN !!
	IndexType sizeH = p.endwrite[0]-p.beginwrite[0]+1;

	double sliceH[sizeH];

	if (red == false) {
		// oberer Nachbar vorhanden
		if (neighbors[2] >= 0) {
			// kopiere oberen rand in array
			for (int i=0; i < sizeH; i++) {
				sliceH[i] = p.GetGridFunction()[p.beginwrite[0]+i][p.endwrite[1]];
			}
			// sende nach rechts
			MPI_Send(&sliceH, sizeH, MPI_DOUBLE, neighbors[2], 0, MPI_COMM_WORLD);

			// warte auf receive von oben (zweig red == true)
			MPI_Recv(&sliceH, sizeH, MPI_DOUBLE, neighbors[2], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// schreibe sliceV in gridfunction
			for (int i = 0; i < sizeH; i++) {
				p.SetGridFunction(p.beginwrite[0]+i,p.endread[1],sliceH[i]);
			}

		}

		// unterer nachbar vorhanden
		if (neighbors[0] >= 0) {
			// schreibe array mit unterem rand
			for (int i=0; i < sizeH; i++) {
				sliceH[i] = p.GetGridFunction()[p.beginwrite[0]+i][p.beginwrite[1]];
			}
			// sende nach unten
			MPI_Send(&sliceH, sizeH, MPI_DOUBLE, neighbors[0], 0, MPI_COMM_WORLD);


			// warte auf receive von unten (zweig red == true)
			MPI_Recv(&sliceH, sizeH, MPI_DOUBLE, neighbors[0], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// schreibe sliceV in gridfunction
			for (int i = 0; i < sizeH; i++) {
				p.SetGridFunction(p.beginwrite[0]+i,p.beginread[1],sliceH[i]);
			}

		}



	}
	else
	{
		// unterer nachbar vorhanden
		if(neighbors[0] >= 0) {
			// warte auf receive von unten (zweig red == false)
			MPI_Recv(&sliceH, sizeH, MPI_DOUBLE, neighbors[0], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// schreibe empfangene Daten in aktuelle gridfunction
			// maybe ToDo gleiche for schleife benutzen
			// schreibe array in gridfunction
			for (int i = 0; i < sizeH; i++) {
				p.SetGridFunction(p.beginwrite[0]+i,p.beginread[1],sliceH[i]);
			}

			// schreibe array fuer kopieren nach unten
			for (int i=0; i < sizeH; i++) {
				sliceH[i] = p.GetGridFunction()[p.beginwrite[0]+i][p.beginwrite[1]];
			}
			// sende array nach links
			MPI_Send(&sliceH, sizeH, MPI_DOUBLE, neighbors[0], 0, MPI_COMM_WORLD);
		}

		// oberer nachbar vorhanden
		if(neighbors[2] >= 0) {
			// warte auf receive von oben (zweig red == false)
			MPI_Recv(&sliceH, sizeH, MPI_DOUBLE, neighbors[2], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// schreibe array in gridfunction
			for (int i = 0; i < sizeH; i++) {
				p.SetGridFunction(p.beginwrite[0]+i,p.endread[1],sliceH[i]);
			}
			// schreibe array fuer kopieren nach oben
			for (int i=0; i < sizeH; i++) {
				sliceH[i] = p.GetGridFunction()[p.beginwrite[0]+i][p.endwrite[1]];
			}
			// sende array nach rechts
			MPI_Send(&sliceH, sizeH, MPI_DOUBLE, neighbors[2], 0, MPI_COMM_WORLD);

		}


	}

	}

