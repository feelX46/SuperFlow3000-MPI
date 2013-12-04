/*
 * gridfunction.cpp
 *
 *  Created on: 06.11.2013
 *      Author: Markus
 */

#include "gridfunction.h"
#include <iostream>
//1
GridFunction::GridFunction(int DimX, int DimY, char indicator){
	 this->InitializeGlobalBoundary(indicator);
	 gridfunction = new RealType*[DimX];
	 for (IndexType i = 0; i < DimX; i++){
		 gridfunction[i] = new RealType [DimY-1];
	 }
	 griddimension[0] = DimX;
	 griddimension[1] = DimY;
	 for (IndexType i=0; i<DimX; i++) {
		 for(IndexType j=0; j<DimY; j++) {
			 this->SetGridFunction(i,j,0);
		 }
	 }
}

//1.1
GridFunction::GridFunction(int DimX, int DimY, RealType value, char indicator){
	 this->InitializeGlobalBoundary(indicator);
	 gridfunction = new RealType*[DimX];
	 for (IndexType i = 0; i < DimX; i++){
		 gridfunction[i] = new RealType [DimY-1];
	 }
	 const MultiIndexType begin(0,0);
	 const MultiIndexType end(DimX-1,DimY-1);
	 SetGridFunction (begin,end,value);
	 griddimension[0] = DimX;
	 griddimension[1] = DimY;
}

//2
GridFunction::GridFunction(const MultiIndexType griddimension_input, char indicator) : griddimension(griddimension_input){
	this->InitializeGlobalBoundary(indicator);
	 gridfunction= new RealType*[griddimension[0]];
	 for (IndexType i = 0; i < griddimension[0]; i++){
		 gridfunction[i] = new RealType[griddimension[1]];
	 }
}

//2.1
GridFunction::GridFunction(const MultiIndexType griddimension_input,RealType value, char indicator) : griddimension(griddimension_input){
	this->InitializeGlobalBoundary(indicator);
	gridfunction= new RealType*[griddimension[0]];
	 for (IndexType i = 0; i < griddimension[0]; i++){
		 gridfunction[i] = new RealType [griddimension[1]];
	 }
	 const MultiIndexType begin(0,0);
	 const MultiIndexType end(griddimension[0]-1,griddimension[1]-1);
	 SetGridFunction (begin,end,value);
}

//3
GridFunction::~GridFunction() {
	for (int i = 0; i < griddimension[0]; i++)
		delete[] gridfunction[i];
	delete[] gridfunction;
}


//4
GridFunctionType GridFunction::GetGridFunction() const{
	return gridfunction;
}
//new
RealType GridFunction::GetGridFunction(const MultiIndexType& index){
	return gridfunction[index[0]][index[1]];
}

//new
MultiIndexType GridFunction::GetGridDimension(){
	return griddimension;
}

//5
void GridFunction::SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType value){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0]; i <=end[0]; i++){
		for (IndexType j = begin[1]; j <=end[1]; j++){
			gridfunction[i][j] = value;
		}
	}
}

//new
void GridFunction::SetGridFunction (const IndexType& i, const IndexType& j,
		RealType value){
		gridfunction[i][j] = value;
}

//6
void GridFunction::SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,RealType factor,
			MultiIndexType& offset){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] = factor * gridfunction[i+offset[0]][j+offset[1]];
		}
	}
}

//8
void GridFunction::SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType factor, GridFunctionType& sourcegridfunction){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] = factor * sourcegridfunction[i][j];
		}
	}
}

//9
void GridFunction::SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType factor, GridFunctionType& sourcegridfunction, MultiIndexType& offset){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] = factor * sourcegridfunction[i+offset[0]][j+offset[1]];
		}
	}
}

//~10
void GridFunction::SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType factor,GridFunctionType& sourcegridfunction, RealType constant){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] = constant + factor * sourcegridfunction[i][j];
		}
	}
}

//~10
void GridFunction::SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType factor,GridFunctionType& sourcegridfunction, MultiIndexType& offset,
		RealType constant){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] = constant + factor * sourcegridfunction[i+offset[0]][j+offset[1]];
		}
	}
}

//10-1
void GridFunction::SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType factor, MultiIndexType& offset, RealType constant){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] = constant + factor * gridfunction[i+offset[0]][j+offset[1]];
		}
	}
}


//7
void GridFunction::ScaleGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType factor){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] = factor * gridfunction[i][j];
		}
	}

}

//11
void GridFunction::AddToGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType factor, GridFunctionType& sourcegridfunction){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] += factor * sourcegridfunction[i][j];
		}
	}
}

//11-1
void GridFunction::AddToGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
		RealType factor, GridFunctionType& sourcegridfunction, MultiIndexType& offset){
	if (CheckInGrid(begin,end)) {exit(0);}
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			gridfunction[i][j] += factor * sourcegridfunction[i+offset[0]][j+offset[1]];
		}
	}
}




//12
RealType GridFunction::MaxValueGridFunction (const MultiIndexType& begin,
		const MultiIndexType& end){
	if (CheckInGrid(begin,end)) {exit(0);}
	RealType max = 0;
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			if (abs(gridfunction[i][j]) > max)
				max = abs(gridfunction[i][j]);
		}
	}
	return max;
}

void GridFunction::PlotGrid(){
	for (IndexType j =griddimension[0]-1; j>=0; j--){
		for (IndexType i = 0; i<griddimension[1]; i++){
			std::cout << gridfunction[i][j] << " ";
		}
		std::cout << "\n";
		}

}

void GridFunction::InitializeGlobalBoundaryPosition(int rank, int mpiSizeH, int mpiSizeV, char indicator){
	if (rank < mpiSizeH){
		// unterer globaler Rand
		globalboundary[0] = true;

		// p und rhs  --> bleibt gleich

		// u und f  --> bleibt gleich

		// v und g
		if (indicator == 'v' || indicator == 'g') {
			beginread[0] = 1; beginread[1] = 1;
			beginwrite[0] = 2; beginwrite[1] = 2;
		}

	}
	if (rank >= (mpiSizeH*mpiSizeV)-mpiSizeH){
		// oberer globaler Rand
		globalboundary[2] = true;

		// p und rhs

		// u und f

		// v und g

		if (indicator == 'v' || indicator == 'g') {
			endread[0] = griddimension[0] - 1; endread[1] = griddimension[1] - 2;
			endwrite[0] = griddimension[0] -2; endwrite[1] = griddimension[1] - 3;
		}

	}
	if (((rank+1) % mpiSizeH) == 0) {
		// rechter globaler Rand
		globalboundary[1] = true;

		// p und rhs --> bleibt gleich

		// u und f

		if (indicator == 'u' || indicator == 'f') {
			endread[0] = griddimension[0] - 2; endread[1] = griddimension[1] - 1;
			endwrite[0] = griddimension[0] - 3; endwrite[1] = griddimension[1] - 2;

		}

		// v und g --> bleibt gleich
	}
	if ((rank%mpiSizeH) == 0) {
		// linker globaler Rand
		globalboundary[3] = true;

		// p und rhs --> bleibt gleich

		// u und f
		if (indicator == 'u' || indicator == 'f') {
			beginread[0] = 1; beginread[1] = 1;
			beginwrite[0] = 2; beginwrite[1] = 2;
		}
		// v und g
	}
}

void GridFunction::InitializeGlobalBoundary(char indicator) {
	// unten 0, dann gegen den Uhrzeigersinn
	globalboundary[0] = false;
	globalboundary[1] = false;
	globalboundary[2] = false;
	globalboundary[3] = false;

	switch (indicator)
	{
	case 'r':
	case 'p': beginread[0] = 1; beginread[1] = 1;
			   endread[0] = griddimension[0]-1; endread[1] = griddimension[1]-1;
			   beginwrite[0] = 2; beginwrite[1] = 2;
			   endwrite[0] = griddimension[0]-2; endwrite[1] = griddimension[1]-2;
			   break;
	case 'u':
	case 'f' : beginread[0] = 0; beginread[1] = 1;
	   	   	   	endread[0] = griddimension[0]-1; endread[1] = griddimension[1]-1;
	   	   	   	beginwrite[0] = 1; beginwrite[1] = 2;
	   	   	   	endwrite[0] = griddimension[0]-2; endwrite[1] = griddimension[1]-2;
	   	   	   	break;

	case 'v':
	case 'g' : beginread[0] = 1; beginread[1] = 0;
	   	   	   	endread[0] = griddimension[0]-1; endread[1] = griddimension[1]-1;
	            beginwrite[0] = 2; beginwrite[1] = 1;
	            endwrite[0] = griddimension[0]-2; endwrite[1] = griddimension[1]-2;
	            break;

	}


	bottomleft[0] = 1;
	bottomleft[1] = 1;

	bottomright[0] = griddimension[0] - 2;
	bottomright[1] = 1;

	upperright[0] = griddimension[0] - 2;
	upperright[1] = griddimension[1] - 2;

	upperleft[0] = 1;
	upperleft[1] = griddimension[1] - 2;
}

/*
RealType* GridFunction::GetSlice(MultiIndexType bb, MultiIndexType ee){
	if (CheckInGrid(bb,ee)) {exit(0);}
	if (ee[0]-bb[0] >0 && ee[1]-bb[1]>0){
		std::cout<< "error: not a slice!\n";
		exit(0);
	}
	IndexType size;

	if (ee[1]-bb[1] > 1) {
		size = ee[1]-bb[1]+1;
		RealType slice[size];
		for (IndexType i = 0; i < size; i++) {
			slice[i] = gridfunction[bb[0]][bb[1]+i];
		}
		return slice;
	}
	else {
		size = ee[0]-bb[0]+1;
		RealType slice[size];
		for (IndexType i = 0; i < size; i++) {
			slice[i] = gridfunction[bb[0]+i][bb[1]];
		}
		return slice;
	}


}
*/


bool GridFunction::CheckInGrid(const MultiIndexType& begin, const MultiIndexType& end){
	bool error = false;
	if (begin[0]>end[0]){
		std::cout<< "error: begin[0]>end[0]\n";
		error = true;}
	if (begin[1]>end[1]){
		std::cout<< "error: begin[1]>end[1]\n";
		error = true;}
	if (begin[0]>griddimension[0]-1){
		std::cout<< "error: begin[0]>griddimension[0]-1\n";
		error = true;}
	if (end[0]>griddimension[0]-1){
		std::cout<< "error: end[0]>griddimension[0]-1\n";
		error = true;}
	if (begin[1]>griddimension[1]-1){
		std::cout<< "error: begin[1]>griddimension[1]-1\n";
		error = true;}
	if (end[1]>griddimension[1]-1){
		std::cout<< "error: end[1]>griddimension[1]-1\n";
		error = true;}
	return error;
}

