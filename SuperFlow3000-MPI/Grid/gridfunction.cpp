/*
 * gridfunction.cpp
 *
 *  Created on: 06.11.2013
 *      Author: Markus
 */

#include "gridfunction.h"
#include <iostream>
//1
GridFunction::GridFunction(int DimX, int DimY){
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
GridFunction::GridFunction(int DimX, int DimY, RealType value){
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
GridFunction::GridFunction(const MultiIndexType griddimension_input) : griddimension(griddimension_input){
	 gridfunction= new RealType*[griddimension[0]];
	 for (IndexType i = 0; i < griddimension[0]; i++){
		 gridfunction[i] = new RealType[griddimension[1]];
	 }
}

//2.1
GridFunction::GridFunction(const MultiIndexType griddimension_input,RealType value) : griddimension(griddimension_input){
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
	IndexType max = 0;
	for (IndexType i = begin[0];i<=end[0]; i++){
		for (IndexType j = begin[1]; j<=end[1]; j++){
			if (gridfunction[i][j] > max)
				max = gridfunction[i][j];
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

