/*
 * stencil.cpp
 *
 *  Created on: Nov 6, 2013
 *      Author: ischifx
 */

#include "stencil.h"
#include <iostream>
#include <math.h>


Stencil::Stencil(int stencilwidth_input, const PointType& h_input) : stencilwidth(stencilwidth_input), h(h_input){
	 stencil = new RealType*[stencilwidth];
	 for (int i = 0; i < stencilwidth; i++){
		 stencil[i] = new RealType [stencilwidth];
	 }
}

Stencil::~Stencil() {
	for (int i = 0; i < stencilwidth; i++)
			delete[] stencil[i];
		delete[] stencil;
}

void Stencil::ApplyStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunction,
		GridFunction& imagegridfunction){
	// Berechne die Ableitungen
	// (0,0) ist bei allen allen drei Matrtzen (Stencil, sourcegrid, imagegrid) oben links
	RealType tmp;
	for (IndexType i=gridwritebegin[0]; i<=gridwriteend[0]; i++){
		for (IndexType j=gridwritebegin[1]; j<=gridwriteend[1]; j++){
			tmp = 0.0;
			for(IndexType k=0; k<3; k++){
				for(IndexType l=0; l<3; l++){
					tmp = tmp + sourcegridfunction[i+k-1][j+l-1]*stencil[k][l];
				}
			}
			imagegridfunction.SetGridFunction(i,j,tmp);
	 	}
	}

}

void Stencil::ApplyFxStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunction,
		GridFunction& imagegridfunction){
	setFxStencil();
	ApplyStencilOperator(gridreadbegin,gridreadend,gridwritebegin,gridwriteend, sourcegridfunction, imagegridfunction);
}

void Stencil::ApplyFyStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunction,
		GridFunction& imagegridfunction){
	setFyStencil();
	ApplyStencilOperator(gridreadbegin,gridreadend,gridwritebegin,gridwriteend, sourcegridfunction, imagegridfunction);
}



void Stencil::ApplyFxxStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunction,
		GridFunction& imagegridfunction){
	setFxxStencil();
	ApplyStencilOperator(gridreadbegin,gridreadend,gridwritebegin,gridwriteend, sourcegridfunction, imagegridfunction);
}

void Stencil::ApplyFyyStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunction,
		GridFunction& imagegridfunction){
	setFyyStencil();
	ApplyStencilOperator(gridreadbegin,gridreadend,gridwritebegin,gridwriteend, sourcegridfunction, imagegridfunction);
}

void Stencil::ApplyPxStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunction,
		GridFunction& imagegridfunction){
	setPxStencil();
	ApplyStencilOperator(gridreadbegin,gridreadend,gridwritebegin,gridwriteend, sourcegridfunction, imagegridfunction);
}

void Stencil::ApplyUSqxStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunction,
		GridFunction& imagegridfunction,
		RealType alpha){
	RealType tmp;
	for (IndexType i=gridwritebegin[0]; i<=gridwriteend[0]; i++){
		for (IndexType j=gridwritebegin[1]; j<=gridwriteend[1]; j++){
			tmp = 0.5*(sourcegridfunction[i][j]*sourcegridfunction[i+1][j]-
					sourcegridfunction[i-1][j]*sourcegridfunction[i][j]+
					sourcegridfunction[i+1][j]*sourcegridfunction[i+1][j]*0.5-
					sourcegridfunction[i-1][j]*sourcegridfunction[i-1][j]*0.5)+
				  alpha*0.25*(abs(sourcegridfunction[i][j]+sourcegridfunction[i+1][j])*
						  (sourcegridfunction[i][j]-sourcegridfunction[i+1][j])-
						      abs(sourcegridfunction[i-1][j]+sourcegridfunction[i][j])*
						      (sourcegridfunction[i-1][j]-sourcegridfunction[i][j])
						  );

			imagegridfunction.SetGridFunction(i,j,tmp/h[0]);

		}
	}
}

void Stencil::ApplyVSqyStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunction,
		GridFunction& imagegridfunction,
		RealType alpha){
	RealType tmp;
	for (IndexType i=gridwritebegin[0]; i<=gridwriteend[0]; i++){
		for (IndexType j=gridwritebegin[1]; j<=gridwriteend[1]; j++){
			tmp = 0.5*(sourcegridfunction[i][j]*sourcegridfunction[i][j+1]-
					sourcegridfunction[i][j-1]*sourcegridfunction[i][j]+
					sourcegridfunction[i][j+1]*sourcegridfunction[i][j+1]*0.5-
					sourcegridfunction[i][j-1]*sourcegridfunction[i][j-1]*0.5)+
				  alpha*0.25*(abs(sourcegridfunction[i][j]+sourcegridfunction[i][j+1])*
						  (sourcegridfunction[i][j]-sourcegridfunction[i][j+1])-
						      abs(sourcegridfunction[i][j-1]+sourcegridfunction[i][j])*
						      (sourcegridfunction[i][j-1]-sourcegridfunction[i][j])
						  );

			imagegridfunction.SetGridFunction(i,j,tmp/h[1]);

		}
	}
}


void Stencil::ApplyUVyStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunctionU,
		const GridFunctionType& sourcegridfunctionV,
		GridFunction& imagegridfunction,
		RealType alpha){
	RealType tmp;
	for (IndexType i=gridwritebegin[0]; i<=gridwriteend[0]; i++){
		for (IndexType j=gridwritebegin[1]; j<=gridwriteend[1]; j++){
			tmp = 0.25*((sourcegridfunctionV[i][j]+sourcegridfunctionV[i+1][j])*
						(sourcegridfunctionU[i][j]+sourcegridfunctionU[i][j+1])-
						(sourcegridfunctionV[i][j-1]+sourcegridfunctionV[i+1][j-1])*
						(sourcegridfunctionU[i][j-1]+sourcegridfunctionU[i][j])
			)+
							  alpha*0.25*(abs(sourcegridfunctionV[i][j]+sourcegridfunctionV[i+1][j])*
									  (sourcegridfunctionU[i][j]-sourcegridfunctionU[i][j+1])-
									      abs(sourcegridfunctionV[i][j-1]+sourcegridfunctionV[i+1][j-1])*
									      (sourcegridfunctionU[i][j-1]-sourcegridfunctionU[i][j])
									  );

			imagegridfunction.SetGridFunction(i,j,tmp/h[1]);
		}
	}
}


void Stencil::ApplyUVxStencilOperator(const MultiIndexType& gridreadbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend,
		const GridFunctionType& sourcegridfunctionU,
		const GridFunctionType& sourcegridfunctionV,
		GridFunction& imagegridfunction,
		RealType alpha){
	RealType tmp;
	for (IndexType i=gridwritebegin[0]; i<=gridwriteend[0]; i++){
		for (IndexType j=gridwritebegin[1]; j<=gridwriteend[1]; j++){
			tmp = 0.25*((sourcegridfunctionU[i][j]+sourcegridfunctionU[i][j+1])*
						(sourcegridfunctionV[i][j]+sourcegridfunctionV[i+1][j])-
						(sourcegridfunctionU[i-1][j]+sourcegridfunctionU[i-1][j+1])*
						(sourcegridfunctionV[i-1][j]+sourcegridfunctionV[i][j])
			)+
							  alpha*0.25*(abs(sourcegridfunctionU[i][j]+sourcegridfunctionU[i][j+1])*
									  (sourcegridfunctionV[i][j]-sourcegridfunctionV[i+1][j])-
									      abs(sourcegridfunctionU[i-1][j]+sourcegridfunctionU[i-1][j+1])*
									      (sourcegridfunctionV[i-1][j]-sourcegridfunctionV[i][j])
									  );

			imagegridfunction.SetGridFunction(i,j,tmp/h[0]);
		}
	}
}



void Stencil::setFxStencil() {
	stencil[0][0] = 0;
	stencil[1][0] = 0;
	stencil[2][0] = 0;
	stencil[0][1] = -1/h[0];
	stencil[1][1] = 1/h[0];
	stencil[2][1] = 0;
	stencil[0][2] = 0;
	stencil[1][2] = 0;
	stencil[2][2] = 0;
}

void Stencil::setFyStencil() {
	stencil[0][0] = 0;
	stencil[1][0] = -1/h[1];
	stencil[0][0] = 0;
	stencil[0][1] = 0;
	stencil[1][1] = 1/h[1];
	stencil[2][1] = 0;
	stencil[0][2] = 0;
	stencil[1][2] = 0;
	stencil[2][2] = 0;

}


void Stencil::setFxxStencil() {
	stencil[0][0] = 0;
	stencil[1][0] = 0;
	stencil[2][0] = 0;
	stencil[0][1] = 1/(h[0]*h[0]);
	stencil[1][1] = -2/(h[0]*h[0]);
	stencil[2][1] = 1/(h[0]*h[0]);
	stencil[0][2] = 0;
	stencil[1][2] = 0;
	stencil[2][2] = 0;
}

void Stencil::setFyyStencil() {
	stencil[0][0] = 0;
	stencil[1][0] = 1/(h[1]*h[1]);
	stencil[0][0] = 0;
	stencil[0][1] = 0;
	stencil[1][1] = -2/(h[1]*h[1]);
	stencil[2][1] = 0;
	stencil[0][2] = 0;
	stencil[1][2] = 1/(h[1]*h[1]);
	stencil[2][2] = 0;
}

void Stencil::setPxStencil() {
	stencil[0][0] = 0;
	stencil[1][0] = 0;
	stencil[0][0] = 0;
	stencil[0][1] = 0;
	stencil[1][1] = -1/h[0];
	stencil[2][1] = 1/h[0];
	stencil[0][2] = 0;
	stencil[1][2] = 0;
	stencil[2][2] = 0;
}

