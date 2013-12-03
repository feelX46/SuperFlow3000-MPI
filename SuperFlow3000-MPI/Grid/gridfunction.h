/*
 * gridfunction.h
 *
 *  Created on: 08.11.2013
 *      Author: Markus
 */

#ifndef GRIDFUNCTION_H_
#define GRIDFUNCTION_H_

#include"../Misc/typedef.h"

class GridFunction {

public:
	/*! Construktor (1)
	   * @param DimX Dimension in X-Direction
	   * @param DimY Dimension in Y-Direction
	   */
	GridFunction(int DimX, int DimY, char indicator);

	/*! Construktor (1.1)
		   * @param DimX Dimension in X-Direction
		   * @param DimY Dimension in Y-Direction
		   * @param value initial value
		   */
	GridFunction(int DimX, int DimY, RealType value, char indicator);

	/*! Construktor (2)
		   * @param griddimension_input grid dimension
		   */
	GridFunction(const MultiIndexType griddimension_input, char indicator);

	/*! Construktor (2.1)
		   * @param griddimension_input grid dimension
		   * @param value inital value
		   */
	GridFunction(const MultiIndexType griddimension_input, RealType value, char indicator);

	//! Destructor (3)
	~GridFunction();

	//! get grid data (pointer to gridfunction) (4)
	GridFunctionType GetGridFunction() const;

	//! get value of index
	RealType GetGridFunction(const MultiIndexType& index);

	//! get dimension
	MultiIndexType GetGridDimension();

	//! set a whole block to value (5)
	void SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,RealType value);

	//! set one value
	void SetGridFunction (const IndexType& i, const IndexType& j, RealType value);

	//! copy values with offset (6)
	void SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,RealType factor,
				MultiIndexType& offset);

    //! copy factor*source in grid (8)
	void SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,RealType factor,
			GridFunctionType& sourcegridfunction);

	//! copy factor*source(with offset) in grid (9)
	void SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end, RealType factor,
			GridFunctionType& sourcegridfunction, MultiIndexType& offset);

	//! copy constant + factor*source(with offset) in grid (~10)
	void SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,RealType factor,
			GridFunctionType& sourcegridfunction, MultiIndexType& offset, RealType constant);

	//! copy constant + factor*source in grid (~10)
	void SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,RealType factor,
			GridFunctionType& sourcegridfunction, RealType constant);

	//!  constant + factor*grid(with offset) (combination of 6 and 10)
	void SetGridFunction (const MultiIndexType& begin, const MultiIndexType& end,
			RealType factor, MultiIndexType& offset, RealType constant);

	//! scale grid (7)
	void ScaleGridFunction (const MultiIndexType& begin, const MultiIndexType& end, RealType factor);

	//! add factor*source to grid (11)
	void AddToGridFunction (const MultiIndexType& begin, const MultiIndexType& end,RealType factor,
			GridFunctionType& sourcegridfunction);
	//! add factor*source to grid (11-1)
	void AddToGridFunction (const MultiIndexType& begin, const MultiIndexType& end,RealType factor,
			GridFunctionType& sourcegridfunction, MultiIndexType& offset);

	//RealType* GetSlice(MultiIndexType bb, MultiIndexType ee);


	//! find maximum (12)
	RealType MaxValueGridFunction (const MultiIndexType& begin, const MultiIndexType& end);

	//! plot grid in console
	void PlotGrid();

	//ToDo wieder private machen
	//! gridfunction Datamatrix
	GridFunctionType gridfunction;
	//! dimension of the grid
	MultiIndexType	griddimension;

	bool globalboundary[4];

	MultiIndexType bottomleft;
	MultiIndexType bottomright;
	MultiIndexType upperright;
	MultiIndexType upperleft;

	MultiIndexType beginread;
	MultiIndexType endread;
	MultiIndexType beginwrite;
	MultiIndexType endwrite;

	void InitializeGlobalBoundaryPosition(int rank, int mpiSizeH, int mpiSizeV, char indicator);

	void InitializeGlobalBoundary(char indicator);

	bool CheckInGrid(const MultiIndexType& begin, const MultiIndexType& end);
};



#endif /* GRIDFUNCTION_H_ */
