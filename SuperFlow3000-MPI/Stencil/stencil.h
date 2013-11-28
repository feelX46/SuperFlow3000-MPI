/*
 * stencil.h
 *
 *  Created on: Nov 6, 2013
 *      Author: ischifx
 */

#ifndef STENCIL_H_
#define STENCIL_H_


#include"../Misc/typedef.h"
#include"../Grid/gridfunction.h"

class Stencil {
public:

	Stencil(int stencilwidth_input, const PointType& h_input);
	~Stencil();

	StencilType stencil;
	int stencilwidth;
	const PointType& h;

	void ApplyStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunction,
			GridFunction& imagegridfunction);

	void ApplyFxStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunction,
			GridFunction& imagegridfunction);

	void ApplyFyStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunction,
			GridFunction& imagegridfunction);

	void ApplyFxxStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunction,
			GridFunction& imagegridfunction);

	void ApplyFyyStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunction,
			GridFunction& imagegridfunction);

	void ApplyPxStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunction,
			GridFunction& imagegridfunction);

	void ApplyUSqxStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunction,
			GridFunction& imagegridfunction,
			RealType alpha);

	void ApplyVSqyStencilOperator(const MultiIndexType& gridreadbegin,
				const MultiIndexType& gridreadend,
				const MultiIndexType& gridwritebegin,
				const MultiIndexType& gridwriteend,
				const GridFunctionType& sourcegridfunction,
				GridFunction& imagegridfunction,
				RealType alpha);

	void ApplyUVyStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunctionU,
			const GridFunctionType& sourcegridfunctionV,
			GridFunction& imagegridfunction,
			RealType alpha);

	void ApplyUVxStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			const GridFunctionType& sourcegridfunctionU,
			const GridFunctionType& sourcegridfunctionV,
			GridFunction& imagegridfunction,
			RealType alpha);

	void setFxStencil();
	void setFyStencil();
	void setFxxStencil();
	void setFyyStencil();
	void setPxStencil();
};

#endif /* STENCIL_H_ */
