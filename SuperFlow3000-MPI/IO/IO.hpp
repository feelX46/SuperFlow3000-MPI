//! The class implements the IO 
/*!
 * @author diehlpk
 * @date 2012
 */

#include "../Misc/typedef.h"

#ifndef IO_HPP_
#define IO_HPP_

#include <iostream>
#include <fstream>
#include "../Structs/simparam.h"
#include "../Grid/gridfunction.h"

class IO
{
public:

  /*! Construktor
   * @param input Path to the file with the simulation parameters.
   * @param outout Path to the directory for the vtk files.
   */
  IO (char *input, char *output);

  //! Destructor
   ~IO ();

  //! Method writes the GridFunctions u,v,p in the vtk data format to the hard disk. 
  //! The files are named in the following convention: field_(step).vts.
  /*!
   * \param griddimension The dimension of the gridfunctions.
   * \param u The gridfunction u.
   * \param v The gridfunction v.
   * \param p The gridfunction p.
   * \param delta The mesh width in x direction and y direction
   * \param step The number of the timestep.
   */
  void writeVTKFile (const MultiIndexType& griddimension,
		     GridFunctionType u, GridFunctionType v,
		     GridFunctionType p, const PointType& delta, int step, int mpiRank);

  // new!!! for MPI
    void writeVTKMasterfile(const IndexType& mpiSizeH, const IndexType& mpiSizeV, int step,
  			int localgriddimensionX, int localgriddimensionY);

    void writeVTKSlavefile (GridFunction& u_gridfunction,
  		  GridFunction& v_gridfunction, GridFunction& p_gridfunction,
  		  const PointType& delta, int mpiSizeH, int mpiSizeV, int step,
  		  int rank);




  //! Method that returns private member variable simparam
  Simparam getSimparam () {
	  return simparam;
  }

private:

//! Path where to write the vtk files.
  char *output;
  Simparam simparam;
 //! Structure with simulation parameters from input file




/*!
   * Methods reads the simulation parameters from the specified input file.
   * 
   * @param filename The name of the file with the simulations paremters
   */
  void readInputfile (char *filename);


  //! Method interpolates the velocity for u in the staggered grid.
  /*!
   * \param x Value of the x coordinate.
   * \param y Value of the y coordinate.
   * \param u The gridfunction u.
   * \param delta The mesh width in x direction and y direction.
   */
  RealType interpolateVelocityU (RealType x, RealType y, GridFunctionType & u,
				 const PointType & delta);

  //! Method interpolates the velocity for v in the staggered grid.
  /*!
   * \param x Value of the x coordinate.
   * \param y Value of the y coordinate.
   * \param u The gridfunction v.
   * \param delta The mesh width in x direction and y direction.
   */
  RealType interpolateVelocityV (RealType x, RealType y, GridFunctionType & v,
				 const PointType & delta);

};

#endif
