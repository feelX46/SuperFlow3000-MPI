#include "IO.hpp"
#include <iostream>
#include <math.h>
#include <stdio.h>
using namespace std;


IO::IO (char *input, char *output) : output(output)
{
  //ToDo Read the file with the simulations parameters
 readInputfile(input);

}

IO::~IO ()
{

}

void
IO::readInputfile (char *filename)
{
	// open input file - char line for input lines
	ifstream file;
	char line[100];
	file.open(filename,ios::in);

	// read input file line by line
	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.xLength);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.yLength);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%d", &simparam.iMax);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%d", &simparam.jMax);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.tEnd);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.deltaT);


	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.tau);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.deltaVec);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%d", &simparam.iterMax);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.eps);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.omg);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.alpha);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.RE);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.GX);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.GY);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.UI);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.VI);

	file.getline(line,sizeof(line));
	sscanf(line, "%*[^0-9]%lf", &simparam.PI);


}

#define Element(field,ic) ((field)[(ic)[0]][(ic)[1]])

RealType
  IO::interpolateVelocityU (RealType x, RealType y, GridFunctionType & u,
			    const PointType & delta)
{

  RealType deltaX = delta[0];
  RealType deltaY = delta[1];

  MultiIndexType index;

  // Computation of u(x,y)
  index[0] = ((int) (x / deltaX)) + 1;
  index[1] = ((int) ((y + (deltaY / 2)) / deltaY)) + 1;

  // The coordinates of the cell corners

  RealType x1 = (index[0] - 1) * deltaX;
  RealType x2 = index[0] * deltaX;
  RealType y1 = ((index[1] - 1) - 0.5) * deltaY;
  RealType y2 = (index[1] - 0.5) * deltaY;

  MultiIndexType offset;

  offset[0] = index[0] - 1;
  offset[1] = index[1] - 1;

  RealType u1 = Element (u, offset);	// datafields->u->getField ()[i - 1][j - 1];

  offset[0] = index[0];
  offset[1] = index[1] - 1;

  RealType u2 = Element (u, offset);	//datafields->u->getField ()[i][j - 1];

  offset[0] = index[0] - 1;
  offset[1] = index[1];

  RealType u3 = Element (u, offset);	//datafields->u->getField ()[i - 1][j];
  RealType u4 = Element (u, index);

  RealType
    uInterploated =
    (1.0 / (deltaX * deltaY)) * ((x2 - x) * (y2 - y) *
				 u1 + (x - x1) * (y2 -
						  y) *
				 u2 + (x2 - x) * (y -
						  y1) *
				 u3 + (x - x1) * (y - y1) * u4);

  return uInterploated;
}


RealType
  IO::interpolateVelocityV (RealType x, RealType y, GridFunctionType & v,
			    const PointType & delta)
{
  RealType deltaX = delta[0];
  RealType deltaY = delta[1];

  // Computation of v(x,y)
  MultiIndexType index;
  index[0] = ((int) ((x + (deltaX / 2)) / deltaX)) + 1;
  index[1] = ((int) (y / deltaY)) + 1;

  // The coordinates of the cell corners

  RealType x1 = ((index[0] - 1) - 0.5) * deltaX;
  RealType x2 = (index[0] - 0.5) * deltaX;
  RealType y1 = (index[1] - 1) * deltaY;
  RealType y2 = index[1] * deltaY;

  MultiIndexType offset;

  offset[0] = index[0] - 1;
  offset[1] = index[1] - 1;

  RealType v1 = Element (v, offset);	//datafields->v->getField ()[i - 1][j - 1];

  offset[0] = index[0];
  offset[1] = index[1] - 1;

  RealType v2 = Element (v, offset);	//datafields->v->getField ()[i][j - 1];

  offset[0] = index[0] - 1;
  offset[1] = index[1];

  RealType v3 = Element (v, offset);	//datafields->v->getField ()[i - 1][j];


  RealType v4 = Element (v, index);	//datafields->v->getField ()[i][j];

  RealType
    vInterpolated =
    (1.0 / (deltaX * deltaY)) * ((x2 - x) * (y2 - y) *
				 v1 + (x - x1) * (y2 -
						  y) *
				 v2 + (x2 - x) * (y -
						  y1) *
				 v3 + (x - x1) * (y - y1) * v4);
  return vInterpolated;
}

void
IO::writeVTKFile (const MultiIndexType& griddimension, GridFunctionType u,
		  GridFunctionType v, GridFunctionType p,
		  const PointType& delta, int step, int mpiRank)
{
  RealType deltaX = delta[0];
  RealType deltaY = delta[1];

  IndexType iMax = griddimension[0] - 1;
  IndexType jMax = griddimension[1] - 1;

  char numstr[21];
  sprintf (numstr, "%d", step);
  char rank[21];
  sprintf(rank, "%d", mpiRank);
  std::string filename;
  filename.append ("./");
  filename.append (output);
  filename.append ("/");
  filename.append (rank);
  filename.append ("/");
  filename.append ("field_");
  filename.append (numstr);
  filename.append (".vts");
  std::cout << filename << std::endl;
  std::filebuf fb;
  fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
  std::ostream os (&fb);
  os << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"StructuredGrid\">" << std::endl
    << "<StructuredGrid WholeExtent=\""
    << "0" << " " << (iMax - 1) << " "
    << "0" << " " << (jMax - 1) << " "
    << "0" << " " << "0" << " "
    << "\" GhostLevel=\"" << "1" << "\">" << std::endl
    << "<Piece Extent=\""
    << "0" << " " << (iMax - 1) << " "
    << "0" << " " << (jMax - 1) << " "
    << "0" << " " << "0" << " "
    << "\">" << std::endl
    << "<Points>" << std::endl
    <<
    "<DataArray type=\"Float64\" format=\"ascii\" NumberOfComponents=\"3\"> "
    << std::endl;
  for (int i = 0; i < iMax; ++i)
    {
      for (int j = 0; j < jMax; ++j)
	{
	  os << std::scientific << i * deltaX << " " << j *
	    deltaY << " " << 0.0 << std::endl;
	}
    }
  os << "</DataArray>" << std::endl
    << "</Points>" << std::endl
    << "<PointData Vectors=\"field\"  Scalars=\"P\">"
    << std::endl <<
    "<DataArray Name=\"field\" NumberOfComponents=\"3\" type=\"Float64\" >" <<
    std::endl;
  for (int i = 0; i < iMax; ++i)
    {
      RealType x = i * deltaX;

      for (int j = 0; j < jMax; ++j)
	{
	  RealType y = j * deltaY;

	  os << std::scientific << interpolateVelocityU (x, y, u,
							 delta) << " " <<
	    interpolateVelocityV (x, y, v, delta) << " " << 0. << std::endl;
	}

    }
  os << "</DataArray>" << std::endl
    << "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">" <<
    std::endl;
  for (int i = 0; i <= iMax; ++i)
    {
      for (int j = 0; j <= jMax; ++j)
	{
	  os << std::scientific << p[i][j] << " ";

	}
      os << std::endl;

    }

  os << "</DataArray>" << std::endl
    << "</PointData>" << std::endl
    << "</Piece>" << std::endl
    << "</StructuredGrid>" << std::endl << "</VTKFile>" << std::endl;
  fb.close ();
  //std::cout<<"saved in "<< filename<<std::endl;
}


void IO::writeVTKMasterfile(const IndexType& mpiSizeH, const IndexType& mpiSizeV, int step,
		int localgriddimensionX, int localgriddimensionY){

	  int x1, x2, x3, x4;
	  int globalgriddimensionX = localgriddimensionX*mpiSizeH;
	  int globalgriddimensionY = localgriddimensionY*mpiSizeV;
	  int processorgridcoordX;
	  int processorgridcoordY;
	  int nprocessor = mpiSizeH*mpiSizeV;
	  char numstr[21];
	  sprintf (numstr, "%d", step);
	  std::string filename;
	  filename.append ("./");
	  filename.append (output);
	  filename.append ("/");
	  filename.append ("sol_");
	  filename.append (numstr);
	  filename.append ("_master.pvtr");

	  std::filebuf fb; // ?
	  fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
	  std::ostream os (&fb);
	  os << "<?xml version=\"1.0\"?>" << std::endl
	    << "<VTKFile type=\"PRectilinearGrid\">" << std::endl
	    << "<PRectilinearGrid WholeExtent=\""
	    << "0" << " " << (globalgriddimensionX - 2) << " "
	    << "0" << " " << (globalgriddimensionY - 2) << " "
	    << "0" << " " << "0" << " "
	    << "\" GhostLevel=\"" << "1" << "\">" << std::endl;
	  os <<"<PCoordinates>"<<std::endl;
	  os <<"<PDataArray type=\"Float64\"/>" << std::endl;
	  os <<"<PDataArray type=\"Float64\"/>" << std::endl;
	  os <<"<PDataArray type=\"Float64\"/>" << std::endl;
	  os <<"</PCoordinates>"<<std::endl;
	  for(int rank=0; rank<nprocessor; rank++){
			processorgridcoordX = rank % mpiSizeH;
			processorgridcoordY = floor(rank / mpiSizeH);
			x1=processorgridcoordX    *localgriddimensionX-processorgridcoordX;
			x2=(processorgridcoordX+1)*localgriddimensionX-processorgridcoordX-1;
			x3=processorgridcoordY    *localgriddimensionY-processorgridcoordY;
			x4=(processorgridcoordY+1)*localgriddimensionY-processorgridcoordY-1;


		  os << "<Piece Extent=\""<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<" 0 0 "
			<<"\" Source=\"sol_"<<numstr<<"_rank"<<rank<<".vtr\"/>"<<std::endl;
	  }

	  os <<"<PPointData>"<<std::endl;
	  os <<"<PDataArray Name=\"field\" NumberOfComponents=\"3\" type=\"Float64\" />"  <<std::endl;
	  os <<"<PDataArray type=\"Float64\" Name=\"p\"/>"<<std::endl
	    <<"</PPointData>"<<std::endl
        << "</PRectilinearGrid>" << std::endl
	    <<"</VTKFile>"<< std::endl;
	  fb.close ();
}

void IO::writeVTKSlavefile(GridFunction& u_gridfunction,
		  GridFunction& v_gridfunction, GridFunction& p_gridfunction,
		  const PointType& delta, int mpiSizeH, int mpiSizeV, int step,
		  int rank){

	double deltaX =delta[0];
	double deltaY =delta[1];
	//ToDo: einkommentieren!
	int ibegin = p_gridfunction.beginwrite[0];
	int iend   = p_gridfunction.endwrite[0];
	int jbegin = p_gridfunction.beginwrite[1];
	int jend   = p_gridfunction.endwrite[1];
	GridFunctionType p = p_gridfunction.GetGridFunction();
	GridFunctionType u = u_gridfunction.GetGridFunction();
	GridFunctionType v = v_gridfunction.GetGridFunction();
	int localgriddimensionX = iend-ibegin+1;
	int localgriddimensionY = jend-jbegin+1;
	int globalgriddimensionX = localgriddimensionX*mpiSizeH;
	int globalgriddimensionY = localgriddimensionY*mpiSizeV;

	int processorgridcoordX = rank % mpiSizeH;
	int processorgridcoordY = floor(rank / mpiSizeH);

	int x1=processorgridcoordX    *localgriddimensionX-processorgridcoordX;
	int x2=(processorgridcoordX+1)*localgriddimensionX-processorgridcoordX-1;
	int x3=processorgridcoordY    *localgriddimensionY-processorgridcoordY;
	int x4=(processorgridcoordY+1)*localgriddimensionY-processorgridcoordY-1;

	  char numstr[21];
	  sprintf (numstr, "%d", step);
	  std::string filename;
	  filename.append ("./");
	  filename.append (output);
	  filename.append ("/");
	  filename.append ("sol_");
	  filename.append (numstr);
	  filename.append ("_rank");
	  sprintf (numstr, "%d", rank);
	  filename.append (numstr);
	  filename.append (".vtr");

	  std::filebuf fb;
	  fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
	  std::ostream os (&fb);
	  os << "<?xml version=\"1.0\"?>" << std::endl
	  //  << "<VTKFile type=\"RectilinearGrid\">" << std::endl
		<< "<VTKFile type=\"RectilinearGrid\">" << std::endl
	    << "<RectilinearGrid WholeExtent=\""
	    << x1 << " " << x2 << " "
	    << x3 << " " << x4 << " "
	    << "0" << " " << "0" << " "
	    << "\" GhostLevel=\"" << "0" << "\">" << std::endl
	    << "<Piece Extent=\""<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<" 0 0 \">" <<std::endl
	    << "<Coordinates>"<<std::endl
	        <<
	        "<DataArray type=\"Float64\" format=\"ascii\"> "
	        << std::endl;
	      //coordinates
	      for (int i = ibegin; i <= iend; ++i)
	      {
	    	  os << std::scientific <<(processorgridcoordX*localgriddimensionX)*deltaX-processorgridcoordX*deltaX+ i * deltaX << " ";
	      }
	      os << std::endl;
	      os << "</DataArray>" << std::endl
	      << "<DataArray type=\"Float64\" format=\"ascii\"> " << std::endl;
	      for (int j = jbegin; j <= jend; ++j)
	      {
	    	  os << std::scientific << (processorgridcoordY*localgriddimensionY)*deltaY-processorgridcoordY*deltaY + j * deltaY << " ";
	      }
	      os << std::endl
	      << "</DataArray>" << std::endl
	      << "<DataArray type=\"Float64\" format=\"ascii\"> " << std::endl
	      << "0 0" << std::endl
	      << "</DataArray>" << std::endl

	        << "</Coordinates>" << std::endl
	        << "<PointData>"
	        << std::endl <<
	        "<DataArray Name=\"field\" NumberOfComponents=\"3\" type=\"Float64\" >" <<
	        std::endl;
	      //velocities
	      for (int j = jbegin; j <= jend; ++j)
	      	    	{
	      	    	  RealType y = j * deltaY;

	      for (int i = ibegin; i <= iend; ++i)
	        {
	          RealType x = i * deltaX;



	    	  os << std::scientific << interpolateVelocityU (x, y, u,
	    							 delta) << " " <<
	    	    interpolateVelocityV (x, y, v, delta) << " " << 0. << std::endl;
	    	}

	        }
	      os << "</DataArray>" << std::endl
	        << "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">" <<
	        std::endl;
	      //pressure
          for (int j = jbegin; j <= jend; ++j)
          {
        	  for (int i = ibegin; i <= iend; ++i)
        	  {
	    	  os << std::scientific << p[i][j] << " ";

	    	}
	          os << std::endl;

	        }

	      os << "</DataArray>" << std::endl
	        << "</PointData>" << std::endl
	        << "</Piece>" << std::endl
	        << "</RectilinearGrid>" << std::endl << "</VTKFile>" << std::endl;
	      fb.close ();

}
