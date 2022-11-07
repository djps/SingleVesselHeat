/**
 @file conductivity.cpp
 */

#include "conductivity.h"


//#define DEBUG_CONDUCTIVITY

/**
 @brief Constructor
 */
conductivity::conductivity()
{

  std::cout << "Beware: the default constructor is a stupid option for this class. ";

  int nx = 1;
  int ny = 1;
  int nz = 1;

  double dx = 0.1;
  double dy = 0.1;
  double dz = 0.1;

  cond_x = new Array3D<double>(M,ny,nz);
  cond_y = new Array3D<double>(nx,M,nz);
  cond_z = new Array3D<double>(nx,ny,M);

  distance_x = new Array3D<double>(M,ny,nz);
  distance_y = new Array3D<double>(nx,M,nz);
  distance_z = new Array3D<double>(nx,ny,M);

  index_x = new Array3D<int>(M,ny,nz);
  index_y = new Array3D<int>(nx,M,nz);
  index_z = new Array3D<int>(nx,ny,M);
}

/**
 @brief Deconstructor of conductivity class
 */
conductivity::~conductivity()
{
  delete index_x;
  delete index_y;
  delete index_z;

  delete cond_x;
  delete cond_y;
  delete cond_z;

  delete distance_x;
  delete distance_y;
  delete distance_z;
}

/**
 @brief Constructor of conductivity class

 Computes the index of the cell in the mesh which is closest to the point query.

 *Arguments*
   nx_in ( int )
          An integer
     ny_in ( int )
         An integer
     nz_in ( int )
         An integer
     dx_in ( int )
         An integer
     dy_in ( int )
         An integer
     dz_in ( int )
         An integer
     number_of_conductivity_parameters ( int )
         An integer
     conductivity_parameters ( double )
         An integer



 */
conductivity::conductivity(int nx_in, int ny_in, int nz_in, double dx_in, double dy_in, \
 double dz_in, int number_of_conductivity_parameters, double *conductivity_parameters)
{

    bool verbose = false;

    if (verbose){
        std::cout << " In constructor" << std::endl;
    }

    nx = nx_in;
	ny = ny_in;
	nz = nz_in;
    dx = dx_in;
	dy = dy_in;
	dz = dz_in;

    cond_x = new Array3D<double>(M,ny,nz);
    cond_y = new Array3D<double>(nx,M,nz);
    cond_z = new Array3D<double>(nx,ny,M);

    if (verbose){
        std::cout << " Got space for cond" << std::endl;
    }

    distance_x = new Array3D<double>(M,ny,nz);
    distance_y = new Array3D<double>(nx,M,nz);
    distance_z = new Array3D<double>(nx,ny,M);

    if (verbose){
        std::cout << " Got space for dist" << std::endl;
    }

    index_x = new Array3D<int>(M,ny,nz);
    index_y = new Array3D<int>(nx,M,nz);
    index_z = new Array3D<int>(nx,ny,M);

    if (verbose){
        std::cout << " Got space for indices" << std::endl;
    }

    GenerateConductivity(number_of_conductivity_parameters, conductivity_parameters);
}



/**
 @brief Subroutine which generates conductivity
 */
void conductivity::GenerateConductivity(int number_of_conductivity_parameters, double *conductivity_parameters )
{
    bool verbose = false;

	// Puts data into the different arrays.
	// OK: The parameters are the following:
	//	1) x and y values of the centre line, in the same units as dx and dy are in, i.e. probably mm.
	//     Note that (0,0,0) is the bottom left hand side of the grid, not the centre.
	//	2) Radius of cylinder
	//	3) Conductivity outside and inside the vessel.
	// Should therefore have 5 parameters.
	// The method is quite simple. In the z direction, either we are in the vessel or we aren't.
	// So conductivity is uniform. In the x and y directions, it is really easy to figure out
	// the positions of the vessel edges, etc.

	if (number_of_conductivity_parameters != 5)
	{
		std::cerr << "Trouble, wrong number of conductivity parameters passed" << std::endl;
		std::exit(2);
	}

    double xC = conductivity_parameters[0];
	double yC = conductivity_parameters[1];
	double radius = conductivity_parameters[2];
	double kapBlood = conductivity_parameters[3];
	double kapTissue = conductivity_parameters[4];

	double tmp1;
	double tmp2;
	int posLeft;
	int posRight;
	// First, go through each of y and z and generate the rows in the x direction.
	for (int y = 0; y < ny; ++y)
	{
		if (((y*dy-yC)*(y*dy-yC)) > radius*radius)
		{
            // Then we are in a fully tissue region.
			for (int z = 0; z < nz; ++z)
			{
				for (int pp = 0; pp < M; pp++)
				{
					cond_x->Input_Value(pp,y,z,kapTissue); // Input the conductivity.
					index_x->Input_Value(pp,y,z,nx-1); // Boundary is right at the end.
					distance_x->Input_Value(pp,y,z,0.01); // No point in having this.  But might avoid a program crash.
				}
			}
		}
		else
		{
			for (int z = 0; z < nz; ++z)
			{
				tmp1= xC/dx;
				tmp2 = sqrt(radius*radius - (y*dy-yC)*(y*dy-yC))/dx;
				posLeft = (floor(tmp1-tmp2) > 0) ? int(floor(tmp1-tmp2)): 0;  // If the test fails, the cylinder juts out of the left-hand side of the picture.
				posRight = (floor(tmp1+tmp2) > nx-1) ? nx-1: int(floor(tmp1+tmp2));
				if (posLeft == 0)
				{
					// then we start in the blood
					cond_x->Input_Value(0,y,z,kapBlood);
					index_x->Input_Value(0,y,z,posRight); // This is the point where we change.

                    double eps = 0.01;
					if ( ((tmp1+tmp2-posRight) > eps) && ((tmp1+tmp2-posRight) < 1-eps) ){
                         distance_x->Input_Value(0,y,z,tmp1+tmp2-posRight); // The distance to the point.
                    }
					else if ( (tmp1+tmp2-posRight) <= eps ) {
                        distance_x->Input_Value(0,y,z,eps);
                    }
					else {
                        distance_x->Input_Value(0,y,z,1-eps);
                    }

					cond_x->Input_Value(1,y,z,kapTissue); // The next value is in tissue, and continues as such.
					index_x->Input_Value(1,y,z,nx-1); // ADDED 10.03.2010. Must put in the final value of index.
				}
				else
				{
					cond_x->Input_Value(0,y,z,kapTissue);
					index_x->Input_Value(0,y,z,posLeft); // This is the point where we change.
					//DELETED 10.03.2010: distance_x->Input_Value(0,y,z,tmp1-tmp2-posLeft); // The distance to the point.

                    double eps = 0.01;
                    if ( ((tmp1-tmp2-posLeft) > eps) && ((tmp1-tmp2-posLeft) < 1-eps) ){
                         distance_x->Input_Value(0,y,z,tmp1-tmp2-posLeft); // The distance to the point.
					}
                    else if ( (tmp1-tmp2-posLeft) <= eps ){
                        distance_x->Input_Value(0,y,z,eps);
                    }
					else{
                        distance_x->Input_Value(0,y,z,1-eps);
                    }

                    cond_x->Input_Value(1,y,z,kapBlood); // The next value is in tissue, and continues as such.
					index_x->Input_Value(1,y,z,posRight); // This is the point where we change.

					if ( ((tmp1+tmp2-posRight) > eps) && ((tmp1+tmp2-posRight) < 1-eps) ){
                        distance_x->Input_Value(1,y,z,tmp1+tmp2-posRight); // The distance to the point.
                    }
					else if ((tmp1+tmp2-posRight) <= eps){
                        distance_x->Input_Value(1,y,z,eps);
                    }
					else {
                        distance_x->Input_Value(1,y,z,1-eps);
                    }

                    cond_x->Input_Value(2,y,z,kapTissue); // The next value is in tissue, and continues as such.
					index_x->Input_Value(2,y,z,nx-1); // ADDED 10.03.2010. Must put in the final value of index.
				}
			}
		}
	}

    if (verbose){
      std::cout << " Done x" << std::endl;
    }

	for (int x = 0; x < nx; ++x)
	{
		if (((x*dx-xC)*(x*dx-xC)) > radius*radius)
		{
			// Then we are in a fully tissue region.
			for (int z = 0; z < nz; ++z)
			{
				for (int pp = 0; pp < M; pp++)
				{
					cond_y->Input_Value(x,pp,z,kapTissue); // Input the conductivity.
					index_y->Input_Value(x,pp,z,ny-1); // Boundary is right at the end.
					distance_y->Input_Value(x,pp,z,0.01); // No point in having this.  But might avoid a program crash.
				}
			}
		}
		else
		{
			for (int z = 0; z<nz; ++z)
			{
				tmp1= yC/dy;
				tmp2 = sqrt(radius*radius - (x*dx-xC)*(x*dx-xC))/dy;
				posLeft = (floor(tmp1-tmp2) > 0) ? int(floor(tmp1-tmp2)) : 0;  // If the test fails, the cylinder juts out of the left-hand side of the picture.
				posRight = (floor(tmp1+tmp2) > ny-1) ? ny-1 : int(floor(tmp1+tmp2));
				if (posLeft == 0)
				{
					// then we start in the blood
					cond_y->Input_Value(x,0,z,kapBlood);
					index_y->Input_Value(x,0,z,posRight); // This is the point where we change.
                    double eps = 0.01;
					if ( ((tmp1+tmp2-posRight) > eps) && ((tmp1+tmp2-posRight) < 1-eps) ){
                        distance_y->Input_Value(x,0,z,tmp1+tmp2-posRight); // The distance to the point.
                    }
					else if ((tmp1+tmp2-posRight) <= eps){
                        distance_y->Input_Value(x,0,z,eps);
                    }
					else{
                        distance_y->Input_Value(x,0,z,1-eps);
                    }
					cond_y->Input_Value(x,1,z,kapTissue); // The next value is in tissue, and continues as such.
					index_y->Input_Value(x,1,z,ny-1); // ADDED 10.03.2010: Include right-hand boundary.
				}
				else
				{
					cond_y->Input_Value(x,0,z,kapTissue);
					index_y->Input_Value(x,0,z,posLeft); // This is the point where we change.
                    double eps = 0.01;
					if ( ((tmp1-tmp2-posLeft) > eps) && ((tmp1-tmp2-posLeft) < 1-eps) ){
                        distance_y->Input_Value(x,0,z,tmp1-tmp2-posLeft); // The distance to the point.
                    }
					else if ( (tmp1-tmp2-posLeft) <= eps ){
                        distance_y->Input_Value(x,0,z,eps);
                    }
					else{
                        distance_y->Input_Value(x,0,z,1-eps);
                    }
					cond_y->Input_Value(x,1,z,kapBlood); // The next value is in tissue, and continues as such.
					index_y->Input_Value(x,1,z,posRight); // This is the point where we change.
					if ( ((tmp1+tmp2-posRight) > eps) && ((tmp1+tmp2-posRight) < 1-eps) ){
                         distance_y->Input_Value(x,1,z,tmp1+tmp2-posRight); // The distance to the point.
                     }
					else if ( (tmp1+tmp2-posRight) <= eps ){
                        distance_y->Input_Value(x,1,z,eps);
                    }
					else{
                        distance_y->Input_Value(x,1,z,1-eps);
                    }
					cond_y->Input_Value(x,2,z,kapTissue); // The next value is in tissue, and continues as such.
					index_y->Input_Value(x,2,z,ny-1); // ADDED 10.03.2010: Include right-hand boundary.
				}
			}
		}
	}

    for (int x = 0; x < nx; ++x)
    {
    	for (int y = 0; y < ny; ++y)
        {
        	// YES, this could be done faster. I don't care. This is a once-off part of the code, better to get right than fast.
        	if ( ((x*dx-xC)*(x*dx-xC) + (y*dy-yC)*(y*dy-yC)) > radius*radius )
        	{
        		// Then we are in a fully tissue region.
        		for (int pp = 0; pp < M; pp++)
        		{
        			cond_z->Input_Value(x,y,pp,kapTissue); // Input the conductivity.
                    index_z->Input_Value(x,y,pp,nz-1); // Boundary is right at the end.
        		}
        	}
        	else
        	{
        		for (int pp = 0; pp < M; pp++)
        		{
        			cond_z->Input_Value(x,y,pp,kapBlood); // Input the conductivity.
                    index_z->Input_Value(x,y,pp,nz-1); // Boundary is right at the end.
        		}
        	}
        }
    }

    if (verbose){
      std::cout << " Done z" << std::endl;
    }

}


/*! Outputs int, index conductivity in specific direction */
int conductivity::OutputIndex(int count, int ind1, int ind2, int option)
{
    if (option==1){
        return index_x->Output_Value(count,ind1,ind2);
    }
    else if (option==2){
        return index_y->Output_Value(ind1,count,ind2);
    }
    else if (option==3){
        return index_z->Output_Value(ind1,ind2,count);
    }
    else{
        return -1;
    }
}


/**
 @brief Outputs double of conductivity in specific direction
 */
double conductivity::OutputCond(int count, int ind1, int ind2, int option)
{
    if (option==1){
        return cond_x->Output_Value(count,ind1,ind2);
    }
    else if (option==2){
        return cond_y->Output_Value(ind1,count,ind2);
    }
    else if (option==3){
        return cond_z->Output_Value(ind1,ind2,count);
    }
    else{
        return -1;
    }
}

/**
 @brief Outputs double of dist conductivity in specific direction
 */
double conductivity::OutputDist(int count, int ind1, int ind2, int option)
{
    if (option==1){
        return distance_x->Output_Value(count,ind1,ind2);
    }
    else if (option==2){
        return distance_y->Output_Value(ind1,count,ind2);
    }
    else if (option==3){
        return distance_z->Output_Value(ind1,ind2,count);
    }
    else{
        return -1;
    }
}

/**
 @brief this is a useful debugging routine.
 */
void conductivity::Write(int ind1, int ind2, int option)
{
    if (option==1){
        for (int p = 0; p < M; ++p)
        {
            std::cout << index_x->Output_Value(p,ind1,ind2) << " ";
        }
    }
    else if (option==2){
        for (int p = 0; p< M; ++p){
            std::cout << index_y->Output_Value(ind1,p,ind2) << " ";
        }
    }
    else if (option==3){
        for (int p = 0; p< M; ++p){
            std::cout << index_z->Output_Value(ind1,ind2,p) << " ";
        }
    }
    std::cout << std::endl;
};
