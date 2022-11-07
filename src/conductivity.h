#ifndef CONDUCTIVITY_H
#define CONDUCTIVITY_H

#include "array3d.cpp"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

/**
 \file conductivity.h
This class contains the structures needed for a reasonably efficient
implementation of the conductivity.  The reason that this is so much weirder
than the class for the velocity is that velocity is assumed smooth, while
conductivity is not. As such, there will be boundaries that must be implemented
and the conductivity contains info regarding these.

This is how I do it. For each pair of (y,z) we generate 3 vectors each with M
elements. These vectors proceed in the direction of increasing x. One of them
stores the values of the conductivity in each region.  One of them stores the x
index at which the boundary occurs (or more precisely the largest value of the
index before the boundary) and the third stores the distance, in units of dx,
away from this final node that the boundary is located. This latter distance
must be a positive number between 0 and 1, but I shall actually force it to be
between 0.01 and 0.99 to ensure that the boundary is not too close to a node.
This structure is repeated also for pairs of y,z and pairs of x,y.

Finally note that if consecutive indices denoting the position of the boundary
are the same then this boundary is ignored completely. This allows an easy way
to make the problem homogeneous.

@author Simon Woodford <simon.woodford@icr.ac.uk>
**/

const int M = 4; // This determines the maximum number of jumps in any direction that the conductivity can have.

const double tolerance = 1.e-10; // sometimes double precision numbers are very close but not quite equal. This allows a certain flexibility in deciding when the conductivity really is uniform and when not.


/**
 @class conductivity
 */
class conductivity{

// The storage of this is rather straightforward.

public:
	// default constructor
    conductivity();

    int nr, nz;
    double dr, dz;

	// note I must use pointers so that I can sensibly instantiate these arrays when I know nx, ny and nz.

    Array3D<double> *cond_x, *distance_x;

    Array3D<int> *index_x;

    Array3D<double> *cond_y, *distance_y;
    Array3D<int> *index_y;

    Array3D<double> *cond_z, *distance_z;
    Array3D<int> *index_z;

	// constructor
    conductivity(int nr, int nz, double dr, double dz, int number_of_conductivity_parameters, double *conductivity_parameters);

	int nx, ny;
	double dx, dy;
	conductivity(int nx, int ny, int nz, double dx, double dy,  double dz, int number_of_conductivity_parameters, double *conductivity_parameters);

	// destructor
    ~conductivity();

    void GenerateConductivity(int number_of_conductivity_parameters, double *conductivity_parameters);

    int OutputIndex(int count, int ind1, int ind2, int option);
    double OutputCond(int count, int ind1, int ind2, int option);
    double OutputDist(int count, int ind1, int ind2, int option);


	void Write(int ind1, int ind2, int option);

};

#endif
