/**
 @file SolveL1.h
 */

#ifndef SOLVE_L1
#define SOLVE_L1

// This file contains the implementation of the "Solve (1 - t/2 L_1)x =b" function.
// It's a long and complicated function, so I've stuck it in a separate header.
// All these headers will be included into run_heat_solver.cpp and that can be compiled into a single
// object file. No need to have too many object files running around.

void solveL1(int &nx, int &ny, int &nz, double &dx, double &dt, conductivity &kappa, conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c, double const&perf)
{
	// This routine steps through y and z, calculating the matrix (1+0.5*dt*L_1), where L_1 is the matrix of the operator
	//	-vd/dx + kappa d^2/dx^2.
    // Note that I already have a routine that calculates 1-0.5dt L_1. So the smartest thing to do is to copy this
	//	routine exactly but change dt to -dt. Fortunately, I introduce a temporary variable HalfTau, so if I just make
	//	this negative, it does all I want AND doesn't overwrite the input variable dt (which would be catastrophic).

    // To ensure that boundaries are correctly accounted for, we do separate sums to the different indices in the conductivity kappa.
    // First: Our working vectors for the tridiagonal matrices.


	double kap1, kap2;
	double diff1, diff2;
	double Halfdx_1 = 0.5/dx; // no point in dividing by 2dx every time!
    double dx_2     = 1.0/(dx*dx); // also useful.
    double HalfTau  = -0.5*dt; // NOTE the minus sign makes this calculate 1-dt/2 L_1.
    int boundary1, boundary2;
    double dist;
    for (int y = 0; y < ny; ++y)
    {
        for (int z = 0; z < nz; ++z)
        {
			// Do the matrix multiply for each value of y and z separately.
            // Recall that there are certain jumps in the conductivity. As such, we need to use a special resolution of spatial derivatives at these points. So instead of simply looping over x, we loop first over the indices between boundaries and then deal with the boundaries. Note that each internal boundary forces 2 points to be dealt with specially.

            boundary1 = -2; // Seems a bit funny but it's used because we deal with the two boundary points (at each internal boundary) at the end of the loop, not at the beginning or symmetrically.
            kap1 = kappa.OutputCond(0,y,z,1); // get the current value of kappa.
            kap2 = kap1; // necessary in case there are no discontinuities and we end up setting kap1=kap2 =undefined.
			diff1 = diff.OutputCond(0,y,z,1); // get the current value of diff.
			diff2 = diff1; // necessary in case there are no discontinuities and we end up setting kap1=kap2 =undefined.
			for (int counter = 0; counter < M-1; counter++) // M is defined in conductivity.h Basically, this loop loops over all discontinuities, so within each loop we need only worry about continuous bits plus the two boundary terms.
            {
                boundary2 = kappa.OutputIndex(counter,y,z,1);
                for (int x = boundary1+2; x < boundary2; x++)
				{
					// loop is over all elements from just beyond the first boundary.
					// Note that  kappa->OutputIndex(counter,y,z,1) yields the highest index before the boundary.
					// So actually, I'm now running over all the elements that will be done without needing to
					// worry about the boundary.

                    //vel = velocity.Output_Value(x,y,z);
					a[x] = HalfTau*diff1*dx_2;
					b[x] = 1.0 -HalfTau*2.0*diff1*dx_2 - HalfTau*perf/3.0;
					c[x] = HalfTau*diff1*dx_2;
                }
				// Now we need to do the boundary. The point before the boundary (indexed by boundary2) is definitely
				// within the domain of computation. The point after this might not be.  So test that we aren't trying to
				// write to a too-large index.
				// First, update boundary1:
                boundary1 = boundary2;
                if (boundary1 < (nx-1)) // then both sides of the boundary are defined.
                {
					// Set up the entries for the point on the left of the boundary.
                    //vel = velocity.Output_Value(boundary1,y,z);
					kap2  = kappa.OutputCond(counter+1,y,z,1);
					diff2 = diff.OutputCond(counter+1,y,z,1);
					dist  = kappa.OutputDist(counter,y,z,1);

                    double tempb = 1.0 + HalfTau*( - diff1*2.0*dx_2*( kap1*(1.0-dist) + kap2*(1.0+dist) )/(kap1*(1.0-dist*dist) + dist*kap2*(1.0+dist)) ) - HalfTau*perf/3.0;

                    double tempc = HalfTau*( diff1*2.0*dx_2*kap2/(kap1*(1.0-dist*dist) + kap2*dist*(1.0+dist)) );

					a[boundary1] = HalfTau*( diff1*2.0*dx_2/(1.0+dist) );

					b[boundary1] = ( (kap1 > 0.0) || (kap2 > 0.0)) ? tempb : 1.0 - HalfTau*perf/3.0;

					c[boundary1] = ( (kap1 > 0.0) || (kap2 > 0.0)) ? tempc : 0.0;

					// Now for the right of the boundary.
					c[boundary1+1] = HalfTau*(diff2*2.0*dx_2/(2.0-dist));

                    tempb = 1.0 + HalfTau*(diff2*2.0*dx_2*(-kap1*(2.0-dist) - kap2*dist ) /((2.0-dist)*(kap1*(1.0-dist) +  kap2*dist))) - HalfTau*perf/3.0;
                    b[boundary1+1] = ( (kap1 > 0.0) || (kap2 > 0.0) )  ?  tempb :  1.0 - HalfTau*perf/3.0;

                    double tempa = HalfTau*(diff2*2.0*dx_2*kap1/((2.0-dist)*(kap1*(1.0- dist) + kap2*dist)));
					a[boundary1+1] = ( (kap1 > 0.0)|| (kap2 > 0.0) )  ?  tempa : 0.0;
                }
                else if (boundary1==(nx-1)) // In this case, we're at the border of our domain. No need for fiddling about,
				// just do it as usual. Implicit here is that the BC is vanishing T on boundary.
                {
                    //vel = velocity.Output_Value(boundary1,y,z);
					a[boundary1] = HalfTau*(diff1*dx_2);
					b[boundary1] = 1.0 - HalfTau*2.0*diff1*dx_2 - HalfTau*perf/3.0;
                    c[boundary1] = 0.0; // This is outside the realm of the tridiagonal matrix.
                }
                else
                {
                    std::cout << "Trouble here: the conductivity extends beyond the correct domain. " << std::endl;
                }

				kap1  = kap2; // finish by updating kappa.
				diff1 = diff2; // finish by updating kappa.
			}

#ifdef ADVECTIVE_BOUNDARIES
			// !!!!!!!!!!!!!!!!!!!!!!!!!Recall that HalfTau has a minus sign in it. This was convenient. Not sure if it still is :-)
            // First, the right boundary. Use the advective-dominated BC iff v > 2*kap/dx.
            a[nx-1] = HalfTau*diff1*dx_2 ;
            b[nx-1] = 1.0 - HalfTau*2.0*diff1*dx_2 - HalfTau*perf/3.0;
            // Now the left boundary
			kap1  = kappa.OutputCond(0,y,z,1);
			diff1 = diff.OutputCond(0,y,z,1);
			c[0] = HalfTau*diff1*dx_2;
			b[0] = 1.0 - HalfTau*2.0*diff1*dx_2 - HalfTau*perf/3.0;

            // The following are essential elements of our representation of the tridiag matrix.
            a[0]    = 0.0;
            c[nx-1] = 0.0;

#elseifdef ADVECTIVE_BOUNDARIES_2
			// I have no idea whether these will be better, but they can't be worse.
			// Basically, we just solve the advection equation using a first order scheme on the boundary.
			// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			// !!!!!!!!!!!!!!!!!!!!! Halftau is negative not positive !!!!!!!!!!!!!!!!!!
			// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			// First, the right boundary.
			a[nx-1] = 0.0;
            b[nx-1] = 1.0 - HalfTau*(perf+RADIATION_FROM_BOUNDARY)/3.0;

			// Now the left boundary
            kap1 = kappa.OutputCond(0,y,z,1);
			c[0] = 0.0;
            b[0] = 1.0 - HalfTau*(perf+RADIATION_FROM_BOUNDARY)/3.0;

            // The following are essential elements of our representation of the tridiag matrix.
            a[0]    = 0.0;
            c[nx-1] = 0.0;

#else
            // If you want other BCs, here's the place for them.

            // I like to have zero derivative boundary conditions. It allows heat flow, which seems pretty reasonable.
            // First, the right boundary
            a[nx-1] = HalfTau*(BOUNDARY_TERM*diff1*dx_2);
			b[nx-1] = 1.0 - HalfTau*2.0*diff1*dx_2 - HalfTau*perf/3.0;

            // Now the left boundary
			diff1 = diff.OutputCond(0,y,z,1);
			diff1 = diff.OutputCond(0,y,z,1);
			c[0] = HalfTau*(BOUNDARY_TERM*diff1*dx_2);
			b[0] = 1.0 - HalfTau*2.0*diff1*dx_2 - HalfTau*perf/3.0;

            // These are essential elements of our representation of the tridiag matrix.
            a[0]    = 0.0;
            c[nx-1] = 0.0;
#endif

            // Now we have the matrix, invert it!
            temperature.TriDiagSolve(a, b, c, 1, nx, y, z); // the 1 refers to multiplication along the x index.
        }
    }
}

#endif
