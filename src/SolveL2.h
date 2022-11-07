#ifndef SOLVE_L2
#define SOLVE_L2


// This file contains the implementation of the "Solve by (1 - t/2 L_2)x=b" function. 
// It's a long and complicated function, so I've stuck it in a separate header. 
// All these headers will be included into run_heat_solver.cpp and that can be compiled into a 
// single object file. No need to have too many object files running around. 

void solveL2(int &nx, int &ny, int &nz, double &dy, double &dt,  conductivity &kappa,  conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c,  double const&perf)
{
	// This routine steps through y and z, calculating the matrix (1-0.5*dt*L_2), where L_2 is the matrix of the operator 
	// -vd/dy + kappa d^2/dy^2. 
    // Note that I already have a routine that calculates 1-0.5dt L_1. So the smartest thing to do is to copy this 
	// routine exactly but change dt to -dt. Fortunately, I introduce a temporary variable HalfTau, so if I just 
	// make this negative, it does all I want AND doesn't overwrite the input variable dt (which would be catastrophic).
    // To ensure that boundaries are correctly accounted for, we do separate sums to the different indices in the 
	// conductivity kappa.    
	
    // First: Our working vectors for the tridiagonal matrices. 
	
    /*
	double *a = new double[ny];  
    double *b = new double[ny];
    double *c = new double[ny];
    double *vel = new double[ny];
	*/
    
	double kap1,kap2;
	double diff1,diff2;
	double Halfdy_1 = 0.5/dy; // no point in dividing by 2dx every time!
    double dy_2 = 1.0/(dy*dy); // also useful. 
    double HalfTau = - 0.5*dt;  // Note the minus sign.
    int boundary1,boundary2;
    double dist;
    for (int x = 0; x<nx; ++x)
    {
        for (int z = 0; z<nz; ++z)
        { 
			// Do the matrix multiply for each value of y and z separately. 
            // Recall that there are certain jumps in the conductivity. As such, we need to use a 
			// special resolution of spatial derivatives at these points. So instead of simply looping over 
			// x, we loop first over the indices between boundaries and then deal with the boundaries. 
			// Note that each internal boundary forces 2 points to be dealt with specially. 
            boundary1 = -2; // Seems a bit funny but it's used because we deal with the two boundary points 
			// (at each internal boundary) at the end of the loop, not at the beginning or symmetrically. 
            kap1 = kappa.OutputCond(0, x ,z,2); // get the current value of kappa.              
            kap2 = kap1; // necessary in case there are no discontinuities and we end up setting kap1=kap2=undefined. 
			diff1 = diff.OutputCond(0,x,z,2); // get the current value of diff.              
			diff2 = diff1; // necessary in case there are no discontinuities and we end up setting kap1=kap2=undefined. 
			for(int counter = 0; counter < M-1;++counter) // M is defined in conductivity.h Basically, this loop 
			// loops over all discontinuities, so within each loop we need only worry about continuous bits plus 
			// the two boundary terms.    
            {
                boundary2 = kappa.OutputIndex( counter,x, z,2);
                for (int y =boundary1+2; y<boundary2; ++y)
				{
					//loop is over all elements from just beyond the first boundary. Note that 
					// kappa->OutputIndex(counter,y,z,1) yields the highest index before the boundary. 
					// So actually, I'm now running over all the elements that will be done without 
					// needing to worry about the boundary.
					//vel = velocity.Output_Value(x,y,z);
					a[y] = HalfTau*( diff1*dy_2);
					b[y] = 1.0 -HalfTau*2.0*diff1*dy_2 - HalfTau*perf*0.333333333333;
					c[y] = HalfTau*( diff1*dy_2);
                }
				// Now we need to do the boundary. The point before the boundary (indexed by boundary2) is definitely within the domain of computation. The point after this might not be.  So test that we aren't trying to write to a too-large index. 
				// First, update boundary1: 
                boundary1 = boundary2;
                if (boundary1<(ny-1)) // then both sides of the boundary are defined. 
                {
					// Set up the entries for the point on the left of the boundary. 
                    //vel = velocity.Output_Value(x,boundary1, z);            
					kap2 = kappa.OutputCond( counter+1,x, z,2);
					diff2 = diff.OutputCond( counter+1,x, z,2);
					dist = kappa.OutputDist( counter,x, z,2);
					a[boundary1] = HalfTau*(diff1*2.0*dy_2/(1.0+dist));
					b[boundary1] = ((kap1>0.0)||(kap2>0.0)) ?  1.0 +HalfTau*(- diff1*2.0*dy_2*(kap1*(1.0-dist)+kap2*(1.0+dist))/(kap1*(1.0-dist*dist) + dist*kap2*(1.0+dist))) - HalfTau*perf*0.333333333333 : 1.0 - HalfTau*perf*0.333333333333;
					c[boundary1] = ((kap1>0.0)||(kap2>0.0)) ?  HalfTau*( diff1*2.0*dy_2*kap2/(kap1*(1.0-dist*dist) + kap2*dist*(1.0+dist))) : 0.0;
                   
                   // Now for the right of the boundary. 
					c[boundary1+1] = HalfTau*(diff2*2.0*dy_2/(2.0-dist));                 
					b[boundary1+1] = ((kap1>0.0)||(kap2>0.0)) ?  1.0 +HalfTau*(diff2*2.0*dy_2*(-kap1*(2.0-dist) - kap2*dist ) /((2.0-dist)*(kap1*(1.0-dist) +  kap2*dist))) - HalfTau*perf*0.333333333333  :  1.0 - HalfTau*perf*0.333333333333;
					a[boundary1+1] = ((kap1>0.0)||(kap2>0.0)) ?  HalfTau*(diff2*2.0*dy_2*kap1/((2.0-dist)*(kap1*(1.0- dist) + kap2*dist))) : 0.0;
                }
                else if (boundary1==(ny-1)) // In this case, we're at the border of our domain. No need for fiddling about, just do it as usual. Implicit here is that the BC is vanishing T on boundary. 
                {
                    //vel = velocity.Output_Value(x,boundary1, z);
					a[boundary1] = HalfTau*(diff1*dy_2);
					b[boundary1] = 1.0 -HalfTau*2.0*diff1*dy_2 - HalfTau*perf*0.333333333333;
                    c[boundary1] = 0.0; // This is outside the realm of the tridiagonal matrix. 
                }
                else std::cout << "Trouble here: the conductivity extends beyond the correct domain. " << std::endl;

				kap1 = kap2; // finish by updating kappa. 
				diff1 = diff2; // finish by updating kappa. 
            
            }
			
#ifdef ADVECTIVE_BOUNDARIES
				// Recall that HalfTau has a minus sign in it. This was convenient. Not sure if it still is.
				// First, the right boundary. Use the advective-dominated BC iff v > 2*kap/dx. 
				a[ny-1] = HalfTau*diff1*dy_2;
				b[ny-1] = 1.0 -HalfTau*2.0*diff1*dy_2 - HalfTau*perf*0.333333333333;
				kap1 =kappa.OutputCond(0,x,z,2);
				diff1 =diff.OutputCond(0,x,z,2);
				c[0] = HalfTau*diff1*dy_2 ;
				b[0] = 1.0 -HalfTau*2.0*diff1*dy_2 - HalfTau*perf*0.333333333333;
				   
				// The following are essential elements of our representation of the tridiag matrix. 
				a[0] = 0.0;
				c[ny-1] = 0.0;
               
               
#elseifdef ADVECTIVE_BOUNDARIES_2            
				// I have no idea whether these will be better, but they can't be worse. 
				// Basically, we just solve the advection equation using a first order scheme on the boundary. 
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// !!!!!!!!!!!!!!!!!!!!! Halftau is negative not positive !!!!!!!!!!!!!!!!!!
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// First, the right boundary.  
                   
				a[ny-1] = 0.0; 
                b[ny-1] = 1.0 - HalfTau*(perf+RADIATION_FROM_BOUNDARY)*0.333333333333;
                
				// Now the left boundary
                kap1 =kappa.OutputCond(0,x,z,2);
				c[0] = 0.0; 
                b[0] = 1.0 - HalfTau*(perf+RADIATION_FROM_BOUNDARY)*0.333333333333;
                
				// The following are essential elements of our representation of the tridiag matrix. 
				a[0] = 0.0;
				c[ny-1] = 0.0;
	       
#else                                      
				// If you want other BCs, here's the place for them. 
                
                // I like to have zero derivative boundary conditions. It allows heat flow, which seems pretty reasonable.
                // First, the right boundary  
				a[ny-1] = HalfTau*(BOUNDARY_TERM*diff1*dy_2);
				b[ny-1] = 1.0 -HalfTau*2.0*diff1*dy_2 - HalfTau*perf*0.333333333333;
                
                // Now the left boundary
				kap1 =kappa.OutputCond(0,x,z,2);
				diff1 =diff.OutputCond(0,x,z,2);
				c[0] = HalfTau*(BOUNDARY_TERM*diff1*dy_2);
				b[0] = 1.0 -HalfTau*2.0*diff1*dy_2 - HalfTau*perf*0.333333333333;
                
				// These are essential elements of our representation of the tridiag matrix. 
				a[0]=0.0;
				c[ny-1]=0.0;
#endif
            // Now we have the matrix, invert it!
            temperature.TriDiagSolve( a,  b,  c, 2, ny, x, z); // the 2 refers to multiplication along the xy index.
        }
    }
	
/*    
	delete[] a;
    delete[] b;
    delete[] c;
    delete[] vel;
	*/
	
}

#endif

