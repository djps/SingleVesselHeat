/***************************************************************************
 *   Copyright (C) 2009 by Simon Woodford   *
 *   simon.woodford@icr.ac.uk   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
 
#ifndef MULT_L3
#define MULT_L3

#include <float.h>

//#define DEBUG_NAN

// This file contains the implementation of the "Multiply by 1 + t/2 L_3" function. 
// It's a long and complicated function, so I've stuck it in a separate header. 
// All these headers will be included into run_heat_solver.cpp and that can be compiled 
// into a single object file. No need to have too many object files running around. 

//EDIT on 02.03.2010: put in explicit checking for kappa = 0.

void multL3(int &nx, int &ny, int &nz, double &dz, double &dt, conductivity &kappa, conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c, double const &perf)
{
	// This routine steps through y and z, calculating the matrix (1+0.5*dt*L_1), 
	// where L_1 is the matrix of the operator -vd/dx + kappa d^2/dx^2. 
    // To ensure that boundaries are correctly accounted for, we do separate sums 
	// to the different indices in the conductivity kappa.
    // First: Our working vectors for the tridiagonal matrices. 
	
    /*
	double *a = new double[nz];  
    double *b = new double[nz];
    double *c = new double[nz];
    double *vel = new double[nz];
	*/
	
    double kap1,kap2;
	double diff1,diff2;
	char comma;
	
#ifdef DEBUG_RUNNER3
    int CheckIfIEnterLoop;
#endif     

    double Halfdz_1 = 0.5/dz; // no point in dividing by 2dx every time!
    double dz_2 = 1.0/(dz*dz); // also useful. 
    double HalfTau = 0.5*dt;
    int boundary1,boundary2;
    double dist;
    for (int y = 0; y<ny; ++y)
    {
        for (int x = 0; x<nx; ++x)
        { 
			// Do the matrix multiply for each value of y and z separately. 
            // Recall that there are certain jumps in the conductivity. As such, we need to use a 
			// special resolution of spatial derivatives at these points. So instead of simply 
			// looping over x, we loop first over the indices between boundaries and then deal 
			// with the boundaries. Note that each internal boundary forces 2 points to be dealt 
			// with specially. 
            boundary1 = -2; // Seems a bit funny but it's used because we deal with the two boundary points (at each internal boundary) at the end of the loop, not at the beginning or symmetrically. 
            kap1 = kappa.OutputCond(0,x,y,3); // get the current value of kappa.              
            kap2 = kap1; // necessary in case there are no discontinuities and we end up setting kap1=kap2 =undefined. 
			diff1 = diff.OutputCond(0,x,y,3); // get the current value of diff.              
			diff2 = diff1; // necessary in case there are no discontinuities and we end up setting kap1=kap2 =undefined. 
			for(int counter = 0; counter < (M-1);++counter) // M is defined in conductivity.h Basically, this loop loops over all discontinuities, so within each loop we need only worry about continuous bits plus the two boundary terms.    
            {

#ifdef DEBUG_RUNNER3
                kappa.Write(x,y,3);
#endif     

                boundary2 = kappa.OutputIndex( counter,x,y,3);

#ifdef DEBUG_RUNNER3
                std::cout << counter << " " << boundary1 << "  " << boundary2 << " " << " Boundaries  " << std::endl;
#endif     

#ifdef DEBUG_RUNNER3
                CheckIfIEnterLoop = 0;
#endif
    
                for (int z =boundary1+2; z<boundary2; ++z)  
				// loop is over all elements from just beyond the first boundary. 
				// Note that kappa->OutputIndex(counter,y,z,1) yields the highest index before 
				// the boundary. So actually, I'm now running over all the elements that will 
				// be done without needing to worry about the boundary.
                {
                    //vel = velocity.Output_Value(x,y,z);
                    a[z] = HalfTau*(diff1*dz_2);
				    b[z] = 1.0 -HalfTau*2.0*diff1*dz_2 - HalfTau*perf*0.333333333333;
                    c[z] = HalfTau*( diff1*dz_2);
					
#ifdef DEBUG_RUNNER3
                    CheckIfIEnterLoop = 1;
                    //std::cout << "Inside Loop " << a[z] << " " <<  b[z] << " " <<  c[z] << std::endl;
#endif     
                }
				
#ifdef DEBUG_RUNNER3
                // std::cout << x << "  " << y << " " << CheckIfIEnterLoop <<"  " << counter << " Mult3, finished loop.  " << std::endl;
#endif
     
                // Now we need to do the boundary. The point before the boundary (indexed by boundary2) 
			    // is definitely within the domain of computation. The point after this might not be.  
			    // So test that we aren't trying to write to a too-large index. 
                // First, update boundary1: 
                boundary1 = boundary2;
                if (boundary1<(nz-1)) 
				// then both sides of the boundary are defined. 
                {
				
#ifdef DEBUG_RUNNER3                        
                    std::cout << "Discontinuity detected! " << std::endl;
#endif     

                    // Set up the entries for the point on the left of the boundary. 
                    //vel = velocity.Output_Value(x,y,boundary1);            
                    kap2 = kappa.OutputCond( counter+1,x,y,3);
                    diff2 = diff.OutputCond( counter+1,x,y,3);
					dist = kappa.OutputDist( counter,x,y,3);
					a[boundary1] = HalfTau*(diff1*2.0*dz_2/(1.0+dist));
					b[boundary1] = ((kap1>0.0)||(kap2>0.0)) ?  1.0 +HalfTau*(- diff1*2.0*dz_2*(kap1*(1.0-dist)+kap2*(1.0+dist))/(kap1*(1.0-dist*dist) + dist*kap2*(1.0+dist))) - HalfTau*perf*0.333333333333 : 1.0 - HalfTau*perf*0.333333333333;
					c[boundary1] = ((kap1>0.0)||(kap2>0.0)) ?  HalfTau*( diff1*2.0*dz_2*kap2/(kap1*(1.0-dist*dist) + kap2*dist*(1.0+dist))) : 0.0;
                   
					// Now for the right of the boundary. 
					c[boundary1+1] = HalfTau*(diff2*2.0*dz_2/(2.0-dist));                 
					b[boundary1+1] = ((kap1>0.0)||(kap2>0.0)) ?  1.0 +HalfTau*( diff2*2.0*dz_2*(-kap1*(2.0-dist) - kap2*dist ) /((2.0-dist)*(kap1*(1.0-dist) +  kap2*dist))) - HalfTau*perf*0.333333333333  :  1.0 - HalfTau*perf*0.333333333333;
					a[boundary1+1] = ((kap1>0.0)||(kap2>0.0)) ?  HalfTau*(diff2*2.0*dz_2*kap1/((2.0-dist)*(kap1*(1.0- dist) + kap2*dist))) : 0.0;
                }
                else if (boundary1==(nz-1)) 
				// In this case, we're at the border of our domain. No need for fiddling about, 
				// just do it as usual. Implicit here is that the BC is vanishing T on boundary. 
                {
				
#ifdef DEBUG_RUNNER3
                    std::cout << "SHOULD SEE THIS SEVERAL TIMES!!" << std::endl;
#endif     

                    //vel = velocity.Output_Value(x,y,boundary1);
                    a[boundary1] = HalfTau*( diff1*dz_2);
                    b[boundary1] = 1.0 -HalfTau*2.0*diff1*dz_2 - HalfTau*perf*0.333333333333;
                    c[boundary1] = 0.0; // This is outside the realm of the tridiagonal matrix. 
#ifdef DEBUG_RUNNER3
                    std::cout << "end values " <<  a[boundary1] << "  "<<  b[boundary1] << "  "<<  c[boundary1] << "  "  << std::endl;
#endif     
                }
                else std::cerr << "Trouble here: the conductivity extends beyond the correct domain. " << std::endl;
                                
#ifdef DEBUG_RUNNER3
                std::cout << x << "  " << y << " " << " Mult3, done boundary.  " << std::endl;
#endif     

                kap1 = kap2; // finish by updating kappa. 
                diff1 = diff2; // finish by updating kappa. 

#ifdef DEBUG_RUNNER3   
                //for (int p = 0; p<nz; ++p) std::cout << "Inside Loop " << a[p] << " " <<  b[p] << " " <<  c[p] << std::endl;
#endif     
            }
                
            
#ifdef ADVECTIVE_BOUNDARIES 
                // First, the right boundary. Use the advective-dominated BC iff v > 2*kap/dx.
                a[nz-1] = HalfTau*diff1*dz_2 ;
                b[nz-1] = 1.0 -HalfTau*2.0*diff1*dz_2 - HalfTau*perf*0.333333333333;
                // Now the left boundary
				kap1 = kappa.OutputCond(0,x,y,3);
				diff1 = diff.OutputCond(0,x,y,3);
				c[0] = HalfTau*diff1*dz_2 ;
				b[0] = 1.0 -HalfTau*2.0*diff1*dz_2 - HalfTau*perf*0.333333333333;
				// The following are essential elements of our representation of the tridiag matrix.
				a[0] = 0.0;
				c[nz-1] = 0.0;
	       
	       
#elseifdef ADVECTIVE_BOUNDARIES_2            
				// I have no idea whether these will be better, but they can't be worse. 
				// Basically, we just solve the advection equation using a first order scheme on the boundary. 
                // First, the right boundary.  
                   
				a[nz-1] = 0.0; 
                b[nz-1] = 1.0 - HalfTau*(perf+RADIATION_FROM_BOUNDARY)*0.333333333333 ;
                
				// Now the left boundary
                kap1 =kappa.OutputCond(0,x,y,3);
				c[0] = 0.0;
                b[0] = 1.0 - HalfTau*(perf+RADIATION_FROM_BOUNDARY)*0.333333333333  ;
                
				// The following are essential elements of our representation of the tridiag matrix. 
				a[0] = 0.0;
				c[nz-1] = 0.0;
	       

            
            
#else            
            // If you want other BCs, here's the place for them. 
                
            // First, the right boundary  
            a[nz-1] = HalfTau*(BOUNDARY_TERM*diff1*dz_2);
	   		b[nz-1] = 1.0 -HalfTau*2.0*diff1*dz_2 - HalfTau*perf*0.333333333333;
                
            // Now the left boundary
			kap1 = kappa.OutputCond(0,x,y,3);
			diff1 = diff.OutputCond(0,x,y,3);
			c[0] = HalfTau*(BOUNDARY_TERM*diff1*dz_2);
			b[0] = 1.0 -HalfTau*2.0*diff1*dz_2 - HalfTau*perf*0.333333333333;
                
            // These are essential elements of our representation of the tridiag matrix. 
            a[0] = 0.0;
            c[nz-1] = 0.0;
#endif
               
#ifdef DEBUG_RUNNER4
            for (int ind = 0; ind < nz; ind++ )  std::cout << "a,b,c  " << ind << ", " <<  a[ind] << "  "<<  b[ind] << "  "<<  c[ind] << "  "  << std::endl;
            std::cin >> comma ;
			for (int ind = 0; ind < nz; ind++ ) std::cout << "x,y,temperature before " << x << ", " << y << ", " << temperature.Output_Value(x,y,ind) << std::endl;   
#endif     
           
            // Now we have the matrix, do the multiplication!
			   
#ifdef DEBUG_NAN
			for (int p = 0 ; p < nz ; p++){
	            if(isnan(a[p])) std::cerr << "Problem in a  at " << p <<std::endl; 
				if(isnan(b[p])) std::cerr << "Problem in b  at " << p <<std::endl; 
				if(isnan(c[p])) std::cerr << "Problem in c  at " << p <<std::endl; 
			}
#endif               
               
            temperature.TriDiagMultiply( a,  b,  c, 3, nz, x, y); // the 3 refers to multiplication along the z index.
			
#ifdef DEBUG_RUNNER4
           for (int ind = 0; ind < nz; ind++ ) std::cout << "x,y,temperature after " << x << ", " << y << ", " << temperature.Output_Value(x,y,ind) << std::endl;   
#endif
     
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
