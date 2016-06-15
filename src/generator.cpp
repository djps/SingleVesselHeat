#include "generator.h"
#include <cstdlib> 
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

//#define DEBUG_GENERATOR

//! subroutine which generates an initial temperature profile from a file
void GenerateTemperature(Array3D<double> &temperature, std::string& restartFileName)
{
	// note: This loads data from a file produced via save_to_file_binary. i.e. 
	// it is in the form of a binary file with three integers nx,ny,nz, 
	// followed by a massive array of doubles.
	std::ifstream tempFile;
	tempFile.open(restartFileName.c_str(), std::ios::in | std::ios::binary);
	if (!(tempFile.is_open()))
	{
		std::cerr << "Couldn't open temperature file " << std::endl;
		exit(3);
	}		
	int x, y,z;
	tempFile.read((char*) &x, sizeof(int));
	tempFile.read((char*) &y, sizeof(int));
	tempFile.read((char*) &z, sizeof(int));
	if ( (x!=temperature.nx) || (y!=temperature.ny) || (z!=temperature.nz) ) 
	{
		std::cerr << "Nonconformant restart file. Sorry " << std::endl;
		exit(4);
	}
	tempFile.read((char*) temperature.array, sizeof(double)*x*y*z); 
	// Just grab the data and stick it into the temperature array. 
	// Assumes the temperature array has been correctly assigned. 
	
}

//! subroutine which generates an initial temperature profile from analytical expression
void GenerateTemperature(int nx, int ny, int nz, double dx, double dy, double dz, Array3D<double> &temperature, \
int no_of_pars, double *pars)
{
	if (no_of_pars == 1)
	{ 
		double tempmax = pars[0]; // uniform temperature.
		double temp;
		// The following is needed if you want to have a Gaussian-type initial distribution. 
		int mid_x = nx/2; // central point of heating in x-direction 
		int mid_y = ny/2;  
		int mid_z = nz/2;  
		double breadth_x = 5 + nx/10; // the width of the heating peak. 
		double breadth_y = 5 + ny/10; // the width of the heating peak. 
		double breadth_z = 5 + nz/10; // the width of the heating peak. 
		for (int x = 0; x<nx; ++x)
		{
			for (int y = 0; y < ny; ++y)
			{
				for (int z = 0; z < nz; ++z)
				{
					temp = tempmax; // *exp(-((x-mid_x)/breadth_x)*((x-mid_x)/breadth_x) -((y-mid_y)/breadth_y)*((y-mid_y)/breadth_y) -((z-mid_z)/breadth_z)*((z-mid_z)/breadth_z)  )  ;// Gaussian heating. 
					temperature.Input_Value(x,y,z,temp);              
				}
			}
		}
	}              
	else if (no_of_pars == 7)
	{ 
		// This is a Gaussian
		double temp;
		double tempMax = pars[0];
		double breadth_x = pars[1]; 
		double breadth_y = pars[2]; 
		double breadth_z = pars[3];
		double mid_x = pars[4]; 
		double mid_y = pars[5]; 
		double mid_z = pars[6]; 
		for (int x=0; x<nx; ++x)
		{
			for (int y=0; y<ny; ++y)
			{
				for (int z=0; z<nz; ++z)
				{
					// Gaussian heating. 
					temp = tempMax*exp( -((x*dx-mid_x)/breadth_x)*((x*dx-mid_x)/breadth_x) \
					- ((y*dy-mid_y)/breadth_y)*((y*dy-mid_y)/breadth_y) \
					- ((z*dz-mid_z)/breadth_z)*((z*dz-mid_z)/breadth_z) );
					temperature.Input_Value(x,y,z,temp);              
				}
			}
		}
	} 	   
	else 
	{ 
	// default is zero.
	memset(temperature.array, 0, nx*ny*nz*sizeof(double));
	}               
}

//! subroutine which generates heat source from file
void GenerateHeating(int nx, int ny, int nz, double dx, double dy, double dz, Array3D<double> &heating, \
int num_heat_pars, double *heat_pars, int num_vel_pars, double *vel_pars, std::string &heatFileName)
{
	/*! 
	In this routine, I load up the heating profile from the file heatFileName 
	(stored as text file, with nx, ny, nz, maxValue as header, then x,y,z,heat, with x,y,z 
	integers >= 0 and < nx,ny,nz). Heating_parameters also has a value that is the maximum 
	heat deposition value to which we must rescale the acoustic profile. 
	We also reduce absorption over the blood vessel. Finally, heating_parameters will also 
	have 4 more parameters, two to determine when it starts switching on and how long until 
	it's  going, and another two to switch it off. But these aren't used in this routine. 
	So num_heat_pars should be 5, but we only use 1. 
	*/
	
  	if (num_heat_pars ==5)
	{
		// First, lets load up the file.
		int x, y, z;
		if (num_vel_pars !=4)
		{
			std::cerr << "Wrong number of velocity parameters " << std::endl;
			exit(3);
		}		 
		double xC = vel_pars[0];
		double yC = vel_pars[1];
		double radius = vel_pars[2];
		// No need to look at the maximum velocity here.
		
		std::ifstream heatFile;
		heatFile.open(heatFileName.c_str());
		if (!(heatFile.is_open()))
		{
			std::cerr << "Couldn't open acoustic file " << std::endl;
			exit(3);
		}		
		double heat, absorb;
		std::string discard;
		char comma;
		double givenHeatMax; // Value given in header of file
		heatFile >> x >> comma >> y >> comma>> z >> comma>>  givenHeatMax ;
		getline(heatFile,discard); // discard the rest of the line.
		if ((x!=nx)||(y!=ny)||(z!=nz))
		{
			std::cerr << "Wrong number of points in the heating file (or at least in the header)... " << std::endl;
			exit(3);
		}
		double heatScale = (givenHeatMax !=  0.0) ? heat_pars[0]/givenHeatMax : 0.0; // desired peak value of heating 		
		double tmpHeatMax = 0.0; // Keeps track of the actual max value that is input. 
		while (!(heatFile.eof()))
		{
			heatFile >> x >> comma >> y >> comma >> z >> comma>> heat ;
#ifdef DEBUG_GENERATOR 
std::cout <<   x << comma << y << comma << z << comma << heat << std::endl ;
#endif
			getline(heatFile,discard); // discard the rest of the line.
			if ( (x>=0) && (x<nx) && (y>=0) && (y<ny) && (z>=0) && (z<nz) )
			{//&&(heat>=0.0)){ // sense check
				if (heat > tmpHeatMax) tmpHeatMax = heat; 
				// Track the largest value of the heating. NB: If this value happens to be in blood, 
				// then it will never actually occur because blood has a lower absorption. However, 
				// this parameter sets the scale for the acoustic parameters, so it makes more sense to 
				// use the non-modified value. Essentially, we are looking for the peak pressure, and 
				// if this happens to be in blood, we don't care. 
				if (((x*dx-xC)*(x*dx-xC) + (y*dy-yC)*(y*dy-yC)) >= radius*radius)  absorb = 1.0; 
				// i.e. if we are outside the radius, then we're in tissue 
				else absorb = ABSORPTION_OF_BLOOD;		
				heating.Input_Value(x,y,z,absorb*heat*heatScale);
			} 
		}
		if (fabs(givenHeatMax-tmpHeatMax) > 1.e-8) 
		{
			std::cerr << "You lied to me! The actual maximum heating that was input was  "<< tmpHeatMax << " and you said it was " << givenHeatMax << std::endl;
			exit(3);
		}
		heatFile.close();		 
	}
	else 
	{ 
		// default is no heating.
		double heat = 0.0;
		//double absorb;
		for (int x=0; x<nx; ++x)
		{
			for (int y=0; y<ny; ++y)
			{
				for (int z=0; z<nz; ++z)
				{
					heating.Input_Value(x,y,z,heat); 
					// Since this is zero, don't need absorption. 
				}        
			}
		}
	}
	
#ifdef DEBUG_GENERATOR 
  std::cerr << "Generator " << heating.Test_NAN() << std::endl;
#endif

}
