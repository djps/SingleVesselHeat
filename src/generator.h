/** 
\file generator.h

\brief Subroutine generating temperature from a given temperate profile. 

These are the functions that generate all the input data for you. In the future, 
I'll probably pass in file names and parameters but for now it's just a pretty 
straightforward setup.  

*/

#ifndef GENERATOR
#define GENERATOR

#include "array3d.h"
#include <cmath> 
#include <iostream>  
#include <fstream>

const double ABSORPTION_OF_BLOOD = 0.4;

/**
	\param heating Array3D double.
	\param restartFileName a constant character.
	\return The test results.
	\sa GenerateHeating(), GenerateTemperature()
*/   
void GenerateTemperature(Array3D<double> &temperature, std::string &restartFileName); 

/** \brief Subroutine generating temperature */
/**
	\param nx an integer argument.
    \param ny an integer argument.
	\param nz an integer argument.
	\param dx an integer argument.
    \param dy an integer argument.
	\param dz an integer argument.
	\param temperature Array3D double
	\param no_of_pars an integer argument.
	\param pars an double pointer.
	\return The test results.
	\sa GenerateHeating()
*/
void GenerateTemperature(int nx, int ny, int nz, double dx, double dy, double dz, Array3D<double> &temperature, \
int no_of_pars, double *pars);

/** \brief Subroutine generating heating profile */
/**
	\param nx an integer argument.
    \param ny an integer argument.
	\param nz an integer argument.
	\param dx an integer argument.
    \param dy an integer argument.
	\param dz an integer argument.
	\param heating Array3D double.
	\param heat_vel_pars an integer argument.
	\param heat_pars an double pointer.
    \param num_vel_pars an integer argument.
	\param vel_pars an double pointer.
	\param heatFileName a constant character
	\return The test results.
	\sa GenerateHeating()
*/
void GenerateHeating(int nx, int ny, int nz, double dx, double dy, double dz, Array3D<double> &heating, \
int num_heat_pars, double *heat_pars, int num_vel_pars, double *vel_pars, std::string &heatFileName);

#endif // ifndef GENERATOR


