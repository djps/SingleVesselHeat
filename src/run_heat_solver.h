// This is a general runner routine.
#ifndef RUN_HEAT_SOLVER
#define RUN_HEAT_SOLVER

#include <cmath>

/*
I'm encapsulating all the running steps in this file.
*/

//#define DEBUG_RUNNER
//#define DEBUG_RUNNER3
//#define DEBUG_RUNNER4

// If BOUNDARY_TERM is 1, then we have zero temperature BCs. If it is 2, then we have insulating BCs.
#define BOUNDARY_TERM  1.0

// MAX_TEMPERATURE is the maximum allowable temperature at any point. Exceeding this leads to a cutoff of the heating.
#define MAX_TEMPERATURE  50.0


//#define ADVECTIVE_BOUNDARIES
// If you have defined this, then slightly more stable BC's are used. Woohoo! They basically use T = 0
// on boundary if there is no flow, inflow, or slow outflow. Otherwise they use an advection-dominated BC.

// NOTE Turned off 25.03.2010. I want insulating BCs as far as the conductivity bit is concerned.

//#define ADVECTIVE_BOUNDARIES_2
// This is a potential improvement, that forces the equations of motion to be T_t + v T_x at the boundaries
// and uses backwards or forward differences to avoid the outside world. Not accurate, but should be pretty good.
// I think. Worth a try.

//#define RADIATION_FROM_BOUNDARY 0.1
// It's a good idea to have extra cooling on the boundary. We've switched off all conduction, which means
// that for a nonflowing boundary, anything that accumulates on the boundary stays forever.

#include "array3d.cpp"
#include "conductivity.h"
#include "generator.h"

#include "velocity.h"

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

// Comment this out in final code.
//#define DEBUG_HEAT_SOLVER


void run_heat_solver(int &nmax, int &nx,int &ny,int &nz, double &dx, double &dy, double &dz, double &dt,double &time, int &iterations, conductivity &kappa,conductivity &diff, Velocity &velocity, Array3D<double> &temperature,Array3D<double> &heating, int num_heat_pars, double *heat_pars, Array3D<double> &therm_dose, double *a, double *b, double *c, double *d, double *centralpoint, int &numoftimesteps, double *timevec, double const&perf);

void multL1(int &nx, int &ny, int &nz, double &dx, double &dt, conductivity &kappa, conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c, double const &perf);

void multL2(int &nx, int &ny, int &nz, double &dy, double &dt, conductivity &kappa, conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c, double const &perf);

void multL3(int &nx, int &ny, int &nz, double &dz, double &dt, conductivity &kappa, conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c, double const &perf);

void solveL1(int &nx, int &ny, int &nz, double &dx, double &dt, conductivity &kappa, conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c, double const &perf);

void solveL2(int &nx, int &ny, int &nz, double &dy, double &dt, conductivity &kappa, conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c, double const &perf);

void solveL3(int &nx, int &ny, int &nz, double &dz, double &dt, conductivity &kappa, conductivity &diff, Array3D<double> &temperature, double *a, double *b, double *c, double const &perf);

void  AddToThermalDose(Array3D<double> const& temperature, Array3D<double> &therm_dose, double step);

void  AddToIntegratedTemp(Array3D<double> const& temperature, Array3D<double> &int_temp, double step);

void  DamagedIfExceeds56Degrees(Array3D<double> const& temperature, Array3D<double> &damaged_tissue, double step);

void  SwitchOffIfTemperatureExceedsThreshold(Array3D<double> const& temperature, Array3D<double> &heating, double cutoff, double time);

//void  SwitchOffIfTemperatureExceedsThreshold(Array3D<double> const& temperature, int numHeatPars, double *heatPars, double cutoff, double time)

#endif
