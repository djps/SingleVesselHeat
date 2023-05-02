/**
 @file fdheat.cpp
 */

/**
 @brief finite difference heat equation

 The program:

 We want to solve the heat equation when a blood vessel is present. So we have to deal with motion and conduction, with different conductivities.

 Motion is tricky, simply because of scales. The HIFU focus is small, so we need good spatial resolution. But the blood travels rapidly, so we need good temporal resolution.

 Enter brilliant plan. Instead of evolving the velocity at the same time as the conduction, we split-step it. This means that we can evolve the velocity part exactly. Now, this is of zero value if the velocity bit is very interesting, but I'm almost certain that very fast moving blood is unlikely to heat up much. So it's not going to make a large error having a large timestep, because basically the blood remains cold. However, it's a pain to do simulations on a split-domain, so it is much better to just have an explicit vessel that happens to be cold than a forced boundary condition.

 To simplify things, I use  a structured grid that is rectilinear with equal spacing, the blood vessel is in the z direction, and I also split-step the conduction bit so that I only need to solve tridiagonal systems at each step. I have not considered using full Crank-Nicolson, so I can't say how much faster this "ADI" is. But it seems pretty quick and stable.

 Setting up the conductivity is a little tricky, so I've changed it a bit from how it used to be. See what you think! Basically, I'm forcing the vessel to be in the z direction, which simplifies everything else.

 I'm also going to modify the heating bit. Now, I load an acoustic map from file, apply an absorption criterion (basically 1 for tissue, 0.4 for blood) and that's my heater. However, I also insert parameters that allow it to ramp on (linear switch on) and ramp off. This means that I'm not regenerating the heat profile at every timestep, which is quite good at saving time. I think...

 Most of the code is self-explanatory (how many times have you heard that?) Basically, I set up the initial data then run the solver. At each time-step, we multiply the temperature array by 3 matrices (tridiagonal, each one is multiplied along a different array index of the temperature corresponding to the x, y and z directions), then invert 3 matrix multiplications.

 Since this involves choosing a definite x,y,z order and I don't know the stability of choosing such an order, I do 6 different steps corresponding to each different permutation of x,y,z. I think it's essential to keep the multiply in the reverse order of the matrix solve, so there are only 6 permutations.

*/


//#define DEBUG_CONDUCTIVITY
//#define DEBUG_MAIN

//#define RUN_TEST_SUITE

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <string>

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include <iomanip>
#include <cstdlib>

#include "HeatConfig.h"

#include "array3d.cpp"
#include "conductivity.h"
#include "generator.h"
#include "run_heat_solver.h"
#include "save_to_file.h"

#include "heatReader.h"

#include "xtime.hpp"

using namespace std;


/*
   PARAMETER FILE LOOKS LIKE:
   nx, ny, nz
   iterations_per_save, number_of_saves, save_option, save_option2
   (save_option tells whether to save temperature profiles, save_option2 is about thermal dose.)
   dx, dy, dz, dt
   resultsDirectoryName
   perfusion_parameter
   num_heat_pars
   heat_pars
   heatFileName
   num_vel_pars
   vel_pars
   num_temp_pars
   temp_pars
   temperatureFileName, preceded by a 1 if we use the file, by a 0 if not.
   num_conductivity_pars
   cond_pars
   num_diffusivity_pars
   diff_pars
*/
int main(int argc, char *argv[])
{

  TTime  xtimer;
  xtimer.start();

  bool verbose = true;

  // First things first, declare and define the parameters.
  int nx;
  int ny;
  int nz;
  int number_of_iterations_per_save;
  int number_of_saves;
  int save_option;
  int save_option2;

  // note that nx is the number of nodes in the x direction, not the number of spaces.
  double dx;
  double dy;
  double dz;
  double dt;

  // DS added 16/06/2016
  fprintf(stdout,"%s Version %d.%d\n", argv[0], Heat_VERSION_MAJOR, Heat_VERSION_MINOR);

  string parameter_file_name, datafilename;
  if (argc > 1)
  {
    parameter_file_name = argv[1];
  }
  else
  {
    parameter_file_name = "DefaultParameterFile.dat";
  }

  // const char *filename = "example.xml";
  // if (argc > 1)
  // {
  //   filename = argv[1];
  //   cout << "data read from file: " << filename << endl;
  // }

  HeatInputData input;

  if (verbose){ fprintf(stdout,"\tReading data from %s\n", parameter_file_name.c_str() ); }

  input = traverse_xml(parameter_file_name);

  string jsonfilename = "heat.json";
  HeatInputData jsoninput;
  jsoninput = traverse_json(jsonfilename);

/*
  parameter_file_name = "Sim2.dat";

  cout << parameter_file_name.c_str() << endl;

  ifstream parameterinput;
  parameterinput.open(parameter_file_name.c_str());
  char Comma;
  string DiscardMe;
  parameterinput >> nx >> Comma >> ny >> Comma >> nz;
  if (verbose){ cout << "nx " << nx << " ny " << ny << " nz " << nz << endl; }
  getline(parameterinput,DiscardMe);
  parameterinput >> number_of_iterations_per_save >> Comma >> number_of_saves >> Comma >> save_option >> Comma >> save_option2; // 19.02.2010, introduce save_option. If this is nonzero then we save temperature profiles. If it is zero, we don't. Similarly, if save_option2 is nonzero, we save the thermal dose.
  double *centralpoint = new double [number_of_iterations_per_save*number_of_saves];
  // This keeps track of the central value of the temperature.
  double *timevector = new double [number_of_iterations_per_save*number_of_saves];
  // This keeps track of the time at each iteration.
  if (verbose){ cout << "iters  " << number_of_iterations_per_save << " saves " << number_of_saves << endl; }

  getline(parameterinput,DiscardMe);
  parameterinput >> dx >> Comma >> dy >> Comma >> dz >> Comma >> dt;
  if (verbose){ cout << "dx,dy,dz,dt " << dx << " " << dy << " " << dz << " " << dt << endl; }

  getline(parameterinput,DiscardMe);
  parameterinput >> datafilename;
  if (verbose){ cout << datafilename << endl; } // sims2/

  //cin >> Comma ;
  getline(parameterinput,DiscardMe);
  if (verbose){ cout << DiscardMe << endl; }
  double perfusion_parameter;
  parameterinput >> perfusion_parameter;
  if (verbose){ cout << "perf " << perfusion_parameter << endl; } //0

  getline(parameterinput,DiscardMe);
  if (verbose){ cout << "Line after perf " << DiscardMe << endl; }

  int number_of_heating_parameters;
  parameterinput >> number_of_heating_parameters;
  if (number_of_heating_parameters != 5)
  {
    cerr << "Need 5 heating parameters! The code is new, old ways won't work!  Number given is " << number_of_heating_parameters <<  endl;
    exit(5);
  }

  if (verbose){ cout << "No: " << number_of_heating_parameters << endl; }
  getline(parameterinput,DiscardMe);
  double *heating_parameters;
  if (number_of_heating_parameters > 0)
  {
    heating_parameters = new double [number_of_heating_parameters];
    // These parameters are useful for determining heating data.
    for (int i = 0; i < number_of_heating_parameters-1; ++i)
    {
      parameterinput >> heating_parameters[i] >> Comma;
      if (verbose){ cout << i << ", " << heating_parameters[i] << endl; }
        }
    parameterinput >> heating_parameters[number_of_heating_parameters-1];
    if (verbose){ cout << number_of_heating_parameters-1<< ", " << heating_parameters[number_of_heating_parameters-1] << endl; }
  }
  else
  {
    heating_parameters = NULL;
  }

  // NEW: HEAT FILE!!!!!!!
  getline(parameterinput,DiscardMe);
  string heatFileName;
  parameterinput >>heatFileName ;
  getline(parameterinput,DiscardMe);
  if (verbose){ cout << heatFileName << endl; }

  int number_of_velocity_parameters;
  parameterinput >> number_of_velocity_parameters;
  if (verbose){ cout <<  "Vel: " << number_of_velocity_parameters << endl; }

  getline(parameterinput,DiscardMe);
  double *velocity_parameters;
  if (number_of_velocity_parameters > 0)
  {
    velocity_parameters = new double [number_of_velocity_parameters];
    // These parameters are useful for determining heating data.
    for (int i = 0; i < number_of_velocity_parameters-1; ++i)
        {
      parameterinput >> velocity_parameters[i] >> Comma;
      if (verbose){ cout << velocity_parameters[i] << endl; }
        }
    parameterinput >>  velocity_parameters[number_of_velocity_parameters-1];
  }
  else
  {
    velocity_parameters = NULL;
  }

  getline(parameterinput,DiscardMe);
  int number_of_temperature_parameters;
  parameterinput >> number_of_temperature_parameters;
  if (verbose){ cout << "No: " << number_of_temperature_parameters << endl; }

  getline(parameterinput,DiscardMe);
  double *temperature_parameters;
  if(number_of_temperature_parameters >0)
  {
    temperature_parameters = new double [number_of_temperature_parameters];
    // These parameters are useful for determining heating data.
    for (int i=0; i<number_of_temperature_parameters-1; ++i)
        {
      parameterinput >>  temperature_parameters[i] >> Comma ;
      if (verbose){ cout << i << ", " << temperature_parameters[i] << endl; }
        }
    parameterinput >>  temperature_parameters[number_of_temperature_parameters-1];
  }
  else
  {
    temperature_parameters = NULL;
  }

  getline(parameterinput,DiscardMe);

  // Insert the ability to restart from a saved temperature profile.
  // The idea is that either you have a 0 or you have a 1 followed by a comma then the
  // filename of the restore file.
  int restart;
  string restartFileName;

  parameterinput >> restart;
  if (restart)
  {
    parameterinput >> Comma >> restartFileName;
    if (verbose){ cout << restartFileName << endl; }
  }
  getline(parameterinput,DiscardMe);

  int number_of_conductivity_parameters;
  parameterinput >> number_of_conductivity_parameters;
    if (verbose){ cout << "No: " << number_of_conductivity_parameters << endl; }

  getline(parameterinput,DiscardMe);
  double *conductivity_parameters;
  if (number_of_conductivity_parameters >0)
  {
    conductivity_parameters = new double [number_of_conductivity_parameters];
    // These parameters are useful for determining heating data.
    for (int i=0; i < number_of_conductivity_parameters-1; ++i)
        {
      parameterinput >> conductivity_parameters[i] >> Comma ;
      if (verbose){ cout << i << ", " << conductivity_parameters[i] << endl; }
        }
    parameterinput >> conductivity_parameters[number_of_conductivity_parameters-1];
    getline(parameterinput,DiscardMe);
  }
  else{
    conductivity_parameters = NULL;
  }

  int number_of_diffusivity_parameters;
  parameterinput >> number_of_diffusivity_parameters;
  if (verbose){ cout << "No: " << number_of_diffusivity_parameters << endl; }

  getline(parameterinput,DiscardMe);
  double *diffusivity_parameters;
  if (number_of_diffusivity_parameters >0)
  {
    diffusivity_parameters = new double [number_of_diffusivity_parameters];
    // These parameters are useful for determining heating data.
    for (int i=0; i < number_of_diffusivity_parameters-1; ++i)
    {
      parameterinput >> diffusivity_parameters[i] >> Comma;
      if (verbose){ cout << i << ", " << diffusivity_parameters[i] << endl; }
    }
    parameterinput >> diffusivity_parameters[number_of_diffusivity_parameters-1];
  }
  else{
    diffusivity_parameters = NULL;
  }

    if (verbose){ cout << nx << ", "<< ny << ", "<< nz << endl; }
  if (verbose){ cout << dx << ", "<< dy << ", "<< dz << endl; }

  parameterinput.close();
*/
  double con[input.number_of_conductivity_parameters];
  std::copy(input.conductivity_parameters.begin(), input.conductivity_parameters.end(), con);

  double dif[input.number_of_diffusivity_parameters];
  std::copy(input.diffusivity_parameters.begin(), input.diffusivity_parameters.end(), dif);

  double vel[input.number_of_velocity_parameters];
  std::copy(input.velocity_parameters.begin(), input.velocity_parameters.end(), vel);

  double hea[input.number_of_heating_parameters];
  std::copy(input.heating_parameters.begin(), input.heating_parameters.end(), hea);

  double tem[input.number_of_temperature_parameters];
  std::copy(input.temperature_parameters.begin(), input.temperature_parameters.end(), tem);

  double *timevector = new double [input.number_of_iterations_per_save * input.number_of_saves];
  double *centralpoint = new double [input.number_of_iterations_per_save * input.number_of_saves];

  // LOADED

  // This is different to how it used to be, because I'm going directly to defining the conductivity. It is very simplistic: there is a single vessel, which is cylindrical, so it's pretty easy

  // The definitions of the conductivity parameters will be found in conductivity.cpp.

/*  conductivity kappa(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, input.number_of_conductivity_parameters, input.conductivity_parameters);

  conductivity diff(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, input.number_of_diffusivity_parameters, input.diffusivity_parameters);
*/
  conductivity kappa(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, input.number_of_conductivity_parameters, con);

  conductivity diff(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, input.number_of_diffusivity_parameters, dif);

  if (verbose){ cout << " Done cond" << std::endl; }

  Array3D<double>  temperature(input.nx, input.ny, input.nz);
  Array3D<double>  heating(input.nx, input.ny, input.nz);
  // Heating  will be ramped on and off, and will be scaled by the acoustic absorption of blood in the vessel voxels.
  Array3D<double> ThermalDose(input.nx, input.ny, input.nz);
  // Keep track of the thermal dose.

  // Now call some generic functions that generate data for these arrays.

  if (verbose){ cout << "Assigned arrays. " << std::endl; }

  if (input.restart)
  {
    GenerateTemperature(temperature, input.restartFileName); // Restart from saved file.
  }
  else
  {
    /*
    GenerateTemperature(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, temperature, input.number_of_temperature_parameters, input.temperature_parameters);
    */
    GenerateTemperature(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, temperature, input.number_of_temperature_parameters, tem);
  }

  if (verbose){ cout << " Done temp" << endl; }
/*      GenerateHeating(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, heating, input.number_of_heating_parameters, input.heating_parameters, input.number_of_velocity_parameters, input.velocity_parameters, input.heatFileName);
*/
  GenerateHeating(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, heating, input.number_of_heating_parameters, hea, input.number_of_velocity_parameters, vel, input.heatFileName);
  // Interestingly, most of the heating parameters aren't used here. All I do here is load up an acoustic
  // array and reduce the heating value for the bits that are in blood. HOWEVER, I do need to rescale
  // everything so that the maximum heating value is equal to the chosen value (of course, the acoustic
  // stuff is just acoustics, and the values aren't necessarily dimensionally useful).
    // Velocity parameters are xCentre, yCentre, radius, maxVelocity.

  if (verbose){ cout << " Done heating" << endl; }


  Velocity velocity;
/*  velocity.initVelocity(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, input.number_of_velocity_parameters, input.velocity_parameters);
*/
  velocity.initVelocity(input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, input.number_of_velocity_parameters, vel);

  // Set it up!

  double time = 0;
  // this is a time variable that keeps track of the time for which the heat equation has run.
  // Very useful for when we have pulsed heating, etc.

  // It's also useful to track the largest of the n's.
  int nmax;
  if ( (input.nx > input.ny) && (input.nx > input.nz) ){
    nmax = nx;
  }
  else if (input.ny > input.nz){
    nmax = input.ny;
  }
  else{
    nmax = input.nz;
  }

  // Now, check that the directory to which you want to save things didn't exist before now and then make it.
  std::string temporary = "mkdir " + input.dataFileName;
  int didthisfail = system( temporary.c_str() ) ;
  if (verbose){ cout << "filename " << input.dataFileName << endl; }
  if (verbose){ cout << "result " << didthisfail << endl; }
  if (didthisfail)
  {
    // I don't give the user the chance to overwrite data. The thing is, I often just
    // press yes without thinking REALLY hard. Or I press enter or something. And I could wipe
    // out my drive by doing this. So don't do it.
    std::cerr << "Directory exists. Don't overwrite old data. " << std::endl;
    return -1;
  }

  /*
  New addition: Let's make some work vectors for use in the tridiagonal routines.
  Then I don't have to constantly delete and reallocate them.
  This is actually a generic idea of mine. Any class that uses work vectors will just hold onto
  them on the assumption that it will need them at some point.
  In principle, my velocity class (has 3 work vectors) could just store pointers to the following
  workspace variables. But this assumes that I'm not going to screw anything up... It's definitely
  better to just encapsulate the vectors within velocity. Much more <class>y... :-)
  It means that range checking is easier (if I really wanted this, I would use vector<double> not double*, but I'm a sucker for the old ways).
  */

  double *worka = new double[nmax];
  double *workb = new double[nmax];
  double *workc = new double[nmax];
  double *workd = new double[nmax];

  int countnumberoftimesteps = 0;
  for (int i = 0; i < input.number_of_saves; ++i)
  {
    if (verbose){ cout << "In main loop, iteration: " << i << std::endl; }

    // Recall the principle of the split-step method: if we do half a velocity step first, then do full conduction
    // steps followed by full velocity steps, then things are more accurate than if we don't do half a velocity step first.
    // i.e. we do: (0.5*L1)  L2 L1 ...  L2 L1  L2 (0.5*L1) T(0), which gives good accuracy.
    // Problem is that our intermediate solutions aren't real solutions. However, they are accurate everywhere except in the
    // moving bit, so they're fine for calculating thermal dose.
    // * My first idea was to put the half-step bit in this loop, so that run_heat_solver only works with full steps.
    // But I don't like it. Basically, I imagine that at some point, someone will say that we can look at the thermal
    // dose of the blood, and that is half a velocity step inaccurate. I know this is barely relevant. But it matters to me.
    // * So I'm going to have the half-steps at the start and end of each loop within run_heat_solver.
    // Since run_heat_solver iterates over 6 steps per loop, we only have an extra step every 6 time steps. Anyway,
    // the velocity is supposed to be cheap! So just do it and stop whining.

    run_heat_solver(nmax, input.nx, input.ny, input.nz, input.dx, input.dy, input.dz, input.dt, time, input.number_of_iterations_per_save, \
     kappa, diff, velocity, temperature, heating, input.number_of_heating_parameters, \
     hea, ThermalDose, worka, workb, workc, workd, centralpoint, \
     countnumberoftimesteps, timevector, input.perfusion_parameter);

    // Pass by reference.
    if (input.save_option)
    {
      savetofile_binary(input.dataFileName, temperature, time);
    }

    // Introduced a save_option, to allow us to avoid saving if we don't want to.
    if (input.save_option2)
    {
      savetofile_integrated_quantities_binary(input.dataFileName, ThermalDose, centralpoint, countnumberoftimesteps, timevector, time);
    }

  }

  if (verbose){ cout << "deleted 1. " << endl; }
  delete[] worka;
  if (verbose){ cout << "deleted 2. " << endl; }
  delete[] workb;
  if (verbose){ cout << "deleted 3. " << endl; }
  delete[] workc;
  if (verbose){ cout << "deleted 4. " << endl; }
  delete[] workd;
  if (verbose){ cout << "deleted 5. " << endl; }

  xtimer.stop();
  double total_time = xtimer.elapsed()/100.0;
  int hours = floor(total_time/3600.0);
  int minutes = floor((total_time-hours*3600.0)/60.0);
  double seconds = total_time - hours*3600.0 - minutes*60.0;
  cout.precision(4);
  cout << "\nElapsed time: " << hours << " hours, " << minutes << " minutes, " << seconds << " seconds\n" << endl;

  return  0;
}
