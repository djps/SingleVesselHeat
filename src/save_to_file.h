/**
 @file save_to_file.h
 */

#ifndef SAVE_TO_FILE
#define SAVE_TO_FILE

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include "array3d.cpp"
#include "conductivity.h"
#include "generator.h"
#include "run_heat_solver.h"

// This handles all the problems of saving to file. At the moment, it's just a
// single function so we only need a header file and we just #include the whole
// thing into main. I guess if we make it a bit more specialized it would be
// better to have it compile to a separate .o file. Should do this all as
// binary files, but at the moment just use text.

void savetofile(std::ofstream &data, Array3D<double> &temperature, double &time);

/**
 @brief ostreaming
 */
void savetofile(std::ofstream &data, Array3D<double> &temperature, double &time)
{
  int nx = temperature.nx;
  int ny = temperature.ny;
  int nz = temperature.nz;
  data << time << std::endl;
  for (int z = 0; z < nz; ++z)
  {
    for (int y = 0; y < ny; ++y)
    {
      for (int x = 0; x < nx; ++x)
      {
        data << x << ", " << y << ", " << z  << ", " << \
       temperature.Output_Value(x,y,z) << "\n";
      }
    }
  }
  data << std::flush;
}

/**
 @brief saves in another way: this breaks the symmetry of the way things are saved to try make it easy to plot.
 */
void savetofile2(std::ofstream &data, Array3D<double> &temperature, double &time)
{
  int nx = temperature.nx;
  int ny = temperature.ny;
  int nz = temperature.nz;

  for (int z = 0; z < nz; z++)
  {
    for (int y = 0; y < ny; y++)
    {
      data << y << ", " << z;
      for (int x = 0; x < nx; x++)
      {
        data  << ", " << temperature.Output_Value(x, y, z);
      }
      data << std::endl;
    }
  }
}


/**
 @brief This method takes in a string, which should be viewed as a directory name. It then concatenates the file name TempXXX.dat where XXX is int(time X 100).
 */
void savetofile(std::string &filename, Array3D<double> &temperature, double &time)
{
  // std::cout << "Now at time " << time << std::endl;
  int nx = temperature.nx;
  int ny = temperature.ny;
  int nz = temperature.nz;

  std::ostringstream pathname ;
  pathname << filename << "/" << "Temp" << int(100*time + 0.5) << ".dat";
// Note: the +0.5 ensures that we do a round off not an integer truncation.

  std::ofstream data;
  std::string fname = pathname.str();
//std::cout << fname << std::endl;
  data.open(fname.c_str());

  for (int z = 0; z < nz; ++z)
  {
    for (int y = 0; y < ny; ++y)
    {
      for (int x = 0; x < nx; ++x)
      {
        data << x << ", " << y << ", " << z  << ", " << temperature.Output_Value(x,y,z) << "\n";
      }
    }
  }
  data<<std::flush;
  data.close();
}


/**
 @brief Saves the accumulated data files.
*/
void savetofile_integrated_quantities(std::string &filename, \
 Array3D<double> &ThermalDose, Array3D<double> &AccumulatedTemperature, \
 double *central, int NNN, double *time)
{
  int nx = ThermalDose.nx;
  int ny = ThermalDose.ny;
  int nz = ThermalDose.nz;

  std::ostringstream pathname1, pathname2, pathname3 ;
  pathname1 << filename << "/" << "ThermalDose.dat";
  pathname2 << filename << "/" << "IntegratedTemp.dat";
  pathname3 << filename << "/" << "CentralTemp.dat";

  std::ofstream data1,data2,data3;
  std::string fname1 = pathname1.str();
  std::string fname2 = pathname2.str();
  std::string fname3 = pathname3.str();

  data1.open(fname1.c_str());
  data2.open(fname2.c_str());
  data3.open(fname3.c_str());

  for (int z = 0; z < nz; ++z)
  {
    for (int y = 0; y < ny; ++y)
    {
      for (int x = 0; x < nx; ++x)
      {
        data1 << x << ", " << y << ", " << z  << ", " << ThermalDose.Output_Value(x,y,z) << "\n";
        data2 << x << ", " << y << ", " << z  << ", " << AccumulatedTemperature.Output_Value(x,y,z) << "\n";
      }
    }
  }

  for (int ppp = 0; ppp < NNN; ppp++)
  {
    data3 << time[ppp] << ",  " << central[ppp] << "\n";
  }
  data1.close();
  data2.close();
  data3.close();
}

/**
 @brief Saves file to binary
 */
void savetofile_binary(std::string &filename, Array3D<double> &temperature, \
 double &time)
{
  // As previously, but this produces a binary file. This binary file has the
  // format: nx, ny, nz (each integers), temperature (which is an array of
  // nx*ny*nz doubles). This method takes in a string, which should be viewed
  // as a directory name. It then concatenates the file name TempXXX.dat where
  // XXX is int(time X 100). I'm hoping this will work.

  // std::cout << "Now at time " << time << std::endl;
  int nx = temperature.nx;
  int ny = temperature.ny;
  int nz = temperature.nz;

  bool verbose = true;

  std::ostringstream pathname;
  pathname << filename << "/" << "Temp" << int(100*time + 0.5) << ".dat";
  // Note that the +0.5 ensures that we do a round off not an integer truncation.

  if (verbose){ std::cout << pathname.str() << std::endl; }

  std::ofstream data;
  std::string fname = pathname.str();
  data.open(fname.c_str(), std::ios::out | std::ios::binary);
  data.write((char*) &nx, sizeof(int)); // with casts to char stream
  data.write((char*) &ny, sizeof(int));
  data.write((char*) &nz, sizeof(int));
  data.write((char*) temperature.array, sizeof(double)*(nx*ny*nz));

  data.close();
}

/**
 @brief save integrated quantities to binary
 */
void savetofile_integrated_quantities_binary(std::string &filename, \
 Array3D<double> &ThermalDose, double *central, int NNN, double *time, \
 double &timeValue)
{
  // As previously, but this produces a binary file. The first 2 binary files
  // have the format: nx, ny, nz (each integers), 3darray (which is an array
  // of nx*ny*nz doubles). The third has the format NNN (int), time
  // (NNN*doubles), centraltemperature (NNN*doubles).

  // Saves the accumulated data files.
  int nx = ThermalDose.nx;
  int ny = ThermalDose.ny;
  int nz = ThermalDose.nz;

  std::ostringstream pathname1, pathname3;
  pathname1 << filename << "/" << "ThermalDose" << int(100*timeValue + 0.5) << ".dat";
  //pathname2 << filename << "/" << "IntegratedTemp.dat";
  pathname3 << filename << "/" << "CentralTemp.dat";

  std::ofstream data1, data3;
  std::string fname1 = pathname1.str();
  //std::string fname2 = pathname2.str();
  std::string fname3 = pathname3.str();

  data1.open(fname1.c_str(), std::ios::out | std::ios::binary);
  //data2.open(fname2.c_str(), std::ios::out | std::ios::binary);
  data3.open(fname3.c_str(), std::ios::out | std::ios::binary);

  data1.write((char*) &nx, sizeof(int) );
  data1.write((char*) &ny, sizeof(int) );
  data1.write((char*) &nz, sizeof(int) );
  data1.write((char*) ThermalDose.array, sizeof(double)*(nx*ny*nz) );

  data3.write((char*) &NNN, sizeof(int) );
  data3.write((char*) time, sizeof(double)*NNN );
  data3.write((char*) central, sizeof(double)*NNN );

  data1.close();
  //data2.close();
  data3.close();
}

#endif
