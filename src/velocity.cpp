
/**
 \file velocity.cpp
*/

#include "velocity.h"


//#define DEBUG_VELOCITY_EVOLVER

/**
 @brief Constructor of velocity class
 */
Velocity::Velocity()
{
  nz = 1;
  numberOfStreamlines = 1;
  dz = 0.1;
  work  = new double[nz];
  work2 = new double[nz];
  work3 = new double[nz];
  work4 = new double[nz];
  streamlineX = new int[numberOfStreamlines];
  streamlineY = new int[numberOfStreamlines];
  velocity = new double [numberOfStreamlines];
}


/**
 @brief Destructor of velocity class
 */
Velocity::~Velocity()
{
  // Simple destructor
  delete[] work;
  delete[] work2;
  delete[] work3;
  delete[] work4;
  delete[] streamlineX;
  delete[] streamlineY;
  delete[] velocity;
}


/**
 @brief Get initial velocity
 */
void Velocity::initVelocity(int nx, int ny, int a_nz,  double dx, double dy, double a_dz, int a_numVelocityPars, double* a_velocityPars)
{
  // Need nx,ny, dx,dy to figure out where the vessel is in space. The dz and nz are just for later use when actually doing the evolution.
  nz = a_nz;
  dz = a_dz;
  if (a_numVelocityPars != 4){
    std::cerr << "Trouble in setting up velocity, wrong number of parameters " << std::endl;
    exit(5);
  }
  double xC = a_velocityPars[0];
  double yC = a_velocityPars[1];
  double radius = a_velocityPars[2];
  double velocityMax = a_velocityPars[3];

  // First, figure out how many streamlines there are.
  numberOfStreamlines = 0;
  if (radius > 0.0)
  {
    // No point doing anything for a vessel that ain't there.
    for (int x = 0; x < nx; x++)
    {
      for (int y=0; y < ny; y++)
      {
        if (((x*dx-xC)*(x*dx-xC) + (y*dy-yC)*(y*dy-yC)) < radius*radius)
        {
          numberOfStreamlines++;
        }
      }
    }
  }

#ifdef DEBUG_VELOCITY_EVOLVER
  std::cout << "Number of streamlines " << numberOfStreamlines << std::endl;
#endif

  // Now scrap the old vectors and grab the correct amount of space.
  delete[] work;
  delete[] work2;
  delete[] work3;
  delete[] work4;
  work = new double [nz];
  work2 = new double [nz];
  work3 = new double [nz];
  work4 = new double [nz];

  if (numberOfStreamlines > 0)
  {
    // Should sense check things! Note, if numberOfStreamlines ==0, then the evolution of
    // the velocity doesn't happen, so no worries.
    delete[] streamlineX;
    delete[] streamlineY;
    delete[] velocity;
    streamlineX = new int[numberOfStreamlines];
    streamlineY = new int[numberOfStreamlines];
    velocity = new double[numberOfStreamlines];

    // Now get insert locations and things.
    double tmpDist2; // Measures the distance from the centre of the circle.
    int count = 0; // Keep track of which streamline we're on.
    if (radius > 0.0)
    {  // No point doing anything for a vessel that ain't there.
      for (int x = 0; x < nx; x++)
      {
        for (int y = 0; y < ny; y++)
        {
          tmpDist2 = (((x*dx-xC)*(x*dx-xC) + (y*dy-yC)*(y*dy-yC)));
          if (tmpDist2 < radius*radius)
          {
            // i.e. if we're in a region where there is flow.
            streamlineX[count] = x; // save location
            streamlineY[count] = y;
            velocity[count] = velocityMax*(1.0 - tmpDist2/(radius*radius)); // save velocity, Poiseuille flow.
            count++;
          }
        }
      }
    }
  }
  // Finished!
}


/**
 @brief Evolve the velocity profile
 */
void Velocity::evolveVelocity(Array3D<double> &temperature, double time)
{
  bool verbose = false;

  // Doesn't seem much point putting this in a separate file.

  int tmpM; // This is a temp value for storing the integer number of dz steps traversed in a single time step.
  double tmpF ; // This is for the fractional number of steps traversed by the velocity.
  for (int n = 0; n < numberOfStreamlines; n++)
  {
    // Step over each streamline.

#ifdef DEBUG_VELOCITY_EVOLVER
    //std::cerr << "Diagnostic info" << std::endl;
    //std::cerr << streamlineX[n] << ",  " <<  streamlineY[n]<< std::endl;
    //std::cerr << work  << ", " << work+nz << ", " << streamlineX <<  std::endl;
#endif

    temperature.Output_Array_Of_Values(streamlineX[n], streamlineY[n], 0, work, 3);
    // NB the 0 is a dummy, the 3 indicates that we move in the z-direction. Now we have a
    // streamline of temperatures.

    tmpM = (velocity[n] > 0) ? int(floor(velocity[n]*time/dz)) : int(ceil(velocity[n]*time/dz));
    // Get the integer number of steps moved. Essentially, I'm truncating the fraction, but I don't
    // think C has a trunc() function. Can't use modf because tmpM must be integer, not double.

    tmpF = velocity[n]*time/dz - tmpM;
    // tmpF has same sign as tmpM by choice. Must test bounds on tmpM because otherwise you can
    // overwrite other data (or rather, I tested bounds on tmpF, overwrote data, generated SIGSEGV
    // and now I'm rewriting...)

#ifdef DEBUG_VELOCITY_EVOLVER
  if (tmpF*tmpM < 0) std::cerr << "ERROR  " << tmpM << ",  " << tmpF<< std::endl;
#endif

    if (tmpM > 0)
    {
      // In this case, everything moves to the right, so we should start at the right to
      // calculate new values. i.e. start downstream, because it doesn't matter if it is overwritten.
      for (int pp = nz-1; pp >= tmpM ; pp--)
      {
        // Shift downwards
        work[pp] = work[pp-tmpM];
      }
      for (int pp = tmpM-1; pp >= 0 ; pp--)
      {
        // zero the ones that are from outside the region.
        work[pp] = 0.0;
      }
      // Moved the integer number of steps. Now interpolate to get the rest.
      cubicSplineInterpolateFOR_EQUIDISTANT_GRID(work, nz, dz, tmpF, work2, work3, work4);
    }
    else if (tmpM < 0)
    {
      // In this case, everything moves to the left, so we should start at the left to
      // calculate new values. i.e. start downstream, because it doesn't matter if it is overwritten.
      for (int pp=0; pp < nz+tmpM; pp++)
      {
        // Shift downwards
        work[pp] = work[pp-tmpM];  // Same formula, but because tmpM is now negative, things move left.
      }
      for (int pp=nz+tmpM; pp < nz; pp++)
      {
        // zero the ones that are from outside the region.
        work[pp] = 0.0;
      }
      // Moved the integer number of steps. Now interpolate to get the rest.
      cubicSplineInterpolateFOR_EQUIDISTANT_GRID(work, nz, dz, tmpF, work2, work3, work4);
    }
    else if ( (tmpM == 0) && (tmpF > 0.0) )
    {
      // move to the right by fractional part.
      cubicSplineInterpolateFOR_EQUIDISTANT_GRID(work, nz, dz, tmpF, work2, work3, work4);
    }
    else
    {
      // This is the case where tmpM = 0 and tmpF <= 0. It should be very rare that both are zero,
      // so just include it here.
      cubicSplineInterpolateFOR_EQUIDISTANT_GRID(work, nz, dz, tmpF, work2, work3, work4);
    }

    // Done streamline. Put it back into temperature.
    temperature.Input_Array_Of_Values(streamlineX[n], streamlineY[n], 0, work, 3);
  }
  // end of loop. i.e. done it for all streamlines.

#ifdef DEBUG_VELOCITY_EVOLVER
  std::cerr << "Diagnostic info" << std::endl;
  for (int n = 0; n < numberOfStreamlines; n++){
    // Step over each streamline.
    std::cerr << streamlineX[n] << ",  " <<  streamlineY[n]<< std::endl;
  }
#endif

} // end of function.



/*
 @brief Solver for tridiagonal system.

Solve a tridiagonal scheme that is a little different to the usual, in that a lot of the values
are equal so we can save memory.
I'm considering the matrix      (a b                          )
                (c d e                        )
                (  c d e                      )
        A =       (    c d e                    )
                (    ....          )
                (             c d e)
                (                          f g)

and solving the problem A*soln = rhs.
This type of matrix is very common on equidistant grids, and the solution should be faster than the
usual method. I've taken the method from Numerical Recipes in C, modified accordingly
**/

void nonGeneralTridiagonalSolver(double a, double b, double c, double d, double e, double f, double g, double *rhs, int n, double *soln, double *work)
{
  int j;
  double bet;
  if (a == 0.0)
  {
    std::cerr << "Error in tridiag solver: a mustn't be zero" << std::endl;
    exit(10);
  }
  soln[0] = rhs[0] / (bet=a);

  // Forward substitution, must be done differently because j=1 means j-1 = 0 which is a special row.
  work[1] = b / bet;
  bet = d - c * work[1];

  if (bet == 0.0)
  {
    std::cerr << "Error 2.1 in tridag" << std::endl;
    exit(10);
  }
  soln[1] = (rhs[1] - c * soln[0]) / bet;

  for (j = 2; j < (n-1); j++)
  {
  // Decomposition and forward substitution.
  // NOTE I've had to chop out j=1 and j=n-1 and treat them separately.
    work[j] = e / bet;
    bet = d - c * work[j];
    if (bet == 0.0)
    {
      std::cerr << "Error 2.2 in tridag" << std::endl;
      exit(10);
    }
    soln[j] = (rhs[j] - c * soln[j-1]) / bet;
  }

  // Finally do it for the j = nz-1
  work[n-1] = e / bet;
  bet = g - f * work[n-1];

  if (bet == 0.0)
  {
    std::cerr << "Error 2.3 in tridag" << std::endl;
    exit(10);
  }
  soln[n-1] = (rhs[n-1] - f * soln[n-2]) / bet;

  // Next, back substitution. This doesn't seem to care about the form of the matrix,
  // so we just use theirs, but shifted down 1 in index as usual
  for ( j = (n-2); j >= 0; j--)
  {
    // Backsubstitution
    soln[j] -= work[j+1] * soln[j+1];
  }
}

/**
 @brief cubic spline interpolation.

 This function uses a cubic spline to evolve the velocity for the small fractional part.

 The main part of it is done directly: that is, I assume that the\f$n\f$ points in \f$y\f$ are equidistantly spaced (spacing \f$dx\f$), and that I'm going to push all of them by an amount \f$f dx\f$, where \f$|f| < 1\f$. Because of these assumptions, this is much faster than the general cubic spline method from Numerical Recipes book.

 If \f$f > 0\f$ then things move to the right and therefore we have to do backwards interpolation (to find out where the new value comes from). \f$f < 0\f$ implies forward interpolation. Basically, we have to get data from upstream.

 There are many ways to solve the cspline problem, basically down to your choice of boundary conditions. These are implemented as preprocessor directives. Perhaps the most typical situation typical one is that the first derivative vanishes at the boundary, but I'm not sure.

 The equidistant criterion greatly simplifies calculating the second derivatives. So a simple implemented a simpler tridiag solver to take advantage of this.
**/
void cubicSplineInterpolateFOR_EQUIDISTANT_GRID(double *y, int n,  double dx, double f, double *workA, double *workB, double *workC)
{
  // Generate some constants that will be used throughout:
  double oneMinusF = 1.0 - fabs(f);
  double coeff1 = (f * f - 1.0) * fabs(f) * dx * dx / 6.0;
  double coeff2 = (oneMinusF*oneMinusF - 1.0) * oneMinusF * dx * dx /6.0;
  double coeff3 = 6.0 / (dx * dx);

  // First step of the cubic spline method is to get the array of second derivatives.
  // These will be stored in workA. To get them requires the inversion of a tridiagonal,
  // so I need to store the RHS of this system in workB. workC is used as a workspace region.
  for (int p = 1; p < (n-1); p++)
  {
    workB[p] = coeff3*(y[p-1] - 2.0 * y[p] + y[p+1]);
  }

#if (VELOCITY_CUBIC_SPLINE_CHOICE == 1)
  // This case corresponds to zero second derivatives at both boundaries.
  workB[0] = 0.0;
  workB[n-1] = 0.0;
  nonGeneralTridiagonalSolver(1.0, 0.0, 1.0, 4.0, 1.0, 0.0, 1.0, workB, n, workA, workC);
  // Note that when second deriv is 0, we have (1  0  0 ...) and (0 0 .... 0 1)  as first and last rows of matrix.

#elif VELOCITY_CUBIC_SPLINE_CHOICE == 2
  // This case corresponds to zero second derivatives at left, zero first at right
  workB[0] = 0.0;
  workB[n-1] = coeff3*(y[n-2] - y[n-1]);
  nonGeneralTridiagonalSolver(1.0, 0.0, 1.0, 4.0, 1.0, 1.0, 2.0, workB, n, workA, workC);
  // Note that when first deriv is 0, we have  (0 0 .... 1 2)  as   last row of matrix.

#elif VELOCITY_CUBIC_SPLINE_CHOICE == 3
  // This case corresponds to zero second derivatives at right, zero first at left.
  workB[0] = coeff3*(y[1]-y[0]);
  workB[n-1] = 0.0;
  nonGeneralTridiagonalSolver(2.0, 1.0, 1.0, 4.0, 1.0, 0.0, 1.0, workB, n, workA, workC);
  // Note that when first deriv is 0, we have (2 1   0 ...)  as first  row of matrix.

#elif VELOCITY_CUBIC_SPLINE_CHOICE == 4
  // This case corresponds to zero first derivatives at both boundaries.
  workB[0] = coeff3*(y[1]-y[0]);
  workB[n-1] = coeff3*(y[n-2] - y[n-1]);
  nonGeneralTridiagonalSolver(2.0, 1.0, 1.0, 4.0, 1.0, 1.0, 2.0, workB, n, workA, workC);
  // Note that when second deriv is 0, we have (1  0  0 ...) and (0 0 .... 0 1)  as first and last rows of matrix.

#else
  // DEFAULT
  // This case corresponds to zero second derivatives at both boundaries.
  workB[0] = 0.0;
  workB[n-1] = 0.0;
  nonGeneralTridiagonalSolver(1.0,0.0, 1.0, 4.0, 1.0, 0.0, 1.0, workB, n, workA, workC);
  // Note that when second deriv is 0, we have (1  0  0 ...) and (0 0 .... 0 1)  as first and last rows of matrix.

#endif

  // OK, now we've generated the second derivatives. Now we need to do the interpolation.
  if (f >= 0.0)
  {
    for (int p = n-1; p > 0; p--)
    {
      y[p] = f * y[p-1] + oneMinusF * y[p] + coeff1 * workA[p-1] + coeff2 * workA[p];
    }
    y[0] = 0.0; // flow from outside.
  }
  else
  {
    for (int p = 0; p < (n-1); p++)
    {
      y[p] = oneMinusF * y[p] + fabs(f) * y[p+1] + coeff2 * workA[p] + coeff1 * workA[p+1];
    }
    y[n-1] = 0.0; // flow from outside.
  }
}
