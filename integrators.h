/**
 @file integrators.h
 */

/**
 @brief Here we implement different integrators. i.e. these take in a 1-d array, a number of points and step size and return the integral.
*/

#ifndef INTEGRATORS
#define INTEGRATORS

#include <iostream>

double simpson_rule(double *y, double dx, int N);

#endif
