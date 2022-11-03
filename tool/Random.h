// Copyright (c) 2013, Wenjun Ying, Department of Mathematics and Institute of 
// Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240, P.R. China.

// The graphics package was developed for the course of numerical analysis and
// scientific computing that the author taught in the Zhiyuan College of 
// Shanghai Jiao Tong University in the Spring of 2013. 

// This file is free software; as a special exception the author gives unlimited
// permission to copy and/or distribute it, with or without modifications, as 
// long as this notice is preserved.

// This file is distributed in the hope that it will be useful, but WITHOUT ANY 
// WARRANTY, to the extent permitted by law; without even the implied warranty 
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#ifndef __Random_h_IS_INCLUDED__
#define __Random_h_IS_INCLUDED__ 

#include <cstdlib> 
#include <cmath> 

#include <time.h>
#include <float.h>
#include <limits.h>

/**
 *****************************************************************************
 * generate a random real number in the interval [a,b]:                      *
 *****************************************************************************
 **/
inline double Random(double a, double b)
{
  static unsigned int seed = 0;
  const int NUMBER_OF_ROTATION_BITS = 20; 
  const unsigned int MAX_UNSIGNED_RANDOM = 1048575;   // 2^20-1
  if (!seed) {
    unsigned int t = time(0); 
    int n = 0;
    unsigned int temp = t; 
    while (temp) {
      temp /= 10; 
      n++;
    }
    int m = 1; 
    for (int i = 0; i < (n * 3) / 4; i++) {
      m *= 10; 
    }
    unsigned int s = t / m;
    unsigned int r = t - s * m;
    seed = r * m + s + r; 
  }

  srand(seed); 
  unsigned int k = rand()%MAX_UNSIGNED_RANDOM;
  double s = static_cast<double>(k) / (MAX_UNSIGNED_RANDOM);
  seed = rand(); 

  double r = a + s * (b - a); 
  return r;
}

/**
 *****************************************************************************
 * generate a random real number in the interval [0,1]:                      *
 *****************************************************************************
 **/
inline double Random(void)
{
  static unsigned int seed = 0;
  const int NUMBER_OF_ROTATION_BITS = 20;
  const unsigned int MAX_UNSIGNED_RANDOM = 1048575;   // 2^20-1
  if (!seed) {
    unsigned int t = time(0);
    int n = 0;
    unsigned int temp = t; 
    while (temp) {
      temp /= 10; 
      n++; 
    }
    int k;
    int m = 1;
    for (int i = 0; i < (n * 3) / 4; i++) {
      m *= 10;
    }
    unsigned int s = t / m;
    unsigned int r = t - s * m; 
    seed = r * m + s + r; 
  }

  srand(seed);
  unsigned int k = rand()%MAX_UNSIGNED_RANDOM;
  double r = static_cast<double>(k) / MAX_UNSIGNED_RANDOM;
  seed = rand();
  return r; 
}

/**
 *****************************************************************************
 * generate a random integer number:                                         *
 *****************************************************************************
 **/
inline int RandInt(int na, int nb) 
{
  if (na < nb) {
    double a = static_cast<double>(na);
    double b = static_cast<double>(nb + 1); 
    int k = static_cast<int>(Random(a, b)); 
    if (k < na) {
      k++;
    }
    if (k > nb) {
      k--;
    }
    return k;
  } else {
    double b = static_cast<double>(nb);
    double a = static_cast<double>(na + 1); 
    int k = static_cast<int>(Random(a, b));
    if (k > na) {
      k--; 
    }
    if (k < nb) {
      k++; 
    }
    return k;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double RandomExponential(void) 
{
  bool done = false; 
  double xi = Random(); 
  while ((xi + DBL_EPSILON) > 1.0) {
    xi = Random(); 
  }

  double t = - log(1.0 - xi);
  return t;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double RandomGaussian(void)
{
  double theta = (M_PI + M_PI) * Random(); 
  double t = RandomExponential();
  double r = sqrt(t + t); 
  double x = r * cos(theta); 
  //double y = r*sin(theta);
  return x;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void RandomGaussian(double &x, double &y)
{
  double theta = (M_PI + M_PI) * Random(); 
  double t = RandomExponential();
  double r = sqrt(t + t); 
  x = r * cos(theta); 
  y = r * sin(theta);
}

#endif 
