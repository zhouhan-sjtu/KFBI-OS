// Copyright (c) Wenjun Ying, Department of Mathematics and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is free software; as a special exception the author gives
// unlimited permission to copy and/or distribute it, with or without
// modifications, as long as this notice is preserved.

// This file is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY, to the extent permitted by law; without even 
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
// PURPOSE.

#ifndef __Tridiag_h_IS_INCLUDED__
#define __Tridiag_h_IS_INCLUDED__ 

#include <cassert>

// Array a[m] is the diagonal of a tridiagonal matrix.
// Array b[m] is the lower sub-diagonal of a tridiagonal matrix.
// Array c[m] is the upper sub-diagonal of a tridiagonal matrix.

/********************************************************************\
 * a[0]  b[0]                                                       *
 * c[0]  a[1]  b[1]                                                 *
 *       c[1]  a[2]  b[2]                                           *
 *             c[2]  a[3]  b[3]                                     *
 *                   ...   ...    ...                               *
 *                         ...    ...    ...                        *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 *    1                                                             *
 * gamma[0]    1                                                    *
 *          gamma[1]   1                                            *
 *                   gamma[2]   1                                   *
 *                             ...    ...                           *
 *                                    ...    ...                    *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * alpha[0]  beta[0]                                                *
 *          alpha[1]  beta[1]                                       *
 *                   alpha[2]  beta[2]                              *
 *                              ...     ...                         *
 *                                      ...    ...                  *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 *  alpha[0] = a[0];                                                *
 *  for (int i = 1; i < n; i++) {                                   *
 *    beta[i - 1] = b[i - 1];                                       *
 *    gamma[i - 1] = c[i - 1] / alpha[i - 1];                       *
 *    alpha[i] = a[i] - gamma[i - 1] * beta[i - 1];                 *
 *  }                                                               *
 *                                                                  *
\********************************************************************/


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void Tridiag(const double *a, const double *b, const double *c, 
                    const double *f, int n, double *x) 
{
  // make LU decomposition for the tri-diagonal matrix.

  double *alpha = new double[n];
  double *beta = new double[n - 1]; 
  double *gamma = new double[n - 1];

  alpha[0] = a[0];

  for (int i = 1; i < n; i++) {
    beta[i - 1] = b[i - 1]; 
    gamma[i - 1] = c[i - 1] / alpha[i - 1];
    alpha[i] = a[i] - gamma[i - 1] * beta[i - 1];
  }

  // forward substitution

  double *y = new double[n]; 

  y[0] = f[0]; 
  for (int i = 1; i < n; i++) {
    y[i] = f[i] - gamma[i - 1] * y[i - 1]; 
  }

  // backward substitution

  int n1 = n - 1; 
  x[n1] = y[n1] / alpha[n1];
  for (int i = n - 2; i >= 0; i--) {
    x[i] = (y[i] - beta[i] * x[i + 1]) / alpha[i];
  }

  delete[] y; 

  delete[] alpha;
  delete[] gamma;
  delete[] beta;
}

/********************************************************************\
 * a[0]  b[0]                                                       *
 * c[0]  a[1]  b[1]                                                 *
 *       c[1]  a[2]  b[2]                                           *
 *             c[2]  a[3]  b[3]                                     *
 *                   ...   ...    ...                               *
 *                         ...    ...    ...                        *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 *    1                                                             *
 * gamma[0]    1                                                    *
 *          gamma[1]   1                                            *
 *                   gamma[2]   1                                   *
 *                             ...    ...                           *
 *                                    ...    ...                    *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * alpha[0]  beta[0]                                                *
 *          alpha[1]  beta[1]                                       *
 *                   alpha[2]  beta[2]                              *
 *                              ...     ...                         *
 *                                      ...    ...                  *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 *  alpha[0] = a[0];                                                *
 *  for (int i = 1; i < n; i++) {                                   *
 *    beta[i - 1] = b[i - 1];                                       *
 *    gamma[i - 1] = c[i - 1] / alpha[i - 1];                       *
 *    alpha[i] = a[i] - gamma[i - 1] * beta[i - 1];                 *
 *  }                                                               *
 *                                                                  *
\********************************************************************/

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void solveByThomasAlgorithm(double a, double b, double c, 
																	 const double *f, int n, double *x)
{
  assert((0 != f) && (n > 0) && (0 != x));

  // make LU decomposition for the tri-diagonal matrix.

  double *alpha = new double[n]; 
  double *beta = new double[n - 1];
  double *gamma = new double[n - 1]; 

  alpha[0] = a; 
#ifdef DEBUG 
  if (fabs(alpha[0]) < DBL_EPSILON) {
    std::cout << "Failed in the tridiagonal solver." << std::endl;
    exit(EXIT_FAILURE); 
  }
#endif 

  for (int i = 1; i < n; i++) {
    beta[i - 1] = b; 
    gamma[i - 1] = c / alpha[i - 1];
    alpha[i] = a - gamma[i - 1] * beta[i - 1]; 
#ifdef DEBUG 
    if (fabs(alpha[i]) < DBL_EPSILON) {
      std::cout << "Failed in the tridiagonal solver." << std::endl; 
      exit(EXIT_FAILURE);
    }
#endif
  }

  // forward substitution

  double *y = new double[n]; 

  y[0] = f[0]; 
  for (int i = 1; i < n; i++) {
    y[i] = f[i] - gamma[i - 1] * y[i - 1]; 
  }

  // backward substitution

  int n1 = n - 1; 
  x[n1] = y[n1] / alpha[n1]; 
  for (int i = n - 2; i >= 0; i--) {
    x[i] = (y[i] - beta[i] * x[i + 1]) / alpha[i];
  }

  delete[] y; 

  delete[] alpha; 
  delete[] gamma;
  delete[] beta; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void solveByThomasAlgorithm(const double *a, const double *b, 
																	 const double *c, const double *f, 
																	 int n, double *x)
{
  assert((0 != a) && (0 != b) && (0 != c) && (0 != f) && (n > 0) && (0 != x));

  // make LU decomposition for the tri-diagonal matrix.

  double *alpha = new double[n]; 
  double *beta = new double[n - 1];
  double *gamma = new double[n - 1]; 

  alpha[0] = a[0]; 
#ifdef DEBUG 
  if (fabs(alpha[0]) < DBL_EPSILON) {
    std::cout << "Failed in the tridiagonal solver." << std::endl;
    exit(EXIT_FAILURE); 
  }
#endif 

  for (int i = 1; i < n; i++) {
    beta[i - 1] = b[i - 1]; 
    gamma[i - 1] = c[i - 1] / alpha[i - 1];
    alpha[i] = a[i] - gamma[i - 1] * beta[i - 1]; 
#ifdef DEBUG 
    if (fabs(alpha[i]) < DBL_EPSILON) {
      std::cout << "Failed in the tridiagonal solver." << std::endl; 
      exit(EXIT_FAILURE);
    }
#endif
  }

  // forward substitution

  double *y = new double[n]; 

  y[0] = f[0]; 
  for (int i = 1; i < n; i++) {
    y[i] = f[i] - gamma[i - 1] * y[i - 1]; 
  }

  // backward substitution

  int n1 = n - 1; 
  x[n1] = y[n1] / alpha[n1]; 
  for (int i = n - 2; i >= 0; i--) {
    x[i] = (y[i] - beta[i] * x[i + 1]) / alpha[i];
  }

  delete[] y; 

  delete[] alpha; 
  delete[] gamma;
  delete[] beta; 
}

/**
 *****************************************************************************
 * Keywords: cyclic, Cyclic tridiagonal matrix solver.                       *
 *                                                                           *
 *          a[0]  b[0]                                  gamma                *
 *          c[0]  a[1]  b[1]                                                 *
 *                c[1]  a[2]  b[2]                                           *
 *                      c[2]  a[3]  b[3]                                     *
 *                            ...   ...    ...                               *
 *          beta                    ...    ...    ...                        *
 *                                                                           *
 *****************************************************************************
 **/
inline void solveByThomasAlgorithm(const double a[], const double b[],
                                   const double c[], double beta, double gamma, 
                                   const double rhs[], int n, double x[])
{
  assert(n > 2);

  const double alpha = - a[0];

  double *diag = new double[n];
  for (int i = 0; i < n; i++) {
    diag[i] = a[i];
  }

  diag[0] -= alpha;
  diag[n - 1] -= beta * gamma / alpha;

  double *u = new double[n]; 
  double *z = new double[n]; 

  for (int i = 0; i < n; i++) {
    u[i] = z[i] = 0; 
  }
  u[0] = alpha;
  u[n - 1] = beta; 

  solveByThomasAlgorithm(diag, b, c, rhs, n, x); 

  solveByThomasAlgorithm(diag, b, c, u, n, z);

  double vy = x[0] + gamma * x[n - 1] / alpha;
  double vz = z[0] + gamma * z[n - 1] / alpha; 
  double r = vy / (1.0 + vz);
  //std::cout << "r = " << r << std::endl;
  //std::cout << std::endl;
  for (int i = 0; i < n; i++) {
    x[i] -= r * z[i];
  }

  delete[] diag; 
  delete[] u; 
  delete[] z; 
}

#endif 
