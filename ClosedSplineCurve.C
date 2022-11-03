// Copyright (c) Wenjun Ying, School of Mathematical Sciences and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is free software; as a special exception the author gives
// unlimited permission to copy and/or distribute it, with or without
// modifications, as long as this notice is preserved.

// This file is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#include <cmath> 
#include <iostream> 

#include "Const.h" 
#include "Random.h" 
#include "Pointers.h" 
#include "MathTools.h" 

#include "ClosedSplineCurve.H" 

#define printBugInfo(expr) \
std::cout << "Bug in file : " << __FILE__ << ", line : " \
					<< __LINE__  << "\n" << expr << std::endl;

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
template <typename T>
void solveByThomasAlgorithm(const T *a, const T *b, const T *c,
                            const T *f, int n, T *x)
{
  int n_1 = n - 1; 

  T *alpha = new T[n];
  T *beta = new T[n_1];
  T *gamma = new T[n_1]; 

  alpha[0] = a[0]; 

  for (int i = 1; i < n; i++) {
    beta[i - 1] = b[i - 1];
    gamma[i - 1] = c[i - 1] / alpha[i - 1];
    alpha[i] = a[i] - gamma[i - 1] * beta[i - 1];
  }

  T *y = new T[n]; 

  y[0] = f[0];
  for (int i = 1; i < n; i++) {
    y[i] = f[i] - gamma[i - 1] * y[i - 1]; 
  }

  x[n_1] = y[n_1] / alpha[n_1];
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

template <typename T>
void solveByThomasAlgorithm(const T a[], const T b[], 
                            const T c[], T beta, T gamma,
                            const T rhs[], int n, T x[]) 
{
  const T alpha = - a[0];

  T *diag = new T[n]; 
  for (int i = 0; i < n; i++) {
    diag[i] = a[i];
  }

  diag[0] -= alpha; 
  diag[n - 1] -= beta * gamma / alpha;

  T *u = new T[n]; 
  T *z = new T[n]; 

  for (int i = 0; i < n; i++) {
    u[i] = z[i] = 0; 
  }
  u[0] = alpha; 
  u[n - 1] = beta; 

  solveByThomasAlgorithm(diag, b, c, rhs, n, x); 

  solveByThomasAlgorithm(diag, b, c, u, n, z);

  T vy = x[0] + gamma * x[n - 1] / alpha; 
  T vz = z[0] + gamma * z[n - 1] / alpha; 
  T r = vy / (1.0 + vz); 
  for (int i = 0; i < n; i++) {
    x[i] -= r * z[i];
  }

  delete[] diag;
  delete[] u; 
  delete[] z;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ClosedSplineCurve::define(const double x[], const double y[], int n)
{
  if ((_n > 0) && (_n != n)) {
    freePointers(); 
    _n = n;
    _n1 = n + 1;
    allocatePointers(n);
  } else if (_n == 0) {
    _n = n; 
    _n1 = n + 1;
    allocatePointers(n);
  }

  _h = 1.0; 

  for (int i = 0; i < _n; i++) {
    _x[i] = x[i];
    _y[i] = y[i]; 
  }
  _x[_n] = _x[0]; 
  _y[_n] = _y[0];

  computeSecondDerivative(_x, _n, _h, _ddx); 
  computeSecondDerivative(_y, _n, _h, _ddy); 

  _ddx[_n] = _ddx[0]; 
  _ddy[_n] = _ddy[0];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ClosedSplineCurve::reset(const double x[], const double y[], int n) 
{
  if ((_n > 0) && (_n != n)) {
    freePointers();
    _n = n; 
    _n1 = n + 1;
    allocatePointers(n);
  } else if (_n == 0) {
    _n = n;
    _n1 = n + 1; 
    allocatePointers(n);
  }

  _h = 1.0; 

  for (int i = 0; i < _n; i++) {
    _x[i] = x[i]; 
    _y[i] = y[i]; 
  }
  _x[_n] = _x[0];
  _y[_n] = _y[0];

  computeSecondDerivative(_x, _n, _h, _ddx);
  computeSecondDerivative(_y, _n, _h, _ddy); 

  _ddx[_n] = _ddx[0]; 
  _ddy[_n] = _ddy[0];
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ClosedSplineCurve::computeSecondDerivative(const double z[], int n, 
                                                double h, double M[]) 
{
  int n_1 = n - 1; 

  double *a = allocateVector(n); 
  double *b = allocateVector(n); 
  double *c = allocateVector(n);

  double one_sixth = 1.0 / 6.0; 
  double two_thirds = 2.0 / 3.0;

  for (int i = 0; i < n; i++) {
    a[i] = two_thirds; 
    b[i] = one_sixth; 
    c[i] = one_sixth;
  }

  double h2 = h * h;
  double r_h2 = 1.0 / h2;

  double *rhs = allocateVector(n);

  for (int i = 1; i < n_1; i++) {
    rhs[i] = (z[i + 1] + z[i - 1] - z[i] - z[i]) * r_h2; 
  }
  rhs[n_1] = (z[0] + z[n - 2] - z[n_1] - z[n_1]) * r_h2; 
  rhs[0] = (z[1] + z[n_1] - z[0] - z[0]) * r_h2; 

  solveByThomasAlgorithm(a, b, c, b[n_1], c[n_1], rhs, n, M); 

  delete[] rhs;

  delete[] a; 
  delete[] b;
  delete[] c; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ClosedSplineCurve::S(const double z[], const double M[],
                            int n, double h, double t) const
{
  double period = n * h; 
  while ((t + EPSILON) < 0.0) {
    t += period; 
  }

  while (t > (period - EPSILON)) {
    t -= period; 
  }
  int i = static_cast<int>(t / h);
  if (i >= n) {
    t -= period; 
    i -= n; 
  }
  if (i < 0) {
    t += period;
    i += n; 
    std::cout << "i = " << i << std::endl;
    std::cout << "t = " << t << std::endl; 
  }

  double ti = i * h; 

  double r2 = t - ti; 
  double r1 = ti + h - t; 
  double s1 = (r1 * r1 * r1 * M[i] + r2 * r2 * r2 * M[i + 1]) / (6.0 * h);
  double s2 = (r1 * z[i] + r2 * z[i + 1]) / h; 
  double s3 = (r1 * M[i] + r2 * M[i + 1]) * h / 6.0; 
  double r = s1 + s2 - s3;

  return r; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ClosedSplineCurve::DS1(const double z[], const double M[], 
                              int n, double h, int i, double xi) const 
{
  double r2 = xi; 
  double r1 = h - xi;

  double s1 = (r2 * r2 * M[i + 1] - r1 * r1 * M[i]) / (2.0 * h);
  double s2 = (z[i + 1] - z[i]) / h;
  double s3 = (M[i + 1] - M[i]) * h / 6.0;
  double r = s1 + s2 - s3; 

  return r; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ClosedSplineCurve::DS1(const double z[], const double M[],
                              int n, double h, double t) const 
{
  double period = n * h;

  while ((t + EPSILON) < 0.0) {
    t += period; 
  }
  while (t > (period - EPSILON)) {
    t -= period;
  }
  int i = static_cast<int>(t / h); 
  if (i >= n) {
    t -= period; 
    i -= n;
  }
  if (i < 0) {
    t += period;
    i += n; 
  }

  double ti = i * h;

  double r2 = t - ti; 
  double r1 = ti + h - t; 

  double s1 = (r2 * r2 * M[i + 1] - r1 * r1 * M[i]) / (2.0 * h);
  double s2 = (z[i + 1] - z[i]) / h;
  double s3 = (M[i + 1] - M[i]) * h / 6.0; 
  double r = s1 + s2 - s3;

  return r;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ClosedSplineCurve::DS1(const double z[], const double M[], 
                              int n, double h, int i) const
{
  double s1 = - 0.5 * M[i] * h;
  double s2 = (z[i + 1] - z[i]) / h;
  double s3 = (M[i + 1] - M[i]) * h / 6.0;
  double r = s1 + s2 - s3;
  return r; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ClosedSplineCurve::DS2(const double M[], int n, double h, double t) const
{
  double period = n * h;

  while ((t + EPSILON) < 0.0) {
    t += period;
  }
  while (t > (period - EPSILON)) {
    t -= period; 
  }
  int i = static_cast<int>(t / h); 
  if (i >= n) {
    t -= period; 
    i -= n; 
  }
  if (i < 0) {
    t += period; 
    i += n;
  }

  double ti = i * h; 

  double r2 = t - ti;
  double r1 = ti + h - t; 
  double r = (r1 * M[i] + r2 * M[i + 1]) / h;

  return r; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ClosedSplineCurve::computeArcLength(void) const
{
  double h = _h; 

  const double m = 0.5 * h; 

  const double quadrature[5][2] = {
    {-0.906179845939, 0.236926885056},
    {-0.538469310106, 0.478628670499},
    {-0.000000000000, 0.568888888889},
    { 0.538469310106, 0.478628670499}, 
    { 0.906179845939, 0.236926885056}, 
  }; 

  double dx = 0.0;
  double dy = 0.0; 
  double len = 0.0; 
  for (unsigned int k = 0; k < _n; k++) {
    for (int i = 0; i < 5; i++) {
      double t = m * (1.0 + quadrature[i][0]); 
      dx = DS1(_x, _ddx, _n, _h, k, t);
      dy = DS1(_y, _ddy, _n, _h, k, t);
      len += quadrature[i][1] * sqrt(dx * dx + dy * dy);
    }
  }
  len *= m;

  return len; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ClosedSplineCurve::getPoint(double theta, double &x, double &y) const
{
  double s = (theta / M_2PI) * (_n * _h);

  x = S(_x, _ddx, _n, _h, s); 
  y = S(_y, _ddy, _n, _h, s);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ClosedSplineCurve::getPoint(double theta, double &x, double &y, 
                                 double &tx, double &ty, double &nx,
                                 double &ny) const
{
  double fac = _n * _h / M_2PI;

  double s = fac * theta;

  x = S(_x, _ddx, _n, _h, s);
  y = S(_y, _ddy, _n, _h, s); 

  double p = DS1(_x, _ddx, _n, _h, s);
  double q = DS1(_y, _ddy, _n, _h, s);

  tx = fac * p;  ty = fac * q; 

  double len = sqrt(tx * tx + ty * ty); 
  nx = ty / len;  ny = - tx / len; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ClosedSplineCurve::getX(double theta, double &x,
                             double &dx, double &ddx) const
{
  double fac = _n * _h / M_2PI;

  double s = fac * theta;

  x = S(_x, _ddx, _n, _h, s);

  double p = DS1(_x, _ddx, _n, _h, s); 

  dx = fac * p;

  p = DS2(_ddx, _n, _h, s); 

  double fac2 = fac * fac; 

  ddx = fac2 * p;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ClosedSplineCurve::getY(double theta, double &y,
                             double &dy, double &ddy) const 
{
  double fac = _n * _h / M_2PI; 

  double s = fac * theta; 

  y = S(_y, _ddy, _n, _h, s);

  double q = DS1(_y, _ddy, _n, _h, s); 

  dy = fac * q;

  q = DS2(_ddy, _n, _h, s);

  double fac2 = fac * fac;

  ddy = fac2 * q; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ClosedSplineCurve::getPoint2(double theta, double &x, double &y,
                                  double &dx, double &dy, double &ddx, 
                                  double &ddy) const 
{
  getX(theta, x, dx, ddx);
  getY(theta, y, dy, ddy); 
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
