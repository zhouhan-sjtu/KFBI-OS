// Copyright (c) Wenjun Ying, School of Mathematical Sciences and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is free software; as a special exception the author gives
// unlimited permission to copy and/or distribute it, with or without
// modifications, as long as this notice is preserved.

// This file is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#ifndef __ClosedSplineCurve_H_IS_INCLUDED__
#define __ClosedSplineCurve_H_IS_INCLUDED__

#include "Variables.h"
#include "ParametricCurve.H"

class ClosedSplineCurve : public ParametricCurve
{
  public :

    ClosedSplineCurve(void); 
    virtual ~ClosedSplineCurve(void); 
    ClosedSplineCurve(const ClosedSplineCurve &);
    ClosedSplineCurve(const double (*p)[2], int n); 
    ClosedSplineCurve(const VectorX2d &p, int n); 
    ClosedSplineCurve(const double *x, const double *y, int n);

    const ClosedSplineCurve &operator = (const ClosedSplineCurve &);

    void getPoint(double theta, double &x, double &y) const;

    void getX(double theta, double &x, double &dx, double &ddx) const;
    void getY(double theta, double &y, double &dy, double &ddy) const; 

    void getPoint(double theta, double &x, double &y, double &tx,
                  double &ty, double &nx, double &ny) const; 
    void getPoint2(double theta, double &x, double &y, double &dx,
                   double &dy, double &ddx, double &ddy) const;

    void initialize(const double x[], const double y[], int n); 
    void define(const double x[], const double y[], int n); 
    void reset(const double x[], const double y[], int n); 

  public :

    double computeArcLength(void) const; 

  private :
    int LU(const double a[], const double b[], const double c[],
           const double w[], int m, double u[]) const; 

    void computeSecondDerivative(const double z[], int n, 
                                 double h, double M[]);

    double S(const double z[], const double M[], 
             int n, double h, double t) const;

    double DS1(const double z[], const double M[],
               int n, double h, double t) const; 
    double DS1(const double z[], const double M[],
               int n, double h, int i) const;
    double DS1(const double z[], const double M[],
               int n, double h, int i, double xi) const;
    double DS2(const double M[], int n, double h, double t) const;

    void allocatePointers(int n);
    void allocatePointers(void); 
    void initDataMembers(void);
    void freePointers(void);

  private :
    int _n;
    int _n1; 

    double _h; 

    double *_x; 
    double *_y;

    double *_ddx; 
    double *_ddy; 
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ClosedSplineCurve::ClosedSplineCurve(void) : ParametricCurve()
{
  initDataMembers();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ClosedSplineCurve::~ClosedSplineCurve(void)
{
  freePointers();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void ClosedSplineCurve::allocatePointers(void) 
{
  _x = new double[_n1]; 
  for (unsigned int i = 0; i < _n1; i++) {
    _x[i] = 0.0;
  }

  _y = new double[_n1]; 
  for (unsigned int i = 0; i < _n1; i++) {
    _y[i] = 0.0;
  }

  _ddx = new double[_n1];
  for (unsigned int i = 0; i < _n1; i++) {
    _ddx[i] = 0.0; 
  }

  _ddy = new double[_n1]; 
  for (unsigned int i = 0; i < _n1; i++) {
    _ddy[i] = 0.0;
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void ClosedSplineCurve::allocatePointers(int n) 
{
  _x = new double[_n1];
  for (unsigned int i = 0; i <= n; i++) {
    _x[i] = 0.0;
  }

  _y = new double[_n1]; 
  for (unsigned int i = 0; i <= n; i++) {
    _y[i] = 0.0; 
  }

  _ddx = new double[_n1]; 
  for (unsigned int i = 0; i <= n; i++) {
    _ddx[i] = 0.0; 
  }

  _ddy = new double[_n1];
  for (unsigned int i = 0; i <= n; i++) {
    _ddy[i] = 0.0;
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ClosedSplineCurve::ClosedSplineCurve(const ClosedSplineCurve &c) : ParametricCurve(c)
{
  _h = c._h;
  _n = c._n; 
  _n1 = c._n1; 

  if (c._x != 0) {

    _x = new double[_n1]; 
    for (unsigned int i = 0; i <= _n; i++) {
      _x[i] = c._x[i]; 
    }

    _y = new double[_n1];
    for (unsigned int i = 0; i <= _n; i++) {
      _y[i] = c._y[i];
    }

    _ddx = new double[_n1]; 
    for (unsigned int i = 0; i <= _n; i++) {
      _ddx[i] = c._ddx[i];
    }

    _ddy = new double[_n1];
    for (unsigned int i = 0; i <= _n; i++) {
      _ddy[i] = c._ddy[i];
    }

  } else {

    _x = 0;
    _y = 0;
    _ddx = 0;
    _ddy = 0;
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ClosedSplineCurve::ClosedSplineCurve(const double (*p)[2], int n) : ParametricCurve()
{
  _n = n; 
  _h = 1.0; 
  _n1 = n + 1;

  allocatePointers(n);

  for (unsigned int i = 0; i < n; i++) {
    _x[i] = p[i][0]; 
    _y[i] = p[i][1]; 
  }
  _x[n] = _x[0];
  _y[n] = _y[0]; 

  _h = 1.0;

  computeSecondDerivative(_x, _n, _h, _ddx); 
  computeSecondDerivative(_y, _n, _h, _ddy); 

  _ddx[_n] = _ddx[0]; 
  _ddy[_n] = _ddy[0]; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ClosedSplineCurve::ClosedSplineCurve(const VectorX2d &p, int n) : ParametricCurve()
{
  _n = n; 
  _h = 1.0; 
  _n1 = n + 1;

  allocatePointers(n);

  for (unsigned int i = 0; i < n; i++) {
    _x[i] = p[i][0]; 
    _y[i] = p[i][1]; 
  }
  _x[n] = _x[0];
  _y[n] = _y[0]; 

  _h = 1.0;

  computeSecondDerivative(_x, _n, _h, _ddx); 
  computeSecondDerivative(_y, _n, _h, _ddy); 

  _ddx[_n] = _ddx[0]; 
  _ddy[_n] = _ddy[0]; 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ClosedSplineCurve::
ClosedSplineCurve(const double *x, const double *y, int n) : ParametricCurve()
{
  _n = n;
  _n1 = n + 1;

  allocatePointers(n);

  for (unsigned int i = 0; i < n; i++) {
    _x[i] = x[i]; 
    _y[i] = y[i];
  }
  _x[n] = _x[0];
  _y[n] = _y[0]; 

  _h = 1.0;

  computeSecondDerivative(_x, _n, _h, _ddx);
  computeSecondDerivative(_y, _n, _h, _ddy);

  _ddx[_n] = _ddx[0];
  _ddy[_n] = _ddy[0]; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const ClosedSplineCurve &ClosedSplineCurve::
operator = (const ClosedSplineCurve &c)
{
  if (this != &c) {
    if (_n == c._n) {
    } else {
      freePointers(); 
      _n = c._n; 
      _n1 = c._n1; 
      allocatePointers(_n);
    }
    for (unsigned int i = 0; i <= _n; i++) {
      _x[i] = c._x[i];  _ddx[i] = c._ddx[i]; 
      _y[i] = c._y[i];  _ddy[i] = c._ddy[i]; 
    }
  }
  return *this;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void ClosedSplineCurve::initDataMembers(void)
{
  _n = 0;
  _x = 0;
  _y = 0; 
  _n1 = 0; 
  _ddx = 0; 
  _ddy = 0; 
  _h = 0.0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void ClosedSplineCurve::freePointers(void)
{
  _n = 0; 
  _n1 = 0; 
  _h = 0.0;
  if (0 != _x) {
    delete[] _x;
    _x = 0;
  }
  if (0 != _y) {
    delete[] _y;
    _y = 0;
  }
  if (0 != _ddx) {
    delete[] _ddx;
    _ddx = 0;
  }
  if (0 != _ddy) {
    delete[] _ddy; 
    _ddy = 0; 
  }
}

#endif
