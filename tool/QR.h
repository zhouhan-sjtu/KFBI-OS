//  Copyright (c) Wenjun Ying, Department of Mathematics and Institute
//  of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

//  This file is free software; as a special exception the author gives
//  unlimited permission to copy and/or distribute it, with or without
//  modifications, as long as this notice is preserved.

//  This file is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY, to the extent permitted by law; without even 
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
//  PURPOSE.

#ifndef __QR_h_IS_INCLUDED__
#define __QR_h_IS_INCLUDED__ 

#include <iostream>
#include <cmath> 

#include "Norm.h" 
#include "Const.h"
#include "Product.h" 
//#include "InnerProduct.h" 

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void makeHouseholderTransform(const double *v, int n, 
                                     const double *x, double *z) 
{
  double s = computeInnerProduct(v, x, n);
  for (int i = 0; i < n; i++) {
    double w = s * v[i];
    z[i] = x[i] - (w + w);
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline bool computeHouseholderVector(const double *x, int m,
                                     double *y, double *v)
{
  for (int i = 0; i < m; i++) {
    y[i] = 0.0; 
  }

  double norm_x = computeVectorNorm(x, m); 
  if (norm_x / sqrt(m) < EPSILON12) {
    for (int i = 0; i < m; i++) {
      y[i] = x[i];
      v[i] = 0.0; 
    }
    return false;
  }

  if (x[0] > DBL_EPSILON) {
    y[0] = - norm_x;
  } else {
    y[0] = norm_x; 
  }

  for (int i = 0; i < m; i++) {
    v[i] = x[i]; 
  }
  v[0] -= y[0];

  double norm_v = sqrt(norm_x * norm_x - x[0] * x[0] + v[0] * v[0]);
  if (norm_v < EPSILON) {
    for (int i = 0; i < m; i++) {
      v[i] = 0.0;
    }
  } else {
    double r_norm_v = 1.0 / norm_v;
    for (int i = 0; i < m; i++) {
      v[i] *= r_norm_v; 
    }
  }

  return true; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <int N>
bool solveByQRdecomposition(double A[N][N], const double rhs[N],
                            double u[N], int n) 
{
  double b[N];
  for (int i = 0; i < n; i++) {
    b[i] = rhs[i]; 
  }

  double v[N]; 
  double x[N];
  double y[N]; 

  double R[N][N];

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      R[i][j] = A[i][j];
    }
  }

  for (int k = 0, m = n; k < n - 1; k++, m--) {
    for (int i = 0; i < m; i++) {
      x[i] = R[i + k][k];
    }
    bool status = computeHouseholderVector(x, m, y, v); 
    if (status) {
      for (int i = 0; i < m; i++) {
        R[i + k][k] = y[i]; 
      }

      for (int j = k + 1; j < n; j++) {
        for (int i = 0; i < m; i++) {
          x[i] = R[i + k][j]; 
        }
        makeHouseholderTransform(v, m, x, y);
        for (int i = 0; i < m; i++) {
          R[i + k][j] = y[i];
        }
      }

      makeHouseholderTransform(v, m, b + k, y);
      for (int i = 0; i < m; i++) {
        b[k + i] = y[i];
      }

    } else {
      return false;
    }
  }

  bool status = true; 
  if (fabs(R[n - 1][n - 1]) > EPSILON12) {
    for (int i = n - 1; i >= 0; i--) {
      double sum = 0.0; 
      for (int j = i + 1; j < n; j++) {
        sum += R[i][j] * u[j]; 
      }
      u[i] = (b[i] - sum) / R[i][i];
    }
  } else {
    status = false;
  }

  return status; 
}

#endif 
