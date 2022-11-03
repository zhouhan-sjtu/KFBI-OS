/*=============================================================================
*   
*   Filename : MathTools.h
*   Creator : Han Zhou
*   Date : 11/02/21
*   Description : 
*
=============================================================================*/

#ifndef _MATHTOOLS_H
#define _MATHTOOLS_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <cstring>

#include "Const.h"
#include "Random.h"
#include "Variables.h"

#ifdef USE_MATHGL
#include "MathGL.h"
#endif

#ifdef USE_OPENGL
#include "NewMathGL.h"
#include "MathGLUT2d.h"
#endif

#define printBugInfo(expr) \
std::cout << "Bug in file : " << __FILE__ << ", line : " \
					<< __LINE__  << "\n" << expr << std::endl;


#ifdef USE_MATHGL
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void mglPlotInteriorNodeDataByIsolines(const VectorXd &x0, const VectorXd &y0,
												const MatrixXb &side, const MatrixXd &v, int m, int n,
												int isovalue_num)
{
	double *x = new double[m+1];
	double *y = new double[n+1];

	for(int i = 0; i <= m; i++){
		x[i] = x0[i];
	}
	for(int i = 0; i <= n; i++){
		y[i] = y0[i];
	}

	double **vect = new double*[m+1];
	bool **interior = new bool*[m+1];

	for(int i = 0; i <= m; i++){
		vect[i] = new double[n+1];
		interior[i] = new bool[n+1];

		for(int j = 0; j <= n; j++){
			vect[i][j] = v[i][j];
			interior[i][j] = side[i][j];
		}
	}

	//mglPlotInteriorNodeDataByIsolines(x, y, interior, vect, m, n, isovalue_num);
	//mglPlotNodeDataByIsolines(x, y, vect, m, n, isovalue_num);
	mglPlotNodeDataByIsoLines(x, y, vect, m, isovalue_num);

	for(int i = 0; i <= m; i++){
		delete[] vect[i];
		delete[] interior[i];
	}
	delete[] vect;
	delete[] interior;
	delete[] x;
	delete[] y;
}
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void intToString(int ind, char *str)
{
  if (ind == 0) {
    str[0] = '0'; 
    str[1] = '\0';
    return; 
  }

  if (ind < 0) {
    str[0] = '-'; 
    intToString(- ind, str + 1); 
    return;
  }

  int width = 0; 
  int cpy = ind;
  while (cpy) {
    width++;
    cpy /= 10;
  }

  for (int i = width - 1; i >= 0; i--) {
    str[i] = '0' + ind%10;
    ind /= 10;
  }

  str[width] = '\0'; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void floatToString(double t, char t_str[]) 
{
  if (t < 0) {
    t_str[0] = '-'; 
    floatToString(-t, t_str + 1); 
    return;
  }

  int p = floor(t + EPSILON12);
  double f = t - p; 

  intToString(p, t_str); 

  if (fabs(f) < EPSILON) {
    return; 
  }

  int len = strlen(t_str); 

  t_str[len] = '.'; 
  len++;

  int d = 0; 
  int max_d = 4;
  do {
    d++;
    f *= 10;
    p = floor(f + EPSILON12); 
    f = f - p;
    t_str[len++] = '0' + p;
  } while ((d < max_d) && (fabs(f) > EPSILON)); 

  t_str[len] = '\0';
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void showCurrentTime2d(const double low[0], const double high[2], double t) 
{

  char cur_time[32]; 
  floatToString(t, cur_time);

  char str[80] = {"t = "};
  strcat(str, cur_time);

	//double word_height = 0.1; 
	double word_height = 0.01 * (high[0] - low[0]); 

  double y = 0.825 + 0.6 * word_height;
  double x = - 0.9;

  double cx = 0.5 * (high[0] + low[0]);
  double cy = 0.5 * (high[1] + low[1]);
  double rx = 0.5 * (high[0] - low[0]);
  double ry = 0.5 * (high[1] - low[1]);

  double xi = cx + rx * x;
  double eta = cy + ry * y; 

#ifdef USE_OPENGL
	mglSetColor(0.1, 0.1, 0.1);
	mglPlotString(xi, eta, str);
#endif
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void splitNeumannBoundaryValue(double tx, double ty, double ut, 
															 				double nx, double ny, double un,
															 				double &ux, double &uy)
{
	double A[2][2], b[2];

	A[0][0] = tx;
	A[0][1] = ty;
	b[0] = ut;

	A[1][0] = nx;
	A[1][1] = ny;
	b[1] = un;

	double det = A[0][0] * A[1][1] - A[1][0] * A[0][1];

	ux = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
	uy = (A[0][0] * b[1] - A[1][0] * b[0]) / det;

}
 
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void computeDrawingArea(const double lo[2], const double hi[2],	
															 double low[2], double high[2], double ratio)
{
	double c[2], r[2];

	for(int i = 0; i < 2; i++){
		c[i] = 0.5 * (lo[i] + hi[i]);
		r[i] = 0.5 * (hi[i] - lo[i]);
	}
	for(int i = 0; i < 2; i++){
		low[i] = c[i] - ratio * r[i];
		high[i] = c[i] + ratio * r[i];
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void rotateToRange(double &t)
{
	int k = static_cast<int>(t / M_2PI);
	t -= k * M_2PI;

	while (t < 0.0) {
		t += M_2PI;
	}
	while (t >= M_2PI) {
		t -= M_2PI;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void normalize(double n[2])
{
	double r_norm = 1.0 / sqrt(n[0] * n[0] + n[1] * n[1]);
	n[0] *= r_norm;
	n[1] *= r_norm;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double distance(double x1, double y1, double x2, double y2) 
{
  double dx = x2 - x1; 
  double dy = y2 - y1; 
	return sqrt(dx * dx + dy * dy);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
T max(T a, T b)
{
	return a > b ? a : b;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
T min(T a, T b)
{
	return a < b ? a : b;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
T Atan(T x, T y) 
{
  T theta = 0.0;
  if (fabs(x) < DBL_EPSILON) {
    if (y > DBL_EPSILON) {
      theta = M_PI_2;
    } else if (y + DBL_EPSILON < 0) {
      theta = M_PI * 1.5; 
    } else {
      theta = 0; 
    }
  } else {
    theta = atan(y / x);
    if (x < 0) {
      theta += M_PI;
    } else {
      if (y + DBL_EPSILON < 0) {
        theta += M_2PI; 
      }
    }
  }
  return(theta);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline int Log2(int m)
{
  int k = 0;
  while (m > 1) {
    m >>= 1;
    k++; 
  }
  return k;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline int Power2(int m)
{
  int k = 1;
  for (int i = 0; i < m; i++) {
    k <<= 1; 
  }

  return k;
}


//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

