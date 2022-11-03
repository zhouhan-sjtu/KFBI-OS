/*=============================================================================
*   
*   Filename : Interpolation.c
*   Creator : Han Zhou
*   Date : 09/06/21
*   Description : 
*
=============================================================================*/
   
#include <iostream>
#include <cmath>
#include <cfloat>
 
#include "MathTools.h"
#include "Interpolation.h"
#include "Const.h"
#include "QR.h"
//#include "CubicFormula.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double makeLinearInterpolation(const double t[2], const double f[2], double xi)
{
	double alpha = (t[1] - xi) / (t[1] - t[0]);
	double beta = 1.0 - alpha;

	return f[0] * alpha + f[1] * beta;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void makeLinearInterpolation(const double t[2], const double f[2], 
														 double xi, double g[2])
{
	double alpha = (t[1] - xi) / (t[1] - t[0]);
	double beta = 1.0 - alpha;

	g[0] = f[0] * alpha + f[1] * beta;
	g[1] = (f[1] - f[0]) / (t[1] - t[0]);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double makeCubicInterpolation(const double t[4], const double f[4], double xi)
{
  double diff[4];

  for (int j = 0; j < 4; j++) {
    diff[j] = f[j];
  }
  for (int j = 1; j < 4; j++) {
    for (int i = 3; i >= j; i--) {
      diff[i] = (diff[i] - diff[i - 1]) / (t[i] - t[i - j]); 
    }
  }

  double p = diff[3];
  for (int j = 2; j >= 0; j--) {
    p = diff[j] + p * (xi - t[j]);
  }

  return p;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void makeCubicInterpolation(const double t[4], const double f[4], 
                            double xi, double g[3])
{
  double diff[4]; 

  for (int j = 0; j < 4; j++) {
    diff[j] = f[j];
  }
  for (int j = 1; j < 4; j++) {
    for (int i = 3; i >= j; i--) {
      diff[i] = (diff[i] - diff[i - 1]) / (t[i] - t[i - j]); 
    }
  }

  double p = diff[3]; 
  for (int j = 2; j >= 0; j--) {
    p = diff[j] + p * (xi - t[j]);
  }

  g[0] = p;

  double c1 = (xi - t[0]) * (xi - t[1]) + (xi - t[1]) * (xi - t[2])
             + (xi - t[2])*(xi - t[0]);

  g[1] = diff[1] + diff[2] * (xi + xi - t[0] - t[1]) + diff[3] * c1;

  double c2 = 6.0 * xi - 2.0 * (t[0] + t[1] + t[2]); 

  g[2] = diff[2] + diff[2] + diff[3] * c2;
}



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double makeQuadraticInterpolation(const double t[3], const double f[3], double xi)
{
  double diff[3];

  for (int j = 0; j < 3; j++) {
    diff[j] = f[j];
  }
  for (int j = 1; j < 3; j++) {
    for (int i = 2; i >= j; i--) {
      diff[i] = (diff[i] - diff[i - 1]) / (t[i] - t[i - j]); 
    }
  }

  double p = diff[2];
  for (int j = 1; j >= 0; j--) {
    p = diff[j] + p * (xi - t[j]);
  }

  return p;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void makeQuadraticInterpolation(const double t[3], const double f[3], 
                            		double xi, double g[2])
{
  double diff[3]; 

  for (int j = 0; j < 3; j++) {
    diff[j] = f[j];
  }
  for (int j = 1; j < 3; j++) {
    for (int i = 2; i >= j; i--) {
      diff[i] = (diff[i] - diff[i - 1]) / (t[i] - t[i - j]); 
    }
  }

  double p = diff[2]; 
  for (int j = 1; j >= 0; j--) {
    p = diff[j] + p * (xi - t[j]);
  }

  g[0] = p;

  g[1] = diff[1] + diff[2] * (xi + xi - t[0] - t[1]);

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void makeQuadraticInterpolation2(const double t[3], const double f[3], 
                            		 double xi, double g[3])
{
	double mat[3][3], c[3];

	for(int r = 0; r < 3; r++){
		mat[r][2] = 1.0;
		mat[r][1] = t[r] - xi;
		mat[r][0] = mat[r][1] * mat[r][1];
	}

	solveByQRdecomposition< 3 > (mat, f, c, 3);

	g[0] = c[2];
	g[1] = c[1];
	g[2] = c[0] + c[0];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*bool findZeroOnQuadraticInterpolant(double &x, const double t[3],
															 			const double f[3])
{
	// normalization

	double x0 = x;

	double dt[3];
	for(int r = 0; r < 3; r++){
		dt[r] = fabs(t[r] - t[(r+1)%3]);
	}

	double h = min(dt[0], dt[1], dt[2]);
	double max_f = max(fabs(f[0]), fabs(f[1]), fabs(f[2]));

	double mat[3][3], c[3], b[3];
	for(int r = 0; r < 3; r++){
		mat[r][2] = 1.0;
  	mat[r][1] = (t[r] - x0) / h;	
		mat[r][0] = mat[r][1] * mat[r][1];	
		b[r] = f[r] / max_f;
	}

	bool status = solveByQRdecomposition< 3 > (mat, b, c, 3);
	if(!status) {
		std::cout << "Pass : " << __FILE__ << ", line : " << __LINE__ << std::endl;
		std::cout << "t = " << t[0] << ", " << t[1] << ", " << t[2] << std::endl;
		std::cout << "f = " << f[0] << ", " << f[1] << ", " << f[2] << std::endl;
		std::cout << "failed to solve by Qr decomposition." << std::endl;
		exit(1);
	}

	double delta = c[1] * c[1] - 4.0 * c[0] * c[2];

	if (delta < DBL_EPSILON) {
		double x1 = - c[1] / (c[0] + c[0]);
		x = x0 + x1 * h;
		return false;
	}

	double delta_r = sqrt(delta);

	double x1, x2;

	if(c[1] >= 0) {

		double tmp = c[1] + sqrt(delta);

		x1 = - tmp / (c[0] + c[0]);
		x2 = - (c[2] + c[2]) / tmp;

	} else {

		double tmp = sqrt(delta) - c[1];

		x1 = (c[2] + c[2]) / tmp;
		x2 = tmp / (c[0] + c[0]);

	}

	if (fabs(x1) < fabs(x2)) {
		x = x0 + x1 * h;
	} else {
		x = x0 + x2 * h;
	}

	return true;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool findZerosOnQuadraticInterpolant(const double t[3], const double f[3], 
																		 double z[2])
{
	// normalization

	double min_t = min(t[0], t[1], t[2]);
	double max_t = max(t[0], t[1], t[2]);
	double x0 = 0.5 * (max_t + min_t);

	double dt[3];

	for(int r = 0; r < 3; r++){
		dt[r] = fabs(t[r] - t[(r+1)%3]);
	}

	double h = min(min(dt[0], dt[1]), dt[2]);
	double max_f = max(fabs(f[0]), fabs(f[1]), fabs(f[2]));

	double mat[3][3], c[3], b[3];
	for(int r = 0; r < 3; r++){
		mat[r][2] = 1.0;
  	mat[r][1] = (t[r] - x0) / h;	
		mat[r][0] = mat[r][1] * mat[r][1];	
		b[r] = f[r] / max_f;
	}

	bool status = solveByQRdecomposition< 3 > (mat, b, c, 3);
	if(!status) {
		std::cout << "Pass : " << __FILE__ << ", line : " << __LINE__ << std::endl;
		std::cout << "t = " << t[0] << ", " << t[1] << ", " << t[2] << std::endl;
		std::cout << "f = " << f[0] << ", " << f[1] << ", " << f[2] << std::endl;
		std::cout << "failed to solve by Qr decomposition." << std::endl;
		exit(1);
	}

	double delta = c[1] * c[1] - 4.0 * c[0] * c[2];

	if (delta < DBL_EPSILON) {
		double x1 = - c[1] / (c[0] + c[0]);
		z[0] = z[1] = x0 + x1 * h;
		return false;
	}

	double delta_r = sqrt(delta);

	double x1, x2;

	if(c[1] >= 0) {

		double tmp = c[1] + sqrt(delta);

		x1 = - tmp / (c[0] + c[0]);
		x2 = - (c[2] + c[2]) / tmp;

	} else {

		double tmp = sqrt(delta) - c[1];

		x1 = (c[2] + c[2]) / tmp;
		x2 = tmp / (c[0] + c[0]);

	}

	z[0] = x0 + x1 * h;
	z[1] = x0 + x2 * h;

	return true;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool findZeroOnCubicInterpolant2(double &x, const double t[4],
															 	const double f[4])
{
	// normalization

	double mat[4][4], c[4], b[4];

	for(int r = 0; r < 4; r++){
		mat[r][3] = 1.0;
		mat[r][2] = t[r] - x;
		mat[r][1] = mat[r][2] * mat[r][2];
		mat[r][0] = mat[r][1] * mat[r][2];
		b[r] = f[r];
	}

	bool status = solveByQRdecomposition< 4 > (mat, b, c, 4);

	if(!status) {
		std::cout << "failed to solve by Qr decomposition." << std::endl;
		exit(1);
	}

	double xi = 0.0;

	xi = findCubicRoot(c[0], c[1], c[2], c[3], xi);

	x += xi;

	return true;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool findZeroOnCubicInterpolant(double &x, const double t[4],
															 	const double f[4])
{
	double x0 = x;

	// normalization

	double dt[6];
	for(int r = 0; r < 4; r++){
		dt[r] = fabs(t[r] - t[(r+1)%4]);
	}
	dt[4] = fabs(t[0] - t[2]);
	dt[5] = fabs(t[1] - t[3]);

	double h = min(min(dt[0], dt[1], dt[2], dt[3]), min(dt[4], dt[5]));
	double max_f = max(fabs(f[0]), fabs(f[1]), fabs(f[2]), fabs(f[3]));

	double mat[4][4], c[4], b[4];

	for(int r = 0; r < 4; r++){
		mat[r][3] = 1.0;
		mat[r][2] = (t[r] - x0) / h;
		mat[r][1] = mat[r][2] * mat[r][2];
		mat[r][0] = mat[r][1] * mat[r][2];
		b[r] = f[r] / max_f;
	}

	bool status = solveByQRdecomposition< 4 > (mat, b, c, 4);
	if(!status) {
		std::cout << "Pass : " << __FILE__ << ", line : " << __LINE__ << std::endl;
		std::cout << "t = " << t[0] << ", " << t[1] << ", " << t[2] << ", " << t[3] << std::endl;
		std::cout << "f = " << f[0] << ", " << f[1] << ", " << f[2] << ", " << f[3] << std::endl;
		std::cout << "failed to solve by Qr decomposition." << std::endl;
		exit(1);
	}

	double xi = 0.0;
	double F = xi * (xi * (xi * c[0] + c[1]) + c[2]) + c[3];

	const double atol = 1.0e-14;
	const double rtol = 1.0e-12;

	if (fabs(F) < atol) {
		return true;
	}

	double tol = atol + fabs(F) * rtol;

	const int max_itr_num = 50;
	int itr_num = 0;

	do {

		double dF = xi * (xi * 3.0 * c[0] + 2.0 * c[1]) + c[2];

		xi -= F / dF;

		F = xi * (xi * (xi * c[0] + c[1]) + c[2]) + c[3];

		itr_num++;

		//std::cout << "itr_num = " << itr_num << ", res_norm = " << fabs(F) << std::endl;

	} while (fabs(F) > tol && itr_num < max_itr_num);

	if(itr_num == max_itr_num) {

		x = x0;
		bool status = findZeroOnQuadraticInterpolant(x, t, f);
		x = x0;
		if (!status) {
			status = findZeroOnQuadraticInterpolant(x, t+1, f+1);
		}
		return status;
  }

	x = x0 + h * xi;

	return true;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double computeDistanceToParabola0(double t_min, double t_max, const double c[3], 
																 double x0, double y0)
{

	const int n = 400;
	double diam = t_max - t_min;
	double ds = 3.0 * diam / n;

	double min_d = DBL_MAX;

	double x = t_min - diam;

	for(int i = 0; i <= n; i++){
		double y = x * (c[0] * x + c[1]) + c[2];

		double dis = distance(x0, y0, x, y);

		if (dis < min_d) {
			min_d = dis;
		}

		x += ds;
	}

	return min_d;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double computeDistanceToParabola(const double t[3], const double f[3], 
																 double x0, double y0)
{

	double min_t = min(t[0], t[1], t[2]);
	double max_t = max(t[0], t[1], t[2]);
	double min_f = min(f[0], f[1], f[2]);
	double max_f = max(f[0], f[1], f[2]);

	double x_ctr = 0.5 * (max_t + min_t);
	double y_ctr = 0.5 * (max_f + min_f);

	double dt[3];
	for(int r = 0; r < 3; r++){
		dt[r] = fabs(t[r] - t[(r+1)%3]);
	}

	double h = min(min(dt[0], dt[1]), dt[2]);

	double mat[3][3], b[3], c[3];
	for(int r = 0; r < 3; r++){
		mat[r][2] = 1.0;
  	mat[r][1] = (t[r] - x_ctr) / h;	
		mat[r][0] = mat[r][1] * mat[r][1];	
		b[r] = (f[r] - y_ctr) / h;
	}

	double x1 = (x0 - x_ctr) / h;
	double y1 = (y0 - y_ctr) / h;

	bool status = solveByQRdecomposition< 3 > (mat, b, c, 3);
	if(!status) {
		std::cout << "failed to solve by Qr decomposition." << std::endl;
		exit(1);
	}

	double x2 = 0.0;
	double y2 = x2 * (c[0] * x2 + c[1]) + c[2];

	double F = 2.0 * ((y2 - y1) * (2.0 * c[0] * x2 + c[1]) + (x2 - x1));

	const double atol = 1.0e-15;
	const double rtol = 1.0e-12;

	if (fabs(F) < atol) {
		double dx = x2 - x1;
		double dy = y2 - y1;
		return sqrt(dx*dx + dy*dy) * h;
	}

	double tol = atol + rtol * fabs(F);

	const int max_itr_num = 50;
	int itr_num = 0;

	do {

		double dF = 2.0*(2.0*c[0]*x2+c[1])*(2.0*c[0]*x2+c[1])+4.0*c[0]*(y2-y1)+2.0;

		x2 -= F / dF;
		y2 = x2 * (c[0] * x2 + c[1]) + c[2];

		F = 2.0 * ((y2 - y1) * (2.0 * c[0] * x2 + c[1]) + (x2 - x1));

		itr_num++;

		//std::cout << "itr_num = " << itr_num << ", res_norm = " << fabs(F) << std::endl;

	} while (fabs(F) > tol && itr_num < max_itr_num);

	if (itr_num == max_itr_num){
		status = false;
#ifdef DEBUG
		std::cout << "Pass : " << __FILE__ << ", line : " << __LINE__ << std::endl;
		std::cout << "Newton's method failed to converge." << std::endl;
#endif
		exit(1);
	}

	assert(status);

	double dx = x2 - x1;
	double dy = y2 - y1;
	double dist = sqrt(dx*dx + dy*dy) * h;

	return dist;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double computeDistanceToCubicInterpolant(const double t[4], const double f[4], 
																 		double x0, double y0)
{
	double mat[4][4], b[4], c[4];
	for(int r = 0; r < 4; r++){
		mat[r][3] = 1.0;
		mat[r][2] = t[r];
		mat[r][1] = t[r] * t[r];
		mat[r][0] = t[r] * mat[r][1];
		b[r] = f[r];
	}

	bool status = solveByQRdecomposition< 4 > (mat, b, c, 4);
	if(!status) {
		std::cout << "failed to solve by Qr decomposition." << std::endl;
		exit(1);
	}

	double min_t = min(min(t[0], t[1]), t[2]);
	double max_t = max(max(t[0], t[1]), t[2]);
	double diam = fabs(max_t - min_t);

	double x = t[1];

	double y = x * (x * (x * c[0] + c[1]) + c[2]) + c[3];
	double yp = x * (x * 3.0 * c[0] + 2.0 * c[1]) + c[2];
	double ypp = 6.0 * c[0] * x + 2.0 * c[1];

	double F = 2.0 * (y - y0) * yp + 2.0 * (x - x0);

	const double atol = 1.0e-15;
	const double rtol = 1.0e-12;

	if (fabs(F) < atol) {
		double dx = x - x0;
		double dy = y - y0;
		return sqrt(dx*dx + dy*dy);
	}

	double tol = atol + fabs(F) * rtol;
	const int max_itr_num = 20;
	int itr_num = 0;

	do {

		double dF = 2.0 * (yp * yp + (y - y0) * ypp) + 2.0;

		x -= F / dF;

		y = x * (x * (x * c[0] + c[1]) + c[2]) + c[3];
		yp = x * (x * 3.0 * c[0] + 2.0 * c[1]) + c[2];
		ypp = 6.0 * c[0] * x + 2.0 * c[1];

		F = 2.0 * (y - y0) * yp + 2.0 * (x - x0);

		itr_num++;

		//std::cout << "itr_num = " << itr_num << ", res_norm = " << fabs(F) << std::endl;

	} while (fabs(F) > tol && itr_num < max_itr_num);

	if (itr_num == max_itr_num || x - max_t > diam || min_t - x > diam){
		status = false;
#ifdef DEBUG
		std::cout << "Pass : " << __FILE__ << ", line : " << __LINE__ << std::endl;
		std::cout << "Newton's method failed to converge." << std::endl;
#endif
	}

	double distc = 0.0;

	if (status) {


		double dx = x - x0;
		double dy = y - y0;
		distc = sqrt(dx*dx + dy*dy);

	} else {
		distc = computeDistanceToParabola0(min_t, max_t, c, x0, y0);
	}

	return distc;
}*/

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
