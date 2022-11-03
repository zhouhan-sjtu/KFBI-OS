// Copyright (c) Wenjun Ying 2018, School of Mathematical Sciences and 
// Institute of Natural Sciences, Shanghai Jiao Tong University, Shanghai 
// 200240.

// This file is free software; as a special exception the author gives
// unlimited permission to copy and/or distribute it, with or without
// modifications, as long as this notice is preserved.

// This file is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY, to the extent permitted by law; without even 
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.

#include <cstdlib>
#include <cstdio> 
#include <cmath>

//#include <eigen3/Eigen/Dense>

#include "QR.h" 
#include "MathTools.h" 
#include "Variables.h" 

/**
 ******************************************************************************
 * solve linear system through LU decomposition with pivoting.                *
 *                                                                            *
 * Factorize PA = LU with pivoting:                                           *
 *   The lower and upper triangular matrices are still stored in the original *
 *   matrix and the permutation matrix "P" is stored in the vector "int *p".  *
 ******************************************************************************
 **/
int LUfactorize6(double A[6][6], int p[6])
{
  int i, j, k; 
  double m, ajj, lij;

  for (j = 0; j < 6; j++) {
    p[j] = j; 
  }

  for (j = 0; j < 6; j++) {

    k = j;
    m = A[p[j]][j]; 
    for (i = j + 1; i < 6; i++) {
      if (fabs(A[p[i]][j]) > fabs(m)) {
        m = A[p[i]][j];
        k = i;
      }
    }

    if (k != j) {
      int temp = p[j]; 
      p[j] = p[k];
      p[k] = temp; 
    }

    ajj = A[p[j]][j]; 
    if (fabs(ajj) < 1.0E-15) {
      printf("failure in the LU factorization.\n"); 
      return (-1); 
    }

    for (i = j + 1; i < 6; i++) {
      lij = A[p[i]][j] / ajj; 
      A[p[i]][j] = lij; 
      for (k = j + 1; k < 6; k++) {
        A[p[i]][k] -= lij * A[p[j]][k];
      }
    }
  }

  return 0;
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
void LUsolve_internal6(double A[6][6], const int p[6], 
                       const double b[6], double x[6]) 
{
  int i, j, n; 
  double rowsum;

  x[0] = b[p[0]]; 
  for (i = 1; i < 6; i++) {
    x[i] = b[p[i]]; 
    rowsum = 0.0;
    for (j = 0; j < i; j++) {
      rowsum += A[p[i]][j] * x[j];
    }
    x[i] -= rowsum; 
  }

  n = 6; 

  x[n - 1] = x[n - 1] / A[p[n - 1]][n - 1]; 
  for (i = n - 2; i >= 0; i--) {
    rowsum = 0.0;
    for (j = n - 1; j > i; j--) {
      rowsum += A[p[i]][j] * x[j]; 
    }
    x[i] = (x[i] - rowsum) / A[p[i]][i]; 
  }
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
bool LUsolve6(double A[6][6], const double b[6], double x[6])
{
  double B[6][6]; 
  int i, j, p[6], status; 

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      B[i][j] = A[i][j];
    }
  }

  status = LUfactorize6(B, p); 

  if (status != 0) {
    printf("Failed in the LU factorization.\n");
    return false;
  }

  LUsolve_internal6(B, p, b, x);

  return true; 
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
double makeQuadraticInterpolation(const double p[2], 
                                  const double coord[6][2], 
                                  const double b[6], bool &status)
{
  int i, j, k;

  double dx, dy; 

  double mat[6][6]; 

  double u[6], c[6], sum;

  c[0] = 1.0; 

  for (i = 0; i < 6; i++) {

    dx = coord[i][0] - p[0]; 
    dy = coord[i][1] - p[1]; 

    c[1] = dx;
    c[2] = dy;

    c[3] = 0.5 * dx * dx;
    c[4] = 0.5 * dy * dy; 

    c[5] = dx * dy; 

    for (j = 0; j < 6; j++) {
      mat[i][j] = c[j];
    }
  }

  status = LUsolve6(mat, b, u); 

  return u[0];
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
bool makeQuadraticInterpolation(const double p[2],
                                const double coord[6][2], 
                                const double b[6], double &u0)
{
  int i, j, k;

  double dx, dy; 

  double mat[6][6]; 

  double u[6], c[6], sum; 

  c[0] = 1.0; 

  for (i = 0; i < 6; i++) {

    dx = coord[i][0] - p[0];
    dy = coord[i][1] - p[1];

    c[1] = dx; 
    c[2] = dy; 

    c[3] = 0.5 * dx * dx;
    c[4] = 0.5 * dy * dy;

    c[5] = dx * dy; 

    for (j = 0; j < 6; j++) {
      mat[i][j] = c[j]; 
    }
  }

  //bool status = LUsolve6(mat, b, u); 
	bool status = solveByQRdecomposition<6>(mat, b, u, 6);
  if (status == false) {
    return false; 
  }

  u0 = u[0]; 

  return status;
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
bool makeQuadraticInterpolation(const double p[2], 
                                const double coord[6][2],
                                const double b[6], double u[6])
{
  int i, j, k;

  double dx, dy; 

  double mat[6][6];

  double c[6], sum;

  c[0] = 1.0; 

  for (i = 0; i < 6; i++) {

    dx = coord[i][0] - p[0]; 
    dy = coord[i][1] - p[1];

    c[1] = dx;
    c[2] = dy; 

    c[3] = 0.5 * dx * dx; 
    c[4] = 0.5 * dy * dy; 

    c[5] = dx * dy;

    for (j = 0; j < 6; j++) {
      mat[i][j] = c[j];
    }
  }

  bool status = LUsolve6(mat, b, u);
  //bool status = solveByQRdecomposition<6>(mat, b, u, 6);
  if (status == false) {
    return false; 
  }

  return status; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*template<int N>
bool makeLeastSquareInterpolation(const double p[2], 
                                	const double coord[N][2],
                                	const double b[N], double u[6])
{
  int i, j, k;

  double dx, dy; 

  double A[N][6];

  double c[6];

  c[0] = 1.0; 

  for (i = 0; i < N; i++) {

    dx = coord[i][0] - p[0]; 
    dy = coord[i][1] - p[1];

    c[1] = dx;
    c[2] = dy; 

    c[3] = 0.5 * dx * dx; 
    c[4] = 0.5 * dy * dy; 

    c[5] = dx * dy;

    for (j = 0; j < 6; j++) {
      A[i][j] = c[j];
    }
  }

	double mat[6][6], rhs[6];

	for(int i = 0; i < 6; i++){
		double s1 = 0.0;
		for(int j = 0; j < N; j++){
			s1 += A[j][i] * b[j];
		}
		rhs[i] = s1;
	}

	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			double s = 0.0;
			for(int k = 0; k < N; k++){
				s += A[k][i] * A[k][j];
			}
			mat[i][j] = s;
		}
	}

  //bool status = LUsolve6(mat, b, u);
  bool status = solveByQRdecomposition<6> (mat, rhs, u, 6);
  if (status == false) {
    return false; 
  }

  return status; 
}*/

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*template<int K>
bool solveByLeastSquareMethod1(const double A[][K], const double b[], 
															 double u[K], int n)
{
	Eigen::MatrixXd M(n, K);
	Eigen::VectorXd f(n), x(K);

	for(int i = 0; i < n; i++){
		f[i] = b[i];
		for(int j = 0; j < K; j++){
			M.coeffRef(i, j) = A[i][j];
		}
	}

	//x = M.colPivHouseholderQr().solve(f);
	//x = M.fullPivHouseholderQr().solve(f);
	//x = M.householderQr().solve(f);
	//x = (M.transpose()*M).ldlt().solve(M.transpose()*f);
	x = M.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(f);

	for(int k = 0; k < K; k++){
		u[k] = x[k];
	}

	//Eigen::MatrixXd M1 = (M.transpose()*M);
	//Eigen::MatrixXd M2 = M1.inverse();
	//std::cout << "LS condition number = " << sqrt(M1.norm()*M2.norm()) << std::endl;

	return true;

}*/

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int K>
bool solveByLeastSquareMethod(const double A[][K], const double b[], 
															 double u[K], int n)
{
	double mat[K][K], rhs[K];

	for(int i = 0; i < K; i++){
		for(int j = i; j < K; j++){
			double s = 0.0;
			for(int k = 0; k < n; k++){
				s += A[k][i] * A[k][j];
			}
			mat[i][j] = s;
		}
	}

	for(int i = 0; i < K; i++){
		for(int j = 0; j < i; j++){
			mat[i][j] = mat[j][i];
		}
	}

	for(int i = 0; i < K; i++){
		double s1 = 0.0;
		for(int j = 0; j < n; j++){
			s1 += A[j][i] * b[j];
		}
		rhs[i] = s1;
	}

  bool status = solveByQRdecomposition<K>(mat, rhs, u, K);
	return status;

	//assert(status);

	//int itr_num = solveByCGMethod<K>(mat, rhs, u);
	//assert(itr_num >= 0);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool makeLeastSquareInterpolation(const double p[2], 
                                	const VectorX2d &coord,
                                	const VectorXd &v, int N, double u[6])
{
  int i, j, k;

  double dx, dy; 

	double (*A)[6] = new double[N][6];
	double *b = new double[N];

  double c[6];

  for (i = 0; i < N; i++) {

    dx = coord[i][0] - p[0]; 
    dy = coord[i][1] - p[1];

		double d2 = dx*dx+dy*dy;
		double d3 = d2*sqrt(d2);

		double w = 1.0 / (d3+1.0e-2);
		//double w = 1.0;

  	c[0] = 1.0; 

    c[1] = dx;
    c[2] = dy; 

    c[3] = dx * dx; 
    c[4] = dy * dy; 

    c[5] = dx * dy;

    for (j = 0; j < 6; j++) {
      A[i][j] = c[j] * w;
    }
		b[i] = v[i] * w;
  }

	bool status = solveByLeastSquareMethod<6>(A, b, u, N);

	c[3] *= 2.0;
	c[4] *= 2.0;

	delete[] A;
	delete[] b;

  return status; 
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
