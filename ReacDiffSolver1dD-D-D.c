/*=============================================================================
*   
*   Filename : ReacDiffSolver1dD-D-D.c
*   Creator : Han Zhou
*   Date : 11/29/21
*   Description : 
*
=============================================================================*/
   
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "Variables.h"
#include "Tridiag.h"
#include "QR.h"

#include "ReacDiffSolver1dD-D-D.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static double makeLinearInterpolation(const double t[2], const double f[2], 
																			double xi)
{
	double alpha = (t[1] - xi) / (t[1] - t[0]);
	double beta = 1.0 - alpha;

	return f[0] * alpha + f[1] * beta;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void makeQuadraticInterpolation(const double t[3], const double f[3], 
                            		 			 double xi, double g[3])
{
	double mat[3][3], c[3];

	for(int r = 0; r < 3; r++){
		mat[r][2] = 1.0;
		mat[r][1] = t[r] - xi;
		mat[r][0] = mat[r][1] * mat[r][1];
	}

	solveByQRdecomposition<3> (mat, f, c, 3);

	g[0] = c[2];
	g[1] = c[1];
	g[2] = c[0] + c[0];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void makeCubicInterpolation(const double t[4], const double f[4], 
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
static void interpolateTwoSideData(const VectorXd &x, const VectorXb &interior,
																	 const VectorXd &F, double a, double &fl,
																	 double &fr)
{
	double h = x[1] - x[0];

	int ia = static_cast<int>((a - x[0]) / h);

	if (interior[ia] == interior[ia+1]) {
		if (interior[ia-1] != interior[ia]) {
			ia--;
		} else if (interior[ia+1] != interior[ia+2]) {
			ia++;
		} else {
			std::cout << "false location." << std::endl;
		}
	}

	{
		double t[3] = {0.0, -1.0, -2.0};
		double f[3] = {F[ia], F[ia-1], F[ia-2]};
		double xi = (a - x[ia]) / h;
		double g[3];

		if (interior[ia] == interior[ia-1]) {

			if (interior[ia] == interior[ia-2]) {
				makeQuadraticInterpolation(t, f, xi, g);
				fl = g[0];
			} else {
				fl = makeLinearInterpolation(t, f, xi);
			}

		} else {
			fl = f[0];
		}
	}

	{
		double t[3] = {0.0, 1.0, 2.0};
		double f[3] = {F[ia+1], F[ia+2], F[ia+3]};
		double xi = (a - x[ia+1]) / h;
		double g[3];

		if (interior[ia+1] == interior[ia+2]) {

			if (interior[ia+1] == interior[ia+3]) {
				makeQuadraticInterpolation(t, f, xi, g);
				fr = g[0];
			} else {
				fr = makeLinearInterpolation(t, f, xi);
			}

		} else {
			fr = f[0];
		}
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void 
callThomasAlgorithm(const VectorXd &f, const VectorXd &eta, int n, 
										double vl, double vr, VectorXd &v)
{
	// Dirichlet Green function

	// 	The Thomas solver solves the 1d modified Helmholtz equantion
	//	- v"(x) + \kappa v(x) = f(x)
	//	v(x_0) = va, v(x_n) = vb;
	// 	with the inputs : f_i = f(x_i) * h^2, eta = \kappa * h^2, 

	int dof = n - 1;
	double *a = new double[dof];
	double *b = new double[dof-1];
	double *c = new double[dof-1];
	double *rhs = new double[dof];
	double *sol = new double[dof];

	for(int i = 1, s = 0; i < n; i++, s++){
		rhs[s] = f[i];
		a[s] = 2.0 + eta[i];
	}

	for(int i = 0; i < dof - 1; i++){
		b[i] = c[i] = - 1.0;
	}

	rhs[0] += vl;
	rhs[dof-1] += vr;

	solveByThomasAlgorithm(a, b, c, rhs, dof, sol);

	for(int i = 1, s = 0; i < n; i++, s++){
		v[i] = sol[s];
	}

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] rhs;
	delete[] sol;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*static void 
callThomasAlgorithm(const VectorXd &f, const VectorXd &eta, int n, 
										double vl, double vr, VectorXd &v)
{
	// Neumann Green function

	// 	The Thomas solver solves the 1d modified Helmholtz equantion
	//	- v"(x) + \kappa v(x) = f(x)
	//	v(x_0) = va, v(x_n) = vb;
	// 	with the inputs : f_i = f(x_i) * h^2, eta = \kappa * h^2, 

	int dof = n + 1;
	double *a = new double[dof];
	double *b = new double[dof-1];
	double *c = new double[dof-1];
	double *rhs = new double[dof];
	double *sol = new double[dof];

	for(int i = 0, s = 0; i <= n; i++, s++){
		rhs[s] = f[i];
		a[s] = 2.0 + eta[i];
	}

	for(int i = 0; i < dof - 1; i++){
		b[i] = c[i] = - 1.0;
	}
	b[0] = c[dof-2] = - 2.0;

	rhs[0] += vl;
	rhs[dof-1] += vr;

	solveByThomasAlgorithm(a, b, c, rhs, dof, sol);

	for(int i = 0, s = 0; i <= n; i++, s++){
		v[i] = sol[s];
	}

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] rhs;
	delete[] sol;
}*/

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void 
makeCorrection(const VectorXd &x, const VectorXb &interior, int n, 
							 const VectorXd &bdry_crd, int bdry_num, 
							 double kap_i, double kap_e, const double *phi, 
							 VectorXd &f)
{
	double h = x[1] - x[0];

	double d, c, jmp0, jmp1, jmp2;

	for(int k = 0; k < bdry_num; k++){

		double a = bdry_crd[k];
		int ia = static_cast<int>((a - x[0]) / h);

		if (interior[ia] == interior[ia+1]) {
			if (interior[ia-1] != interior[ia]) {
				ia--;
			} else if (interior[ia+1] != interior[ia+2]) {
				ia++;
			} else {
				std::cout << "false location." << std::endl;
			}
		}

		if (interior[ia]) {

			assert(!interior[ia+1]);

			jmp0 = phi[k];
			jmp1 = 0.0;
			jmp2 = 0.0;

			d = x[ia + 1] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;
			f[ia] += c;

			d = x[ia] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;
			f[ia + 1] -= c;

		} else {

			assert(!interior[ia]);

			jmp0 = phi[k];
			jmp1 = 0.0;
			jmp2 = 0.0;

			d = x[ia + 1] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;
			f[ia] -= c;

			d = x[ia] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;
			f[ia + 1] += c;

		}
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void 
extractDirichletData(const VectorXd &x, const VectorXb &interior, int n, 
										 const VectorXd &bdry_crd, int bdry_num, 
										 double kap_i, double kap_e, const double *phi,
										 const VectorXd &v, VectorX6d &v_bdry_data)
{
	double h = x[1] - x[0];
	double h2 = h * h;

	double jmp0, jmp1, jmp2, c, d;

	double t[3] = {-1.0, 0.0, 1.0};
	double f[3], xi, g[3];

	for(int k = 0; k < bdry_num; k++){

		double a = bdry_crd[k];
		int ia = static_cast<int>((a - x[0]) / h);

		if (interior[ia] == interior[ia+1]) {
			if (interior[ia-1] != interior[ia]) {
				ia--;
			} else if (interior[ia+1] != interior[ia+2]) {
				ia++;
			} else {
				std::cout << "false location." << std::endl;
			}
		}

		if (interior[ia]) {

			assert(!interior[ia+1]);

			jmp0 = phi[k];
			jmp1 = 0.0;
			jmp2 = 0.0;

			d = x[ia+1] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;

			f[0] = v[ia - 1];
			f[1] = v[ia];
			f[2] = v[ia + 1] + c;

			xi = (a - x[ia]) / h;

		} else {

			assert(!interior[ia]);

			jmp0 = phi[k];
			jmp1 = 0.0;
			jmp2 = 0.0;

			d = x[ia] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;

			f[0] = v[ia] + c;
			f[1] = v[ia + 1];
			f[2] = v[ia + 2];

			xi = (a - x[ia + 1]) / h;

		}

		makeQuadraticInterpolation(t, f, xi, g);

		v_bdry_data[k][0] = g[0];
		v_bdry_data[k][1] = g[1] / h;
		v_bdry_data[k][2] = g[2] / h2;
		v_bdry_data[k][3] = v_bdry_data[k][0] - jmp0;
		v_bdry_data[k][4] = v_bdry_data[k][1] - jmp1;
		v_bdry_data[k][5] = v_bdry_data[k][2] - jmp2;

	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void 
makeCorrectionV(const VectorXd &x, const VectorXb &interior, int n, 
								const VectorXd &bdry_crd, int bdry_num, 
								double kap_i, double kap_e, const VectorXd &F, 
								const VectorXd &u_bdry, VectorXd &f)
{
	double h = x[1] - x[0];

	double d, c, jmp0, jmp1, jmp2, f_jmp;

	for(int k = 0; k < bdry_num; k++){

		double a = bdry_crd[k];
		int ia = static_cast<int>((a - x[0]) / h);

		if (interior[ia] == interior[ia+1]) {
			if (interior[ia-1] != interior[ia]) {
				ia--;
			} else if (interior[ia+1] != interior[ia+2]) {
				ia++;
			} else {
				std::cout << "false location." << std::endl;
			}
		}

		double bdry_val = u_bdry[k];

		double fl = F[ia];
		double fr = F[ia + 1];

		interpolateTwoSideData(x, interior, F, a, fl, fr);

		if (interior[ia]) {

			assert(!interior[ia+1]);

			f_jmp = fl - fr;

			jmp0 = 0.0;
			jmp1 = 0.0;
			jmp2 = bdry_val * (kap_i - kap_e) - f_jmp;

			d = x[ia + 1] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;
			f[ia] += c;

			d = x[ia] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;
			f[ia + 1] -= c;

		} else {

			assert(!interior[ia]);

			f_jmp = fr - fl;

			jmp0 = 0.0;
			jmp1 = 0.0;
			jmp2 = bdry_val * (kap_i - kap_e) - f_jmp;

			d = x[ia + 1] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;
			f[ia] -= c;

			d = x[ia] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;
			f[ia + 1] += c;

		}
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void 
extractDirichletDataV(const VectorXd &x, const VectorXb &interior, int n, 
											const VectorXd &bdry_crd, int bdry_num, 
											double kap_i, double kap_e, 
											const VectorXd &F, const VectorXd &u_bdry,
											const VectorXd &v, VectorX6d &v_bdry_data)
{
	double h = x[1] - x[0];
	double h2 = h * h;

	double jmp0, jmp1, jmp2, f_jmp, c, d;

	double t[3] = {-1.0, 0.0, 1.0};
	double f[3], xi, g[3];

	for(int k = 0; k < bdry_num; k++){

		double a = bdry_crd[k];
		int ia = static_cast<int>((a - x[0]) / h);

		if (interior[ia] == interior[ia+1]) {
			if (interior[ia-1] != interior[ia]) {
				ia--;
			} else if (interior[ia+1] != interior[ia+2]) {
				ia++;
			} else {
				std::cout << "false location." << std::endl;
			}
		}

		double bdry_val = u_bdry[k];

		double fl = F[ia];
		double fr = F[ia + 1];

		interpolateTwoSideData(x, interior, F, a, fl, fr);

		if (interior[ia]) {

			assert(!interior[ia+1]);

			f_jmp = fl - fr;

			jmp0 = 0.0;
			jmp1 = 0.0;
			jmp2 = bdry_val * (kap_i - kap_e) - f_jmp;

			d = x[ia+1] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;

			f[0] = v[ia - 1];
			f[1] = v[ia];
			f[2] = v[ia + 1] + c;

			xi = (a - x[ia]) / h;


		} else {

			assert(!interior[ia]);

			f_jmp = fr - fl;

			jmp0 = 0.0;
			jmp1 = 0.0;
			jmp2 = bdry_val * (kap_i - kap_e) - f_jmp;

			d = x[ia] - a;
			c = jmp0 + d * jmp1 + 0.5 * d * d * jmp2;

			f[0] = v[ia] + c;
			f[1] = v[ia + 1];
			f[2] = v[ia + 2];

			xi = (a - x[ia + 1]) / h;

		}

		makeQuadraticInterpolation(t, f, xi, g);

		v_bdry_data[k][0] = g[0];
		v_bdry_data[k][1] = g[1] / h;
		v_bdry_data[k][2] = g[2] / h2;
		v_bdry_data[k][3] = v_bdry_data[k][0] - jmp0;
		v_bdry_data[k][4] = v_bdry_data[k][1] - jmp1;
		v_bdry_data[k][5] = v_bdry_data[k][2] - jmp2;

	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void 
computeMatrixVectorProduct(const VectorXd &x, const VectorXb &interior, int n,
													 const VectorXd &bdry_crd, int bdry_num,
													 double kap_i, double kap_e,
													 const double *phi, double *w_bdry)
{
	int n1 = n + 1;
	double h = x[1] - x[0];
	double h2 = h * h;

	VectorXd f(n1), w(n1);;

	VectorXd eta(n1);

	for(int i = 0; i <= n; i++){
		if (interior[i]) {
			eta[i] = kap_i * h2;
		} else {
			eta[i] = kap_e * h2;
		}
	}

	f.fill(0.0);
	w.fill(0.0);

	makeCorrection(x, interior, n, bdry_crd, bdry_num, kap_i, kap_e, phi, f);

	w[0] = 0.0;
	w[n] = 0.0;

	callThomasAlgorithm(f, eta, n, w[0], w[n], w);

	VectorX6d w_bdry_data(bdry_num);

	extractDirichletData(x, interior, n, bdry_crd, bdry_num, kap_i, kap_e,
											 phi,  w, w_bdry_data);

	for(int k = 0; k < bdry_num; k++){
		w_bdry[k] = w_bdry_data[k][0];
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static double makeQRdecomposition(double **H, double **R, double *cs,
                           				double *sn, double *b, int m, double atol)
{
  int m1 = m + 1;

  int ell = m - 1;
  for (int j = 0; j < m1; j++) {
    R[j][ell] = H[j][ell]; 
  }

  for (int k = 0; k < ell; k++) {
    double t = cs[k] * R[k][ell] + sn[k] * R[k + 1][ell]; 
    R[k + 1][ell] = - sn[k] * R[k][ell] + cs[k] * R[k + 1][ell]; 
    R[k][ell] = t; 
  }

  if (fabs(R[ell + 1][ell]) > atol) {
    double x = R[ell][ell];
    double y = R[ell + 1][ell]; 
    double r = sqrt(x * x + y * y);
    cs[ell] = x / r; 
    sn[ell] = y / r;

    int k = ell;

    double t = cs[k] * R[k][ell] + sn[k] * R[k + 1][ell];
    R[k + 1][ell] = - sn[k] * R[k][ell] + cs[k] * R[k + 1][ell]; 
    R[k][ell] = t; 

    t = cs[k] * b[k] + sn[k] * b[k + 1];
    b[k + 1] = - sn[k] * b[k] + cs[k] * b[k + 1];
    b[k] = t;
  }

  double res = b[m];

  return res; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static bool GMRES(const VectorXd &x, const VectorXb &interior, int grid_n,
									const VectorXd &bdry_crd, int bdry_num, 
									double kap_i, double kap_e, 
									const double *b, double *u, int n, int max_m, 
           				int max_itr_num, double rtol, int &itr_num)
{
  double atol = 1.0E-15; 

  itr_num = 0; 

  double *r = new double[n]; 
  double *w = new double[n]; 
  for (int i = 0; i < n; i++) {
    r[i] = w[i] = 0.0; 
  }

  computeMatrixVectorProduct(x, interior, grid_n, bdry_crd, bdry_num, kap_i, kap_e, u, w); 

  for (int i = 0; i < n; i++) {
    r[i] = b[i] - w[i];
  }

  double norm_r0 = computeVectorNorm(r, n);
  if ((norm_r0 / sqrt(n)) < atol) {
    delete[] r; r = 0; 
    delete[] w; w = 0;
    return true; 
  }

  int max_m1 = max_m + 1; 
  int max_m2 = max_m + 2;

  double **V = new double*[max_m2]; 
  double **H = new double*[max_m2];
  double **R = new double*[max_m2];

  for (int i = 0; i < max_m2; i++) {
    V[i] = new double[n];
    H[i] = new double[max_m1]; 
    R[i] = new double[max_m1];
  }

  for (int i = 0; i < max_m2; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = 0.0;
    }
    for (int j = 0; j < max_m1; j++) {
      H[i][j] = 0.0;
      R[i][j] = 0.0;
    }
  }

  double *cs = new double[max_m1]; 
  double *sn = new double[max_m1];
  for (int i = 0; i < max_m1; i++) {
    cs[i] = 1.0;
    sn[i] = 0.0; 
  }

  double *z = new double[max_m1];
  for (int i = 0; i < max_m1; i++) {
    z[i] = 0.0; 
  }

  int m = 0; 
  double beta = norm_r0;
  double r_beta = 1.0 / beta;
  for (int j = 0; j < n; j++) {
    V[m][j] = r[j] * r_beta;
  }

  double *c = new double [max_m2]; 
  for (int i = 1; i < max_m2; i++) {
    c[i] = 0.0;
  }
  c[0] = beta; 

  int done = 0; 
  double tol = rtol * norm_r0; 
  while ((itr_num < max_itr_num) && (!done)) {

  	computeMatrixVectorProduct(x, interior, grid_n, bdry_crd, bdry_num, kap_i, kap_e, V[m], w); 

    for (int j = 0; j <= m; j++) {
      H[j][m] = computeInnerProduct(V[j], w, n);
      for (int i = 0; i < n; i++) {
        w[i] -= V[j][i] * H[j][m]; 
      }
    }
    H[m + 1][m] = computeVectorNorm(w, n);

    if (fabs(H[m + 1][m]) > atol) {
      for (int i = 0; i < n; i++) {
        V[m + 1][i] = w[i] / H[m + 1][m];
      }
    } else {
      done = 1; 
    }

    double rz = makeQRdecomposition(H, R, cs, sn, c,  m + 1, atol); 

    itr_num++;

    if (fabs(rz) < tol) {

      for (int i = m; i >= 0; i--) {
        double s = 0.0; 
        for (int j = m; j > i; j--) {
          s += R[i][j] * z[j];
        }
        z[i] = (c[i] - s) / R[i][i];
      }

      for (int i = 0; i < n; i++) {
        double sum = 0.0; 
        for (int j = 0; j <= m; j++) {
          sum += V[j][i] * z[j];
        }
        u[i] += sum;
      }

      done = 1; 

    } else {

      if (m == max_m) {

        for (int i = m; i >= 0; i--) {
          double s = 0.0; 
          for (int j = m; j > i; j--) {
            s += R[i][j] * z[j];
          }
          z[i] = (c[i] - s) / R[i][i]; 
        }

        for (int i = 0; i < n; i++) {
          double sum = 0.0;
          for (int j = 0; j <= m; j++) {
            sum += V[j][i] * z[j];
          }
          u[i] += sum;
        }

  			computeMatrixVectorProduct(x, interior, grid_n, bdry_crd, bdry_num, kap_i, kap_e, u, w); 
        for (int i = 0; i < n; i++) {
          r[i] = b[i] - w[i];
        }
        norm_r0 = computeVectorNorm(r, n); 
        if ((norm_r0 / sqrt(n)) < atol) {
          done = 1; 
        } else {
          m = 0;
          beta = norm_r0; 
          tol = rtol * norm_r0;
          double r_beta = 1.0 / beta;
          for (int j = 0; j < n; j++) {
            V[m][j] = r[j] * r_beta; 
          }
          for (int i = 1; i < max_m2; i++) {
            c[i] = 0.0;
          }
          c[0] = beta; 
        }
      } else {
        m++;
      }
    }
  }


  if (itr_num >= max_itr_num) {
    std::cout << "GMRES iteration failed." << std::endl;
    exit(EXIT_FAILURE);
    return false; 
  }

  delete[] cs; cs = 0;
  delete[] sn; sn = 0; 

  for (int i = 0; i < max_m2; i++) {
    delete[] V[i]; V[i] = 0;
    delete[] H[i]; H[i] = 0;
    delete[] R[i]; R[i] = 0;
  }
  delete[] H; H = 0; 
  delete[] R; R = 0;
  delete[] V; V = 0;

  delete[] r; r = 0;
  delete[] w; w = 0; 

  delete[] z; z = 0;
  delete[] c; c = 0;


  return true;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ReacDiffSolver1d_DDD(const VectorXd &x, const VectorXb &interior, int n,
											 const VectorXd &bdry_crd, int bdry_num, double kap_i,
											 double kap_e, const VectorXd &F, 
											 const VectorXd &u_bdry, double ul, double ur,
											 VectorXd &u, VectorX6d &u_bdry_data)

{
	// this version solves two point BVP subject to Dirichlet BC,
	// double layer potential is used (Fredholm IE of the 2nd kind)
	// Dirichlet BC for Green's function

	int n1 = n + 1;
	double h = x[1] - x[0];
	double h2 = h * h;

	VectorXd f(n1), v(n1), w(n1);
	VectorX6d w_bdry_data(bdry_num), v_bdry_data(bdry_num);

	VectorXd eta(n1);

	for(int i = 0; i <= n; i++){
		if (interior[i]) {
			eta[i] = kap_i * h2;
		} else {
			eta[i] = kap_e * h2;
		}

		f[i] = F[i] * h2;
	}

	makeCorrectionV(x, interior, n, bdry_crd, bdry_num, 
									kap_i, kap_e, F, u_bdry, f);

	v[0] = ul;
	v[n] = ur;

	callThomasAlgorithm(f, eta, n, v[0], v[n], v);

	extractDirichletDataV(x, interior, n, bdry_crd, bdry_num, 
												kap_i, kap_e, F, u_bdry, v, v_bdry_data);

	double *phi = new double[bdry_num];
	double *g = new double[bdry_num];

	for(int k = 0; k < bdry_num; k++){
		phi[k] = 0.0;
		g[k] = u_bdry[k] - v_bdry_data[k][0];
	}

	int max_itr_num = 100;
	int max_m = 20;
	if (max_m > bdry_num) {
		max_m = bdry_num;
	}

	int itr_num = 0;
	bool status = GMRES(x, interior, n, bdry_crd, bdry_num, kap_i, kap_e, 
											g, phi, bdry_num, max_m, max_itr_num, 1.0E-8, itr_num);


	f.fill(0.0);

	makeCorrection(x, interior, n, bdry_crd, bdry_num, kap_i, kap_e, phi, f);

	w[0] = 0.0;
	w[n] = 0.0;

	callThomasAlgorithm(f, eta, n, w[0], w[n], w);

	extractDirichletData(x, interior, n, bdry_crd, bdry_num, 
											 kap_i, kap_e, phi, w, w_bdry_data);

	for(int i = 0; i <= n; i++){
		u[i] = w[i] + v[i];
	}

	for(int k = 0; k < bdry_num; k++){
		for(int s = 0; s < 6; s++){
			u_bdry_data[k][s] = v_bdry_data[k][s] + w_bdry_data[k][s];
		}
	}

	delete[] phi;
	delete[] g;
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
