/*=============================================================================
*   
*   Filename : ReacDiffSolver1d-011.c
*   Creator : Han Zhou
*   Date : 01/11/22
*   Description : 
*
=============================================================================*/
   
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

//#include "Interpolation.h"   
#include "Stencil.h"
#include "Variables.h"
#include "Tridiag.h"
#include "QR.h"

#include "DiffusionSolver1d-001-O4.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int N> // degree N-1, node number N
static void makePolynomialInterpolation1d(const double t[N], 
																					const double f[N], 
																	 				double xi, 
																					double g[N])
{
	double mat[N][N];

	for(int i = 0; i < N; i++){
		double d = t[i] - xi;
		mat[i][0] = 1.0;
		for(int j = 1; j < N; j++){
			mat[i][j] = mat[i][j-1] * d;
		}
	}

	solveByQRdecomposition<N>(mat, f, g, N);

	int fac = 1;
	for(int i = 2; i < N; i++){
		fac *= i;
		g[i] *= fac;
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void interpolateTwoSideData(const VectorXd &x, 
																	 const VectorXb &interior,
											 						 const VectorXd &F, 
																	 double a, double &f_jmp,
											 						 double &fx_jmp, 
																	 double &fxx_jmp)
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

	double fr[4] = {0.0, 0.0, 0.0, 0.0};
	double fl[4] = {0.0, 0.0, 0.0, 0.0};

	{
		double t[4] = {0.0, -1.0, -2.0, -3.0};
		double f[4] = {F[ia], F[ia-1], F[ia-2], F[ia-3]};
		double xi = (a - x[ia]) / h;

		int count = 0;
		for(int m = 0, k = ia; m < 4; m++, k--){
			if (interior[k] == interior[ia]) {
				count++;
			} else {
				break;
			}
		}

		if (count == 0) {
			std::cout << "invalid. " << std::endl;
			exit(1);
		} else if (count == 1) {
				makePolynomialInterpolation1d<1>(t, f, xi, fl);
		} else if (count == 2) {
				makePolynomialInterpolation1d<2>(t, f, xi, fl);
		} else if (count == 3) {
				makePolynomialInterpolation1d<3>(t, f, xi, fl);
		} else if (count == 4) {
				makePolynomialInterpolation1d<4>(t, f, xi, fl);
		}

	}

	{
		double t[4] = {0.0, 1.0, 2.0, 3.0};
		double f[4] = {F[ia+1], F[ia+2], F[ia+3], F[ia+4]};
		double xi = (a - x[ia+1]) / h;

		int count = 0;
		for(int m = 0, k = ia+1; m < 4; m++, k++){
			if (interior[k] == interior[ia+1]) {
				count++;
			} else {
				break;
			}
		}

		if (count == 0) {
			std::cout << "invalid. " << std::endl;
			exit(1);
		} else if (count == 1) {
				makePolynomialInterpolation1d<1>(t, f, xi, fr);
		} else if (count == 2) {
				makePolynomialInterpolation1d<2>(t, f, xi, fr);
		} else if (count == 3) {
				makePolynomialInterpolation1d<3>(t, f, xi, fr);
		} else if (count == 4) {
				makePolynomialInterpolation1d<4>(t, f, xi, fr);
		}
	}

	if (interior[ia]) {
		f_jmp = fl[0] - fr[0];
		fx_jmp = fl[1] - fr[1];
		fxx_jmp = fl[2] - fr[2];
	} else {
		f_jmp = fr[0] - fl[0];
		fx_jmp = fr[1] - fl[1];
		fxx_jmp = fr[2] - fl[2];
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void computeRHSJumps(const VectorXd &x, 
														const VectorXb &interior,
														const VectorXd &F, 
														const VectorXd &bdry_crd,
														int bdry_num, 
														double *f_jmp, 
														double *fx_jmp,
														double *fxx_jmp)
{
	for(int k = 0; k < bdry_num; k++){
		interpolateTwoSideData(x, interior, F, bdry_crd[k], 
													 f_jmp[k], fx_jmp[k], fxx_jmp[k]);
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void callThomasAlgorithm(const VectorXd &f, 
																int n, 
																double eta, 
																double vl, 
																double vr, 
																VectorXd &v)
{
	// Dirichlet Green function

	// 	The Thomas solver solves the 1d modified Helmholtz equantion (compact)
	//	- v"(x) + \kappa v(x) = f(x)
	//	v(x_0) = va, v(x_n) = vb;
	// 	with the inputs : f_i = f(x_i) * h^2, eta = \kappa * h^2 /deno, 

	int dof = n - 1;
	double *a = new double[dof];
	double *b = new double[dof-1];
	double *c = new double[dof-1];
	double *rhs = new double[dof];
	double *sol = new double[dof];

	double diag = 2.0 + eta;

	for(int i = 1, s = 0; i < n; i++, s++){
		a[s] = diag;
		rhs[s] = f[i];
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
static void makeCorrection(const VectorXd &x, 
													 const VectorXb &interior, int n, 
							 						 const VectorXd &bdry_crd, int bdry_num, 
													 double kappa,
							 						 const double *phi, 
													 const double *psi, 
													 const double *f_jmp,
							 						 const double *fx_jmp, 
													 const double *fxx_jmp, 
													 VectorXd &f)
{
	double h = x[1] - x[0];

	double d, c, jmp0, jmp1, jmp2, jmp3, jmp4;

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
			jmp1 = psi[k];
			jmp2 = kappa * jmp0 - f_jmp[k];
			jmp3 = kappa * jmp1 - fx_jmp[k];
			jmp4 = kappa * jmp2 - fxx_jmp[k];

			d = x[ia + 1] - a;
			c = (((jmp4 / 24.0 * d + jmp3 / 6.0) * d + jmp2 / 2.0) * d + jmp1) * d + jmp0;
			f[ia] += c;

			d = x[ia] - a;
			c = (((jmp4 / 24.0 * d + jmp3 / 6.0) * d + jmp2 / 2.0) * d + jmp1) * d + jmp0;
			f[ia + 1] -= c;

		} else {

			assert(!interior[ia]);

			jmp0 = phi[k];
			jmp1 = - psi[k];
			jmp2 = kappa * jmp0 - f_jmp[k];
			jmp3 = kappa * jmp1 - fx_jmp[k];
			jmp4 = kappa * jmp2 - fxx_jmp[k];

			d = x[ia + 1] - a;
			c = (((jmp4 / 24.0 * d + jmp3 / 6.0) * d + jmp2 / 2.0) * d + jmp1) * d + jmp0;
			f[ia] -= c;

			d = x[ia] - a;
			c = (((jmp4 / 24.0 * d + jmp3 / 6.0) * d + jmp2 / 2.0) * d + jmp1) * d + jmp0;
			f[ia + 1] += c;

		}
	}

}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void 
extractBoundaryData(const VectorXd &x, const VectorXb &interior, int n, 
										 const VectorXd &bdry_crd, int bdry_num, double kappa, 
										 const double *phi, const double *psi, const double *f_jmp,
										 const double *fx_jmp, const double *fxx_jmp,
										 const VectorXd &v, VectorX10d &v_bdry_data)
{
	double h = x[1] - x[0];
	double h_2 = 0.5 * h;

	double hp[5];
	hp[0] = 1.0;
	for(int m = 1; m < 5; m++){
		hp[m] = hp[m-1] * h;
	}

	for(int k = 0; k < bdry_num; k++){

		double a = bdry_crd[k];

		int ia = static_cast<int>((a - x[0])/h + 0.5);
		int ia2 = static_cast<int>((a - x[0])/h);

		if (interior[ia2] == interior[ia2+1]) {
			if (interior[ia2-1] != interior[ia2]) {
				ia2--;
			} else if (interior[ia2+1] != interior[ia2+2]) {
				ia2++;
			} else {
				std::cout << "false location." << std::endl;
			}
			assert(interior[ia2] != interior[ia2+1]);
		}

		int dir;
		if (a > x[ia]) {
			dir = 1;
		} else {
			dir = -1;
		}

		int stc[5];
		get1DStencil_5(dir, stc);

		int count = 0;
		for(int m = 0; m < 5; m++){

			int i0 = ia + stc[m];
			bool tag1 = false, tag2 = false;

			if (fabs(x[ia] - a) < h) {
				tag1 = tag2 = true;
			} else if (x[i0] > a+h_2 && interior[i0] == interior[ia2+1]) {
				tag1 = true;
			} else if (x[i0] < a-h_2 && interior[i0] == interior[ia2]) {
				tag2 = true;
			}

			if (tag1 || tag2) {
				count++;
			} else {
				break;
			}
		}

		double jmp[5];

		jmp[0] = phi[k];
		
		if (interior[ia2]) {
			jmp[1] = psi[k];
		} else {
			jmp[1] = -psi[k];
		}

		jmp[2] = kappa * jmp[0] - f_jmp[k];
		jmp[3] = kappa * jmp[1] - fx_jmp[k];
		jmp[4] = kappa * jmp[2] - fxx_jmp[k];

		double t[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
		double f[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
		double g[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
		double xi = (a - x[ia]) / h;

		for(int m = 0; m < count; m++){
			t[m] = stc[m];
			int i0 = ia + stc[m];
			f[m] = v[i0];
			if (!interior[i0]) {
				double d = x[i0] - a;
				double c = (((jmp[4]/24.0*d + jmp[3]/6.0)*d + jmp[2]/2.0) * d + jmp[1]) * d + jmp[0];
				f[m] += c;
			}
		}


		if (count == 0) {
			std::cout << "a = " << a << std::endl;
			std::cout << "x0 = " << x[ia + stc[0]] << std::endl;
			std::cout << "invalid." << std::endl;
			std::cout << "Pass : " << __FILE__ << ", line : " << __LINE__ << std::endl;
			exit(1);
		} else if (count == 1) {
			makePolynomialInterpolation1d<1>(t, f, xi, g);
		} else if (count == 2) {
			makePolynomialInterpolation1d<2>(t, f, xi, g);
		} else if (count == 3) {
			makePolynomialInterpolation1d<3>(t, f, xi, g);
		} else if (count == 4) {
			makePolynomialInterpolation1d<4>(t, f, xi, g);
		} else if (count == 5) {
			makePolynomialInterpolation1d<5>(t, f, xi, g);
		}

		for(int m = 0; m < 5; m++){
			double val = g[m] / hp[m];
			v_bdry_data[k][m] = val;
			v_bdry_data[k][m+5] = val - jmp[m];
		}

	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
static void computeMatrixVectorProduct(const VectorXd &x, 
																			 const VectorXb &interior, int n,
													 						 const VectorXd &bdry_crd, int bdry_num,
													 						 double kappa, 
																			 const double *phi, double *w_bdry)
																			 //const double *psi, double *w_bdry)
{
	int n1 = n + 1;
	double h = x[1] - x[0];
	double h2 = h * h;

	double deno = 1.0 - kappa * h2 / 12.0;
	double eta = kappa * h2 / deno;

	double *psi = new double[bdry_num];
	//double *phi = new double[bdry_num];
	double *f_jmp = new double[bdry_num];
	double *fx_jmp = new double[bdry_num];
	double *fxx_jmp = new double[bdry_num];

	for(int k = 0; k < bdry_num; k++){
		psi[k] = f_jmp[k] = fx_jmp[k] = fxx_jmp[k] = 0.0;
		//phi[k] = f_jmp[k] = fx_jmp[k] = fxx_jmp[k] = 0.0;
	}

	VectorXd rhs(n1), w(n1);;
	rhs.fill(0.0);
	w.fill(0.0);

	makeCorrection(x, interior, n, bdry_crd, bdry_num, kappa, 
								 phi, psi, f_jmp, fx_jmp, fxx_jmp, rhs);

	w[0] = 0.0;
	w[n] = 0.0;
	callThomasAlgorithm(rhs, n, eta, w[0], w[n], w);

	VectorX10d w_bdry_data(bdry_num);

	extractBoundaryData(x, interior, n, bdry_crd, bdry_num, kappa, 
											phi, psi, f_jmp, fx_jmp, fxx_jmp,	w, w_bdry_data);

	for(int k = 0; k < bdry_num; k++){
		w_bdry[k] = w_bdry_data[k][0];
	}

	delete[] psi;
	//delete[] phi;
	delete[] f_jmp;
	delete[] fx_jmp;
	delete[] fxx_jmp;

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
static bool GMRES(const VectorXd &x, 
									const VectorXb &interior, int grid_n,
									const VectorXd &bdry_crd, int bdry_num, 
									double kappa, 
									const double *b, double *u, int n, 
									int max_m, int max_itr_num, 
									double rtol, int &itr_num)
{
  double atol = 1.0E-15; 

  itr_num = 0; 

  double *r = new double[n]; 
  double *w = new double[n]; 
  for (int i = 0; i < n; i++) {
    r[i] = w[i] = 0.0; 
  }

  computeMatrixVectorProduct(x, interior, grid_n, bdry_crd, bdry_num, kappa, u, w); 

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

  	computeMatrixVectorProduct(x, interior, grid_n, bdry_crd, bdry_num, kappa, V[m], w); 

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

  			computeMatrixVectorProduct(x, interior, grid_n, bdry_crd, bdry_num, kappa, u, w); 
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
void computeRightHandSide(const VectorXb &interior, 
													const VectorXd &F, 
													VectorXd &rhs, int n, 
													double kappa,	double h2)
{
	double deno = 1.0-kappa*h2/12.0;

	for(int i = 1; i < n; i++){

		double f = F[i];
		double lf = 0.0;

		if (interior[i] == interior[i+1] && interior[i] == interior[i-1]) {

			lf = F[i+1] + F[i-1] - 2.0 * F[i];

		} else if (interior[i] != interior[i-1]) {

			double t[4] = {0.0, 1.0, 2.0, 3.0};
			double w[4] = {F[i], F[i+1], F[i+2], F[i+3]};
			double g[4] = {0.0, 0.0, 0.0, 0.0};

			int count = 0;
			for(int m = 0, k = i; m < 4; m++, k++){
				if (interior[k] == interior[i]) {
					count++;
				} else {
					break;
				}
			}

			if (count == 3) {
				makePolynomialInterpolation1d<3>(t, w, 0.0, g);
			} else if (count == 4) {
				makePolynomialInterpolation1d<4>(t, w, 0.0, g);
			}

			lf = g[2];

		} else if (interior[i] != interior[i+1]) {

			double t[4] = {0.0, -1.0, -2.0, -3.0};
			double w[4] = {F[i], F[i-1], F[i-2], F[i-3]};
			double g[4] = {0.0, 0.0, 0.0, 0.0};

			int count = 0;
			for(int m = 0, k = i; m < 4; m++, k--){
				if (interior[k] == interior[i]) {
					count++;
				} else {
					break;
				}
			}

			if (count == 3) {
				makePolynomialInterpolation1d<3>(t, w, 0.0, g);
			} else if (count == 4) {
				makePolynomialInterpolation1d<4>(t, w, 0.0, g);
			}

			lf = g[2];
		}

		rhs[i] = (f + lf / 12.0) * h2 / deno;
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DiffusionSolver1d_001_O4(const VectorXd &x, 
															const VectorXb &interior, int n,
											 				const VectorXd &bdry_crd, 
															int bdry_num, double kappa,
											 				const VectorXd &F, 
															const VectorXd &u_bdry, 
											 				double ul, double ur,	
															VectorXd &u, 
															VectorX6d &u_bdry_data)
															//VectorX10d &u_bdry_data)

{
	// this version solves two point BVP subject to Dirichlet BC,
	// single layer potential is used (Fredholm IE of the 1st kind)
	// Dirichlet BC for Green's function

	int n1 = n + 1;
	double h = x[1] - x[0];
	double h2 = h * h;

	double deno = 1.0 - kappa * h2 / 12.0;
	double eta = kappa * h2 / deno;

	if (bdry_num == 0) {

		VectorXd rhs(n1);

		for(int i = 1; i < n; i++){
			rhs[i] = (F[i]+(F[i+1]+F[i-1]-2.0*F[i])/12.0)*h2/deno;
		}

		u[0] = ul;	
		u[n] = ur;
		callThomasAlgorithm(rhs, n, eta, u[0], u[n], u);

		return;
	}

	// compute volume integral

	VectorXd rhs(n1);
	rhs.fill(0.0);

	computeRightHandSide(interior, F, rhs, n, kappa, h2);

	double *f_jmp = new double[bdry_num];
	double *fx_jmp = new double[bdry_num];
	double *fxx_jmp = new double[bdry_num];

	computeRHSJumps(x, interior, F, bdry_crd, bdry_num, f_jmp, fx_jmp, fxx_jmp);

	double *phi = new double[bdry_num];
	double *psi = new double[bdry_num];

	for(int k = 0; k < bdry_num; k++){
		phi[k] = psi[k] = 0.0;
	}

	makeCorrection(x, interior, n, bdry_crd, bdry_num, kappa, 
								 phi, psi, f_jmp, fx_jmp, fxx_jmp, rhs);


	VectorXd v(n1);

	v[0] = ul;	
	v[n] = ur;
	callThomasAlgorithm(rhs, n, eta, v[0], v[n], v);

	VectorX10d v_bdry_data(bdry_num);

	extractBoundaryData(x, interior, n, bdry_crd, bdry_num, kappa, 
											phi, psi, f_jmp, fx_jmp, fxx_jmp,	v, v_bdry_data);

	// solve boundary integral equation

	double *g = new double[bdry_num];
	for(int k = 0; k < bdry_num; k++){
		g[k] = u_bdry[k] - v_bdry_data[k][0];
	}

	int max_m = 20;
	if (max_m > bdry_num) {
		max_m = bdry_num;
	}

	int itr_num = 0;
	bool status = GMRES(x, interior, n, bdry_crd, bdry_num, kappa,
											g, phi, bdry_num, max_m, 100, 1.0E-14, itr_num);
	//bool status = GMRES(x, interior, n, bdry_crd, bdry_num, kappa,
	//										g, psi, bdry_num, max_m, 100, 1.0E-14, itr_num);
	if (!status) {
		std::cout << "GMRES failed." << std::endl;
	}

	// evaluate solution

	for(int k = 0; k < bdry_num; k++){
		psi[k] = f_jmp[k] = fx_jmp[k] = fxx_jmp[k] = 0.0;
		//phi[k] = f_jmp[k] = fx_jmp[k] = fxx_jmp[k] = 0.0;
	}

	rhs.fill(0.0);

	makeCorrection(x, interior, n, bdry_crd, bdry_num, kappa, 
								 phi, psi, f_jmp, fx_jmp, fxx_jmp, rhs);

	VectorXd w(n1);
	w[0] = w[n] = 0.0;
	callThomasAlgorithm(rhs, n, eta, w[0], w[n], w);

	VectorX10d w_bdry_data(bdry_num);

	extractBoundaryData(x, interior, n, bdry_crd, bdry_num, kappa, 
											phi, psi, f_jmp, fx_jmp, fxx_jmp,	w, w_bdry_data);


	for(int i = 0; i < n1; i++){
		u[i] = v[i] + w[i];
	}

	//for(int k = 0; k < bdry_num; k++){
	//	for(int s = 0; s < 10; s++){
	//		u_bdry_data[k][s] = v_bdry_data[k][s] + w_bdry_data[k][s];
	//	}
	//}

	delete[] g;
	delete[] phi;
	delete[] psi;
	delete[] f_jmp;
	delete[] fx_jmp;
	delete[] fxx_jmp;

}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
