/*=============================================================================
*   
*   Filename : GrayScottEqn.C
*   Creator : Han Zhou
*   Date : 11/20/21
*   Description : 
*
=============================================================================*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
 
#include "GrayScottEqn.H"

#include "ReacDiffSolver1d3-2.h"
#include "ReacDiffSolver1dN-D.h"
#include "MathTools.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::initialize(double t0, double T)
{
	_t = t0;
	_T = T;

	//double deno = _I/32*10;
	//double deno = _I/16*10;
	double deno = _I/16*5;

	//_dt = (_T - _t) / deno;
	_dt = 1.0 / deno;

	_check_interval = static_cast<int>(_print_info_interval/_dt + 0.5);
	_plot_check_interval = static_cast<int>(_plot_interval/_dt + 0.5);
	_max_step_num = static_cast<int>((_T - _t) / _dt + 0.5);
	_step_count = 0;


	_u.reallocate(_I1, _J1);
	_uxx.reallocate(_I1, _J1);

	_v.reallocate(_I1, _J1);
	_vxx.reallocate(_I1, _J1);

	_work_u.reallocate(_I1, _J1);
	_work_v.reallocate(_I1, _J1);

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){

			bool side = _interior[i][j];
			if (side) {

				_u[i][j] = U(_x[i], _y[j], _t, side);
				_uxx[i][j] = Uxx(_x[i], _y[j], _t, side);

				_v[i][j] = V(_x[i], _y[j], _t, side);
				_vxx[i][j] = Vxx(_x[i], _y[j], _t, side);

			} else {
				_u[i][j] = _v[i][j] = _uxx[i][j] = _v[i][j] = 0.0;
			}
		}
	}


	_u_tau_past1.reallocate(_total_intersect_node_n);
	_u_tau_past2.reallocate(_total_intersect_node_n);

	_v_tau_past1.reallocate(_total_intersect_node_n);
	_v_tau_past2.reallocate(_total_intersect_node_n);

	extractTangentU(_u, _u_tau_past1);
	extractTangentU(_v, _v_tau_past1);

	_u_tau_past2 = _u_tau_past1;
	_v_tau_past2 = _v_tau_past1;

	//std::cout << "Solver is ready." << std::endl;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::F1(double u, double v) const
{
	return (_gama * (1.0 - u) - u * v * v) / _eps0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::F1_u(double u, double v) const
{
	return - (_gama + v * v) / _eps0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::F1_v(double u, double v) const
{
	return (- 2.0 * u * v) / _eps0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::F2(double u, double v) const
{
	return (u * v * v - (_gama + _kappa) * v) / _eps0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::F2_u(double u, double v) const
{
	return (v * v) / _eps0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::F2_v(double u, double v) const
{
	return (2.0 * u * v - (_gama + _kappa)) / _eps0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::U(double x, double y, double t, bool side) const
{
	double u;

	if (side) {
		u = _Ui(x, y, t);
	} else {
		u = _Ue(x, y, t);
	}

	return u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::Ux(double x, double y, double t, bool side) const
{
	const double delta = 1.0e-5;
	const double delta2 = delta + delta;

	double u0, u1;

	if (side) {

		u0 = _Ui(x - delta, y, t);
		u1 = _Ui(x + delta, y, t);

	} else {

		u0 = _Ue(x - delta, y, t);
		u1 = _Ue(x + delta, y, t);
	}

	double ux = (u1 - u0) / delta2;

	return ux;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::Uy(double x, double y, double t, bool side) const
{
	const double delta = 1.0e-5;
	const double delta2 = delta + delta;

	double u0, u1;

	if (side) {

		u0 = _Ui(x, y - delta, t);
		u1 = _Ui(x, y + delta, t);

	} else {

		u0 = _Ue(x, y - delta, t);
		u1 = _Ue(x, y + delta, t);
	}

	double uy = (u1 - u0) / delta2;

	return uy;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::Uxx(double x, double y, double t, bool side) const
{
	const double delta = 1.0e-5;
	const double delta2 = delta * delta;

	double u0, u1, u2;

	if (side) {

		u0 = _Ui(x - delta, y, t);
		u1 = _Ui(x, y, t);
		u2 = _Ui(x + delta, y, t);

	} else {

		u0 = _Ue(x - delta, y, t);
		u1 = _Ue(x, y, t);
		u2 = _Ue(x + delta, y, t);

	}

	double uxx = (u0 + u2 - u1 - u1) / delta2;
	return uxx;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::Uyy(double x, double y, double t, bool side) const
{
	const double delta = 1.0e-5;
	const double delta2 = delta * delta;

	double u0, u1, u2;

	if (side) {

		u0 = _Ui(x, y - delta, t);
		u1 = _Ui(x, y, t);
		u2 = _Ui(x, y + delta, t);

	} else {

		u0 = _Ue(x, y - delta, t);
		u1 = _Ue(x, y, t);
		u2 = _Ue(x, y + delta, t);
	}

	double uyy = (u0 + u2 - u1 - u1) / delta2;
	return uyy;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::V(double x, double y, double t, bool side) const
{
	double u;

	if (side) {
		u = _Vi(x, y, t);
	} else {
		u = _Ve(x, y, t);
	}

	return u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::Vx(double x, double y, double t, bool side) const
{
	const double delta = 1.0e-5;
	const double delta2 = delta + delta;

	double u0, u1;

	if (side) {

		u0 = _Vi(x - delta, y, t);
		u1 = _Vi(x + delta, y, t);

	} else {

		u0 = _Ve(x - delta, y, t);
		u1 = _Ve(x + delta, y, t);
	}

	double ux = (u1 - u0) / delta2;

	return ux;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::Vy(double x, double y, double t, bool side) const
{
	const double delta = 1.0e-5;
	const double delta2 = delta + delta;

	double u0, u1;

	if (side) {

		u0 = _Vi(x, y - delta, t);
		u1 = _Vi(x, y + delta, t);

	} else {

		u0 = _Ve(x, y - delta, t);
		u1 = _Ve(x, y + delta, t);
	}

	double uy = (u1 - u0) / delta2;

	return uy;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::Vxx(double x, double y, double t, bool side) const
{
	const double delta = 1.0e-5;
	const double delta2 = delta * delta;

	double u0, u1, u2;

	if (side) {

		u0 = _Vi(x - delta, y, t);
		u1 = _Vi(x, y, t);
		u2 = _Vi(x + delta, y, t);

	} else {

		u0 = _Ve(x - delta, y, t);
		u1 = _Ve(x, y, t);
		u2 = _Ve(x + delta, y, t);

	}

	double uxx = (u0 + u2 - u1 - u1) / delta2;
	return uxx;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double GrayScottEqn::Vyy(double x, double y, double t, bool side) const
{
	const double delta = 1.0e-5;
	const double delta2 = delta * delta;

	double u0, u1, u2;

	if (side) {

		u0 = _Vi(x, y - delta, t);
		u1 = _Vi(x, y, t);
		u2 = _Vi(x, y + delta, t);

	} else {

		u0 = _Ve(x, y - delta, t);
		u1 = _Ve(x, y, t);
		u2 = _Ve(x, y + delta, t);
	}

	double uyy = (u0 + u2 - u1 - u1) / delta2;
	return uyy;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::
solveByNewtonMethod_BE(double u_old, double v_old, double dt_2, 
											 double &u_new, double &v_new) const // BE method
{
	const double abs_tol = 1.0E-15;
	const double rel_tol = 1.0E-12;

	double rhs[2], sol[2], f[2];

	rhs[0] = u_old;
	rhs[1] = v_old;

	sol[0] = u_old + dt_2 * F1(u_old, v_old);
	sol[1] = v_old + dt_2 * F2(u_old, v_old);

	f[0] = sol[0] - dt_2 * F1(sol[0], sol[1]) - rhs[0];
	f[1] = sol[1] - dt_2 * F2(sol[0], sol[1]) - rhs[1];

	double res_norm = max(fabs(f[0]), fabs(f[1]));

	if(res_norm < abs_tol) {
		u_new = sol[0];
		v_new = sol[1];
	}

	double tol = abs_tol + res_norm * rel_tol;

	const int max_itr_num = 50;
	int itr_num = 0;

	double JF[2][2], d[2];

	do {

		JF[0][0] = 1.0 - dt_2 * F1_u(sol[0], sol[1]);
		JF[0][1] = - dt_2 * F1_v(sol[0], sol[1]);

		JF[1][0] = - dt_2 * F2_u(sol[0], sol[1]);
		JF[1][1] = 1.0 - dt_2 * F2_v(sol[0], sol[1]);

		double det = JF[0][0] * JF[1][1] - JF[0][1] * JF[1][0];

		d[0] = (JF[1][1] * f[0] - JF[0][1] * f[1]) / det;
		d[1] = (JF[0][0] * f[1] - JF[1][0] * f[0]) / det;

		sol[0] -= d[0];
		sol[1] -= d[1];

		f[0] = sol[0] - dt_2 * F1(sol[0], sol[1]) - rhs[0];
		f[1] = sol[1] - dt_2 * F2(sol[0], sol[1]) - rhs[1];

		res_norm = max(fabs(f[0]), fabs(f[1]));

		itr_num++;

		//std::cout << "itr_num = " << itr_num << ", res_norm = " << res_norm << std::endl;

	} while (res_norm > tol && itr_num < max_itr_num);

	if(itr_num == max_itr_num) {
		std::cout << "Newton's method failed to converge." << std::endl;
		exit(1);
	}

	u_new = sol[0];
	v_new = sol[1];

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::
solveByNewtonMethod_CN(double u_old, double v_old, double dt_2, 
											 double &u_new, double &v_new) const // CN method
{
	double dt_4 = 0.5 * dt_2;
	double u_mid = u_old + dt_4 * F1(u_old, v_old);
	double v_mid = v_old + dt_4 * F2(u_old, v_old);
	solveByNewtonMethod_BE(u_mid, v_mid, dt_4, u_new, v_new);

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::
getSplittedNeumannBC(int i_nd_idx, double t, double &ux, double &uy,
										 double &vx, double &vy) const
{
	double ut_2 = _u_tau_past2[i_nd_idx];
	double ut_1 = _u_tau_past1[i_nd_idx];
	double ut = ut_1 + ut_1 - ut_2;

	double vt_2 = _v_tau_past2[i_nd_idx];
	double vt_1 = _v_tau_past1[i_nd_idx];
	double vt = vt_1 + vt_1 - vt_2;

	double x = _intersect_node_crd[i_nd_idx][0];
	double y = _intersect_node_crd[i_nd_idx][1];

	double nx = _intersect_node_nml[i_nd_idx][0];
	double ny = _intersect_node_nml[i_nd_idx][1];

	double tx = _intersect_node_tan[i_nd_idx][0];
	double ty = _intersect_node_tan[i_nd_idx][1];

	double norm = sqrt(tx * tx + ty * ty);

	tx /= norm;
	ty /= norm;

	double un = 0.0;
	double vn = 0.0;

	splitNeumannBoundaryValue(tx, ty, ut, nx, ny, un, ux, uy);
	splitNeumannBoundaryValue(tx, ty, vt, nx, ny, vn, vx, vy);

	//std::cout << "tan = " << tx << ", " << ty << ", error = "
	//					<< fabs(ux - ux0) << ", " << fabs(uy - uy0) << std::endl;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::advanceByStrang_N(double t, double dt) // x-y-x
{
	double dt_2 = 0.5 * dt;
	double t_new = t + dt;
	double t_mid = t + dt_2;

	double nu[2] = {1.0 / (_eps1 * dt), 1.0 / (_eps2 * dt)};

	// x - forward Euler

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {

				_u[i][j] += dt_2 * _eps1 * _uxx[i][j];
				_v[i][j] += dt_2 * _eps2 * _vxx[i][j];

				double u0 = _u[i][j];
				double v0 = _v[i][j];

				_u[i][j] += dt_2 * F1(u0, v0);
				_v[i][j] += dt_2 * F2(u0, v0);
				//solveByNewtonMethod_CN(u0, v0, dt_2, _u[i][j], _v[i][j]);

			}
		}
	}


		// y - mid point

#pragma omp parallel for
	for(int i = 0; i < _I1; i++){

		if (_y_irreglar_line[i]) {

			double kappa[2];
			kappa[0] = nu[0] + nu[0];
			kappa[1] = nu[1] + nu[1];

			int bdry_num = _y_line_intersect_num[i];

			VectorXd u_bdry(bdry_num), uf(_J1), tmp_u(_J1);
			VectorXd v_bdry(bdry_num), vf(_J1), tmp_v(_J1);

			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					uf[j] = kappa[0] * _u[i][j];
					vf[j] = kappa[1] * _v[i][j];
				} else {
					uf[j] = vf[j] = 0.0;
				}
			}


			for(int j = 0, k = 0; j < _J; j++){
				if (_y_irregular_edge[i][j]) {

					int idx = _y_intersect_idx[i][j];

					assert(idx >= 0);

					double ux, uy, vx, vy;

					getSplittedNeumannBC(idx, t_new, ux, uy, vx, vy);

					u_bdry[k] = uy;
					v_bdry[k] = vy;
					k++;
				}
			}

			VectorX6d u_bdry_data(bdry_num), v_bdry_data(bdry_num);

			ReacDiffSolver1dN_D(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
													bdry_num, kappa[0], kappa[0], uf,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			ReacDiffSolver1dN_D(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
													bdry_num, kappa[1], kappa[1], vf,
													v_bdry, 0.0, 0.0, tmp_v, v_bdry_data);


			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					_work_u[i][j] = tmp_u[j];
					_work_v[i][j] = tmp_v[j];
					_u[i][j] = tmp_u[j] + tmp_u[j] - _u[i][j];
					_v[i][j] = tmp_v[j] + tmp_v[j] - _v[i][j];
				}
			}

		}

	}

	_u_tau_past2 = _u_tau_past1;
	extractTangentU(_work_u, _u_tau_past1);

	_v_tau_past2 = _v_tau_past1;
	extractTangentU(_work_v, _v_tau_past1);

		// reaction - backward Euler

#pragma omp parallel for collapse(2)
		for(int i = 0; i < _I1; i++){
			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					double u0 = _u[i][j];
					double v0 = _v[i][j];
					solveByNewtonMethod_BE(u0, v0, dt_2, _u[i][j], _v[i][j]);
					//solveByNewtonMethod_CN(u0, v0, dt_2, _u[i][j], _v[i][j]);
				}
			}
		}



	// x - backward Euler

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			double kappa[2];
			kappa[0] = nu[0] + nu[0];
			kappa[1] = nu[1] + nu[1];

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), uf(_I1), tmp_u(_I1);
			VectorXd v_bdry(bdry_num), vf(_I1), tmp_v(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					uf[i] = kappa[0] * _u[i][j];
					vf[i] = kappa[1] * _v[i][j];
				} else {
					vf[i] = uf[i] = 0.0;
				}
			}

			for(int i = 0, k = 0; i < _I; i++){
				if (_x_irregular_edge[i][j]) {

					int idx = _x_intersect_idx[i][j];
					assert(idx >= 0);
					double ux, uy, vx, vy;

					getSplittedNeumannBC(idx, t_mid, ux, uy, vx, vy);

					u_bdry[k] = ux;
					v_bdry[k] = vx;
					k++;
				}
			}

			VectorX6d u_bdry_data(bdry_num), v_bdry_data(bdry_num);

			ReacDiffSolver1dN_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
													bdry_num, kappa[0], kappa[0], uf,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			ReacDiffSolver1dN_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
													bdry_num, kappa[1], kappa[1], vf,
													v_bdry, 0.0, 0.0, tmp_v, v_bdry_data);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					_uxx[i][j] = (tmp_u[i] - _u[i][j]) / (dt_2 * _eps1);
					_vxx[i][j] = (tmp_v[i] - _v[i][j]) / (dt_2 * _eps2);

					_u[i][j] = tmp_u[i];
					_v[i][j] = tmp_v[i];
				}
			}

		}

	}

	_u_tau_past2 = _u_tau_past1;
	extractTangentU(_u, _u_tau_past1);

	_v_tau_past2 = _v_tau_past1;
	extractTangentU(_v, _v_tau_past1);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*void GrayScottEqn::advanceByADI(double t, double dt)
{
	double dt_2 = 0.5 * dt;
	double t_new = t + dt;
	double t_mid = t + dt_2;

	double nu = 1.0 / (_epsilon * dt);

	double kappa = nu + nu;

		// x - sweep

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), f(_I1), tmp_u(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					f[i] = kappa * _u[i][j] + _uyy[i][j];
				} else {
					f[i] = 0.0;
				}
			}

			for(int k = 0; k < bdry_num; k++){
				double x = _x_line_intersec_crd[j][k];

				u_bdry[k] = _Ui(x, _y[j], t_mid);
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1d3_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
													bdry_num, kappa, kappa, f,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			for(int i = 0; i < _I1; i++){
				//if (_interior[i][j]) {
					_u[i][j] = tmp_u[i];
					_uxx[i][j] = kappa * tmp_u[i] - f[i];
				//}
			}

		}

	}


	// y - sweep

#pragma omp parallel for
	for(int i = 0; i < _I1; i++){

		if (_y_irreglar_line[i]) {

			int bdry_num = _y_line_intersect_num[i];

			VectorXd u_bdry(bdry_num), f(_J1), tmp_u(_J1);

			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					f[j] = kappa * _u[i][j] + _uxx[i][j];
				} else {
					f[j] = 0.0;
				}
			}

			for(int k = 0; k < bdry_num; k++){
				double y = _y_line_intersec_crd[i][k];

				u_bdry[k] = _Ui(_x[i], y, t_new);
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1d3_D(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
													bdry_num, kappa, kappa, f,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);


			for(int j = 0; j < _J1; j++){
				//if (_interior[i][j]) {
					_u[i][j] = tmp_u[j];
					_uyy[i][j] = kappa * tmp_u[j] - f[j];
				//}
			}

		}

	}

}*/

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::
computeInteriorError(const MatrixXd &u, double &max_err, double &l2_err) const
{
	max_err = 0.0;
	l2_err = 0.0;

	int int_count= 0;

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {

				double ue = _Ui(_x[i], _y[j], _t);
				double e = fabs(u[i][j] - ue);
				max_err = max(max_err, e);
				l2_err += e * e;
				int_count++;

			}
		}
	}

	l2_err = sqrt(l2_err / int_count);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::advance(double t, double dt)
{
	advanceByStrang_N(_t, _dt);
	//advanceByADI(_t, _dt);

	_t += _dt;

	if ((++_step_count)%_check_interval == 0) {

		std::cout << "T = " << _T << ", t = " << _t 
							<< ", dt = " << _dt	<< std::endl;
	}


}

static int count0 = 0;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::solve(void)
{
	for(int k = 0; k < _max_step_num; k++){

#ifdef USE_MATHGL
		if (k%_plot_check_interval == 0) {
			plotToEPS();
			//std::cout << "saved." << std::endl;
		}
#endif

		advance(_t, _dt);

	}

#ifdef USE_MATHGL
	plotToEPS();
#endif
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef USE_MATHGL
void GrayScottEqn::plotToEPS(void) const
{
		mglInitBackgroundColor(1.0, 1.0, 1.0);

		mglNewPage(_t);

		mglColor(MGL_BLACK);
		mglSetTimeStamp(_t);

		mglColor(MGL_BLACK);
		mglLineWidth(2.0);
		mglPlot2dRectangle(_low[0], _low[1], _high[0], _high[1]);

		mglColor(MGL_BLACK);
		mglLineWidth(2.0);
		_curve->plotCurve();

		mglLineWidth(0.5);
		mglPlotInteriorNodeDataByIsolines(_x, _y, _interior, _v, _I, _J, 50);

		mglFlush(0.0);

}
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#ifdef USE_OPENGL

void GrayScottEqn::setPlotOption(int i)
{
	_plot_opt = i;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GrayScottEqn::plotSolution(void) const
{
	if (_plot_opt == 0) {

		//plotNodeDataByColormap(_u);
		plotInteriorNodeDataByColormap(_v);
		//plotExteriorNodeDataByColormap(_v);

	} else if (_plot_opt == 1) {

		//plotNodeDataByIsolines(_u);
		plotInteriorNodeDataByIsolines(_v);
		//plotExteriorNodeDataByIsolines(_v);

	}
}

#endif

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

