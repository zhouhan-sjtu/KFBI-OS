/*=============================================================================
*   
*   Filename : WaveEqn.C
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
 
#include "WaveEqn.H"

#include "ReacDiffSolver1d3-2.h"
#include "ReacDiffSolver1dN-D.h"
#include "DiffusionSolver1d-001-O4.h"
#include "MathTools.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void WaveEqn::initialize(double t0, double T)
{
	_t = t0;
	_T = T;

	double deno = _I/16*5;
	//double deno = _I/16*10;

	//_dt = (_T - _t) / deno;
	_dt = 1.0 / deno;


	_check_interval = static_cast<int>(_print_info_interval/_dt + 0.5);
	_max_step_num = static_cast<int>((_T - _t) / _dt + 0.5);
	_step_count = 0;


	_u.reallocate(_I1, _J1);
	_u_past.reallocate(_I1, _J1);

	_work1.reallocate(_I1, _J1);
	_work2.reallocate(_I1, _J1);

	_uxx.reallocate(_I1, _J1);
	_uxx_past.reallocate(_I1, _J1);
	_uxx_work.reallocate(_I1, _J1);

	_uyy.reallocate(_I1, _J1);
	_uyy_past.reallocate(_I1, _J1);
	_uyy_work.reallocate(_I1, _J1);

	_f_work.reallocate(_I1, _J1);

	double t_past = _t - _dt;

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){

			bool side = _interior[i][j];

			if (side) {
				_u[i][j] = U(_x[i], _y[j], _t, side);
				_uxx[i][j] = Uxx(_x[i], _y[j], _t, side);
				_uyy[i][j] = Uyy(_x[i], _y[j], _t, side);

				_u_past[i][j] = U(_x[i], _y[j], t_past, side);
				_uxx_past[i][j] = Uxx(_x[i], _y[j], t_past, side);
				_uyy_past[i][j] = Uyy(_x[i], _y[j], t_past, side);

			} else {

				_u[i][j] = 0.0;
				_uxx[i][j] = 0.0;
				_uyy[i][j] = 0.0;

				_u_past[i][j] = 0.0;
				_uxx_past[i][j] = 0.0;
				_uyy_past[i][j] = 0.0;

			}

		}
	}

	_u_tau_past1.reallocate(_total_intersect_node_n);
	_u_tau_past2.reallocate(_total_intersect_node_n);

	//for(int k = 0; k < _total_intersect_node_n; k++){

	//	double x = _intersect_node_crd[k][0];
	//	double y = _intersect_node_crd[k][1];

	//	double tx = _intersect_node_tan[k][0];
	//	double ty = _intersect_node_tan[k][1];

	//	double norm = sqrt(tx * tx + ty * ty);
	//	tx /= norm;
	//	ty /= norm;
	//	double ux = Ux(x, y, _t, true);
	//	double uy = Uy(x, y, _t, true);

	//	_u_tau_past2[k] = tx * ux + ty * uy;
	//}

	extractTangentU(_u, _u_tau_past1);

	double err1 = _u_tau_past1.max_error(_u_tau_past2);
	double err2 = _u_tau_past1.l2_error_scaled(_u_tau_past2);

	//std::cout << "ut error = " << err1 << ", " << err2 << std::endl;

	_u_tau_past2 = _u_tau_past1;

	//std::cout << "Solver is ready." << std::endl;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double WaveEqn::U(double x, double y, double t, bool side) const
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
double WaveEqn::Ux(double x, double y, double t, bool side) const
{
	const double h = 1.0e-5;
	const double h2 = h + h;

	double u0, u1;

	if (side) {

		u0 = _Ui(x - h, y, t);
		u1 = _Ui(x + h, y, t);

	} else {

		u0 = _Ue(x - h, y, t);
		u1 = _Ue(x + h, y, t);
	}

	double ux = (u1 - u0) / h2;

	return ux;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double WaveEqn::Uy(double x, double y, double t, bool side) const
{
	const double h = 1.0e-5;
	const double h2 = h + h;

	double u0, u1;

	if (side) {

		u0 = _Ui(x, y - h, t);
		u1 = _Ui(x, y + h, t);

	} else {

		u0 = _Ue(x, y - h, t);
		u1 = _Ue(x, y + h, t);
	}

	double uy = (u1 - u0) / h2;

	return uy;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double WaveEqn::Uxx(double x, double y, double t, bool side) const
{
	const double h = 1.0e-5;
	const double h2 = h * h;

	double u0, u1, u2;

	if (side) {

		u0 = _Ui(x - h, y, t);
		u1 = _Ui(x, y, t);
		u2 = _Ui(x + h, y, t);

	} else {

		u0 = _Ue(x - h, y, t);
		u1 = _Ue(x, y, t);
		u2 = _Ue(x + h, y, t);

	}

	double uxx = (u0 + u2 - u1 - u1) / h2;
	return uxx;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double WaveEqn::Uyy(double x, double y, double t, bool side) const
{
	const double h = 1.0e-5;
	const double h2 = h * h;

	double u0, u1, u2;

	if (side) {

		u0 = _Ui(x, y - h, t);
		u1 = _Ui(x, y, t);
		u2 = _Ui(x, y + h, t);

	} else {

		u0 = _Ue(x, y - h, t);
		u1 = _Ue(x, y, t);
		u2 = _Ue(x, y + h, t);
	}

	double uyy = (u0 + u2 - u1 - u1) / h2;
	return uyy;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void WaveEqn::getSplittedNeumannBC(int i_nd_idx, double t, 
																				double &ux, double &uy) const
{
	double ut_2 = _u_tau_past2[i_nd_idx];
	double ut_1 = _u_tau_past1[i_nd_idx];

	double ut = ut_1 + ut_1 - ut_2;
	//double ut = ut_1;

	double x = _intersect_node_crd[i_nd_idx][0];
	double y = _intersect_node_crd[i_nd_idx][1];

	double nx = _intersect_node_nml[i_nd_idx][0];
	double ny = _intersect_node_nml[i_nd_idx][1];

	double tx = _intersect_node_tan[i_nd_idx][0];
	double ty = _intersect_node_tan[i_nd_idx][1];

	double norm = sqrt(tx * tx + ty * ty);

	tx /= norm;
	ty /= norm;

	double ux0 = Ux(x, y, t, true);
	double uy0 = Uy(x, y, t, true);

	double un = nx * ux0 + ny * uy0;

	splitNeumannBoundaryValue(tx, ty, ut, nx, ny, un, ux, uy);

	//std::cout << "tan = " << tx << ", " << ty << ", error = "
	//					<< fabs(ux - ux0) << ", " << fabs(uy - uy0) << std::endl;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*void WaveEqn::advanceByStrang(double t, double dt)
{
	double dt_2 = 0.5 * dt;
	double t_new = t + dt;
	double t_mid = t + dt_2;

	double nu = 1.0 / (_epsilon * dt);

	double kappa = nu + nu;

	// y - forward Euler

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {
				_u[i][j] += dt_2 * _epsilon * _uyy[i][j];
			}
		}
	}


		// x - backward Euler

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), f(_I1), tmp_u(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					f[i] = kappa * _u[i][j];
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
				if (_interior[i][j]) {
					_u[i][j] = tmp_u[i] + tmp_u[i] - _u[i][j];
				}
			}

		}

	}


	// y - backward Euler

#pragma omp parallel for
	for(int i = 0; i < _I1; i++){

		if (_y_irreglar_line[i]) {

			int bdry_num = _y_line_intersect_num[i];

			VectorXd u_bdry(bdry_num), f(_J1), tmp_u(_J1);

			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					f[j] = kappa * _u[i][j];
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
				if (_interior[i][j]) {
					_uyy[i][j] = (tmp_u[j] - _u[i][j]) / (dt_2 * _epsilon);
					_u[i][j] = tmp_u[j];
				}
			}

		}

	}

}*/

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void WaveEqn::advanceByADI0(double t, double dt)
{
	//std::cout << "Pass : " << __FILE__ << ", line : " << __LINE__ << std::endl;

	double dt_2 = 0.5 * dt;
	double t_new = t + dt;
	double t_mid = t + dt_2;

	double dt2 = dt * dt;

	double kappa = 1.0 / (dt2 * _eta * _epsilon);
	double one_eta2 = 1.0 - (_eta+_eta);

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {

				double s0 = _F(_x[i], _y[j], t-dt);
				double s = _F(_x[i], _y[j], t);
				double s1 = _F(_x[i], _y[j], t+dt);

				_f_work[i][j] = _eta * s0 + one_eta2 * s + _eta * s1;

			}
		}
	}


	// x - sweep

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), f(_I1), tmp_u(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {

					double f0 = dt2 * _epsilon * (one_eta2 * _uxx[i][j] + _eta * _uxx_past[i][j]);
					double f1 = _u[i][j] + _u[i][j] - _u_past[i][j];
					double f2 = dt2 * _epsilon * _uyy[i][j] + dt2 * _f_work[i][j];
					f[i] = kappa * (f0 + f1 + f2);

				} else {
					f[i] = 0.0;
				}
			}

			for(int k = 0; k < bdry_num; k++){
				double x = _x_line_intersec_crd[j][k];

				u_bdry[k] = _Ui(x, _y[j], t_new);
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1d3_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
													bdry_num, kappa, kappa, f,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					//_work1[i][j] = tmp_u[i];
					double tmp  =	tmp_u[i] - 2.0 * _u[i][j] + _u_past[i][j] 
											- dt2 * _epsilon * _uyy[i][j] - dt2 * _f_work[i][j];
					_work1[i][j] = tmp;
				}
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
					double f0 = dt2 * _epsilon * (one_eta2 * _uyy[i][j] + _eta * _uyy_past[i][j]);
					double f1 = _work1[i][j] + (_u[i][j] + _u[i][j] - _u_past[i][j]) + dt2 * _f_work[i][j];
					f[j] = kappa * (f0 + f1);
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
				if (_interior[i][j]) {
					_work1[i][j] = tmp_u[j];
					_uyy_work[i][j] = kappa * tmp_u[j] - f[j];
				}
			}

		}

	}

	///////////////////////////////////////////////////////////////////// 

	// y - sweep

#pragma omp parallel for
	for(int i = 0; i < _I1; i++){

		if (_y_irreglar_line[i]) {

			int bdry_num = _y_line_intersect_num[i];

			VectorXd u_bdry(bdry_num), f(_J1), tmp_u(_J1);

			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					//f[j] = kappa * _work2[i][j] + (_uyy_past[i][j] - 2.0 * _uyy[i][j]);
					double f0 = dt2 * _epsilon * (one_eta2 * _uyy[i][j] + _eta * _uyy_past[i][j]);
					double f1 = _u[i][j] + _u[i][j] - _u_past[i][j];
					double f2 = dt2 * _epsilon * _uxx[i][j] + dt2 * _f_work[i][j];
					f[j] = kappa * (f0 + f1 + f2);
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
				if (_interior[i][j]) {
					//_work2[i][j] = tmp_u[j];
					double tmp  =	tmp_u[j] - 2.0 * _u[i][j] + _u_past[i][j] 
											- dt2 * _epsilon * _uxx[i][j] - dt2 * _f_work[i][j];
					_work2[i][j] = tmp;
				}
			}

		}

	}

		// x - sweep

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), f(_I1), tmp_u(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					//f[i] = kappa * _work2[i][j] + (_uxx_past[i][j] - 2.0 * _uxx[i][j]);
					double f0 = dt2 * _epsilon * (one_eta2 * _uxx[i][j] + _eta * _uxx_past[i][j]);
					double f1 = _work2[i][j] + (_u[i][j] + _u[i][j] - _u_past[i][j]) + dt2 * _f_work[i][j];
					f[i] = kappa * (f0 + f1);
				} else {
					f[i] = 0.0;
				}
			}

			for(int k = 0; k < bdry_num; k++){
				double x = _x_line_intersec_crd[j][k];

				u_bdry[k] = _Ui(x, _y[j], t_new);
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1d3_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
													bdry_num, kappa, kappa, f,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					_work2[i][j] = tmp_u[i];
					_uxx_work[i][j] = kappa * tmp_u[i] - f[i];
				}
			}

		}

	}

	///////////////////////////////////////////////////////////////////// 

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {

			_u_past[i][j] = _u[i][j];
			_u[i][j] = 0.5 * (_work1[i][j] + _work2[i][j]);

			_uxx_past[i][j] = _uxx[i][j];
			_uxx[i][j] = _uxx_work[i][j];

			_uyy_past[i][j] = _uyy[i][j];
			_uyy[i][j] = _uyy_work[i][j];

			}
		}
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void WaveEqn::advanceByADI(double t, double dt)
{
	double dt_2 = 0.5 * dt;
	double t_new = t + dt;
	double t_mid = t + dt_2;

	double dt2 = dt * dt;

	double kappa = 1.0 / (dt2 * _eta * _epsilon);

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {

				double tmp = dt2 * _epsilon * (_uxx[i][j] + _uyy[i][j]);
				double tmp_u = _u[i][j] + _u[i][j] - _u_past[i][j] + tmp;

				double s = _F(_x[i], _y[j], t);
				double s0 = _F(_x[i], _y[j], t-dt);
				double s1 = _F(_x[i], _y[j], t+dt);
				double stt = (s1 + s0 - s - s) / dt2;
				//double rhs = dt2 * (s + _eta * dt2);
				double rhs = dt2 * (s + _eta * dt2 * stt);

				_work1[i][j] = _work2[i][j] = tmp_u + rhs;

			}
		}
	}


	// x - sweep

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), f(_I1), tmp_u(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					f[i] = kappa * _work1[i][j] + (_uxx_past[i][j] - 2.0 * _uxx[i][j]);
				} else {
					f[i] = 0.0;
				}
			}

			for(int k = 0; k < bdry_num; k++){
				double x = _x_line_intersec_crd[j][k];

				u_bdry[k] = _Ui(x, _y[j], t_new);
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1d3_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
													bdry_num, kappa, kappa, f,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);
			//DiffusionSolver1d_001_O4(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
			//										bdry_num, kappa, f,
			//										u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					_work1[i][j] = tmp_u[i];
				}
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
					f[j] = kappa * _work1[i][j] + (_uyy_past[i][j] - 2.0 * _uyy[i][j]);
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
			//DiffusionSolver1d_001_O4(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
			//										bdry_num, kappa, f,
			//										u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);


			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					_work1[i][j] = tmp_u[j];
					_uyy_work[i][j] = kappa * tmp_u[j] - f[j];
				}
			}

		}

	}

	///////////////////////////////////////////////////////////////////// 

	// y - sweep

#pragma omp parallel for
	for(int i = 0; i < _I1; i++){

		if (_y_irreglar_line[i]) {

			int bdry_num = _y_line_intersect_num[i];

			VectorXd u_bdry(bdry_num), f(_J1), tmp_u(_J1);

			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					f[j] = kappa * _work2[i][j] + (_uyy_past[i][j] - 2.0 * _uyy[i][j]);
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
			//DiffusionSolver1d_001_O4(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
			//										bdry_num, kappa, f,
			//										u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);


			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					_work2[i][j] = tmp_u[j];
				}
			}

		}

	}

		// x - sweep

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), f(_I1), tmp_u(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					f[i] = kappa * _work2[i][j] + (_uxx_past[i][j] - 2.0 * _uxx[i][j]);
				} else {
					f[i] = 0.0;
				}
			}

			for(int k = 0; k < bdry_num; k++){
				double x = _x_line_intersec_crd[j][k];

				u_bdry[k] = _Ui(x, _y[j], t_new);
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1d3_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
													bdry_num, kappa, kappa, f,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);
			//DiffusionSolver1d_001_O4(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
			//										bdry_num, kappa, f,
			//										u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					_work2[i][j] = tmp_u[i];
					_uxx_work[i][j] = kappa * tmp_u[i] - f[i];
				}
			}

		}

	}

	///////////////////////////////////////////////////////////////////// 

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {

			_u_past[i][j] = _u[i][j];
			_u[i][j] = 0.5 * (_work1[i][j] + _work2[i][j]);

			_uxx_past[i][j] = _uxx[i][j];
			_uxx[i][j] = _uxx_work[i][j];

			_uyy_past[i][j] = _uyy[i][j];
			_uyy[i][j] = _uyy_work[i][j];
			}
		}
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*void WaveEqn::advanceByStrang_N(double t, double dt)
{
	double dt_2 = 0.5 * dt;
	double t_new = t + dt;
	double t_mid = t + dt_2;

	double nu = 1.0 / (_epsilon * dt);

	double kappa = nu + nu;

	// y - forward Euler

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {
				_u[i][j] += dt_2 * _epsilon * _uyy[i][j];
			}
		}
	}


		// x - backward Euler

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), f(_I1), tmp_u(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					f[i] = kappa * _u[i][j];
				} else {
					f[i] = 0.0;
				}
			}

			//for(int k = 0; k < bdry_num; k++){
			//	double x = _x_line_intersec_crd[j][k];

			//	u_bdry[k] = _Ui(x, _y[j], t_mid);
			//}
			for(int i = 0, k = 0; i < _I; i++){
				if (_x_irregular_edge[i][j]) {

					int idx = _x_intersect_idx[i][j];
					assert(idx >= 0);
					double ux, uy;

					getSplittedNeumannBC(idx, t_mid, ux, uy);

					u_bdry[k++] = ux;
				}
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1dN_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
													bdry_num, kappa, kappa, f,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					_work[i][j] = tmp_u[i];
					_u[i][j] = tmp_u[i] + tmp_u[i] - _u[i][j];
				}
			}

		}

	}

	_u_tau_past2 = _u_tau_past1;
	extractTangentU(_work, _u_tau_past1);


	// y - backward Euler

#pragma omp parallel for
	for(int i = 0; i < _I1; i++){

		if (_y_irreglar_line[i]) {

			int bdry_num = _y_line_intersect_num[i];

			VectorXd u_bdry(bdry_num), f(_J1), tmp_u(_J1);

			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					f[j] = kappa * _u[i][j];
				} else {
					f[j] = 0.0;
				}
			}

			//for(int k = 0; k < bdry_num; k++){
			//	double y = _y_line_intersec_crd[i][k];

			//	u_bdry[k] = _Ui(_x[i], y, t_new);
			//}
			for(int j = 0, k = 0; j < _J; j++){
				if (_y_irregular_edge[i][j]) {

					int idx = _y_intersect_idx[i][j];
					assert(idx >= 0);
					double ux, uy;

					getSplittedNeumannBC(idx, t_new, ux, uy);

					u_bdry[k++] = uy;
				}
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1dN_D(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
													bdry_num, kappa, kappa, f,
													u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);


			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					_uyy[i][j] = (tmp_u[j] - _u[i][j]) / (dt_2 * _epsilon);
					_u[i][j] = tmp_u[j];
				}
			}

		}

	}

	_u_tau_past2 = _u_tau_past1;
	extractTangentU(_u, _u_tau_past1);

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void WaveEqn::advanceByADI_N(double t, double dt)//Neumann BC
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

			//for(int k = 0; k < bdry_num; k++){
			//	double x = _x_line_intersec_crd[j][k];

			//	//u_bdry[k] = _Ui(x, _y[j], t_mid);
			//	u_bdry[k] = Ux(x, _y[j], t_mid, true);
			//}

			for(int i = 0, k = 0; i < _I; i++){
				if (_x_irregular_edge[i][j]) {

					int idx = _x_intersect_idx[i][j];
					assert(idx >= 0);
					double ux, uy;

					getSplittedNeumannBC(idx, t_mid, ux, uy);

					u_bdry[k++] = ux;
				}
			}


			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1dN_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
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

	_u_tau_past2 = _u_tau_past1;
	extractTangentU(_u, _u_tau_past1);

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

			//for(int k = 0; k < bdry_num; k++){
			//	double y = _y_line_intersec_crd[i][k];

			//	//u_bdry[k] = _Ui(_x[i], y, t_new);
			//	u_bdry[k] = Uy(_x[i], y, t_new, true);
			//}

			for(int j = 0, k = 0; j < _J; j++){
				if (_y_irregular_edge[i][j]) {

					int idx = _y_intersect_idx[i][j];
					assert(idx >= 0);
					double ux, uy;

					getSplittedNeumannBC(idx, t_new, ux, uy);

					u_bdry[k++] = uy;
				}
			}

			VectorX6d u_bdry_data(bdry_num);

			ReacDiffSolver1dN_D(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
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

	_u_tau_past2 = _u_tau_past1;
	extractTangentU(_u, _u_tau_past1);

}*/

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void WaveEqn::
computeInteriorError(const MatrixXd &u, double &max_err, double &l2_err) const
{
	max_err = 0.0;
	l2_err = 0.0;

	int int_count = 0;

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
void WaveEqn::advance(double t, double dt)
{
	//advanceByStrang(_t, _dt);
	advanceByADI0(_t, _dt); // ADI
	//advanceByADI(_t, _dt); // LOD

	//advanceByStrang_N(_t, _dt);
	//advanceByADI_N(_t, _dt);

	_t += _dt;

	if ((++_step_count)%_check_interval == 0) {

		double max_err = 0.0, l2_err = 0.0;
		computeInteriorError(_u, max_err, l2_err);

		std::cout << "T = " << _T << ", t = " << _t << ", dt = " << _dt	
							<< ", max-err = " << max_err << ", l2-err = "
							<< l2_err << std::endl;
	}


}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void WaveEqn::solve(void)
{
	for(int k = 0; k < _max_step_num; k++){
		advance(_t, _dt);
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef USE_OPENGL

void WaveEqn::setPlotOption(int i)
{
	_plot_opt = i;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void WaveEqn::plotSolution(void) const
{
	if (_plot_opt == 0) {

		//plotNodeDataByColormap(_u);
		plotInteriorNodeDataByColormap(_u);
		//plotExteriorNodeDataByColormap(_u);

	} else if (_plot_opt == 1) {

		//plotNodeDataByIsolines(_u);
		plotInteriorNodeDataByIsolines(_u);
		//plotExteriorNodeDataByIsolines(_u);

	}
}

#endif

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

