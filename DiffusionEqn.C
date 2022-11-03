/*=============================================================================
*   
*   Filename : DiffusionEqn.C
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
 
#include "DiffusionEqn.H"

#include "DiffusionSolver1d-001-O4.h"
#include "ReacDiffSolver1d3-2.h"
#include "ReacDiffSolver1dN-D.h"
#include "ReacDiffSolver1dD-D-D.h"
#include "MathTools.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DiffusionEqn::initialize(double t0, double T)
{
	_t = t0;
	_T = T;


	//double deno = _I*2;
	double deno = _I/16*5;
	//double deno = _I/16*80;
	//double deno_2 = _I/64*10;
	//double deno = deno_2 * deno_2;

	//_dt = (_T - _t) / deno;
	_dt = 1.0 / deno;

	double dt0 = 0.5 * (_h * _h) / _epsilon;

	//std::cout << "dt / dx = " << _dt / _h << std::endl;
	//std::cout << "dt / dt_{ex} ratio = " << _dt / dt0 << std::endl;

	_check_interval = static_cast<int>(_print_info_interval/_dt + 0.5);
	_max_step_num = static_cast<int>((_T - _t) / _dt + 0.5);
	_step_count = 0;

	_u.reallocate(_I1, _J1);
	_ux.reallocate(_I1, _J1);
	_uy.reallocate(_I1, _J1);
	_uxx.reallocate(_I1, _J1);
	_uyy.reallocate(_I1, _J1);

	_work.reallocate(_I1, _J1);

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){

			bool side = _interior[i][j];

			_u[i][j] = U(_x[i], _y[j], _t, side);
			_ux[i][j] = Ux(_x[i], _y[j], _t, side);
			_uy[i][j] = Uy(_x[i], _y[j], _t, side);
			_uxx[i][j] = Uxx(_x[i], _y[j], _t, side);
			_uyy[i][j] = Uyy(_x[i], _y[j], _t, side);
		}
	}

	_u_tau_past1.reallocate(_total_intersect_node_n);
	_u_tau_past2.reallocate(_total_intersect_node_n);

	extractTangentU(_u, _u_tau_past1);
	_u_tau_past2 = _u_tau_past1;


	_u_x_past1.reallocate(_total_intersect_node_n);
	_u_x_past2.reallocate(_total_intersect_node_n);
	_u_y_past1.reallocate(_total_intersect_node_n);
	_u_y_past2.reallocate(_total_intersect_node_n);

	extractBoundryData(_u, _u_x_past1, _u_y_past1);
	_u_x_past2 = _u_x_past1;
	_u_y_past2 = _u_y_past1;


	//std::cout << "Solver is ready." << std::endl;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double DiffusionEqn::U(double x, double y, double t, bool side) const
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
double DiffusionEqn::Ux(double x, double y, double t, bool side) const
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
double DiffusionEqn::Uy(double x, double y, double t, bool side) const
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
double DiffusionEqn::Uxx(double x, double y, double t, bool side) const
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
double DiffusionEqn::Uyy(double x, double y, double t, bool side) const
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
void DiffusionEqn::getSplittedNeumannBC(int i_nd_idx, double t, 
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
void DiffusionEqn::getSplittedNeumannBC2(int i_nd_idx, double t, 
																				double &ux, double &uy) const
{

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

	const double thresh = M_PI / 3.0;

	if (i_nd_idx < _intersect_node_n[0]) { // x intersection node

		if (fabs(nx) > cos(thresh)) {

			double uy1 = _u_y_past1[i_nd_idx] + _u_y_past1[i_nd_idx] - _u_y_past2[i_nd_idx];
			ux = (un - ny * uy1) / nx;

		} else {
			ux = _u_x_past1[i_nd_idx] + _u_x_past1[i_nd_idx] - _u_x_past2[i_nd_idx];
		}

	} else {

		if (fabs(ny) > cos(thresh)) {

			double ux1 = _u_x_past1[i_nd_idx] + _u_x_past1[i_nd_idx] - _u_x_past2[i_nd_idx];
			uy = (un - nx * ux1) / ny;

		} else {
			uy = _u_y_past1[i_nd_idx] + _u_y_past1[i_nd_idx] - _u_y_past2[i_nd_idx];
		}
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DiffusionEqn::advanceByStrang(double t, double dt)
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

			//ReacDiffSolver1d3_D(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
			//										bdry_num, kappa, kappa, f,
			//										u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);
			ReacDiffSolver1d_DDD(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
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

			//ReacDiffSolver1d3_D(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
			//										bdry_num, kappa, kappa, f,
			//										u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);
			ReacDiffSolver1d_DDD(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
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

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DiffusionEqn::advanceByADI(double t, double dt)
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
			//DiffusionSolver1d_001_O4(_x, _x_line_side[j], _I, _x_line_intersec_crd[j],
			//												 bdry_num, kappa, f,
			//												 u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					_u[i][j] = tmp_u[i];
					_uxx[i][j] = kappa * tmp_u[i] - f[i];
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
			//DiffusionSolver1d_001_O4(_y, _y_line_side[i], _J, _y_line_intersec_crd[i],
			//												 bdry_num, kappa, f,
			//												 u_bdry, 0.0, 0.0, tmp_u, u_bdry_data);


			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					_u[i][j] = tmp_u[j];
					_uyy[i][j] = kappa * tmp_u[j] - f[j];
				}
			}

		}

	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DiffusionEqn::advanceByStrang_N(double t, double dt)
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
void DiffusionEqn::advanceByADI_N(double t, double dt)//Neumann BC
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

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DiffusionEqn::advanceByADI_N_2(double t, double dt)//Neumann BC
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

					getSplittedNeumannBC2(idx, t_mid, ux, uy);

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

	_u_x_past2 = _u_x_past1;
	_u_y_past2 = _u_y_past1;
	extractBoundryData(_u, _u_x_past1, _u_y_past1);

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

					getSplittedNeumannBC2(idx, t_new, ux, uy);

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

	_u_x_past2 = _u_x_past1;
	_u_y_past2 = _u_y_past1;
	extractBoundryData(_u, _u_x_past1, _u_y_past1);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DiffusionEqn::
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
void DiffusionEqn::advance(double t, double dt)
{
	//advanceByADI(_t, _dt);
	advanceByStrang(_t, _dt);

	//advanceByStrang_N(_t, _dt);
	//advanceByADI_N(_t, _dt);
	//advanceByADI_N_2(_t, _dt);

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
void DiffusionEqn::solve(void)
{
	for(int k = 0; k < _max_step_num; k++){
		advance(_t, _dt);
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef USE_OPENGL

void DiffusionEqn::setPlotOption(int i)
{
	_plot_opt = i;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DiffusionEqn::plotSolution(void) const
{
	if (_plot_opt == 0) {

		//plotNodeDataByColormap(_u);
		plotInteriorNodeDataByColormap(_u);
		plotExteriorNodeDataByColormap(_u);

	} else if (_plot_opt == 1) {

		//plotNodeDataByIsolines(_u);
		plotInteriorNodeDataByIsolines(_u);
		plotExteriorNodeDataByIsolines(_u);

	}
}

#endif

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

