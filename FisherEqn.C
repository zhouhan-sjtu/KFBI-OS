/*=============================================================================
*   
*   Filename : FisherEqn.C
*   Creator : Han Zhou
*   Date : 11/20/21
*   Description : 
*
=============================================================================*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
 
#include "FisherEqn.H"

#include "ReacDiffSolver1d3-2.h"
#include "MathTools.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void FisherEqn::initialize(double t0, double T)
{
	_t = t0;
	_T = T;

	double deno = _I/16*5;
	//double deno = _I/16*20;

	//_dt = (_T - _t) / deno;
	_dt = 1.0 / deno;

	_check_interval = static_cast<int>(_print_info_interval/_dt + 0.5);
	_max_step_num = static_cast<int>((_T - _t) / _dt + 0.5);
	_step_count = 0;

	_u.reallocate(_I1, _J1);
	//_ux.reallocate(_I1, _J1);
	//_uy.reallocate(_I1, _J1);
	_uxx.reallocate(_I1, _J1);
	_uyy.reallocate(_I1, _J1);

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){

			bool side = _interior[i][j];

			_u[i][j] = U(_x[i], _y[j], _t, side);
			//_ux[i][j] = Ux(_x[i], _y[j], _t, side);
			//_uy[i][j] = Uy(_x[i], _y[j], _t, side);
			_uxx[i][j] = Uxx(_x[i], _y[j], _t, side);
			_uyy[i][j] = Uyy(_x[i], _y[j], _t, side);
		}
	}

	//std::cout << "Solver is ready." << std::endl;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double FisherEqn::F(double u) const
{
	return 2.0 * (1.0 - u) * u * u / _epsilon;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double FisherEqn::U(double x, double y, double t, bool side) const
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
double FisherEqn::Ux(double x, double y, double t, bool side) const
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
double FisherEqn::Uy(double x, double y, double t, bool side) const
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
double FisherEqn::Uxx(double x, double y, double t, bool side) const
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
double FisherEqn::Uyy(double x, double y, double t, bool side) const
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
double FisherEqn::solveByNewtonMethod_CN(double u_old, double dt_2) const 
{
	double dt_4 = 0.5 * dt_2;
	double u_mid = u_old + dt_4 * F(u_old);
	double u_new = solveByNewtonMethod_BE(u_mid, dt_4);
	return u_new;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double FisherEqn::solveByNewtonMethod_BE(double u_old, double dt_2) const 
{
	const double abs_tol = 1.0E-15;
	const double rel_tol = 1.0E-12;

	double nu = dt_2 / _epsilon;
	double nu2 = nu + nu;

	double x = u_old + dt_2 * F(u_old);

	double f = nu2 * (x - 1.0) * x * x + x - u_old;

	if(fabs(f) < abs_tol) {
		return x;
	}

	double tol = abs_tol + fabs(f) * rel_tol;

	const int max_itr_num = 50;
	int itr_num = 0;

	do {

		double df = nu * (6.0 * x * x - 4.0 * x) + 1.0;

		x += - f / df;

		f = nu2 * (x - 1.0) * x * x + x - u_old;

		itr_num++;

		//std::cout << itr_num << ", " << fabs(f) << std::endl;

	} while (fabs(f) > tol && itr_num < max_itr_num);

	if(itr_num == max_itr_num) {
		std::cout << "Newton's method failed to converge." << std::endl;
		exit(1);
		return u_old;
	}

	return x;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void FisherEqn::advanceByStrang(double t, double dt) // x-y-x
{
	double dt_2 = 0.5 * dt;
	double t_new = t + dt;
	double t_mid = t + dt_2;

	double nu = 1.0 / (_epsilon * dt);

	double kappa = nu + nu;

	// y - forward Euler +  reaction - forward Euler

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {

				_u[i][j] += dt_2 * _epsilon * _uyy[i][j];
				double f = F(_u[i][j]);
				_u[i][j] += dt_2 * f;
				//double u0 = _u[i][j];
				//_u[i][j] = solveByNewtonMethod_CN(u0, dt_2);

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

	// reaction - backward Euler

#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {
				double u0 = _u[i][j];
				_u[i][j] = solveByNewtonMethod_BE(u0, dt_2);
				//_u[i][j] = solveByNewtonMethod_CN(u0, dt_2);
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

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*void FisherEqn::advanceByStrang(double t, double dt) // x-y-x
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
void FisherEqn::advanceByADI(double t, double dt)
{
	double dt_2 = 0.5 * dt;
	double t_new = t + dt;
	double t_mid = t + dt_2;

	double nu = 1.0 / (_epsilon * dt);

	double kappa = nu + nu;

	MatrixXd work(_I1, _J1), work2(_I1, _J1);

	MatrixXd uxx_new(_I1, _J1);
	MatrixXd uyy_new(_I1, _J1);


#pragma omp parallel for collapse(2)
	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {

				double tmp = _u[i][j] + dt * _epsilon * (_uxx[i][j] + _uyy[i][j]) 
									 + dt_2 * F(_u[i][j]);
				_u[i][j] = solveByNewtonMethod_BE(tmp, dt_2);

			}
		}
	}

	work = _u;

		// x - sweep

#pragma omp parallel for
	for(int j = 0; j < _J1; j++){

		if (_x_irreglar_line[j]) {

			int bdry_num = _x_line_intersect_num[j];

			VectorXd u_bdry(bdry_num), f(_I1), tmp_u(_I1);

			for(int i = 0; i < _I1; i++){
				if (_interior[i][j]) {
					f[i] = kappa * _u[i][j] - _uxx[i][j];
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
					_u[i][j] = tmp_u[i];
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
					f[j] = kappa * _u[i][j] - _uyy[i][j];
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
					_u[i][j] = tmp_u[j];
					uyy_new[i][j] = kappa * tmp_u[j] - f[j];
				}
			}

		}

	}

	///////////////////////////////////////////////////////////////////// 

	_u = work;

	// y - sweep

#pragma omp parallel for
	for(int i = 0; i < _I1; i++){

		if (_y_irreglar_line[i]) {

			int bdry_num = _y_line_intersect_num[i];

			VectorXd u_bdry(bdry_num), f(_J1), tmp_u(_J1);

			for(int j = 0; j < _J1; j++){
				if (_interior[i][j]) {
					f[j] = kappa * _u[i][j] - _uyy[i][j];
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
					_u[i][j] = tmp_u[j];
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
					f[i] = kappa * _u[i][j] - _uxx[i][j];
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
					_u[i][j] = tmp_u[i];
					uxx_new[i][j] = kappa * tmp_u[i] - f[i];
				}
			}

		}

	}

	_uxx = uxx_new;
	_uyy = uyy_new;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void FisherEqn::
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
void FisherEqn::advance(double t, double dt)
{
	//advanceByADI(_t, _dt);
	advanceByStrang(_t, _dt);

	_t += _dt;

	if ((++_step_count)%_check_interval == 0) {

		double max_err = 0.0, l2_err = 0.0;

		computeInteriorError(_u, max_err, l2_err);

		std::cout << "T = " << _T << ", t = " << _t << ", dt = " << _dt	
							<< ", max-err = " << max_err << " l2-err = "
							<< l2_err << std::endl;
	}


}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void FisherEqn::solve(void)
{
	for(int k = 0; k < _max_step_num; k++){
		advance(_t, _dt);
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef USE_OPENGL

void FisherEqn::setPlotOption(int i)
{
	_plot_opt = i;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void FisherEqn::plotSolution(void) const
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

