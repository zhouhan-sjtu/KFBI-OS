/*=============================================================================
*   
*   Filename : GrayScottEqn.H
*   Creator : Han Zhou
*   Date : 11/20/21
*   Description : 
*
=============================================================================*/
 
 
#ifndef _GRAYSCOTTEQN_H
#define _GRAYSCOTTEQN_H


#include "CartesianGrid.H"
#include "Variables.h"

class GrayScottEqn : public CartesianGrid
{
	public :
		GrayScottEqn(void);
		virtual ~GrayScottEqn(void);

		GrayScottEqn(const CartesianGrid &g);

	public :

#ifdef USE_OPENGL
		void setPlotOption(int);
		void plotSolution(void) const;
#endif

#ifdef USE_MATHGL
		void plotToEPS(void) const;
#endif

	public :

		void setCoefficient(double, double, double, double, double);

		void setExactSolution(double (*U)(double, double, double), 
													double (*V)(double, double, double));

		double F1(double u, double v) const;
		double F1_u(double u, double v) const;
		double F1_v(double u, double v) const;
		double F2(double u, double v) const;
		double F2_u(double u, double v) const;
		double F2_v(double u, double v) const;

		double U(double, double, double, bool) const;
		double Ux(double, double, double, bool) const;
		double Uy(double, double, double, bool) const;
		double Uxx(double, double, double, bool) const;
		double Uyy(double, double, double, bool) const;

		double V(double, double, double, bool) const;
		double Vx(double, double, double, bool) const;
		double Vy(double, double, double, bool) const;
		double Vxx(double, double, double, bool) const;
		double Vyy(double, double, double, bool) const;

		void solveByNewtonMethod_BE(double u_old, double v_old, double dt_2, 
																double &u_new, double &v_new) const;
		void solveByNewtonMethod_CN(double u_old, double v_old, double dt_2, 
																double &u_new, double &v_new) const;

	public :

		void setGrid(const CartesianGrid &g);


		void setPrintInterval(double);
		void setPlotInterval(double);

		void initialize(double t0, double T);

		double getTimeStep(void) const;
		int getTimeStepNumber(void) const;

		void getSplittedNeumannBC(int i, double, double&, double&, double&, double&) const;


		void computeInteriorError(const MatrixXd &u,
															double &max_err, double &l2_err) const;

		void advanceByStrang(double t, double dt);
		void advanceByADI(double t, double dt);

		void advanceByStrang_N(double t, double dt);
		void advanceByADI_N(double t, double dt);

		void advance(double t, double dt);

		void solve(void);

	private :

#ifdef USE_OPENGL
		int _plot_opt;
#endif

		double _eps0;
		double _eps1;
		double _eps2;
		double _gama; // feed rate
		double _kappa; // death rate

		//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		double _T;
		double _t;
		double _dt;

		double _print_info_interval;
		double _plot_interval;

		int _max_step_num;
		int _step_count;
		int _check_interval;
		int _plot_check_interval;
		
		double (*_Ui)(double, double, double);
		double (*_Ue)(double, double, double);

		double (*_Vi)(double, double, double);
		double (*_Ve)(double, double, double);

		MatrixXd _u;
		MatrixXd _work_u;

		MatrixXd _ux;
		MatrixXd _uy;
		MatrixXd _uxx;
		MatrixXd _uyy;

		MatrixXd _v;
		MatrixXd _work_v;

		MatrixXd _vx;
		MatrixXd _vy;
		MatrixXd _vxx;
		MatrixXd _vyy;

		VectorXd _u_tau_past1;
		VectorXd _u_tau_past2;
		VectorXd _v_tau_past1;
		VectorXd _v_tau_past2;
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline GrayScottEqn::GrayScottEqn(void) : CartesianGrid()
{
	_eps0 = 0.01;
	_eps1 = 0.008;
	_eps2 = 0.004;

	//_gama = 0.037; // feed rate
	//_kappa = 0.06; // death rate

	_gama = 0.024;
	_kappa = 0.06;

	_T = 1.0;
	_t = 0.0;
	_dt = 0.01;

	_print_info_interval = 0.1;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline GrayScottEqn::~GrayScottEqn(void)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline GrayScottEqn::
GrayScottEqn(const CartesianGrid &g) : CartesianGrid(g)
{
	_eps0 = 0.01;
	_eps1 = 0.008;
	_eps2 = 0.004;

	//_gama = 0.037; // feed rate
	//_kappa = 0.06; // death rate

	_gama = 0.024;
	_kappa = 0.06;

	_T = 1.0;
	_t = 0.0;
	_dt = 0.01;

	_print_info_interval = 0.1;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void GrayScottEqn::
setExactSolution(double (*U)(double, double, double), 
								 double (*V)(double, double, double))
{
	_Ui = U;
	_Vi = V;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void GrayScottEqn::setCoefficient(double eps0, double eps1, double eps2,
																				 double gama, double kappa)
{
	_eps0 = eps0;
	_eps1 = eps1;
	_eps2 = eps2;
	_gama = gama;
	_kappa = kappa;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void GrayScottEqn::setPrintInterval(double d)
{
	_print_info_interval = d;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void GrayScottEqn::setPlotInterval(double d)
{
	_plot_interval = d;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double GrayScottEqn::getTimeStep(void) const
{
	return _dt;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline int GrayScottEqn::getTimeStepNumber(void) const
{
	return _max_step_num;
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

