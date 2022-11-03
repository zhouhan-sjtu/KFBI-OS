/*=============================================================================
*   
*   Filename : FisherEqn.H
*   Creator : Han Zhou
*   Date : 11/20/21
*   Description : 
*
=============================================================================*/
 
 
#ifndef _FISHEREQN_H
#define _FISHEREQN_H


#include "CartesianGrid.H"
#include "Variables.h"

class FisherEqn : public CartesianGrid
{
	public :
		FisherEqn(void);
		virtual ~FisherEqn(void);

		FisherEqn(const CartesianGrid &g);

	public :

#ifdef USE_OPENGL
		void setPlotOption(int);
		void plotSolution(void) const;
#endif


	public :

		void setGrid(const CartesianGrid &g);

		void setExactSolution(double (*Ui)(double, double, double), 
													double (*Ue)(double, double, double));

		void setCoefficient(double);
		void setPrintInterval(double);

		double F(double u) const;

		double U(double, double, double, bool) const;
		double Ux(double, double, double, bool) const;
		double Uy(double, double, double, bool) const;
		double Uxx(double, double, double, bool) const;
		double Uyy(double, double, double, bool) const;

		void initialize(double t0, double T);

		double getTimeStep(void) const;

		int getTimeStepNumber(void) const;

		double solveByNewtonMethod_BE(double u_old, double dt_2) const;
		double solveByNewtonMethod_CN(double u_old, double dt_2) const;

		void computeInteriorError(const MatrixXd &u,
															double &max_err, double &l2_err) const;

		void advanceByStrang(double t, double dt);
		void advanceByADI(double t, double dt);

		void advance(double t, double dt);

		void solve(void);

	private :

#ifdef USE_OPENGL
		int _plot_opt;
#endif

		double _epsilon;

		double _T;
		double _t;
		double _dt;

		double _print_info_interval;

		int _max_step_num;
		int _step_count;
		int _check_interval;
		
		double (*_Ui)(double, double, double);
		double (*_Ue)(double, double, double);

		MatrixXd _u;

		MatrixXd _ux;
		MatrixXd _uy;
		MatrixXd _uxx;
		MatrixXd _uyy;
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline FisherEqn::FisherEqn(void) : CartesianGrid()
{
	_T = 1.0;
	_t = 0.0;
	_dt = 0.01;

	_epsilon = 1.0;

	_print_info_interval = 0.1;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline FisherEqn::~FisherEqn(void)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline FisherEqn::
FisherEqn(const CartesianGrid &g) : CartesianGrid(g)
{
	_T = 1.0;
	_t = 0.0;
	_dt = 0.01;

	_epsilon = 1.0;

	_print_info_interval = 0.1;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void FisherEqn::
setExactSolution(double (*Ui)(double, double, double), 
								 double (*Ue)(double, double, double))
{
	_Ui = Ui;
	_Ue = Ue;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void FisherEqn::setCoefficient(double k)
{
	_epsilon = k;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void FisherEqn::setPrintInterval(double d)
{
	_print_info_interval = d;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double FisherEqn::getTimeStep(void) const
{
	return _dt;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline int FisherEqn::getTimeStepNumber(void) const
{
	return _max_step_num;
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

