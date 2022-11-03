/*=============================================================================
*   
*   Filename : ParametricCurve.H
*   Creator : Han Zhou
*   Date : 11/18/21
*   Description : 
*
=============================================================================*/
 
#ifndef _PARAMETRICCURVE_H
#define _PARAMETRICCURVE_H

#include "Variables.h"

class ParametricCurve
{

	public:
		ParametricCurve(void);
		virtual ~ParametricCurve(void);

		ParametricCurve(const ParametricCurve &);
		const ParametricCurve &operator = (const ParametricCurve &g);

	public:

		void plotCurve(void) const;

	public:
		
		virtual void getX(double t, double &x, double &dx, double &ddx) const = 0;

		virtual void getY(double t, double &y, double &dy, double &ddy) const = 0;

		virtual void getPoint(double t, double &x, double &y) const = 0;

		virtual void getPoint(double t, double &x, double &y, double &tx, 
													double &ty, double &nx, double &ny) const;

		virtual void getPoint2(double t, double &x, double &y, double &dx, 
														double &dy, double &ddx, double &ddy) const;

		virtual int findZero(void (*F)(double, double &, double &, double &), 
												 double p, double &t0) const;

		virtual int findZero2(void (*F)(double, double &, double &, double &), 
												 double p, double ta, double tb, double &t0) const;

		virtual int findXIntersection(double &t, double &x, double y) const;
		virtual int findYIntersection(double &t, double x, double &y) const;

		virtual int findXIntersection2(double &t, double &x, double y, 
																	 double ta, double tb) const;
		virtual int findYIntersection2(double &t, double x, double &y, 
																	 double ta, double tb) const;

		virtual int findClosestPointOnCurve(double p, double q, double &t0, 
																				double &x0, double &y0) const;

		virtual int findClosestPointOnCurve2(double p, double q, double &t0, 
																				double &x0, double &y0, 
																				double ta, double tb) const;

		virtual double getCurveLength(void) const;
		virtual double getCurvature(double) const;

		virtual double computeArcLength(void) const;

		virtual void setupCurveOnCartesianGrid(const VectorXd &x, const VectorXd &y, 
															int I, int J,	MatrixXb &interior, 
															MatrixXb &irr_edge_x, MatrixXb &irr_edge_y,
															MatrixXd &i_theta_x, MatrixXd &i_theta_y) const;

};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ParametricCurve::ParametricCurve(void)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ParametricCurve::~ParametricCurve(void)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline ParametricCurve::ParametricCurve(const ParametricCurve &s)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const ParametricCurve &ParametricCurve::
							operator = (const ParametricCurve &s)
{
  if (this != &s) {
  }
  return (*this);
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif
