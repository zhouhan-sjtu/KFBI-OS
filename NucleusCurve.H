/*=============================================================================
*   
*   Filename : NucleusCurve.H
*   Creator : Han Zhou
*   Date : 11/30/21
*   Description : 
*
=============================================================================*/

 
#ifndef _NUCLEUSCURVE_H
#define _NUCLEUSCURVE_H

#include <iostream>
#include "ParametricCurve.H"

class NucleusCurve : public ParametricCurve
{

	public :

		NucleusCurve(void);
		NucleusCurve(double, double, double, double);
		virtual ~NucleusCurve(void);

		NucleusCurve(const NucleusCurve &);
		const NucleusCurve &operator = (const NucleusCurve &g);

	public:

		void setShape(double, double, double, double);
		
		virtual void getX(double t, double &x, double &dx, double &ddx) const;

		virtual void getY(double t, double &y, double &dy, double &ddy) const;

		virtual void getPoint(double t, double &x, double &y) const;

		void printCurveInfo(void) const;


	private :
		
		double _r0;
		double _m1;
		double _m2;
		double _eps;

};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline NucleusCurve::~NucleusCurve(void)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline NucleusCurve::NucleusCurve(void) : ParametricCurve()
{
	_r0 = 1.0;
	_m1 = 2.0;
	_m2 = 3.0;
	_eps = 0.1;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline NucleusCurve::NucleusCurve(double r0, double eps, 
																	double m1, double m2) : ParametricCurve()
{
	_r0 = r0;
	_m1 = m1;
	_m2 = m2;
	_eps = eps;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline NucleusCurve::NucleusCurve(const NucleusCurve &c) : ParametricCurve(c)
{
	_r0 = c._r0;
	_m1 = c._m1;
	_m2 = c._m2;
	_eps = c._eps;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const NucleusCurve &NucleusCurve::operator = (const NucleusCurve &c)
{
	if (this != &c) {
		_r0 = c._r0;
		_m1 = c._m1;
		_m2 = c._m2;
		_eps = c._eps;
	}
	return (*this);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void NucleusCurve::setShape(double r0, double eps, double m1, double m2)
{
	_r0 = r0;
	_m1 = m1;
	_m2 = m2;
	_eps = eps;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void NucleusCurve::getX(double t, double &x, double &dx, double &ddx) const
{
	double snm1 = sin(_m1 * t);
	double csm2 = cos(_m2 * t);

	double r = _r0 + _eps * (snm1 + csm2);
	double rp = _eps * (_m1 * cos(_m1 * t) - _m2 * sin(_m2 * t));
	double rpp = - _eps * (_m1 * _m1 * snm1 + _m2 * _m2 * csm2);

	double cs = cos(t);
	double sn = sin(t);

	x = r * cs;
	dx = rp * cs - r * sn;
	ddx = (rpp - r) * cs - 2.0 * rp * sn;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void NucleusCurve::getY(double t, double &y, double &dy, double &ddy) const
{
	double snm1 = sin(_m1 * t);
	double csm2 = cos(_m2 * t);

	double r = _r0 + _eps * (snm1 + csm2);
	double rp = _eps * (_m1 * cos(_m1 * t) - _m2 * sin(_m2 * t));
	double rpp = - _eps * (_m1 * _m1 * snm1 + _m2 * _m2 * csm2);

	double cs = cos(t);
	double sn = sin(t);

	y = r * sn;
	dy = rp * sn + r * cs;
	ddy = (rpp - r) * sn + 2.0*rp*cs;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void NucleusCurve::getPoint(double t, double &x, double &y) const
{
	double r = _r0 + _eps * (sin(_m1 * t) + cos(_m2 * t));

	x = r * cos(t);
	y = r * sin(t);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void NucleusCurve::printCurveInfo(void) const
{
	std::cout << "star shaped curve." << std::endl;
	std::cout << "r0 = " << _r0 << ", m1 = " << _m1 
						<< ", m2 = " << _m2 << ", eps = "
						<< _eps << std::endl;
}

#endif
