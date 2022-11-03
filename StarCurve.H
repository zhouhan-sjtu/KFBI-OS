/*=============================================================================
*   
*   Filename : StarCurve.H
*   Creator : Han Zhou
*   Date : 11/18/21
*   Description : 
*
=============================================================================*/
   
 
#ifndef _STARCURVE_H
#define _STARCURVE_H

#include <iostream>
#include "ParametricCurve.H"

class StarCurve : public ParametricCurve
{

	public :

		StarCurve(void);
		StarCurve(double, double, double, double, double);
		virtual ~StarCurve(void);

		StarCurve(const StarCurve &);
		const StarCurve &operator = (const StarCurve &g);

	public:

		void setShape(double, double, double, double, double);
		
		virtual void getX(double t, double &x, double &dx, double &ddx) const;

		virtual void getY(double t, double &y, double &dy, double &ddy) const;

		virtual void getPoint(double t, double &x, double &y) const;

		void printCurveInfo(void) const;


	private :
		
		double _ra;
		double _rb;
		double _eps;
		double _alpha;
		double _m;

};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline StarCurve::~StarCurve(void)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline StarCurve::StarCurve(void) : ParametricCurve()
{
	_ra = 1.0;
	_rb = 1.0;

	_eps = 0.0;
	_alpha = 0.0;
	_m = 1.0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline StarCurve::StarCurve(double ra, double rb, double e, 
														double t, double k) : ParametricCurve()
{
	_ra = ra;
	_rb = rb;
	_eps = e;
	_alpha = t;
	_m = k;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline StarCurve::StarCurve(const StarCurve &c) : ParametricCurve(c)
{
	_ra = c._ra;
	_rb = c._rb;
	_eps = c._eps;
	_alpha = c._alpha;
	_m = c._m;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const StarCurve &StarCurve::operator = (const StarCurve &c)
{
	if (this != &c) {
		_ra = c._ra;
		_rb = c._rb;
		_eps = c._eps;
		_alpha = c._alpha;
		_m = c._m;
	}
	return (*this);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void StarCurve::setShape(double ra, double rb, double e, 
																double t, double k)
{
	_ra = ra;
	_rb = rb;
	_eps = e;
	_alpha = t;
	_m = k;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void StarCurve::getX(double t, double &x, double &dx, double &ddx) const
{
	double s = t - _alpha;

	double cs1 = cos(t);
	double sn1 = sin(t);

	double cs2 = cos(_m * s);
	double sn2 = sin(_m * s);

	double r = _ra + _eps * cs2;

	double meps = _m * _eps;

	x = r * cs1;
	dx = - r * sn1 - meps * cs1 * sn2;
	ddx = 2.0 * meps * sn1 * sn2 - _m * meps * cs1 * cs2 - cs1 * r;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void StarCurve::getY(double t, double &y, double &dy, double &ddy) const
{
	double s = t - _alpha;

	double cs1 = cos(t);
	double sn1 = sin(t);

	double cs2 = cos(_m * s);
	double sn2 = sin(_m * s);

	double r = _rb + _eps * cs2;

	double meps = _m * _eps;

	y = r * sn1;
	dy = r * cs1 - meps * sn1 * sn2;
	ddy = - 2.0 * meps * cs1 * sn2 - _m * meps * cs2 * sn1 - sn1 * r;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void StarCurve::getPoint(double t, double &x, double &y) const
{
	double s = t - _alpha;

	double cs2 = cos(_m * s);

	x = (_ra + _eps * cs2) * cos(t);
	y = (_rb + _eps * cs2) * sin(t);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void StarCurve::printCurveInfo(void) const
{
	std::cout << "star shaped curve." << std::endl;
	std::cout << "ra = " << _ra << ", rb = " << _rb 
						<< ", c = " << _eps << ", alpha = " << _alpha
						<< ", m = " << _m << std::endl;
}

#endif
