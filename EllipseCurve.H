/*=============================================================================
*   
*   Filename : EllipseCurve.H
*   Creator : Han Zhou
*   Date : 11/19/21
*   Description : 
*
=============================================================================*/
   
#ifndef _ELLIPSECURVE_H
#define _ELLIPSECURVE_H

#include <cmath>

#include "ParametricCurve.H"

class EllipseCurve : public ParametricCurve
{

	public :

		EllipseCurve(void);
		EllipseCurve(double, double, double, double, double);
		virtual ~EllipseCurve(void);

		EllipseCurve(const EllipseCurve &);
		const EllipseCurve &operator = (const EllipseCurve &g);

	public:

		void setShape(double, double, double, double, double);
		
		virtual void getX(double t, double &x, double &dx, double &ddx) const;

		virtual void getY(double t, double &y, double &dy, double &ddy) const;

		virtual void getPoint(double t, double &x, double &y) const;

		void printCurveInfo(void) const;

	private :
		
		double _cx;
		double _cy;
		double _ra;
		double _rb;
		double _cs;
		double _sn;
		double _alpha;

};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline EllipseCurve::~EllipseCurve(void)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline EllipseCurve::EllipseCurve(void) : ParametricCurve()
{
	_cx = 0.0;
	_cy = 0.0;
	_ra = 1.0;
	_rb = 1.0;
	_alpha = 0.0;
	_cs = 1.0;
	_sn = 0.0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline EllipseCurve::EllipseCurve(double cx, double cy, double ra, double rb,
																	double t) : ParametricCurve()
{
	_cx = cx;
	_cy = cy;
	_ra = ra;
	_rb = rb;
	_alpha = t;
	_cs = cos(_alpha);
	_sn = sin(_alpha);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline EllipseCurve::EllipseCurve(const EllipseCurve &c) : ParametricCurve(c)
{
	_cx = c._cx;
	_cy = c._cy;
	_ra = c._ra;
	_rb = c._rb;
	_alpha = c._alpha;
	_cs = c._cs;
	_sn = c._sn;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const EllipseCurve &EllipseCurve::operator = (const EllipseCurve &c)
{
	if (this != &c) {
		_cx = c._cx;
		_cy = c._cy;
		_ra = c._ra;
		_rb = c._rb;
		_alpha = c._alpha;
		_cs = c._cs;
		_sn = c._sn;
	}
	return (*this);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void EllipseCurve::setShape(double cx, double cy, 
																	 double ra, double rb, double t)
{
	_cx = cx;
	_cy = cy;
	_ra = ra;
	_rb = rb;
	_alpha = t;
	_cs = cos(_alpha);
	_sn = sin(_alpha);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void EllipseCurve::
getX(double t, double &x, double &dx, double &ddx) const
{
  double y, dy, ddy; 

  double cs = cos(t); 
  double sn = sin(t); 

  double p = _ra * cs; 
  double q = _rb * sn; 

  double dp = - _ra * sn;
  double dq = _rb * cs; 

  double ddp = - _ra * cs; 
  double ddq = - _rb * sn;

  x = _cx + _cs * p - _sn * q;
  y = _cy + _sn * p + _cs * q; 

  dx = _cs * dp - _sn * dq;
  dy = _sn * dp + _cs * dq; 

  ddx = _cs * ddp - _sn * ddq; 
  ddy = _sn * ddp + _cs * ddq;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void EllipseCurve::
getY(double t, double &y, double &dy, double &ddy) const
{
  double x, dx, ddx;

  double cs = cos(t); 
  double sn = sin(t);

  double p = _ra * cs; 
  double q = _rb * sn;

  double dp = - _ra * sn;
  double dq = _rb * cs; 

  double ddp = - _ra * cs;
  double ddq = - _rb * sn; 

  x = _cx + _cs * p - _sn * q;
  y = _cy + _sn * p + _cs * q;

  dx = _cs * dp - _sn * dq;
  dy = _sn * dp + _cs * dq;

  ddx = _cs * ddp - _sn * ddq; 
  ddy = _sn * ddp + _cs * ddq;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void EllipseCurve::getPoint(double t, double &x, double &y) const
{
  double p = _ra * cos(t);
  double q = _rb * sin(t);
  x = _cx + _cs * p - _sn * q;
  y = _cy + _sn * p + _cs * q; 
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void EllipseCurve::printCurveInfo(void) const
{
	std::cout << "ellipse shaped curve." << std::endl;
	std::cout << "ra = " << _ra << ", rb = " << _rb 
						<< ", cx = " << _cx << ", cy = " << _cy
						<< ", alpha = " << _alpha << std::endl;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#endif
