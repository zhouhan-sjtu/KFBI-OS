/*=============================================================================
*   
*   Filename : CartesianGrid.H
*   Creator : Han Zhou
*   Date : 11/19/21
*   Description : 
*
=============================================================================*/
 
#ifndef _CARTESIANGRID_H
#define _CARTESIANGRID_H

#include "ParametricCurve.H"
#include "Variables.h"

class CartesianGrid
{
	public :

		CartesianGrid(void);
		virtual ~CartesianGrid(void);

		CartesianGrid(const double low[2], const double high[2], int n);
		CartesianGrid(const double low[2], const double high[2], int I, int J);

		CartesianGrid(const CartesianGrid &);
		const CartesianGrid &operator = (const CartesianGrid &g);

	public :

		virtual double getMeshParameter(void) const;

		virtual bool isInterior(int i, int j) const;

		virtual void defineGrid(const double low[2], const double high[2], int n);
		virtual void defineGrid(const double low[2], const double high[2], int I, int J);

		virtual void setCurve(const ParametricCurve *c);
		virtual void setCurve(const ParametricCurve &c);

		virtual void extractTangentU(const MatrixXd &u, VectorXd &u_tau) const;
		virtual void extractBoundryData(const MatrixXd &u, 
																		VectorXd &nd_ux, VectorXd &nd_uy) const;

#ifdef USE_OPENGL
		virtual void plotGridLines(void) const;
		virtual void plotInteriorNodes(void) const;
		virtual void plotControlPoints(void) const;
		virtual void plotIntersectPoints(void) const;
		virtual void plotNormalsAndTangents(void) const;
		virtual void plotNodeDataByColormap(const MatrixXd&) const;
		virtual void plotNodeDataByIsolines(const MatrixXd&) const;
		virtual void plotInteriorNodeDataByColormap(const MatrixXd&) const;
		virtual void plotInteriorNodeDataByIsolines(const MatrixXd&) const;
		virtual void plotExteriorNodeDataByColormap(const MatrixXd&) const;
		virtual void plotExteriorNodeDataByIsolines(const MatrixXd&) const;
#endif

	public :
		
		static void initializeStaticData(void);

	private :

		void initDataMembers(void);
		void initCurveEntities(void);

	protected :

		const ParametricCurve *_curve;

		double _low[2];
		double _high[2];
		
		int _I;
		int _J;
		int _I1;
		int _J1;

		double _dx;
		double _dy;
		double _h;

		VectorXd _x;
		VectorXd _y;

		MatrixXb _interior;

		//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		VectorXb _x_irreglar_line;
		VectorXb _y_irreglar_line;

		VectorXi _x_line_intersect_num;
		VectorXi _y_line_intersect_num;

		Vector< VectorXd > _x_line_intersec_crd;
		Vector< VectorXd > _y_line_intersec_crd;

		Vector< VectorXb > _x_line_side;
		Vector< VectorXb > _y_line_side;

		//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		MatrixXb _x_irregular_edge;
		MatrixXb _y_irregular_edge;

		MatrixXd _x_intersect_theta;
		MatrixXd _y_intersect_theta;

		MatrixXi _x_intersect_idx;
		MatrixXi _y_intersect_idx;

		int _total_intersect_node_n;
		int _intersect_node_n[2];

		VectorX2d _intersect_node_crd;
		VectorX2d _intersect_node_nml;
		VectorX2d _intersect_node_tan;

		//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		int _ctrl_node_n;

		VectorX2d _ctrl_node_crd;
		VectorX2d _ctrl_node_nml;
		VectorX2d _ctrl_node_tan;

		static int _total_stc_n;
		static int _stc_bak[49][2];

};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline CartesianGrid::CartesianGrid(void)
{
	for(int i = 0; i < 2; i++){
		_low[i] = - 1.0;
		_high[i] = 1.0;
	}
	_I = _J = 32;

	initDataMembers();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline CartesianGrid::CartesianGrid(const double low[2], const double high[2],
																		int n)
{
	for(int i = 0; i < 2; i++){
		_low[i] = low[i];
		_high[i] = high[i];
	}
	_I = _J = n;

	initDataMembers();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline CartesianGrid::CartesianGrid(const double low[2], const double high[2],
																		int I, int J)
{
	for(int i = 0; i < 2; i++){
		_low[i] = low[i];
		_high[i] = high[i];
	}
	_I = I;
	_J = J;

	initDataMembers();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void CartesianGrid::defineGrid(const double low[2], const double high[2],
																	int n)
{
	for(int i = 0; i < 2; i++){
		_low[i] = low[i];
		_high[i] = high[i];
	}
	_I = _J = n;

	initDataMembers();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void CartesianGrid::defineGrid(const double low[2], const double high[2],
																  int I, int J)
{
	for(int i = 0; i < 2; i++){
		_low[i] = low[i];
		_high[i] = high[i];
	}
	_I = I;	_J = J;

	initDataMembers();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline CartesianGrid::~CartesianGrid(void)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void CartesianGrid::setCurve(const ParametricCurve *c)
{
	_curve = c;
	initCurveEntities();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void CartesianGrid::setCurve(const ParametricCurve &c)
{
	_curve = &c;
	initCurveEntities();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline bool CartesianGrid::isInterior(int i, int j) const
{
	return _interior[i][j];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double CartesianGrid::getMeshParameter(void) const
{
	return _h;
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

