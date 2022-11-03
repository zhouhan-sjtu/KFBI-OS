/*=============================================================================
*   
*   Filename : CartesianGrid.C
*   Creator : Han Zhou
*   Date : 11/19/21
*   Description : 
*
=============================================================================*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
 
#include "CartesianGrid.H"
#include "MathTools.h"
#include "Const.h"

int CartesianGrid::_total_stc_n = 49;
int CartesianGrid::_stc_bak[49][2] = {0};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::initializeStaticData(void)
{
	for(int i = -3, k = 0; i <= 3; i++){
		for(int j = -3; j <= 3; j++, k++){
			_stc_bak[k][0] = i;
			_stc_bak[k][1] = j;
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::initDataMembers(void) 
{
	_I1 = _I + 1;
	_J1 = _J + 1;

	_dx = (_high[0] - _low[0]) / _I;
	_dy = (_high[1] - _low[1]) / _J;
	_h = 0.5 * (_dx + _dy);

	_x.setAxis(_low[0], _high[0], _I);
	_y.setAxis(_low[1], _high[1], _J);

	initializeStaticData();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::initCurveEntities(void)
{
	_interior.reallocate(_I1, _J1);

	_x_irregular_edge.reallocate(_I, _J1);
	_x_intersect_theta.reallocate(_I, _J1);
	_x_intersect_idx.reallocate(_I, _J1);

	_y_irregular_edge.reallocate(_I1, _J);
	_y_intersect_theta.reallocate(_I1, _J);
	_y_intersect_idx.reallocate(_I1, _J);

	_curve->setupCurveOnCartesianGrid(_x, _y, _I, _J, _interior,
																		_x_irregular_edge, _y_irregular_edge,
																		_x_intersect_theta, _y_intersect_theta);

	int M0 = 0;
	int M1 = 0;

	for(int i = 0; i < _I; i++){
		for(int j = 0; j < _J1; j++){
			if (_x_irregular_edge[i][j]) {
				M0++;
			}
		}
	}

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J; j++){
			if (_y_irregular_edge[i][j]) {
				M1++;
			}
		}
	}

	_intersect_node_n[0] = M0;
	_intersect_node_n[1] = M1;
	_total_intersect_node_n = M0 + M1;

	_intersect_node_crd.reallocate(_total_intersect_node_n);
	_intersect_node_nml.reallocate(_total_intersect_node_n);
	_intersect_node_tan.reallocate(_total_intersect_node_n);

	double t, x, y, tx, ty, nx, ny;

	int count = 0;
	for(int i = 0; i < _I; i++){
		for(int j = 0; j < _J1; j++){
			if (_x_irregular_edge[i][j]) {

				t = _x_intersect_theta[i][j];
				_curve->getPoint(t, x, y, tx, ty, nx, ny);

				_intersect_node_crd[count][0] = x;
				_intersect_node_crd[count][1] = y;

				_intersect_node_nml[count][0] = nx;
				_intersect_node_nml[count][1] = ny;

				_intersect_node_tan[count][0] = tx;
				_intersect_node_tan[count][1] = ty;

				_x_intersect_idx[i][j] = count;

				count++;
			}
		}
	}

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J; j++){
			if (_y_irregular_edge[i][j]) {

				t = _y_intersect_theta[i][j];
				_curve->getPoint(t, x, y, tx, ty, nx, ny);

				_intersect_node_crd[count][0] = x;
				_intersect_node_crd[count][1] = y;

				_intersect_node_nml[count][0] = nx;
				_intersect_node_nml[count][1] = ny;

				_intersect_node_tan[count][0] = tx;
				_intersect_node_tan[count][1] = ty;

				_y_intersect_idx[i][j] = count;

				count++;
			}
		}
	}

	double length = _curve->getCurveLength();

	_ctrl_node_n = static_cast<int>(length / _h + 0.5);

	_ctrl_node_crd.reallocate(_ctrl_node_n);
	_ctrl_node_nml.reallocate(_ctrl_node_n);
	_ctrl_node_tan.reallocate(_ctrl_node_n);

	double ds = M_2PI / _ctrl_node_n;
	double s = 0.0;
	
	for(int k = 0; k < _ctrl_node_n; k++){

		_curve->getPoint(s, x, y, tx, ty, nx, ny);

		_ctrl_node_crd[k][0] = x;
		_ctrl_node_crd[k][1] = y;

		_ctrl_node_nml[k][0] = nx;
		_ctrl_node_nml[k][1] = ny;

		_ctrl_node_tan[k][0] = tx;
		_ctrl_node_tan[k][1] = ty;

		s += ds;
	}


	_x_line_side.reallocate(_J1);
	_y_line_side.reallocate(_I1);

	for(int j = 0; j < _J1; j++){
		_x_line_side[j].reallocate(_I1);

		for(int i = 0; i < _I1; i++){
			_x_line_side[j][i] = _interior[i][j];
		}
	}

	for(int i = 0; i < _I1; i++){
		_y_line_side[i].reallocate(_J1);

		for(int j = 0; j < _J1; j++){
			_y_line_side[i][j] = _interior[i][j];
		}
	}


	///////////////////////////////////////////////////////////////////// 

	_x_irreglar_line.reallocate(_J1);
	_y_irreglar_line.reallocate(_I1);

	_x_line_intersect_num.reallocate(_J1);
	_y_line_intersect_num.reallocate(_I1);

	_x_line_intersec_crd.reallocate(_J1);
	_y_line_intersec_crd.reallocate(_I1);

	for(int j = 0; j < _J1; j++){

		int count = 0;

		for(int i = 0; i < _I; i++){
			if (_x_irregular_edge[i][j]) {
				count++;
			}
		}

		_x_line_intersect_num[j] = count;

		if (count > 0) {

			_x_irreglar_line[j] = true;

			_x_line_intersec_crd[j].reallocate(count);

			count = 0;

			for(int i = 0; i < _I; i++){
				if (_x_irregular_edge[i][j]) {

					int k = _x_intersect_idx[i][j];

					assert(k >= 0);

					_x_line_intersec_crd[j][count] = _intersect_node_crd[k][0];
					count++;
				}
			}

		} else {
			_x_irreglar_line[j] = false;
		}

	}

	for(int i = 0; i < _I1; i++){

		int count = 0;

		for(int j = 0; j < _J; j++){
			if (_y_irregular_edge[i][j]) {
				count++;
			}
		}

		_y_line_intersect_num[i] = count;

		if (count > 0) {

			_y_irreglar_line[i] = true;

			_y_line_intersec_crd[i].reallocate(count);

			count = 0;

			for(int j = 0; j < _J; j++){

				if (_y_irregular_edge[i][j]) {

					int k = _y_intersect_idx[i][j];

					assert(k >= 0);

					_y_line_intersec_crd[i][count] = _intersect_node_crd[k][1];
					count++;
				}
			}

		} else {
			_y_irreglar_line[i] = false;
		}

	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::extractTangentU(const MatrixXd &u, VectorXd &u_tau) const
{
	extern bool makeLeastSquareInterpolation(const double p[2], 
	                                 	const VectorX2d &coord,
	                                 	const VectorXd &b, int N, double u[6]);

	const int used_n = 49;

#pragma omp parallel for 
	for(int l = 0; l < _total_intersect_node_n; l++){

		double g[6];
		VectorXd b(_total_stc_n);
		VectorX2d crd(_total_stc_n);

		double x = _intersect_node_crd[l][0];
		double y = _intersect_node_crd[l][1];

		double tx = _intersect_node_tan[l][0];
		double ty = _intersect_node_tan[l][1];

		double norm = sqrt(tx * tx + ty * ty);
		tx /= norm;
		ty /= norm;

		int i = static_cast<int>((x - _low[0]) / _dx + 0.5);
		int j = static_cast<int>((y - _low[1]) / _dy + 0.5);

		int count = 0;

		for(int m = 0; m < _total_stc_n; m++){

			int r = i + _stc_bak[m][0];
			int s = j + _stc_bak[m][1];

			if (_interior[r][s] && count < used_n) {

				crd[count][0] = _stc_bak[m][0];
				crd[count][1] = _stc_bak[m][1];
				b[count++] = u[r][s];

			}
		}

		double p[2] = {(x - _x[i]) / _dx, (y - _y[j]) / _dy};
		bool status = makeLeastSquareInterpolation(p, crd, b, count, g);
		if (!status) {
			printBugInfo("failed to makeLeastSquareInterpolation .");
			exit(1);
		}

		u_tau[l] = (tx * g[1] + ty * g[2]) / _h;

	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::extractBoundryData(const MatrixXd &u, 
																		VectorXd &nd_ux, VectorXd &nd_uy) const
{
	extern bool makeLeastSquareInterpolation(const double p[2], 
	                                 	const VectorX2d &coord,
	                                 	const VectorXd &b, int N, double u[6]);

	const int used_n = 49;

#pragma omp parallel for 
	for(int l = 0; l < _total_intersect_node_n; l++){

		double g[6];
		VectorXd b(_total_stc_n);
		VectorX2d crd(_total_stc_n);

		double x = _intersect_node_crd[l][0];
		double y = _intersect_node_crd[l][1];

		double tx = _intersect_node_tan[l][0];
		double ty = _intersect_node_tan[l][1];

		double norm = sqrt(tx * tx + ty * ty);
		tx /= norm;
		ty /= norm;

		int i = static_cast<int>((x - _low[0]) / _dx + 0.5);
		int j = static_cast<int>((y - _low[1]) / _dy + 0.5);

		int count = 0;

		for(int m = 0; m < _total_stc_n; m++){

			int r = i + _stc_bak[m][0];
			int s = j + _stc_bak[m][1];

			if (_interior[r][s] && count < used_n) {

				crd[count][0] = _stc_bak[m][0];
				crd[count][1] = _stc_bak[m][1];
				b[count++] = u[r][s];

			}
		}

		double p[2] = {(x - _x[i]) / _dx, (y - _y[j]) / _dy};
		bool status = makeLeastSquareInterpolation(p, crd, b, count, g);
		if (!status) {
			printBugInfo("failed to makeLeastSquareInterpolation .");
			exit(1);
		}

		nd_ux[l] = g[1] / _h;
		nd_uy[l] = g[2] / _h;

	}

}

#ifdef USE_OPENGL

#include "MathGLUT2d.h"
#include "NewMathGL.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotGridLines(void) const
{
	mglSetColor(0.5, 0.5, 0.5);
	mglSetLineWidth(1.0);

	for(int i = 0; i < _I1; i++){
		mglPlotLine(_x[i], _y[0], _x[i], _y[_J]);
	}

	for(int j = 0; j < _J1; j++){
		mglPlotLine(_x[0], _y[j], _x[_J], _y[j]);
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotInteriorNodes(void) const
{
	mglSetColor(0, 0, 0);
	mglSetPointSize(5.0);

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {
				mglPlotPoint(_x[i], _y[j]);
			}
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotIntersectPoints(void) const
{
	mglSetPointSize(5.0);

	int K1 = _intersect_node_n[0];
	int K2 = _total_intersect_node_n;

	mglSetColor(1, 0, 0);
	for(int k = 0; k < K1; k++){
		mglPlotPoint(_intersect_node_crd[k][0], _intersect_node_crd[k][1]);
	}

	mglSetColor(0, 0, 1);
	for(int k = K1; k < K2; k++){
		mglPlotPoint(_intersect_node_crd[k][0], _intersect_node_crd[k][1]);
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotControlPoints(void) const
{
	mglSetPointSize(5.0);
	mglSetColor(0, 0, 0);

	for(int k = 0; k < _ctrl_node_n; k++){
		mglPlotPoint(_ctrl_node_crd[k][0], _ctrl_node_crd[k][1]);
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotNodeDataByColormap(const MatrixXd &u) const
{
	double max_u = u[0][0];
	double min_u = u[0][0];

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			max_u = max(max_u, u[i][j]);
			min_u = min(min_u, u[i][j]);
		}
	}

	mglSetColormap(min_u, max_u);

	double coord[3][2], v[3];

	for(int i = 0; i < _I; i++){
		for(int j = 0; j < _J; j++){

			coord[0][0] = _x[i];
			coord[0][1] = _y[j];
			v[0] = u[i][j];

			coord[1][0] = _x[i + 1];
			coord[1][1] = _y[j];
			v[1] = u[i + 1][j];

			coord[2][0] = _x[i];
			coord[2][1] = _y[j + 1];
			v[2] = u[i][j + 1];

			mglPlotSolidTriangle(coord, v);

			coord[0][0] = _x[i + 1];
			coord[0][1] = _y[j + 1];
			v[0] = u[i + 1][j + 1];

			mglPlotSolidTriangle(coord, v);
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotNodeDataByIsolines(const MatrixXd &u) const
{
	double max_u = u[0][0];
	double min_u = u[0][0];

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			max_u = max(max_u, u[i][j]);
			min_u = min(min_u, u[i][j]);
		}
	}

	double coord[3][2], v[3];

	for(int i = 0; i < _I; i++){
		for(int j = 0; j < _J; j++){

			coord[0][0] = _x[i];
			coord[0][1] = _y[j];
			v[0] = u[i][j];

			coord[1][0] = _x[i + 1];
			coord[1][1] = _y[j];
			v[1] = u[i + 1][j];

			coord[2][0] = _x[i];
			coord[2][1] = _y[j + 1];
			v[2] = u[i][j + 1];

			mglPlotIsolines(coord, v, min_u, max_u, 50);

			coord[0][0] = _x[i + 1];
			coord[0][1] = _y[j + 1];
			v[0] = u[i + 1][j + 1];

			mglPlotIsolines(coord, v, min_u, max_u, 50);
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotInteriorNodeDataByColormap(const MatrixXd &u) const
{
	double max_u = -DBL_MAX;
	double min_u = DBL_MAX;

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {
				max_u = max(max_u, u[i][j]);
				min_u = min(min_u, u[i][j]);
			}
		}
	}

	mglSetColormap(min_u, max_u);

	double coord[3][2], v[3];

	for(int i = 0; i < _I; i++){
		for(int j = 0; j < _J; j++){

			coord[0][0] = _x[i];
			coord[0][1] = _y[j];
			v[0] = u[i][j];

			coord[1][0] = _x[i + 1];
			coord[1][1] = _y[j];
			v[1] = u[i + 1][j];

			coord[2][0] = _x[i];
			coord[2][1] = _y[j + 1];
			v[2] = u[i][j + 1];

			if (_interior[i][j] && _interior[i+1][j] && _interior[i][j+1]) {
				mglPlotSolidTriangle(coord, v);
			}

			coord[0][0] = _x[i + 1];
			coord[0][1] = _y[j + 1];
			v[0] = u[i + 1][j + 1];

			if (_interior[i+1][j+1] && _interior[i+1][j] && _interior[i][j+1]) {
				mglPlotSolidTriangle(coord, v);
			}
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotInteriorNodeDataByIsolines(const MatrixXd &u) const
{
	double max_u = -DBL_MAX;
	double min_u = DBL_MAX;

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (_interior[i][j]) {
				max_u = max(max_u, u[i][j]);
				min_u = min(min_u, u[i][j]);
			}
		}
	}

	double coord[3][2], v[3];

	for(int i = 0; i < _I; i++){
		for(int j = 0; j < _J; j++){

			coord[0][0] = _x[i];
			coord[0][1] = _y[j];
			v[0] = u[i][j];

			coord[1][0] = _x[i + 1];
			coord[1][1] = _y[j];
			v[1] = u[i + 1][j];

			coord[2][0] = _x[i];
			coord[2][1] = _y[j + 1];
			v[2] = u[i][j + 1];

			if (_interior[i][j] && _interior[i+1][j] && _interior[i][j+1]) {
				mglPlotIsolines(coord, v, min_u, max_u, 50);
			}

			coord[0][0] = _x[i + 1];
			coord[0][1] = _y[j + 1];
			v[0] = u[i + 1][j + 1];

			if (_interior[i+1][j+1] && _interior[i+1][j] && _interior[i][j+1]) {
				mglPlotIsolines(coord, v, min_u, max_u, 50);
			}
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotExteriorNodeDataByColormap(const MatrixXd &u) const
{
	double max_u = -DBL_MAX;
	double min_u = DBL_MAX;

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (!_interior[i][j]) {
				max_u = max(max_u, u[i][j]);
				min_u = min(min_u, u[i][j]);
			}
		}
	}
	//std::cout << max_u << ", " << min_u << std::endl;

	mglSetColormap(min_u, max_u);

	double coord[3][2], v[3];

	for(int i = 0; i < _I; i++){
		for(int j = 0; j < _J; j++){

			coord[0][0] = _x[i];
			coord[0][1] = _y[j];
			v[0] = u[i][j];

			coord[1][0] = _x[i + 1];
			coord[1][1] = _y[j];
			v[1] = u[i + 1][j];

			coord[2][0] = _x[i];
			coord[2][1] = _y[j + 1];
			v[2] = u[i][j + 1];

			if (!(_interior[i][j] || _interior[i+1][j] || _interior[i][j+1])) {
				mglPlotSolidTriangle(coord, v);
			}

			coord[0][0] = _x[i + 1];
			coord[0][1] = _y[j + 1];
			v[0] = u[i + 1][j + 1];

			if (!(_interior[i+1][j+1] || _interior[i+1][j] || _interior[i][j+1])) {
				mglPlotSolidTriangle(coord, v);
			}
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotExteriorNodeDataByIsolines(const MatrixXd &u) const
{
	double max_u = -DBL_MAX;
	double min_u = DBL_MAX;

	for(int i = 0; i < _I1; i++){
		for(int j = 0; j < _J1; j++){
			if (!_interior[i][j]) {
				max_u = max(max_u, u[i][j]);
				min_u = min(min_u, u[i][j]);
			}
		}
	}
	//std::cout << max_u << ", " << min_u << std::endl;

	double coord[3][2], v[3];

	for(int i = 0; i < _I; i++){
		for(int j = 0; j < _J; j++){

			coord[0][0] = _x[i];
			coord[0][1] = _y[j];
			v[0] = u[i][j];

			coord[1][0] = _x[i + 1];
			coord[1][1] = _y[j];
			v[1] = u[i + 1][j];

			coord[2][0] = _x[i];
			coord[2][1] = _y[j + 1];
			v[2] = u[i][j + 1];

			if (!(_interior[i][j] || _interior[i+1][j] || _interior[i][j+1])) {
				mglPlotIsolines(coord, v, min_u, max_u, 50);
			}

			coord[0][0] = _x[i + 1];
			coord[0][1] = _y[j + 1];
			v[0] = u[i + 1][j + 1];

			if (!(_interior[i+1][j+1] || _interior[i+1][j] || _interior[i][j+1])) {
				mglPlotIsolines(coord, v, min_u, max_u, 50);
			}
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CartesianGrid::plotNormalsAndTangents(void) const
{
	for(int k = 0; k < _total_intersect_node_n; k++){

		double x = _intersect_node_crd[k][0];
		double y = _intersect_node_crd[k][1];

		double nx = _intersect_node_nml[k][0];
		double ny = _intersect_node_nml[k][1];

		double tx = _intersect_node_tan[k][0];
		double ty = _intersect_node_tan[k][1];

		double norm = sqrt(tx * tx + ty * ty);
		tx /= norm;
		ty /= norm;

		mglSetColor(1.0, 0.0, 0.0);
		mglPlotLine(x, y, x + 0.1 * nx, y + 0.1 * ny);

		mglSetColor(0.0, 0.0, 1.0);
		mglPlotLine(x, y, x + 0.1 * tx, y + 0.1 * ty);


	}
}

#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

