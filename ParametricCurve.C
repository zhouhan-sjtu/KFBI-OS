/*=============================================================================
*   
*   Filename : ParametricCurve.C
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
#include <vector>
#include <queue>
 
#include "ParametricCurve.H"
#include "MathTools.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ParametricCurve::getPoint(double t, double &x, double &y, 
															 double &tx, double &ty, 
										 					 double &nx, double &ny) const 
{
	// normalized

	double dx, dy, ddx, ddy;

	getX(t, x, dx, ddx);
	getY(t, y, dy, ddy);

	double norm = sqrt(dx * dx + dy * dy);

#ifdef DEBUG
	if (norm < 1.0E-10) {
		std::cout << "zero tangent vector." << std::endl;
		exit(1);
	}
#endif

	double r_norm = 1.0 / norm;

	//tx = dx * r_norm;
	//ty = dy * r_norm;
	//nx = ty;
	//ny = - tx;

	tx = dx;
	ty = dy;
	nx = dy * r_norm;
	ny = - dx * r_norm;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ParametricCurve::getPoint2(double t, double &x, double &y, 
																double &dx, double &dy, 
																double &ddx, double &ddy) const
{
	getX(t, x, dx, ddx);
	getY(t, y, dy, ddy);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ParametricCurve::getCurvature(double t) const
{
	double x, y, dx, dy, ddx, ddy;
	getPoint2(t, x, y, dx, dy, ddx, ddy);

	double d2 = dx * dx + dy * dy;
	double d = sqrt(d2);
	double numer = dx * ddy - dy * ddx;
	double kappa = numer / (d * d2);

	return kappa;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ParametricCurve::getCurveLength(void) const
{
	int M = 400;

	double ds = (M_PI + M_PI) / M;
	double s = 0.0;

	double x, dx, ddx;
	double y, dy, ddy;

	double sum = 0.0;

	for(int i = 0; i < M; i++){

		getX(s, x, dx, ddx);
		getY(s, y, dy, ddy);

		double d = sqrt(dx * dx + dy * dy);

		sum += d;

		s += ds;
	}
	sum *= ds;

	return sum;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ParametricCurve::computeArcLength(void) const
{
	return getCurveLength();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int ParametricCurve::
findZero(void (*F)(double, double &, double &, double &),
									 double p, double &t) const
{
	const double atol = 1.0E-14;
	const double rtol = 1.0E-12;

	double x, dx, ddx;

	F(t, x, dx, ddx);

	double f = x - p;
	double res_norm = fabs(f);

	if (res_norm < atol) {
		return 0;
	}

	double tol = atol + rtol * res_norm;

	const int max_itr_num = 50;
	int itr_num = 0;

	do {

		t -= f / dx;

		F(t, x, dx, ddx);

		f = x - p;
		res_norm = fabs(f);

		itr_num++;

	} while (itr_num < max_itr_num && res_norm > tol);

	if (itr_num == max_itr_num) {
		std::cout << "Newton's method failed to converge." << std::endl;
		return -1;
	}

	return itr_num;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int ParametricCurve::
findZero2(void (*F)(double, double &, double &, double &),
										double p, double t1, double t2, double &t) const
{
	const double tol = 1.0E-15;

	double ta, tb;
	if (t1 < t2) {
		ta = t1;
		tb = t2;
	} else {
		ta = t1;
		tb = t2;
	}

	double xa, dxa, ddxa;
	double xb, dxb, ddxb;

	F(ta, xa, dxa, ddxa);
	F(tb, xb, dxb, ddxb);

	double ga = xa - p;
	double gb = xb - p;

	if (fabs(ga + gb) > fabs(ga - gb)) {
		return -1;
	}

	const int max_itr_num = 100;
	int itr_num = 0;

	double error = fabs(ta - tb);
	double tc = 0.5 * (ta + tb);

	double xc, dxc, ddxc;

	while (error > tol && itr_num < max_itr_num) {

		F(tc, xc, dxc, ddxc);
		double gc = xc - p;

		if (fabs(ga + gc) > fabs(ga - gc)) {
			ta = tc;	ga = gc;
		} else {
			tb = tc;	gb = gc;
		}

		tc = 0.5 * (ta + tb);
		error *= 0.5;

		itr_num++;
	};

	if (itr_num == max_itr_num) {
		std::cout << "The bisection method failed to converge." << std::endl;
		return -1;
	}

	t = tc;

	return itr_num;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int ParametricCurve::findXIntersection(double &s, double &x, double y) const
{
	const double atol = 1.0E-14;
	const double rtol = 1.0E-12;

	double xi, eta;
	getPoint(s, xi, eta);

	double f = eta - y;
	double res_norm = fabs(f);

	if(res_norm < atol) {
		x = xi;
		return 0;
	}

	const int max_itr_num = 50;
	int itr_num = 0;

	double tol = atol + rtol * res_norm;

	double dy, ddy;

	do {

		getY(s, eta, dy, ddy);

		if (fabs(dy) < 1.0E-10) {
			printBugInfo("In function findXIntersection(), Zero Jacobian.");
			return -1;
		}

		f = eta - y;
		s -= f / dy;
			
	  res_norm = fabs(f);

		itr_num++;

		//std::cout << "itr_num = " << itr_num 
		//					<< ", res_norm = " << res_norm << std::endl;

	} while((res_norm > tol) && (itr_num < max_itr_num));

	getPoint(s, xi, eta);
	x = xi;

	if(itr_num >= max_itr_num) {
		//printBugInfo("failed to find X intersection by Newton's method.");
		return -1;
	}

	return itr_num;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int ParametricCurve::findXIntersection2(double &s, double &x, double y, 
																					double ta, double tb) const
{
	const double tol = 1.0E-15;

	double tl, tr;
	if (ta < tb) {
		tl = ta;
		tr = tb;
	} else {
		tl = tb;
		tr = ta;
	}

	double pl, ql, pr, qr;

	getPoint(tl, pl, ql);
	getPoint(tr, pr, qr);

	double gl = ql - y;
	double gr = qr - y;

	if (fabs(gl + gr) > fabs(gl - gr)) {
		return -1;
	}

	const int max_itr_num = 100;
	int itr_num = 0;

	double error = fabs(tr - tl);
	double tc = 0.5 * (tl + tr);

	double pc, qc;

	while (error > tol && itr_num < max_itr_num) {

		itr_num++;

		getPoint(tc, pc, qc);
		double gc = qc - y;

		if (fabs(gl + gc) > fabs(gl - gc)) {
			tl = tc;	gl = gc;
		} else {
			tr = tc;	gr = gc;
		}

		tc = 0.5 * (tr + tl);
		error *= 0.5;

		//std::cout << "bisection itr_n = " << itr_num 
		//					<< ", error = " << error << std::endl;
	}

	if (itr_num == max_itr_num) {
		std::cout << "The bisection method failed to converge." << std::endl;
		return -1;
	}

	s = tc;
	getPoint(s, x, qc);

	return itr_num;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int ParametricCurve::findYIntersection(double &s, double x, double &y) const
{
	const double atol = 1.0E-14;
	const double rtol = 1.0E-12;

	double xi, eta;
	getPoint(s, xi, eta);

	double f = xi - x;
	double res_norm = fabs(f);

	if(res_norm < atol) {
		y = eta;
		return 0;
	}

	const int max_itr_num = 50;
	int itr_num = 0;

	double tol = atol + rtol * res_norm;

	double dx, ddx;

	do {

		getX(s, xi, dx, ddx);

		if (fabs(dx) < 1.0E-10) {
			printBugInfo("In function findYIntersection(), Zero Jacobian.");
			return -1;
		}

		f = xi - x;
		s -= f / dx;
			
	  res_norm = fabs(f);

		itr_num++;

		//std::cout << "itr_num = " << itr_num 
		//					<< ", res_norm = " << res_norm << std::endl;

	} while((res_norm > tol) && (itr_num < max_itr_num));

	getPoint(s, xi, eta);
	y = eta;

	if(itr_num >= max_itr_num) {
		//printBugInfo("failed to find Y intersection by Newton's method.");
		return -1;
	}

	return itr_num;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int ParametricCurve::findYIntersection2(double &s, double x, double &y, 
																					double ta, double tb) const
{
	const double tol = 1.0E-15;

	double tl, tr;
	if (ta < tb) {
		tl = ta;
		tr = tb;
	} else {
		tl = tb;
		tr = ta;
	}

	double pl, ql, pr, qr;

	getPoint(tl, pl, ql);
	getPoint(tr, pr, qr);

	double gl = pl - x;
	double gr = pr - x;

	if (fabs(gl + gr) > fabs(gl - gr)) {
		return -1;
	}

	const int max_itr_num = 100;
	int itr_num = 0;

	double error = fabs(tr - tl);
	double tc = 0.5 * (tl + tr);

	double pc, qc;

	while (error > tol && itr_num < max_itr_num) {

		itr_num++;

		getPoint(tc, pc, qc);
		double gc = pc - x;

		if (fabs(gl + gc) > fabs(gl - gc)) {
			tl = tc;	gl = gc;
		} else {
			tr = tc;	gr = gc;
		}
		tc = 0.5 * (tr + tl);
		error *= 0.5;

		//std::cout << "bisection itr_n = " << itr_num 
		//					<< ", error = " << error << std::endl;
	}

	if (itr_num == max_itr_num) {
		std::cout << "The bisection method failed to converge." << std::endl;
		return -1;
	}

	s = tc;
	getPoint(s, pc, y);

	return itr_num;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int ParametricCurve::findClosestPointOnCurve(double p, double q, double &t0, 
																							 double &x0, double &y0) const
{
	const double atol = 1.0E-14;
	const double rtol = 1.0E-12;

	double t = t0;

	double x, y, dx, dy, ddx, ddy;

	getPoint2(t, x, y, dx, dy, ddx, ddy);

	double dis0 = distance(x, y, p, q);

	double f = dx * (x - p) + dy * (y - q);

	double res_norm = fabs(f);

	if (res_norm < atol) {
		t0 = t;
		getPoint(t0, x0, y0);
		return 0;
	}

	double tol = atol + rtol * res_norm;

	const int max_itr_num = 50;
	int itr_num = 0;

	double df = ddx * (x - p) + dx * dx + ddy * (y - q) + dy * dy;

	do {

		if (fabs(df) < 1.0E-10) {
			std::cout << "In findClosestPointOnCurve(), Zero Jacobian." << std::endl;
			return -2;
		}

		t -= f / df;

		getPoint2(t, x, y, dx, dy, ddx, ddy);

		f = dx * (x - p) + dy * (y - q);
		df = ddx * (x - p) + dx * dx + ddy * (y - q) + dy * dy;

		if (distance(x, y, p, q) > dis0 + dis0) {
			double per = M_PI_2 / 180.0;
			t = t0 + Random(-per, per);
		}

		res_norm = fabs(f);

		itr_num++;

		//std::cout << "itr_num = " << itr_num 
		//					<< ", res_norm = " << res_norm << std::endl;

	} while (res_norm > tol && itr_num < max_itr_num);

	t0 = t;
	getPoint(t0, x0, y0);

	if (itr_num == max_itr_num) {
		//std::cout << "failed to find closest point." << std::endl;
		return -1;
	}

	return itr_num;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int ParametricCurve::
findClosestPointOnCurve2(double p, double q, double &t0, 
												 double &x0, double &y0, double ta, double tb) const
{
	const double tol = 1.0E-15;

  double xa, ya, dxa, dya, ddxa, ddya;
  double xb, yb, dxb, dyb, ddxb, ddyb;

  getPoint2(ta, xa, ya, dxa, dya, ddxa, ddya); 
  getPoint2(tb, xb, yb, dxb, dyb, ddxb, ddyb);

  double ga = (xa - p) * dxa + (ya - q) * dya; 
  double gb = (xb - p) * dxb + (yb - q) * dyb;

  if (fabs(ga + gb) > fabs(ga - gb)) {
    return (-1); 
  }

  double xc, yc, dxc, dyc, ddxc, ddyc; 

  int itr_count = 0; 

  int max_itr_num = 100;

  double tc = 0.5 * (ta + tb); 

  double error = fabs(tb - ta);

  while ((error > tol) && (itr_count < max_itr_num)) {

    itr_count++;
    getPoint2(tc, xc, yc, dxc, dyc, ddxc, ddyc);
    double gc = (xc - p) * dxc + (yc - q) * dyc; 

    if (fabs(ga + gc) > fabs(ga - gc)) {
      ta = tc;  ga = gc; 
    } else {
      tb = tc;  gb = gc; 
    }

    tc = 0.5 * (ta + tb); 
    error *= 0.5;

  }

  if (itr_count >= max_itr_num) {
    std::cout << "The bisection method failed to converge." << std::endl; 
		return -1;
  }

	t0 = tc;
	getPoint(t0, x0, y0);

  return itr_count; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ParametricCurve::
setupCurveOnCartesianGrid(const VectorXd &x, const VectorXd &y, 
													int I, int J,	MatrixXb &interior, 
													MatrixXb &irr_edge_x, MatrixXb &irr_edge_y,
													MatrixXd &i_theta_x, MatrixXd &i_theta_y) const
{
	int I1 = I + 1;
	int J1 = J + 1;

	double low[2] = {x[0], y[0]};
	double high[2] = {x[I], y[J]};

	double hx = (high[0] - low[0]) / I;
	double hy = (high[1] - low[1]) / J;
	double h = 0.5 * (hx + hy);

#ifdef DEBUG
	assert(interior.rows() == I1 && interior.cols() == J1);
	assert(irr_edge_x.rows() == I && irr_edge_x.cols() == J1);
	assert(irr_edge_y.rows() == I1 && irr_edge_y.cols() == J);
	assert(i_theta_x.rows() == I && i_theta_x.cols() == J1);
	assert(i_theta_y.rows() == I1 && i_theta_y.cols() == J);
#endif

	MatrixXd dist_func(I1, J1), closest_theta(I1, J1);
	MatrixXb node_tag(I1, J1);

#pragma omp parallel for collapse(2)
	for(int i = 0; i < I1; i++){
		for(int j = 0; j < J1; j++){
			dist_func[i][j] = DBL_MAX;
			closest_theta[i][j] = DBL_MAX;
			node_tag[i][j] = false;
			interior[i][j] = false;
		}
	}

#pragma omp parallel for collapse(2)
	for(int i = 0; i < I; i++){
		for(int j = 0; j < J1; j++){
			irr_edge_x[i][j] = false;
			i_theta_x[i][j] = DBL_MAX;
		}
	}

#pragma omp parallel for collapse(2)
	for(int i = 0; i < I1; i++){
		for(int j = 0; j < J; j++){
			irr_edge_y[i][j] = false;
			i_theta_y[i][j] = DBL_MAX;
		}
	}

	std::queue<Array2i> tube_nd_id0; 

	double curve_length = computeArcLength();

	int M0 = static_cast<int>(curve_length / h + 0.5);
	int M = M0 << 2;

	double dtheta0 = M_2PI / M0;
	double dtheta = M_2PI / M;

	const int width = 2;

	double theta = 0.0;

	for(int l = 0; l < M; l++){

		double xi, eta;
		getPoint(theta, xi, eta);

		int i = static_cast<int>((xi - low[0]) / h + 0.5);
		int j = static_cast<int>((eta - low[1]) / h + 0.5);

		int li = max(i - width, 0);
		int hi = min(i + width, I);
		int lj = max(j - width, 0);
		int hj = min(j + width, J);

		for(int r = li; r <= hi; r++){
			for(int s = lj; s <= hj; s++){

				double dis = distance(x[r], y[s], xi, eta);

				if (dis < dist_func[r][s]) {
					dist_func[r][s] = dis;
					closest_theta[r][s] = theta;
				}

				if (!node_tag[r][s]) {

					Array2i id;
					id[0] = r;	
					id[1] = s;

					tube_nd_id0.push(id);
					node_tag[r][s] = true;
				}
				
			}
		}

		theta += dtheta;
	}

	std::vector<Array2i> tube_nd_id(tube_nd_id0.size()); 
	{
		int k = 0;
		while (!tube_nd_id0.empty()){
			tube_nd_id[k++] = tube_nd_id0.front();
			tube_nd_id0.pop();
		}
	}

	std::queue<Array2i> sch_list;

	for(int l = 0; l < tube_nd_id.size(); l++){

		int i = tube_nd_id[l][0];
		int j = tube_nd_id[l][1];

		bool done = false;

		double t = closest_theta[i][j];

		double xi, eta;
		getPoint(t, xi, eta);

		double dis0 = distance(x[i], y[j], xi, eta);

		int itr_num = findClosestPointOnCurve(x[i], y[j], t, xi, eta);

		double dis1 = distance(x[i], y[j], xi, eta);

		if (itr_num >= 0 && dis1 < dis0 + 1.0E-8) {
			done = true;
		}

		int max_bis_count = 6;

		int bis_count = 0;
		double error = dtheta0;

		while (!done && (bis_count++) < max_bis_count) {

			double ta = closest_theta[i][j] - error;
			double tb = closest_theta[i][j];
			itr_num = findClosestPointOnCurve2(x[i], y[j], t, xi, eta, ta, tb);

			if (itr_num >= 0) {
				dis1 = distance(x[i], y[j], xi, eta);
				if (dis1 < dis0 + 1.0E-8) {
					done = true;
				}
			}

			error *= 0.5;
		}

		bis_count = 0;
		error = dtheta0;

		while (!done && (bis_count++) < max_bis_count) {

			double ta = closest_theta[i][j];
			double tb = closest_theta[i][j] + error;
			itr_num = findClosestPointOnCurve2(x[i], y[j], t, xi, eta, ta, tb);

			if (itr_num >= 0) {
				dis1 = distance(x[i], y[j], xi, eta);
				if (dis1 < dis0 + 1.0E-8) {
					done = true;
				}
			}

			error *= 0.5;
		}

		if (done) {
			rotateToRange(t);
			closest_theta[i][j] = t;
		} else {
			//std::cout << "failed to find cloeset point." << std::endl;
		}

		double tx, ty, nx, ny;

		t = closest_theta[i][j];
		getPoint(t, xi, eta, tx, ty, nx, ny);

		double dx = x[i] - xi;
		double dy = y[j] - eta;

		double prd = dx * nx + dy * ny;

		double dis = sqrt(dx * dx + dy * dy);

		if (prd >= 0.0) {

			dist_func[i][j] = dis;
			interior[i][j] = false;

		} else {

			dist_func[i][j] = - dis;
			interior[i][j] = true;

			Array2i idx;

			idx[0] = i;	
			idx[1] = j;
			sch_list.push(idx);
		}
	}

	const int off_set[4][2] = {{0, 1}, {1, 0}, {-1, 0}, {0, -1}};

	unsigned int sch_count = 0;
	while (!sch_list.empty()) {

		Array2i nd_id = sch_list.front();
		sch_list.pop();

		int i = nd_id[0];
		int j = nd_id[1];

		for(int k = 0; k < 4; k++){

			int i0 = i + off_set[k][0];
			int j0 = j + off_set[k][1];

			if (i0 < 0 || i0 > I || j0 < 0 || j0 > J) {
				printBugInfo("failed to identify interior nodes.");
				exit(1);
			}

			if (!node_tag[i0][j0] && !interior[i0][j0]) {

				Array2i idx;
				idx[0] = i0;	
				idx[1] = j0;
				sch_list.push(idx);

				interior[i0][j0] = true;
			}

			node_tag[i0][j0] = true;
		}

		if ((sch_count++) > I*J) {
			printBugInfo("failed to identify interior nodes.");
			exit(1);
		}
	}


#pragma omp parallel for collapse(2)
	for(int i = 0; i < I; i++){
		for(int j = 0; j < J1; j++){
			if (interior[i][j] != interior[i + 1][j]) {

				irr_edge_x[i][j] = true;

				double t1 = closest_theta[i][j];
				double t2 = closest_theta[i + 1][j];

#ifdef DEBUG
				if (!(fabs(t1 - t2) < M_2PI + M_2PI)) {
					std::cout << "invalid t - x :" <<  t1 << ", " << t2 << std::endl;
					exit(1);
				}
#endif

				if (fabs(t1 - t2) > M_PI_2) {
					int k1 = static_cast<int>(t1 / M_2PI + 0.5);
					int k2 = static_cast<int>(t2 / M_2PI + 0.5);
					t1 -= k1 * M_2PI;
					t2 -= k2 * M_2PI;
				}

#ifdef DEBUG
				if (!(fabs(t1 - t2) < M_PI_4)) {
					std::cout << t1 << ", " << t2 << std::endl;
				}
				assert(fabs(t1 - t2) < M_PI_4);
#endif

				double ta = min(t1, t2);
				double tb = max(t1, t2);

				double x1, y1, x2, y2;

				getPoint(t1, x1, y1);
				getPoint(t2, x2, y2);

				double alpha = (y[j] - y1) / (y2 - y1 + 1.0E-12);
				double beta = 1.0 - alpha;

				double t0 = alpha * t2 + beta * t1;

				double t, xi, eta;

				t = t0;
				getPoint(t, xi, eta);

				const int max_newton_count = 10;
				int newton_count = 0;

				bool done = false;

				do {

					int itr_num = findXIntersection(t, xi, y[j]);

					if (itr_num >= 0 && xi >= x[i] && xi <= x[i+1]) {
						done = true;
					} else {
						t = Random(ta, tb);
						getPoint(t, xi, eta);
					}

					newton_count++;

				} while (!done && newton_count < max_newton_count);

				if (!done) {

					const int max_bis_count = 6;

					int bis_count = 0;
					double diam = fabs(t0 - ta);

					while (!done && (bis_count++) < max_bis_count) {
						double ta1 = t0 - diam;
						double tb1 = t0;
						int itr_num = findXIntersection2(t, xi, y[j], ta1, tb1);
						if (itr_num >= 0 && xi >= x[i] && xi <= x[i+1]) {
							done = true;
						}
						diam *= 0.5;
					}

					bis_count = 0;
					diam = fabs(tb - t0);

					while (!done && (bis_count++) < max_bis_count) {
						double ta1 = t0;
						double tb1 = t0 + diam;
						int itr_num = findXIntersection2(t, xi, y[j], ta1, tb1);
						if (itr_num >= 0 && xi >= x[i] && xi <= x[i+1]) {
							done = true;
						}
						diam *= 0.5;
					}

					if (!done) {
						std::cout << "yj = " << y[j] << std::endl;
						std::cout << "x1, y1 = " << x1 << ", " << y1 << std::endl;
						std::cout << "x2, y2 = " << x2 << ", " << y2 << std::endl;
						std::cout << "failed to find intersection." << std::endl;
						t = 0.5 * (t1 + t2);
						//exit(1);
					}
				}

				//assert(x[i] <= xi && x[i+1] >= xi);

				rotateToRange(t);
				i_theta_x[i][j] = t;

			}
		}
	}


#pragma omp parallel for collapse(2)
	for(int i = 0; i < I1; i++){
		for(int j = 0; j < J; j++){
			if (interior[i][j] != interior[i][j + 1]) {

				irr_edge_y[i][j] = true;

				double t1 = closest_theta[i][j];
				double t2 = closest_theta[i][j+1];

#ifdef DEBUG
				if (!(fabs(t1 - t2) < M_2PI + M_2PI)) {
					std::cout << "invalid t - y :" <<  t1 << ", " << t2 << std::endl;
					exit(1);
				}
#endif

				if (fabs(t1 - t2) > M_PI_2) {
					int k1 = static_cast<int>(t1 / M_2PI + 0.5);
					int k2 = static_cast<int>(t2 / M_2PI + 0.5);
					t1 -= k1 * M_2PI;
					t2 -= k2 * M_2PI;
				}

#ifdef DEBUG
				if (!(fabs(t1 - t2) < M_PI_4)) {
					std::cout << t1 << ", " << t2 << std::endl;
				}
				assert(fabs(t1 - t2) < M_PI_4);
#endif

				double ta = min(t1, t2);
				double tb = max(t1, t2);

				double x1, y1, x2, y2;
				getPoint(t1, x1, y1);
				getPoint(t2, x2, y2);

				double alpha = (x[i] - x1) / (x2 - x1 + 1.0E-12);
				double beta = 1.0 - alpha;

				double t0 = alpha * t2 + beta * t1;

				double t, xi, eta;
				t = t0;
				getPoint(t, xi, eta);

				const int max_newton_count = 10;
				int newton_count = 0;

				bool done = false;

				do {

					int itr_num = findYIntersection(t, x[i], eta);

					if (itr_num >= 0 && eta >= y[j] && eta <= y[j + 1]) {
						done = true;
					} else {
						t = Random(ta, tb);
						getPoint(t, xi, eta);
					}

					newton_count++;

				} while (!done && newton_count < max_newton_count);

				if (!done) {

					const int max_bis_count = 6;
					int bis_count = 0;
					double diam = fabs(ta - t0);

					while (!done && (bis_count++) < max_bis_count) {
						double ta1 = t0 - diam;
						double tb1 = t0;
						int itr_num = findYIntersection2(t, x[i], eta, ta1, tb1);

						if (itr_num >= 0 && eta >= y[j] && eta <= y[j + 1]) {
							done = true;
						}
						diam *= 0.5;
					}

					bis_count = 0;
					diam = fabs(tb - t0);

					while (!done && (bis_count++) < max_bis_count) {
						double ta1 = t0;
						double tb1 = t0 + diam;
						int itr_num = findYIntersection2(t, x[i], eta, ta1, tb1);

						if (itr_num >= 0 && eta >= y[j] && eta <= y[j + 1]) {
							done = true;
						}
						diam *= 0.5;
					}

					if(!done) {
						std::cout << "xi = " << x[i] << std::endl;
						std::cout << "x1, y1 = " << x1 << ", " << y1 << std::endl;
						std::cout << "x2, y2 = " << x2 << ", " << y2 << std::endl;
						std::cout << "failed to find intersection." << std::endl;

						t = 0.5 * (t1 + t2);
						//exit(1);
					}
				}

				//assert(y[j] <= eta && y[j+1] >= eta);

				rotateToRange(t);
				i_theta_y[i][j] = t;

			}
		}
	}
}

#ifdef USE_OPENGL
#include "NewMathGL.h"
#include "MathGLUT2d.h"
#endif

#ifdef USE_MATHGL
#include "MathGL.h"
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ParametricCurve::plotCurve(void) const
{
	int K = 2000;
	int K1 = K + 1;

	double *x = new double[K1];
	double *y = new double[K1];

	double ds = M_2PI / K;
	double s = 0.0;
	
	for(int k = 0; k < K1; k++){
		getPoint(s, x[k], y[k]);
		s += ds;
	}

#ifdef USE_OPENGL
	mglSetColor(0, 0, 0);
	mglSetLineWidth(2.0);
	mglPlotLines(x, y, K1);
#endif

#ifdef USE_MATHGL
	mglColor(MGL_BLACK);
	mglPlot2dLines(x, y, K1);
#endif

	delete[] x;
	delete[] y;
}


//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


