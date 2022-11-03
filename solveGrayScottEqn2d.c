/*=============================================================================
*   
*   Filename : solveGrayScottEqn2d.c
*   Creator : Han Zhou
*   Date : 11/22/21
*   Description : 
*
=============================================================================*/
   
#include <iostream>
#include <fstream> 
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include <queue>
#include <unistd.h> 
#include <sys/stat.h>
#include <sys/wait.h>

#include <GL/glut.h>

#ifdef USE_GL2PS
#include "gl2ps.h"
#endif

#ifdef USE_OPENGL
#include "MathGLUT2d.h"
#include "NewMathGL.h"
#endif

#include "Variables.h"
#include "MathTools.h"

#include "CartesianGrid.H"
#include "ParametricCurve.H"

#include "StarCurve.H"
#include "NucleusCurve.H"
#include "EllipseCurve.H"
#include "ClosedSplineCurve.H"

#include "MyTimer.h"

#include "GrayScottEqn.H"

static int I = 32;
static int J = 32;

static double box_radius = 2.2;

static double low[2] = {-box_radius, -box_radius};
static double high[2] = {box_radius, box_radius};

static bool is_solving_pde = false;

static bool plot_grid = false;
static bool plot_curve = true;
static bool plot_bdry_node = false;
static bool plot_grid_data = true;
static bool plot_interior_node = false;
static bool plot_intersect_node = false;
static bool plot_with_colormap = false;

static int gl2ps_plot_count = 0;
static int gl2ps_check_interval = 0;
static double gl2ps_plot_interval = 0.2;

static int step_count = 0;

static int bdry_node_n = 64;
static VectorX2d bdry_node_crd;

static GrayScottEqn solver;

//static NucleusCurve curve;
//static StarCurve curve;
static EllipseCurve curve;

//static ClosedSplineCurve curve;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

static double T = 0.5;
static double t = 0.0;
static double dt = 0.01;

static const double _gama = 0.037; // feed rate
static const double _kappa = 0.06; // death rate

//static const double _gama = 0.024;
//static const double _kappa = 0.06;

static const double _eps0 = 0.01;
static const double _eps1 = 0.008;
static const double _eps2 = 0.004;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double U(double x, double y, double t)
{
	extern double V(double, double, double);

	return (1.0 - 2.0 * V(x, y, 0.0));
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double V(double x, double y, double t)
{
	//double r = sqrt(x * x + y * y);
	//double xi = 20.0 * (r - 0.05);
	//return 1.0 / (1.0 + exp(xi));

	if (fabs(x) < 0.25 && fabs(y) < 0.25) {
		double sn = sin(4.0 * M_PI * x) * sin(4.0 * M_PI * y);
		return 0.25 * sn * sn;
	}

	return 0.0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void savePlot(int count)
{
#ifdef USE_GL2PS
  int max_n = 10000;

  if (count >= max_n) {
    return;
  }

  char filename[64];

  int k = access("./EPS", F_OK & W_OK);
  if (k != 0) {
    mkdir("./EPS", 0700);
  }

  //sprintf(filename, "./EPS/%d.eps.gz", 1000 + count);
  sprintf(filename, "./EPS/%d.eps", 1000 + count);

  FILE *fp = 0;
  int buffsize = 0;
  int state = GL2PS_OVERFLOW;

  //extern void display(void); 

  fp = fopen(filename, "wb");
  while (state == GL2PS_OVERFLOW){
    buffsize += 1024 * 1024;
    gl2psBeginPage("test", "gl2psTestSimple", NULL, GL2PS_EPS, GL2PS_BSP_SORT,
       GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT | GL2PS_COMPRESS, 
       GL_RGBA, 0, NULL, 0, 0, 0,  buffsize, fp, filename); 
    mgl::display();
    state = gl2psEndPage(); 
  }
  fclose(fp); 

	std::cout << filename << " saved." << std::endl;
#endif
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void plotAll(void)
{
#ifdef USE_OPENGL
	mglNewPage();

	if (plot_grid_data) {
		solver.plotSolution();
	}

	if (plot_grid) {
		solver.plotGridLines();
	}

	if (plot_interior_node) {
		solver.plotInteriorNodes();
	}

	if (plot_intersect_node) {
		solver.plotIntersectPoints();
	}

	if (plot_curve) {
		curve.plotCurve();
	}

	//solver.plotNormalsAndTangents();

	if (plot_bdry_node) {
		solver.plotControlPoints();
	}

	showCurrentTime2d(low, high, t);

	mglSetColor(0, 0, 0);
	mglSetLineWidth(2.0);
	mglPlotRectangle(low[0], low[1], high[0], high[1]);

	mglFlush();
#endif
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void solvePDE(void)
{
	if (t < T && is_solving_pde) {

		solver.advance(t, dt);

		t += dt;

		plotAll();

		if ((++step_count)%gl2ps_check_interval == 0) {
			savePlot(gl2ps_plot_count++);
		}

	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void startSolvingPDE(void)
{
  const double ctrl_node[28][2] = {
    {1.000e+00, 2.763e-01},
    {9.639e-01, 5.482e-01},
    {8.613e-01, 7.787e-01},
    {7.076e-01, 9.328e-01}, 
    {5.263e-01, 9.868e-01}, 
    {3.070e-01, 9.118e-01}, 
    {8.772e-02, 7.724e-01}, 
    {-1.316e-01, 6.974e-01}, 
    {-3.553e-01, 7.500e-01}, 
    {-5.789e-01, 8.026e-01},
    {-7.895e-01, 7.462e-01},
    {-9.436e-01, 5.921e-01},
    {-1.000e+00, 3.816e-01},
    {-9.808e-01, 1.440e-01}, 
    {-9.238e-01, - 8.645e-02}, 
    {-8.308e-01, - 3.026e-01}, 
    {-7.045e-01, - 4.980e-01}, 
    {-5.488e-01, - 6.667e-01}, 
    {-3.684e-01, - 8.035e-01},
    {-1.689e-01, - 9.043e-01},
    {4.381e-02, - 9.661e-01}, 
    {2.632e-01, - 9.868e-01},
    {4.271e-01, - 9.552e-01},
    {5.829e-01, - 8.618e-01},
    {7.226e-01, - 7.113e-01}, 
    {8.392e-01, - 5.113e-01}, 
    {9.270e-01, - 2.717e-01}, 
    {9.815e-01, - 4.763e-03},
  };

	double x[28], y[28];
	for(int k = 0; k < 28; k++){
		x[k] = ctrl_node[k][0];
		y[k] = ctrl_node[k][1];
	}

	//curve.define(x, y, 28);

	//curve.setShape(0.0, 0.0, 1.0, 0.8, 0.0);
	curve.setShape(0.0, 0.0, 2.0, 1.6, 0.0);

	//curve.setShape(0.7, 0.5, 0.2, M_PI/3.0, 7.0);
	//curve.setShape(0.8, 0.8, 0.2, 0.0, 4.0);
	//curve.setShape(0.8, 0.8, 0.2, M_PI/4.0, 5.0);

	//curve.setShape(1.0, 0.1, 2.0, 3.0); // Nucleus

	solver.defineGrid(low, high, I, J);
	solver.setCurve(curve);

	solver.setExactSolution(U, V);
	solver.setCoefficient(_eps0, _eps1, _eps2, _gama, _kappa);


	t = 0.0;
	//T = 40.0;
	T = 60.0;

	solver.setPrintInterval(1.0);
	solver.setPlotInterval(1.0);

	solver.initialize(t, T);

	dt = solver.getTimeStep();

	std::cout << "\ngrid size : " << I << " X " << J 
						<< ", dx = dy = " << (high[0] - low[0]) / I
						<< std::endl << std::endl;;

#ifdef USE_OPENGL
	solver.setPlotOption(1);
	plotAll();
#endif

#ifdef USE_GL2PS
	step_count = 0;
	gl2ps_check_interval = static_cast<int> (gl2ps_plot_interval / dt + 0.5);
	savePlot(gl2ps_plot_count++);
#endif
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void solveIBVP(void)
{
	startSolvingPDE();

	is_solving_pde = true;

	Timer timer;
	timer.start();

	solver.solve();


	timer.stop();
	double total_t = timer.getUsedTime();
	int num = solver.getTimeStepNumber();

	std::cout << "CPU time = " << total_t << " secs"
					  << ", averaged cpu time = " << total_t / num 
						<< ", secs. " << std::endl;

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef USE_OPENGL
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int key_func(unsigned char ch) 
{
	switch (ch) {

		case 'q':
		case 'Q':
		case 27:

			exit(0);
			return 1;
			break;

		case 13:

			startSolvingPDE();

			return 1;
			break;

		case ' ':
			is_solving_pde = !is_solving_pde;
			glutIdleFunc(solvePDE);
			return 1;
			break;

		case '+':
		case '=':

			if (I < 2048 && J < 2048) {
				I <<= 1;
				J <<= 1;
			}

			if (is_solving_pde) {
				startSolvingPDE();
			}

			plotAll();
			return 1;
			break;

		case '-':
		case '_':

			if (I > 32 && J > 32) {
				I >>= 1;
				J >>= 1;
			}

			if (is_solving_pde) {
				startSolvingPDE();
			}

			plotAll();
			return 1;
			break;

		case 's':
		case 'S':
			savePlot(gl2ps_plot_count++);
			return 1;
			break;

		case 'c':
		case 'C':

			plot_curve = !plot_curve;
			plotAll();
			return 1;
			break;

		case 'n':
		case 'N':

			plot_intersect_node = !plot_intersect_node;
			plotAll();
			return 1;
			break;

		case 'b':
		case 'B':

			plot_bdry_node = !plot_bdry_node;
			plotAll();
			return 1;
			break;

		case 'g':
		case 'G':

			plot_grid = !plot_grid;
			plotAll();
			return 1;
			break;

		case 'd':
		case 'D':

			plot_grid_data = !plot_grid_data;
			plotAll();
			return 1;
			break;

		case 'i':
		case 'I':

			plot_interior_node = !plot_interior_node;
			plotAll();
			return 1;
			break;

		default: 
			return 1;
	}

	return 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int special_func(int key) 
{
  if (key == GLUT_KEY_F1) {

		solver.setPlotOption(0);
		plotAll();
		return 1;

	} else if (key == GLUT_KEY_F2) {

		solver.setPlotOption(1);
		plotAll();
		return 1;

	} else if (key == GLUT_KEY_F3) {
		
		is_solving_pde = true;
		solvePDE();
		plotAll();
		return 1;

	} else if (key == GLUT_KEY_F4) {

		plotAll();
		return 1;

	} else if (key == GLUT_KEY_F5) {

		startSolvingPDE();
		return 1;

	}

	return 0;
}
#endif

//=============================================================================
int main (int argc, char* argv[])
{
	int n = 64;
	if (argc > 1) {
		n = atoi(argv[1]);
	}

	I = J = n;

  //std::cout.setf(std::ios::uppercase);
  std::cout.setf(std::ios::scientific);
  std::cout.setf(std::ios::showpoint); 
  std::cout.setf(std::ios::right);
  std::cout.precision(2);

#ifndef USE_OPENGL

#ifdef USE_MATHGL
	//mglInitWindowSize(800, 800);
	mglInitDrawingArea(low[0], low[1], high[0], high[1]); 
	mglInitialization(argc, argv);
#endif

	solveIBVP();

#else 

	double lo[2], hi[2];
	computeDrawingArea(low, high, lo, hi, 1.1);

  mgl::initialize(argc, argv);
  mgl::createWindow(800, 800, "Gray Scott model"); 
  mgl::setDrawingArea(lo[0], hi[0], lo[1], hi[1]); 
  mgl::setSpecialFunc(special_func);
  mgl::setKeyboardFunc(key_func); 
  mgl::enterMainLoop(); 

#endif

	return EXIT_SUCCESS;
}
