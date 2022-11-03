/*=============================================================================
*   
*   Filename : solveDiffusionEqn2d2.c
*   Creator : Han Zhou
*   Date : 11/21/21
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

#include <fstream> 
#include <sstream>
#include <unistd.h> 

#include <GL/glut.h>

#include <sys/stat.h>
#include <sys/wait.h>

#ifdef USE_GL2PS
#include "gl2ps.h"
#endif

#include "MathGLUT2d.h"
#include "NewMathGL.h"

#include "Variables.h"
#include "MathTools.h"

#include "CartesianGrid.H"
#include "ParametricCurve.H"
#include "StarCurve.H"
#include "EllipseCurve.H"

#include "MyTimer.h"

#include "DiffusionEqn.H"

static int I = 32;
static int J = 32;

static double box_radius = 1.2;
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

static int bdry_node_n = 64;
static VectorX2d bdry_node_crd;

static DiffusionEqn solver;

static StarCurve curve;
//static EllipseCurve curve;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

static double T = 0.5;
static double t = 0.0;
static double dt = 0.01;
static double _epsilon = 1.0;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double Ui(double x, double y, double t)
{
	double theta = M_PI / 6.0;

	double cs = cos(theta);
	double sn = sin(theta);
	double xi = (cs * x + sn * y) / sqrt(_epsilon);
	return exp(-t) * sin(xi);

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double Ue(double x, double y, double t)
{
	return 0.0;

	double theta = M_PI / 3.0;
	double cs = cos(theta);
	double sn = sin(theta);
	double xi = (cs * x + sn * y) / sqrt(_epsilon);

	return exp(-t) * sin(xi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

static int gl2ps_plot_count = 0;
static int gl2ps_check_interval = 0;
static double gl2ps_plot_interval = 0.2;

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

	if (plot_bdry_node) {
		solver.plotControlPoints();
	}

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

	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void startSolvingPDE(void)
{
	//curve.setShape(0.7, 0.5, 0.2, M_PI/3.0, 7.0);
	curve.setShape(0.8, 0.8, 0.2, M_PI/4.0, 5.0);

	//curve.setShape(-0.1, 0.2, 0.7, 0.5, M_PI/3.0);
	//curve.setShape(0.0, 0.0, 1.0, 0.5, 0.0);
	//curve.setShape(0.0, 0.0, 1.0, 1.0, 0.0);


	solver.defineGrid(low, high, I, J);
	solver.setCurve(curve);

	solver.setExactSolution(Ui, Ue);
	solver.setCoefficient(_epsilon);

	t = 0.0;
	T = 1.0;

	solver.setPrintInterval(0.1);

	solver.initialize(t, T);

	dt = solver.getTimeStep();

	double h = solver.getMeshParameter();

	std::cout << "\ngrid size : " << I << " X " << J 
						<< ", dx = dy = " << h
						<< ", epsilon = " << _epsilon 
						<< std::endl;
	curve.printCurveInfo();

#ifdef USE_OPENGL
	plotAll();
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

	for(int k = 0; k < 5; k++){
		solveIBVP();
		n <<= 1;
		I = J = n;
	}

#else 

	double lo[2], hi[2];
	computeDrawingArea(low, high, lo, hi, 1.1);

  mgl::initialize(argc, argv);
  mgl::createWindow(800, 800, "Heat equation"); 
  mgl::setDrawingArea(lo[0], hi[0], lo[1], hi[1]); 
  mgl::setSpecialFunc(special_func);
  mgl::setKeyboardFunc(key_func); 
  mgl::enterMainLoop(); 

#endif

	return EXIT_SUCCESS;
}
