/*=============================================================================
*   
*   Filename : Interpolation.h
*   Creator : Han Zhou
*   Date : 09/06/21
*   Description : 
*
=============================================================================*/
   

extern double makeLinearInterpolation(const double t[2], const double f[2], 
																			double xi);

extern void makeLinearInterpolation(const double t[2], const double f[2], 
  										 							double xi, double g[2]);

extern double makeQuadraticInterpolation(const double t[3], const double f[3], 
  																 			 double xi);
 
extern void makeQuadraticInterpolation(const double t[3], const double f[3], 
                       					 			 double xi, double g[2]);

extern void makeQuadraticInterpolation2(const double t[3], const double f[3], 
                       					 			  double xi, double g[3]);
 
extern bool findZeroOnQuadraticInterpolant(double &x, const double t[3],
  											 						 			 const double f[3]);
 
extern bool findZerosOnQuadraticInterpolant(const double t[3],
															 						  const double f[3], double z[2]);

extern double makeCubicInterpolation(const double t[4], const double f[4], 
																		 double xi);
 
extern void makeCubicInterpolation(const double t[4], const double f[4], 
                       			 			 double xi, double g[3]);
 
extern bool findZeroOnCubicInterpolant(double &x, const double t[4],
  											 				 			 const double f[4]);

extern double computeDistanceToParabola(const double t[3], const double f[3], 
																				double, double);

extern double computeDistanceToCubicInterpolant(const double t[4], const double f[4], 
																								double, double);
