// Copyright (c) Wenjun Ying, School of Mathematical Sciences and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is free software; as a special exception the author gives
// unlimited permission to copy and/or distribute it, with or without
// modifications, as long as this notice is preserved.

// This file is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

/**
 *************************************************************************** 
 * The OpenGL-based C/C++ codes were written for students by Wenjun Ying,  *
 * School of Mathematical Sciences and Institute of Natural Sciences,      *
 * Shanghai Jiao Tong University, Shanghai 200240.                         *
 ***************************************************************************
 **/

#ifndef __NewMathGL_h_IS_INCLUDED__
#define __NewMathGL_h_IS_INCLUDED__

extern void mglPlotAll(void); 

extern void mglClear(void); 
extern void mglFlush(void);
extern void mglNewpage(void); 
extern void mglNewPage(void);

extern void mglSetColor(double c); 

extern void mglSetPointSize(double sz); 
extern void mglSetLineWidth(double wd);
extern void mglSetColor(double r, double g, 
                        double b);

extern void mglSetColormap(const double c[2]); 
extern void mglResetColormap(const double c[2]);

extern void mglSetColormap(double lo, double hi);
extern void mglResetColormap(double lo, double hi);

extern void mglPlotPoint(double x, double y);
extern void mglPlotLine(double x1, double y1, 
                        double x2, double y2); 
extern void mglPlotDashLine(double x1, double y1,
                            double x2, double y2);

extern void mglPlotCircle(double x, double y, double r); 
extern void mglPlotDashCircle(double x, double y, double r);

extern void mglPlotTriangle(double x1, double y1,
                            double x2, double y2,
                            double x3, double y3); 
extern void mglPlotRectangle(double lx, double ly,
                             double hx, double hy); 

extern void mglPlotQuadrilateral(double x1, double y1, 
                                 double x2, double y2, 
                                 double x3, double y3, 
                                 double x4, double y4);

extern void mglPlotSolidCircle(double x, double y, double r);

extern void mglPlotSolidTriangle(double x1, double y1,
                                 double x2, double y2, 
                                 double x3, double y3); 

extern void mglPlotSolidRectangle(double lx, double ly, 
                                  double hx, double hy); 

extern void mglPlotSolidQuadrilateral(double x1, double y1, 
                                      double x2, double y2, 
                                      double x3, double y3,
                                      double x4, double y4); 

extern void mglPlotSolidTriangle(const double coord[3][2],
                                 const double u[3]); 

extern void mglPlotSolidTriangle(double x1, double y1, double v1,
                                 double x2, double y2, double v2, 
                                 double x3, double y3, double v3);

extern void mglPlotPolygon(const double x[], const double y[], int n);
extern void mglPlotSolidPolygon(const double x[], const double y[], int n); 
extern void mglPlotPoints(const double x[], const double y[], int n);
extern void mglPlotLines(const double x[], const double y[], int n); 
extern void mglPlotString(double x, double y, const char *str);
extern void mglPlotString(double x, double y, double value); 
extern void mglPlotString(double x, double y, int idx); 

extern void mglPlotIsolines(const double coord[3][2], const double u[3],
                            double min_u, double max_u, int isoval_n);

#endif
