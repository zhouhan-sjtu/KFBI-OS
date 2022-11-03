    /**********************************************************************\
     * Copyright (c) 2003-2005, Wenjun Ying, Duke University, Math Dept.  *
     *                                                                    *
     *                       ALL RIGHTS RESERVED                          *
     *                                                                    *
     * Permission to use, copy, modify, and distribute this software for  *
     * any purpose and without fee is hereby granted, provided that the   *
     * above copyright notice appear in all copies and that both the      *
     * copyright notice and this permission notice appear in supporting   *
     * documentation, and that the name of MathGL (Mathemtical Graphics   *
     * Library) not be used in advertising or publicity pertaining to     *
     * distribution of the software without specific, written prior       *
     * permission.                                                        *
     *                                                                    *
     * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"  *
     * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,   *
     * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR   *
     * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL Y'S MATHGL BE *
     * LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT, SPECIAL, INCIDENTAL,  *
     * INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND, OR ANY DAMAGES      *
     * WHATSOEVER, INCLUDING WITHOUT LIMITATION, LOSS OF PROFIT, LOSS OF  *
     * USE, SAVINGS OR REVENUE, OR THE CLAIMS OF THIRD PARTIES, WHETHER   *
     * OR NOT Y'S MATHGL HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH      *
     * LOSS, HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, ARISING OUT   *
     * OF OR IN CONNECTION WITH THE POSSESSION, USE OR PERFORMANCE OF     *
     * THIS SOFTWARE.                                                     *
     *   __      __      __     _______  _      _    _______    _         *
     *  /  \    /  \    /  \   |__   __|| |    | |  / ______|  | |        *
     *  | | \  / | |   / /\ \     | |   | |____| | / /   ____  | |        *
     *  | |\ \/ /| |  / /__\ \    | |   |  ____  | | |  |_  _| | |    __  *
     *  | | \__/ | | /  ____  \   | |   | |    | | | \___/  |  | |___/ /  *
     *  |_|      |_|/__/    \__\  |_|   |_|    |_|  \____/|_|  |______/   *
     *                                                                    *
     * WENJUN YING SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT  *
     * NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND      *
     * FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER *
     * IS ON AN "AS IS" BASIS, AND WENJUN YING HAS NO OBLIGATION TO       *
     * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR            *
     * MODIFICATIONS.                                                     *
     *                                                                    *
    \**********************************************************************/

/**
 * @file $RCSfile: MathGL.h,v $
 *
 * Copyright (C) 2004 Wenjun Ying, Mathematics Department, Duke University
 *
 * $Id: MathGL.h,v 1.9 2005/03/25 08:00:43 ying Exp $
 *
 * Description:
 *   This file contains ...
 **/

#ifndef __MathGL_h_IS_INCLUDED__
#define __MathGL_h_IS_INCLUDED__ // Mon Apr 19 22:14:37 EDT 2004

#include "X11Colors.h"

// Creation Date: Fri May  3 19:17:20 EDT 2002
// Copyright (c) Wenjun Ying, 2001 (Duke University)

// This program is freely distributable without licensing fees and is
// provided without guarantee or warrantee expressed or implied. This
// program is -not- in the public domain.

// MathGL revision history:

// MATHGL_API_VERSION is updated to reflect incompatible MathGL API
// changes (interface changes, semantic changes, deletions, or
// additions).

// MATHGL_API_VERSION = 1 First public release of MathGL. 10/01/2001
// MATHGL_API_VERSION = 2 First public release of MathGL. 05/01/2002
// MATHGL_API_VERSION = 3 First public release of MathGL. 03/22/2003

#ifndef MATHGL_API_VERSION    // Allow this to be overriden.
#define MATHGL_API_VERSION 3
#endif 

// MathGL implementation revision history:

// MATHGL_OPENGL_IMPLEMENTATION is updated to reflect both MathGL
// API revisions and implementation revisions (i.e., bug fixes).

// MATHGL_OPENGL_IMPLEMENTATION = 1 Ying's first public release of
// MathGL OpenGL-based implementation.   10/01/2001

// MATHGL_OPENGL_IMPLEMENTATION = 2 Ying's first public release of
// MathGL OpenGL-based implementation.   05/01/2002

#ifndef MATHGL_OPENGL_IMPLEMENTATION  // Allow this to be overriden.
#define MATHGL_OPENGL_IMPLEMENTATION 2
#endif 

// Parent process status.

#define MGL_RUNNING         0
#define MGL_FLUSHING        1 
#define MGL_WAITING         2
#define MGL_DUMPING         3

#define MGL_BLACK           0
#define MGL_BLUE            1
#define MGL_RED             2
#define MGL_GREEN           3
#define MGL_YELLOW          4
#define MGL_MAGENTA         5
#define MGL_SKYBLUE         6
#define MGL_GOLDEN          7 
#define MGL_PURPLE          8 
#define MGL_ORANGE          9 
#define MGL_BROWN          10 
#define MGL_CYAN           11 
#define MGL_GRAY           12
#define MGL_GREY           12
#define MGL_WHITE          13
#define MGL_DARKRED        14 
#define MGL_DARKGREEN      15 
#define MGL_LIMEGREEN      16
#define MGL_TURQUOISE      17 
#define MGL_VIOLET         18 
#define MGL_PINK           19
#define MGL_NAVY           20 

#define MGL_REFERENCE_GRID  0
#define MGL_COORD_FRAME_BOX 1
#define MGL_COORD_WIRE_FRAME_BOX 1 
#define MGL_COORD_OPAQUE_FRAME_BOX 2
#define MGL_BLACK_WHITE_MODE 3 
#define MGL_COLOR_MODE 4 

typedef void (*MGL_WORK_PROC)(void);

// MathGL Initialization.

// If ((0 < w <= 1.0) && (0 < h <= 1.0)) {
//   each of the window dimensions (w/h) should be a fraction of the screen.
// } else if ((1.0 < w < 72.0) && (1.0 < h < 72.0)) {
//   the units of the window dimensions should be "inches".
// } else {
//   the units of the window dimensions should be "pixels".
// }
extern void mglInitDrawingAreaSize(double w, double h); 
extern void mglInitWindowSize(double w, double h); 

extern void mglRemapOutput(void); 
extern void mglEnableOutRemap(void); 
extern void mglDisableOutRemap(void); 
extern void mglRemapStandardOut(void);

extern void mglInitSpaceDimension(int dim);

extern void mglInitWindowSize(int w, int h);

extern void mglRemapOutput(const char []); 
extern void mglRemapStandardOut(const char []);

extern void mglInitGraphNote(const char *note);
extern void mglInitGraphTitle(const char *title);
extern void mglInitWindowTitle(const char *title);
extern void mglInitSliderTitle(const char title[]);

extern void mglInitXLabel(const char xlabel[]);
extern void mglInitYLabel(const char xlabel[]);

extern void mglInitXaxisLabel(const char xlabel[]); 
extern void mglInitYaxisLabel(const char xlabel[]);
extern void mglInitXAxisLabel(const char xlabel[]); 
extern void mglInitYAxisLabel(const char xlabel[]);

extern void mglInitXYLabels(const char s1[], 
                            const char s2[]);
extern void mglInitAxisLabels(const char s1[], 
                              const char s2[]);

extern void mglInitBackgroundColor(int); 
extern void mglInitForegroundColor(int); 

extern void mglRunGuiInBackground(void); 
extern void mglRunGuiInForeground(void); 

extern void mglInitBackgroundColor(double r, double g, double b); 
extern void mglInitForegroundColor(double r, double g, double b);

extern void mglInitDrawingArea(const double low[3], const double high[3]); 
extern void mglInitDrawingArea(const double low[3], const double high[3], int);

extern void mglInitDrawingArea2d(const double low[2], const double high[2]);
extern void mglInitDrawingArea(double, double, double, double);

extern void mglInitColormap(int type, double cold, double hot); 
extern void mglInitColormap(double cold, double hot);

extern void mglInitEmissionRate(float rate); 
extern void mglInitDiffuseRate(float rate);

extern int  mglInitialization(int argc, char **argv);

extern void mglInitMemoryCleaner(void (*c)(void)); 

extern void mglEnableGui(int argc, char *argv[]);
extern void mglEnableGui(void);

extern void mglStartGui(int argc, char *argv[]);
extern void mglRunGui(int argc, char *argv[]); 

extern int mglGetMainGuiProcessId(void); 

extern void mglCloseGuiProcess(void);
extern void mglTurnOffGuiProcess(void);
extern void mglTerminateGuiProcess(void);

extern void mglInitMainWork(MGL_WORK_PROC);

extern void mglInitEnable(int state);
extern void mglInitDisable(int state); 

extern void mglInitEnableStretching(void);
extern void mglInitDisableStretching(void); 

extern void mglEnableAutomaticAdjustment(void);
extern void mglDisableAutomaticAdjustment(void);

extern void mglDisableDrawingAreaAdjusting(void);
extern void mglEnableDrawingAreaAdjusting(void);

extern void mglDisableGraphTitleAndAxes(void);
extern void mglEnableGraphTitleAndAxes(void);

extern void mglDisableTailoring(void); 
extern void mglEnableTailoring(void); 

extern void mglDisableZooming(void); 
extern void mglEnableZooming(void); 

extern void mglDisableSlider(void);
extern void mglEnableSlider(void); 

extern void mglHint(const char str[]); 
extern void mglMessage(const char str[]); 
extern void mglWarning(const char str[]); 

// MathGL Setup after mglInitialization(argc, argv).

extern void mglAddTextToEpsFile(const char []); 

extern void mglEnable(int state);
extern void mglDisable(int state);

extern void mglSetTimeTag(float t); 
extern void mglSetTimeTag(double t); 
extern void mglSetTimeStamp(float t); 
extern void mglSetTimeStamp(double t);

extern void mglSetGraphNote(const char *note); 
extern void mglSetGraphTitle(const char *title);
extern void mglSetSliderTitle(const char *title);
extern void mglSetWindowTitle(const char *title); 

extern void mglResetGraphNote(const char note[]); 
extern void mglResetGraphTitle(const char title[]);

extern void mglResetWindowTitle(const char title[]);
extern void mglResetSliderTitle(const char title[]);

extern void mglSetXYLabels(const char *s1, const char *s2); 
extern void mglSetAxisLabels(const char *s1, const char *s2);

extern void mglResetXYLabels(const char s1[], const char s2[]); 
extern void mglResetAxisLabels(const char s1[], const char s2[]); 

extern void mglSetBackgroundColor(double r, double g, double b);
extern void mglSetForegroundColor(double r, double g, double b);

extern void mglSetDrawingArea(double, double, double, double); 

extern void mglSetDrawingArea(const double low[3], const double high[3]); 
extern void mglSetDrawingArea(const double low[3], const double high[3], int); 

extern void mglResetBackgroundColor(double r, double g, double b);
extern void mglResetForegroundColor(double r, double g, double b); 

extern void mglResetDrawingArea(const double low[3], const double high[3], int);
extern void mglResetDrawingArea(const double low[3], const double high[3]); 

extern void mglReset2dDrawingArea(const double low[2], const double high[2]);
extern void mglReset3dDrawingArea(const double low[2], const double high[2]); 
extern void mglResetDrawingArea2d(const double low[2], const double high[2]);
extern void mglResetDrawingArea3d(const double low[3], const double high[3]);

extern void mglResetDrawingArea(double, double, double, double); 
extern void mglResetWindowSize(double width, double height);
extern void mglResetWindowSize(int width, int height); 

extern void mglResetColormap(int type, double cold, double hot);
extern void mglResetColormap(double cold, double hot);

extern void mglEnableColormap(int type); 

extern void mglResetDimension(int dim); 

extern void mglFindBox(const double *x, const double *y, int n, 
                       double low[2], double high[2]);

extern void mglComputeBoundingBox(const double *x, const double *y, int n, 
                                  double low[2], double high[2]); 
extern void mglComputeBoundingBox(double (*f)(double), 
                                  double a, double b, int n,
                                  double low[2], double high[2]); 

extern void mglAddStartingComments(const char *str, int, double);
extern void mglAddEndingComments(const char *str, int, double);
extern void mglAddComments(const char *str, int, double);

extern void mglAddStartingComments(const char *str, double, int); 
extern void mglAddEndingComments(const char *str, double, int); 
extern void mglAddComments(const char *str, double, int); 

extern void mglAddStartingComments(const char *str, double); 
extern void mglAddEndingComments(const char *str, double); 
extern void mglAddComments(const char *str, double); 

extern void mglAddStartingComments(const char *str1, int,
                                   const char *str2, double);
extern void mglAddEndingComments(const char *str1, int, 
                                 const char *str2, double);
extern void mglAddComments(const char *str1, int,
                           const char *str2, double); 

extern void mglAddStartingComments(const char *str1, double,
                                   const char *str2, int); 
extern void mglAddEndingComments(const char *str1, double,
                                 const char *str2, int); 
extern void mglAddComments(const char *str1, double, 
                           const char *str2, int); 

extern void mglAddStartingComments(const char *str, int);
extern void mglAddEndingComments(const char *str, int);
extern void mglAddComments(const char *str, int);

extern void mglAddStartingComments(const char *str); 
extern void mglAddEndingComments(const char *str);
extern void mglAddComments(const char *str); 

// MathGL Facility Functions.

extern void mglWait(void);
extern void mglDump(void); 
extern void mglPause(void);
extern void mglFlush(void);
extern void mglClear(void);
extern void mglSleep(float);
extern void mglFlush(double);

extern void mglNewpage(void);
extern void mglNewPage(void);
extern void mglNewpage(double); 
extern void mglNewPage(double); 

extern void mglNewpage(const char *);
extern void mglNewPage(const char *); 

extern void mglSetTimer(double, double);
extern void mglSetSlider(double, double); 

extern void mglResetTimer(double, double); 
extern void mglResetSlider(double, double); 

extern void mglRegisterTimer(double, double); 
extern void mglRegisterSlider(double, double);

extern void mglReadWindowCoord(double &, double &);

extern void mglInput(int &type, double &x, double &y, double &z); 
extern void mglInput(double &x, double &y, double &z); 
extern void mglInput(int &type, double &x, double &y);
extern void mglInput(double &x, double &y); 

extern int  mglReadGuiInput(const char message[80]); 

extern int  mglSetWindowSignature(const char *);

extern void mglMakeCurrentWindow(int id = 0);
extern void mglMakeCurrentWindow(const char *);

extern int  mglWindowNumber(void);

extern int  mglNewWindow(const char*); 
extern int  mglNewWindow(void); 

// MathGL mglPlot Setting.

extern void mglColor(int);
extern void mglColor(char);
extern void mglColor(double);
extern void mglClearColor(int);
extern void mglClearColor(char); 
extern void mglRandomColor(void); 
extern void mglColor(int n, int i);
extern void mglColor(double, double); 

extern double mglLookupColormap(double v); 

extern void mglMaterialColor(double);
extern void mglMaterialColor(double, double);
extern void mglColor(double red, double green, 
                     double blue, double alpha);
extern void mglClearColor(double red, double green,
                          double blue, double alpha); 
extern void mglDrawingColor(double red, double green, 
                            double blue, double alpha);
extern void mglMaterialColor(double red, double green, 
                             double blue, double alpha);

extern void mglColor(double red, double green, double blue); 
extern void mglClearColor(double red, double green, double blue);
extern void mglDrawingColor(double red, double green, double blue);
extern void mglMaterialColor(double red, double green, double blue);

extern void mglBackMaterialColor(double, double);
extern void mglBackMaterialColor(double red, double green,
                                 double blue, double alpha); 
extern void mglBackMaterialColor(double red, double green, double blue); 

extern void mglFrontMaterialColor(double); 
extern void mglFrontMaterialColor(double, double); 
extern void mglFrontMaterialColor(double red, double green,
                                  double blue, double alpha);
extern void mglFrontMaterialColor(double red, double green, double blue);

extern void mglBackMaterialColor(double); 
extern void mglAllocateColor(double low, double high, double value,
                             double &r, double &g, double &b);
extern void mglAllocateColor(double c, double &r, double &g, double &b); 

extern void mglTranslate(double u, double v, double w);
extern void mglNormal(double u, double v, double w);
extern void mglTranslate(double u, double v);
extern void mglRotate(double angle, double x,
                      double y, double z); 

extern void mglPlotSingleColorOpaqueBoxFrame(void);
extern void mglPlotSingleColorBoxFrame(void);
extern void mglPlotGrayOpaqueBoxFrame(void); 
extern void mglPlotGrayBoxFrame(void);
extern void mglPlotOpaqueBoxFrame(void); 
extern void mglPlotBoundingBox(void);
extern void mglPlotBoxFrame(void); 

// extern void mglTranslate(const SpaceCoord&);
// extern void mglNormal(const SpaceCoord&);

extern void mglLineWidth(double wd);

// MathGL mglPlotting.

extern void mglPlot2dLine(double x1, double y1, 
                          double x2, double y2);

extern void mglPlot2dDashLine(double x1, double y1, 
                              double x2, double y2); 
extern void mglPlot2dDashLine(double x1, double y1, 
                              double x2, double y2, 
                              double width); 

extern void mglPlot2dLine(const double a[2], const double b[2]); 

extern void mglPlotLine(const double a[], const double b[], int n);

extern void mglPlot2dDashLine(const double a[2], const double b[2]);
extern void mglPlot2dDashLine(const double a[2], const double b[2], 
                              double width);

extern void mglPlot2dArrowLine(double x1, double y1, 
                               double x2, double y2);
extern void mglPlot2dArrowHead(double x1, double y1, 
                               double x2, double y2); 

extern void mglPlot2dArrowLine(const double a[2], const double b[2]);
extern void mglPlot2dArrowHead(const double a[2], const double b[2]);

extern void mglPlot2dArrowLine(const double a[2], const double b[2], 
                               double width);
extern void mglPlot2dArrowHead(const double a[2], const double b[2],
                               double size); 
extern void mglPlot2dArrowLine2(const double a[2], const double b[2],
                                double arrow_size); 

extern void mglPlot2dArrowLine(double x1, double y1,
                               double x2, double y2, 
                               double width); 
extern void mglPlot2dArrowHead(double x1, double y1, 
                               double x2, double y2, 
                               double size); 
extern void mglPlot2dArrowLine2(double x1, double y1,
                               double x2, double y2, 
                               double arrow_size);

extern void mglPlot2dPoint(double, double);
extern void mglPlot2dLines(const double[], int);
extern void mglPlot2dPoints(const double[], int);
extern void mglPlot2dLines(const double[], const double[], int); 
extern void mglPlot2dCurve(const double[], const double[], int);
extern void mglPlot2dPoints(const double[], const double[], int);
extern void mglPlot2dLines(double a, double b, const double[], int); 
extern void mglPlot2dPoints(double a, double b, const double[], int);
extern void mglPlot2dClosedCurve(const double[], const double[], int); 
extern void mglPlot2dClosedCurve(const double[], const double[],
                                 int, double width);
extern void mglPlot2dCurve(const double[], const double[], 
                           int, double width);

extern void mglPlot2dLine(double x1, double y1,
                          double x2, double y2, double width);
extern void mglPlot2dLine(const double a[2], const double b[2],
                          double width); 
extern void mglPlot2dLines(const double[], int, double width);
extern void mglPlot2dLines(double a, double b, const double[], 
                           int, double width); 
extern void mglPlot2dLines(const double[], const double[], 
                           int, double width); 

extern void mglPlot2dLines(const double (*)[2], int n);
extern void mglPlot2dCurve(const double (*)[2], int n);
extern void mglPlot2dPoints(const double (*)[2], int n);
extern void mglPlot2dClosedCurve(const double (*)[2], int n); 

extern void mglPlot2dLines(const double (*)[2], int n, double width);
extern void mglPlot2dCurve(const double (*)[2], int n, double width); 
extern void mglPlot2dPoints(const double (*)[2], int n, double width); 
extern void mglPlot2dClosedCurve(const double (*)[2], int n, double width);

extern void mglPlot2dDashLine(double x1, double y1, 
                              double x2, double y2, double width); 
extern void mglPlot2dDashLine(const double a[2], const double b[2], 
                              double width);

extern void mglPlot2dDashLines(const double[], 
                               const double[],
                               int, double width);
extern void mglPlot2dDashLines(const double[],
                               const double[], int);
extern void mglPlot2dDashLines(double a, double b,
                               const double[], int); 
extern void mglPlot2dDashLines(double a, double b, 
                               const double[], int,
                               double width);
extern void mglPlot2dDashLines(const double (*)[2], int n, 
                               double width); 
extern void mglPlot2dDashLines(const double (*)[2], int n);

extern void mglPlot2dLineWithinBoundingBox(double u1, double v1,
                                           double u2, double v2, 
                                           const double low[2], 
                                           const double high[2]);
extern void mglPlot2dLineWithinBoundingBox(double u1, double v1, 
                                           double u2, double v2,
                                           double lx, double ly,
                                           double hx, double hy); 

extern void mglPlot2dDashLineWithinBoundingBox(double u1, double v1, 
                                               double u2, double v2,
                                               const double low[2],
                                               const double high[2]);
extern void mglPlot2dDashLineWithinBoundingBox(double u1, double v1, 
                                               double u2, double v2, 
                                               double lx, double ly,
                                               double hx, double hy);

extern void mglPlot2dArrowArc(const double O[2], const double A[2],
                              double theta);
extern void mglPlot2dArrowArc(const double O[2], const double A[2], 
                              const double B[2]); 
extern void mglPlot2dArrowArc(double Ox, double Oy, double Ax,
                              double Ay, double Bx, double By);
extern void mglPlot2dArrowArc(double Ox, double Oy, double Ax,
                              double Ay, double theta);

extern void mglPlot2dArrowArc(double Ox, double Oy, double Ax, 
                              double Ay, double Bx, double By,
                              double width);
extern void mglPlot2dArrowArc2(double Ox, double Oy, double Ax,
                              double Ay, double theta,
                              double width); 
extern void mglPlot2dArrowArc(const double O[2], const double A[2], 
                              double theta, double width); 
extern void mglPlot2dArrowArc(const double O[2], const double A[2],
                              const double B[2], double width); 

extern void mglPlot2dArrowArc3(double Ox, double Oy, double Ax,
                               double Ay, double Bx, double By,
                               double size);
extern void mglPlot2dArrowArc3(double Ox, double Oy, double Ax,
                               double Ay, double theta,
                               double size); 
extern void mglPlot2dArrowArc3(const double O[2], const double A[2], 
                               double theta, double size); 
extern void mglPlot2dArrowArc3(const double O[2], const double A[2], 
                               const double B[2], double size); 

extern void mglPlot2dArc(const double O[2], const double A[2], double theta);
extern void mglPlot2dArc(const double O[2], const double A[2],
                         const double B[2]); 
extern void mglPlot2dArc(double Ox, double Oy, double Ax, 
                         double Ay, double Bx, double By);
extern void mglPlot2dArc(double Ox, double Oy, double Ax, 
                         double Ay, double theta); 
extern void mglPlot2dArc(double Ox, double Oy, double Ax, 
                         double Ay, double Bx, double By, 
                         double width);
extern void mglPlot2dArc2(double Ox, double Oy, double Ax, 
                         double Ay, double theta,
                         double width);
extern void mglPlot2dArc(double O[2], double A[2], 
                         double theta, double width); 
extern void mglPlot2dArc(double O[2], double A[2], 
                         double B[2], double width);

extern void mglPlotCircle(double x0 = 0, double y0 = 0,
                          double r = 1.0); 
extern void mglPlotCircle(double x0, double y0, double r,
                          double width);

extern void mglPlotSolidCircle(double x0 = 0, double y0 = 0,
                               double r = 1.0);

extern void mglPlotEllipse(double x0 = 0, double y0 = 0, 
                           double ra = 1.0, double rb = 1.0);
extern void mglPlotEllipse(double x0, double y0, double ra,
                           double rb, double width);
extern void mglPlotEllipse(double x0, double y0, double ra,
                           double rb, double theta, double width);

extern void mglPlotSolidEllipse(double x0, double y0, 
                                double ra, double rb); 
extern void mglPlotSolidEllipse(double x0, double y0, 
                                double ra, double rb, double alpha);

extern void mglPlotEllipseArc(double cx, double cy, double rx,
                              double rb, double theta1, double theta2); 
extern void mglPlotEllipseArc(double cx, double cy, double rx, 
                              double ry, double x1, double y1, 
                              double delta); 
extern void mglPlotEllipseArc(double cx, double cy, double rx, 
                              double ry, double x1, double y1,
                              double x2, double y2); 

extern void mglPlotEllipseArc2(double cx, double cy, double rx, 
                               double ry, double x1, double y1, 
                               double delta, double width); 
extern void mglPlotEllipseArc2(double cx, double cy, double rx,
                               double ry, double theta1, double theta2, 
                               double width);
extern void mglPlotEllipseArc2(double cx, double cy, double rx,
                               double ry, double x1, double y1, 
                               double x2, double y2,
                               double width); 

extern void mglPlotDashCircle(double x0 = 0, double y0 = 0,
                              double r = 1.0);
extern void mglPlotDashCircle(double x0, double y0, double r,
                              double width);

extern void mglPlotDashEllipse(double x0 = 0, double y0 = 0, 
                               double ra = 1.0, double rb = 1.0); 
extern void mglPlotDashEllipse(double x0, double y0, double ra,
                               double rb, double width);
extern void mglPlotDashEllipse(double x0, double y0, double ra, 
                               double rb, double theta, double width);

extern void mglPlotDashEllipseArc(double x0, double y0, double ra, 
                                  double rb, double x1, double y1, 
                                  double delta); 
extern void mglPlotDashEllipseArc(double x0, double y0, double ra,
                                  double rb, double theta1, double theta2); 
extern void mglPlotDashEllipseArc(double cx0, double cy0, double ra, 
                                  double rb, double x1, double y1, 
                                  double x2, double y2);

extern void mglPlotDashEllipseArc2(double x0, double y0, double ra,
                                   double rb, double x1, double y1, 
                                   double delta, double width); 
extern void mglPlotDashEllipseArc2(double x0, double y0, double ra, 
                                   double rb, double theta1, double theta2, 
                                   double width); 
extern void mglPlotDashEllipseArc2(double cx0, double cy0, double ra,
                                   double rb, double x1, double y1,
                                   double x2, double y2, 
                                   double width);

extern void mglPlot2dTriangle(const double coord[3][2], double width); 
extern void mglPlot2dTriangle(const double coord[6], double width); 
extern void mglPlot2dTriangle(const double coord[3][2]); 
extern void mglPlot2dTriangle(const double coord[6]); 

extern void mglPlot2dTriangle(double, double, double, double,
                              double, double); 
extern void mglPlot2dTriangle(double, double, double, double,
                              double, double, double width);

extern void mglPlot2dSolidTriangle(double, double, double,
                                   double, double, double); 

extern void mglPlot2dSolidTriangle(const double coord[3][2]); 
extern void mglPlot2dSolidTriangle(const double coord[6]); 

extern void mglPlot2dTriangle(const double a[2], const double b[2],
                              const double c[2]); 
extern void mglPlot2dTriangle(const double a[2], const double b[2],
                              const double c[2], double width); 
extern void mglPlot2dSolidTriangle(const double a[2], const double b[2], 
                                   const double c[2]);

extern void mglPlot2dRectangle(double, double, double, double); 
extern void mglPlot2dRectangle(double, double, double, double, double wd); 
extern void mglPlot2dSolidRectangle(double, double, double, double); 

extern void mglPlot2dRectangle(double, double, double, double, 
                               double, double, double, double);
extern void mglPlot2dRectangle(double, double, double, double,
                               double, double, double, double, double wd);

extern void mglPlot2dRectangle(const double lo[2], const double hi[2], double); 
extern void mglPlot2dRectangle(const double lo[2], const double hi[2]); 

extern void mglPlot2dRectangle2(const double c[2], const double r[2], double); 
extern void mglPlot2dRectangle2(const double c[2], const double r[2]); 

extern void mglPlot2dRectangle(double lx, double ly, double hx, double hy, 
                               const double data[4]); 
extern void mglPlot2dRectangle(const double lo[2], const double hi[2],
                               const double data[4]); 

extern void mglPlot2dDashRectangle(double lx, double ly,
                                   double hx, double hy,
                                   double line_width);
extern void mglPlot2dDashRectangle(double lx, double ly,
                                   double hx, double hy);

extern void mglPlot2dDashRectangle(const double low[2], 
                                   const double high[2], 
                                   double line_width);
extern void mglPlot2dDashRectangle(const double low[2], 
                                   const double high[2]); 

extern void mglPlot2dDashRectangle2(const double center[2],
                                    const double radius[2], 
                                    double line_width); 
extern void mglPlot2dDashRectangle2(const double center[2],
                                    const double radius[2]); 

extern void mglPlot2dQuadrilateral(const double coord[8]);
extern void mglPlot2dQuadrilateral(const double coord[4][2]); 

extern void mglPlot2dQuadrilateral(const double coord[8], double wd); 
extern void mglPlot2dQuadrilateral(const double coord[4][2], double wd);

extern void mglPlot2dQuadrilateral(double, double, double, double,
                                   double, double, double, double);
extern void mglPlot2dQuadrilateral(double, double, double, double,
                                   double, double, double, double, double wd);

extern void mglPlot2dQuadrilateral(const double x[4], const double y[4]); 
extern void mglPlot2dQuadrilateral(const double x[4], const double y[4],
                                   double width); 

extern void mglPlot2dQuadrilateral(const double a[2], const double b[2],
                                   const double c[2], const double d[2]); 
extern void mglPlot2dQuadrilateral(const double a[2], const double b[2], 
                                   const double c[2], const double d[2], 
                                   double width);

extern void mglPlot2dSolidRectangle(double, double, double, double,
                                    double, double, double, double); 
extern void mglPlot2dSolidRectangle(const double x[2], const double y[2]);

extern void mglPlot2dSolidQuadrilateral(double, double, double, double,
                                        double, double, double, double); 
extern void mglPlot2dSolidQuadrilateral(const double x[4], const double y[4]);
extern void mglPlot2dSolidQuadrilateral(const double a[2], const double b[2],
                                        const double c[2], const double d[2]); 

extern void mglPlot2dSolidPolygon(const double x[], const double y[], int n);
extern void mglPlot2dPolygon(const double x[], const double y[], int n);
extern void mglPlot2dPolygon(const double x[], const double y[], int n, 
                             double width); 

extern void mglPlotShadedRectangle(double lx, double ly,
                                   double hx, double hy); 
extern void mglPlotShadedRectangle(double lx, double ly,
                                   double hx, double hy, int n);

extern void mglPlotShadedRectangle(const double low[2], 
                                   const double high[2]); 
extern void mglPlotShadedRectangle(const double low[2],
                                   const double high[2], int n);

extern void mglPlot2dLatticeData(const double *x, int m,
                                 const double *y, int n,
                                 double **const u);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglScaleObject(const double s[3]);
extern void mglTranslateObject(const double d[3]);
extern void mglScaleObject(double sx, double sy, double sz); 
extern void mglTranslateObject(double dx, double dy, double dz); 

extern void mglRotateObject(const double dir[3], double theta); 
extern void mglRotateObject(double dx, double dy, double dz, double theta);

extern void mglRotateObject(double cx, double cy, double cz,
                            double dx, double dy, double dz, 
                            double theta);
extern void mglRotateObject(const double center[3],
                            const double dir[3], 
                            double theta); 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPlot3d(void);
extern void mglPlot3dLines(void);

extern void mglPlot3dPoint(const double[3]);
extern void mglPlot3dLine(const double[3], const double[3]); 
extern void mglPlot3dPoint(double, double, double);
extern void mglPlot3dDashLine(const double[3], const double[3]);
extern void mglPlot3dPoint(const double [3], double size);
extern void mglPlot3dArrowLine(const double a[3], const double b[3]);
extern void mglPlot3dLine(const double[3], const double[3], double lw);
extern void mglPlot3dLines(const double[], const double[],
                           const double[], int);
extern void mglPlot3dPoint(double x, double y, double z, double size); 
extern void mglPlot3dLine(const double a[3], const double b[3], 
                          const double normal[3]);
extern void mglPlot3dLine2(const double a[3], const double b[3], 
                           const double normal[3]);
extern void mglPlot3dLine3(const double a[3], const double b[3],
                           const double normal[3]); 
extern void mglPlot3dDashLine(const double[3], const double[3], double lw); 
extern void mglPlot3dLines(const double[], const double[], const double[],
                           int, double lw); 

extern void mglPlot3dLine(double, double, double, double, double, double); 
extern void mglPlot3dDashLine(double, double, double, double, double, double); 
extern void mglPlot3dArrowLine(double, double, double, double, double, double); 

extern void mglPlot3dLine(double ax, double ay, double az,
                          double bx, double by, double bz, 
                          double nx, double ny, double nz);
extern void mglPlot3dLine2(double ax, double ay, double az,
                           double bx, double by, double bz,
                           double nx, double ny, double nz);
extern void mglPlot3dLine3(double ax, double ay, double az,
                           double bx, double by, double bz,
                           double nx, double ny, double nz); 

extern void mglPlot3dDashLine(double, double, double,
                              double, double, double, 
                              double lw); 
extern void mglPlot3dLine(double, double, double, 
                          double, double, double, 
                          double lw); 

extern void mglPlotCube(const double low[3],
                        const double high[3]);
extern void mglPlotCube2(const double center[3], 
                         const double radius[3]); 
extern void mglPlotCube(const double ctr[3], double r);
extern void mglPlotCube(double lx, double ly, double lz,
                        double hx, double hy, double hz);
extern void mglPlotCube(double cx, double cy, double cz,
                        double radius); 
extern void mglPlotCube(double radius);

extern void mglPlot3dCube(const double low[3], const double high[3]);
extern void mglPlot3dCube2(const double ctr[3], const double radius[3]); 
extern void mglPlot3dWireCube(const double low[3], const double high[3]);
extern void mglPlot3dWireCube2(const double ctr[3], const double radius[3]);

extern void mglPlotWireCube(const double ctr[3], double r);
extern void mglPlotWireCube(const double lo[3], const double hi[3]);
extern void mglPlotWireCube2(const double ctr[3], const double r[3]); 
extern void mglPlotWireCube(double lx, double ly, double lz,
                            double hx, double hy, double hz); 
extern void mglPlotWireCube(double cx, double cy, double cz, 
                            double radius); 
extern void mglPlotWireCube(double radius); 

extern void mglPlotOpaqueWireCube(const double ctr[3], double r); 
extern void mglPlotOpaqueWireCube(const double lo[3], const double hi[3]);
extern void mglPlotOpaqueWireCube(double lx, double ly, double lz,
                                  double hx, double hy, double hz);
extern void mglPlotOpaqueWireCube(double cx, double cy, double cz,
                                  double radius);
extern void mglPlotOpaqueWireCube(double radius); 

extern void mglPlotSolidCube(const double ctr[3], double r); 
extern void mglPlotSolidCube(const double lo[3], const double hi[3]);
extern void mglPlotSolidCube2(const double ctr[3], const double r[3]);
extern void mglPlotSolidCube(double lx, double ly, double lz, 
                             double hx, double hy, double hz);
extern void mglPlotSolidCube(double cx, double cy, double cz,
                             double radius);
extern void mglPlotSolidCube(double radius);

extern void mglPlot3dCurve(const double[], const double[],
                           const double[], int);
extern void mglPlot3dPoints(const double[], const double[], 
                            const double[], int);

extern void mglPlot3dLines(const double x[], const double y[], 
                           const double z[], int); 
extern void mglPlot3dLines(double t, const double x[],
                           const double y[], int); 
extern void mglPlot3dLines(const double x[], const double y[],
                           double t, int); 
extern void mglPlot3dLines(const double x[], double t, 
                           const double y[], int); 
extern void mglPlot3dLines(const double x[], const double y[],
                           const double z[], int, double wd);

extern void mglPlot3dLines(const double x[][3], int n, double wd); 
extern void mglPlot3dLines(const double x[][3], int n);

extern void mglPlot3dPoints(const double x[][3], int n, double wd); 
extern void mglPlot3dPoints(const double x[][3], int n); 

extern void mglPlot3dWireTriangle(const double a[3], const double b[3], 
                                  const double c[3]); 

extern void mglPlot3dWireTriangle(double ax, double ay, double az,
                                  double bx, double by, double bz,
                                  double cx, double cy, double cz);

extern void mglPlot3dWireRectangle(const double x[4], const double y[4],
                                   const double z[4]); 

extern void mglPlot3dWireRectangle(double ax, double ay, double az,
                                   double bx, double by, double bz, 
                                   double cx, double cy, double cz,
                                   double dx, double dy, double dz); 

extern void mglPlot3dDashRectangle(double ax, double ay, double az, 
                                   double bx, double by, double bz, 
                                   double cx, double cy, double cz,
                                   double dx, double dy, double dz); 

extern void mglPlot3dDashQuadrilateral(double ax, double ay, double az, 
                                       double bx, double by, double bz,
                                       double cx, double cy, double cz, 
                                       double dx, double dy, double dz); 

extern void mglPlot3dWireQuadrilateral(double ax, double ay, double az,
                                       double bx, double by, double bz, 
                                       double cx, double cy, double cz,
                                       double dx, double dy, double dz); 

extern void mglPlot3dWireQuadrilateral(const double a[3], const double b[3],
                                       const double c[3], const double d[3]);
extern void mglPlot3dDashQuadrilateral(const double a[3], const double b[3],
                                       const double c[3], const double d[3]);

extern void mglPlot3dWireRectangle(const double a[3], const double b[3],
                                   const double c[3], const double d[3]);
extern void mglPlot3dDashRectangle(const double a[3], const double b[3],
                                   const double c[3], const double d[3]); 
extern void mglPlotWireTetrahedron(const double a[3], const double b[3],
                                   const double c[3], const double d[3]); 
extern void mglPlotDashTetrahedron(const double a[3], const double b[3], 
                                   const double c[3], const double d[3]); 

extern void mglPlot3dRectangle(const double x[4], const double y[4], 
                               const double z[4]); 

extern void mglPlot3dRectangle(const double low[3], const double high[3]);

extern void mglPlot3dTriangle(const double a[3], const double b[3], 
                              const double c[3]); 

extern void mglPlot3dTriangle(double ax, double ay, double az,
                              double bx, double by, double bz,
                              double cx, double cy, double cz); 

extern void mglPlot3dRectangle(double ax, double ay, double az,
                               double bx, double by, double bz, 
                               double cx, double cy, double cz, 
                               double dx, double dy, double dz); 

extern void mglPlot3dQuadrilateral(double ax, double ay, double az, 
                                   double bx, double by, double bz,
                                   double cx, double cy, double cz,
                                   double dx, double dy, double dz);

extern void mglPlot3dRectangle(const double a[3], const double b[3], 
                               const double c[3], const double d[3]);

extern void mglPlot3dQuadrilateral(const double a[3], const double b[3],
                                   const double c[3], const double d[3]);

extern void mglPlot3dRectangle(const double low[2], const double high[2], 
                               double value); 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPlot3dTriangles(const double (*a)[3], const double (*b)[3],
                               const double (*c)[3], int n); 

extern void mglPlot3dWireTriangles(const double (*a)[3], const double (*b)[3],
                                   const double (*c)[3], int n);

extern void mglPlot3dWireTriangles(const double (*a)[3], const double (*b)[3],
                                   const double (*c)[3], int n, double width);

extern void mglPlot3dRectangles(const double (*a)[3], const double (*b)[3],
                                const double (*c)[3], const double (*d)[3],
                                int n); 

extern void mglPlot3dWireRectangles(const double (*a)[3], const double (*b)[3], 
                                    const double (*c)[3], const double (*d)[3], 
                                    int n); 

extern void mglPlot3dwireRectangles(const double (*a)[3], const double (*b)[3],
                                    const double (*c)[3], const double (*d)[3], 
                                    int n, double width);

extern void mglPlot3dQuadrilaterals(const double (*a)[3], const double (*b)[3], 
                                    const double (*c)[3], const double (*d)[3], 
                                    int n); 

extern void mglPlot3dWireQuadrilaterals(const double (*a)[3],
                                        const double (*b)[3],
                                        const double (*c)[3], 
                                        const double (*d)[3],
                                        int n); 

extern void mglPlot3dWireQuadrilaterals(const double (*a)[3],
                                        const double (*b)[3],
                                        const double (*c)[3], 
                                        const double (*d)[3],
                                        int n, double width);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPlotDataOnTetrahedron(const double a[3], const double b[3],
                                     const double c[3], const double d[3], 
                                     const double data[4]);
extern void mglPlotDataOnTetrahedron2(const double a[3], const double b[3],
                                      const double c[3], const double d[3],
                                      const double data[4]); 

extern void mglPlotDataOnQuadrilateral(const double a[3], const double b[3], 
                                       const double c[3], const double d[3], 
                                       const double data[4]); 
// One side is colormap and another side shows isocontours.
extern void mglPlotDataOnQuadrilateral2(const double a[3], const double b[3], 
                                        const double c[3], const double d[3],
                                        const double data[4]);

extern void mglPlotDataOnRectangle(const double a[3], const double b[3], 
                                   const double c[3], const double d[3],
                                   const double data[4]); 
// One side is colormap and another side shows isocontours.
extern void mglPlotDataOnRectangle2(const double a[3], const double b[3], 
                                    const double c[3], const double d[3], 
                                    const double data[4]); 

extern void mglPlotIsolinesOn3dRectangle(const double a[3], const double b[3], 
                                         const double c[3], const double d[3], 
                                         const double data[4]); 

extern void mglPlot3dWireRectangle(const double a[3], const double b[3], 
                                   const double c[3], const double d[3], 
                                   const double data[4]);

extern void mglPlotDataOnTriangle(const double a[3], const double b[3], 
                                  const double c[3], const double data[3]); 
// One side is colormap and another side shows isocontours.
extern void mglPlotDataOnTriangle2(const double a[3], const double b[3], 
                                   const double c[3], const double data[3]);

extern void mglPlotIsolinesOn3dTriangle(const double a[3], const double b[3], 
                                        const double c[3], const double v[3]); 

extern void mglPlotDataOnTriangles(const double (*a)[3], const double (*b)[3],
                                   const double (*c)[3], const double (*u)[3],
                                   int n); 

extern void mglPlotDataOnRectangles(const double (*a)[3], const double (*b)[3],
                                    const double (*c)[3], const double (*d)[3],
                                    const double (*u)[4], int n);

extern void mglPlotDataOnQuadrilaterals(const double (*a)[3], 
                                        const double (*b)[3],
                                        const double (*c)[3], 
                                        const double (*d)[3],
                                        const double (*u)[4], 
                                        int n); 

extern void mglPlotTetrahedron(const double a[3], const double b[3],
                               const double c[3], const double d[3]);

extern void mglPlotSolidHexahedron(const double coord[8][3]); 
extern void mglPlotHexahedron(const double coord[8][3]); 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPlotStarMarkers(const double x[], 
                               const double y[], 
                               int n, double size);
extern void mglPlotPlusMarkers(const double x[], 
                               const double y[],
                               int n, double size); 
extern void mglPlotRectMarkers(const double x[], 
                               const double y[], 
                               int n, double size);
extern void mglPlotCrossMarkers(const double x[], 
                                const double y[],
                                int n, double size); 

extern void mglPlotCircleMarkers(const double x[],
                                 const double y[], 
                                 int n, double size); 
extern void mglPlotDiamondMarkers(const double x[], 
                                  const double y[],
                                  int n, double size);
extern void mglPlotTriangleMarkers(const double x[], 
                                   const double y[], 
                                   int n, double size);
extern void mglPlotRectangleMarkers(const double x[], 
                                    const double y[],
                                    int n, double size);

extern void mglPlotSolidCircleMarkers(const double x[], 
                                      const double y[],
                                      int n, double size);
extern void mglPlotSolidDiamondMarkers(const double x[],
                                       const double y[], 
                                       int n, double size);
extern void mglPlotSolidTriangleMarkers(const double x[], 
                                        const double y[], 
                                        int n, double size); 
extern void mglPlotSolidRectangleMarkers(const double x[], 
                                         const double y[], 
                                         int n, double size);

extern void mglPlotStarMarkers(const double x[], 
                               const double y[], int n); 
extern void mglPlotPlusMarkers(const double x[], 
                               const double y[], int n);
extern void mglPlotRectMarkers(const double x[],
                               const double y[], int n);
extern void mglPlotCrossMarkers(const double x[], 
                                const double y[], int n); 

extern void mglPlotCircleMarkers(const double x[], 
                                 const double y[], int n);
extern void mglPlotDiamondMarkers(const double x[], 
                                  const double y[], int n); 
extern void mglPlotTriangleMarkers(const double x[],
                                   const double y[], int n);
extern void mglPlotRectangleMarkers(const double x[], 
                                    const double y[], int n);

extern void mglPlotSolidCircleMarkers(const double x[],
                                      const double y[], int n);
extern void mglPlotSolidDiamondMarkers(const double x[], 
                                       const double y[], int n);
extern void mglPlotSolidTriangleMarkers(const double x[], 
                                        const double y[], int n); 
extern void mglPlotSolidRectangleMarkers(const double x[],
                                         const double y[], int n); 

extern void mglPlotStarMarkers(const double (*x)[2], int n); 
extern void mglPlotPlusMarkers(const double (*x)[2], int n);
extern void mglPlotRectMarkers(const double (*x)[2], int n); 
extern void mglPlotCrossMarkers(const double (*x)[2], int n);
extern void mglPlotCircleMarkers(const double (*x)[2], int n); 
extern void mglPlotDiamondMarkers(const double (*x)[2], int n);
extern void mglPlotTriangleMarkers(const double (*x)[2], int n);
extern void mglPlotRectangleMarkers(const double (*x)[2], int n); 
extern void mglPlotSolidCircleMarkers(const double (*x)[2], int n); 
extern void mglPlotSolidDiamondMarkers(const double (*x)[2], int n);
extern void mglPlotSolidTriangleMarkers(const double (*x)[2], int n);
extern void mglPlotSolidRectangleMarkers(const double (*x)[2], int n);

extern void mglPlotStarMarkers(const double (*x)[2], int n, double);
extern void mglPlotPlusMarkers(const double (*x)[2], int n, double);
extern void mglPlotRectMarkers(const double (*x)[2], int n, double);
extern void mglPlotCrossMarkers(const double (*x)[2], int n, double); 
extern void mglPlotCircleMarkers(const double (*x)[2], int n, double);
extern void mglPlotDiamondMarkers(const double (*x)[2], int n, double);
extern void mglPlotTriangleMarkers(const double (*x)[2], int n, double);
extern void mglPlotRectangleMarkers(const double (*x)[2], int n, double); 
extern void mglPlotSolidCircleMarkers(const double (*x)[2], int n, double); 
extern void mglPlotSolidDiamondMarkers(const double (*x)[2], int n, double);
extern void mglPlotSolidTriangleMarkers(const double (*x)[2], int n, double); 
extern void mglPlotSolidRectangleMarkers(const double (*x)[2], int n, double);

extern void mglPlotStarMarker(const double x[2]); 
extern void mglPlotPlusMarker(const double x[2]); 
extern void mglPlotRectMarker(const double x[2]); 
extern void mglPlotCrossMarker(const double x[2]);

extern void mglPlotCircleMarker(const double x[2]);
extern void mglPlotDiamondMarker(const double x[2]);
extern void mglPlotTriangleMarker(const double x[2]);
extern void mglPlotRectangleMarker(const double x[2]); 

extern void mglPlotSolidCircleMarker(const double x[2]); 
extern void mglPlotSolidDiamondMarker(const double x[2]); 
extern void mglPlotSolidTriangleMarker(const double x[2]);
extern void mglPlotSolidRectangleMarker(const double x[2]);

extern void mglPlotStarMarker(double x, double y);
extern void mglPlotPlusMarker(double x, double y);
extern void mglPlotRectMarker(double x, double y);
extern void mglPlotCrossMarker(double x, double y);

extern void mglPlotCircleMarker(double x, double y); 
extern void mglPlotDiamondMarker(double x, double y); 
extern void mglPlotTriangleMarker(double x, double y); 
extern void mglPlotRectangleMarker(double x, double y); 

extern void mglPlotSolidCircleMarker(double x, double y);
extern void mglPlotSolidDiamondMarker(double x, double y);
extern void mglPlotSolidTriangleMarker(double x, double y);
extern void mglPlotSolidRectangleMarker(double x, double y);

extern void mglPlotStarMarker(const double x[2], double size); 
extern void mglPlotPlusMarker(const double x[2], double size); 
extern void mglPlotRectMarker(const double x[2], double size); 
extern void mglPlotCrossMarker(const double x[2], double size);

extern void mglPlotCircleMarker(const double x[2], double size); 
extern void mglPlotDiamondMarker(const double x[2], double size);
extern void mglPlotTriangleMarker(const double x[2], double size); 
extern void mglPlotRectangleMarker(const double x[2], double size); 

extern void mglPlotSolidCircleMarker(const double x[2], double size);
extern void mglPlotSolidDiamondMarker(const double x[2], double size);
extern void mglPlotSolidTriangleMarker(const double x[2], double size);
extern void mglPlotSolidRectangleMarker(const double x[2], double size);

extern void mglPlotStarMarker(double x, double y, double size); 
extern void mglPlotPlusMarker(double x, double y, double size);
extern void mglPlotRectMarker(double x, double y, double size); 
extern void mglPlotCrossMarker(double x, double y, double size);

extern void mglPlotCircleMarker(double x, double y, double size);
extern void mglPlotDiamondMarker(double x, double y, double size); 
extern void mglPlotTriangleMarker(double x, double y, double size);
extern void mglPlotRectangleMarker(double x, double y, double size); 

extern void mglPlotSolidCircleMarker(double x, double y, double size);
extern void mglPlotSolidDiamondMarker(double x, double y, double size);
extern void mglPlotSolidTriangleMarker(double x, double y, double size); 
extern void mglPlotSolidRectangleMarker(double x, double y, double size);

extern void mglPlot3dSolidCubeMarker(double, double, double, double size = 1); 
extern void mglPlot3dWireCubeMarker(double, double, double, double size = 1); 
extern void mglPlot3dSphereMarker(double, double, double, double size = 1); 
extern void mglPlot3dCrossMarker(double, double, double, double size = 1); 
extern void mglPlot3dCubeMarker(double, double, double, double size = 1); 
extern void mglPlot3dBallMarker(double, double, double, double size = 1);

extern void mglPlot3dSolidCubeMarker(const double a[3], double size = 1); 
extern void mglPlot3dWireCubeMarker(const double a[3], double size = 1); 
extern void mglPlot3dSphereMarker(const double a[3], double size = 1); 
extern void mglPlot3dCrossMarker(const double a[3], double size = 1); 
extern void mglPlot3dCubeMarker(const double a[3], double size = 1); 
extern void mglPlot3dBallMarker(const double a[3], double size = 1);

extern void mglPlotSolidCubeMarker(double, double, double, double size = 1); 
extern void mglPlotWireCubeMarker(double, double, double, double size = 1);
extern void mglPlotSphereMarker(double, double, double, double size = 1); 
extern void mglPlotCubeMarker(double, double, double, double size = 1); 

extern void mglPlotSolidCubeMarker(const double a[3], double size = 1); 
extern void mglPlotWireCubeMarker(const double a[3], double size = 1);
extern void mglPlotSphereMarker(const double a[3], double size = 1); 
extern void mglPlotCubeMarker(const double a[3], double size = 1);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPrintLeftAlignedWindowString(int x, int y, const char []); 
extern void mglPrintRightAlignedWindowString(int x, int y, const char []); 
extern void mglPrintCenteredWindowString(int x, int y, const char []); 

extern void mglPlotLeftAlignedWindowString(int x, int y, const char []); 
extern void mglPlotRightAlignedWindowString(int x, int y, const char []); 
extern void mglPlotCenteredWindowString(int x, int y, const char []);

extern void mglPrintLeftAlignedWindowText(int x, int y, const char []);
extern void mglPrintRightAlignedWindowText(int x, int y, const char []); 
extern void mglPrintCenteredWindowText(int x, int y, const char []);

extern void mglPlotLeftAlignedWindowText(int x, int y, const char []);
extern void mglPlotRightAlignedWindowText(int x, int y, const char []); 
extern void mglPlotCenteredWindowText(int x, int y, const char []); 

extern void mglPrintLeftAlignedString(double x, double y, const char []);
extern void mglPrintRightAlignedString(double x, double y, const char []); 
extern void mglPrintCenteredString(double x, double y, const char []); 

extern void mglPlotLeftAlignedString(double x, double y, const char []);
extern void mglPlotRightAlignedString(double x, double y, const char []);
extern void mglPlotCenteredString(double x, double y, const char []);

extern void mglPrintLeftAlignedText(double x, double y, const char []); 
extern void mglPrintRightAlignedText(double x, double y, const char []); 
extern void mglPrintCenteredText(double x, double y, const char []); 

extern void mglPlotLeftAlignedText(double x, double y, const char []);
extern void mglPlotRightAlignedText(double x, double y, const char []);
extern void mglPlotCenteredText(double x, double y, const char []); 

extern void mglPlotString(double x, double y, double z, const char []); 
extern void mglPlotString(double x, double y, double z, double); 
extern void mglPlotString(double x, double y, double z, int);
extern void mglPlotString(double x, double y, const char []);
extern void mglPlotString(double x, double y, double);
extern void mglPlotString(double x, double y, int);

extern void mglPlotString(const double x[], int dim, int);
extern void mglPlotString(const double x[], int dim, double); 
extern void mglPlotString(const double x[], int dim, const char []);

extern void mglPrintString(double x, double y, double z, const char []); 
extern void mglPrintString(double x, double y, double z, double); 
extern void mglPrintString(double x, double y, double z, int); 
extern void mglPrintString(double x, double y, const char []); 
extern void mglPrintString(double x, double y, double); 
extern void mglPrintString(double x, double y, int); 

extern void mglPrintText(double x, double y, double z, const char []);
extern void mglPrintText(double x, double y, double z, double); 
extern void mglPrintText(double x, double y, double z, int);
extern void mglPrintText(double x, double y, const char []);
extern void mglPrintText(double x, double y, double);
extern void mglPrintText(double x, double y, int); 

extern void mglPlotText(double x, double y, double z, const char []); 
extern void mglPlotText(double x, double y, double z, double);
extern void mglPlotText(double x, double y, double z, int); 
extern void mglPlotText(double x, double y, const char []); 
extern void mglPlotText(double x, double y, double); 
extern void mglPlotText(double x, double y, int);

extern void mglPrintWindowString(int x, int y, const char []);
extern void mglPrintWindowString(int x, int y, int);

extern void mglPrintWindowText(int x, int y, const char []);
extern void mglPrintWindowText(int x, int y, int);

extern void mglPlotWindowString(int x, int y, const char []);
extern void mglPlotWindowString(int x, int y, int); 

extern void mglPlotWindowText(int x, int y, const char []);
extern void mglPlotWindowText(int x, int y, int); 

// MathGL mglPlot Extension.

extern void mglPlotSphere(const double [3], double);
extern void mglPlot3dSphere(const double [3], double); 
extern void mglPlotSphere(double x, double y, double z, double r);
extern void mglPlot3dSphere(double x, double y, double z, double r);

extern void mglPlot2dFunction(double (*F)(double x), double a, double b, int n);

extern void mglPlot2dCurve(double (*f)(double), double a, double b, int n);
extern void mglPlot2dCurve(double (*f)(double), double a, double b, 
                           int n, double width);
extern void mglPlot2dCurve(double (*x)(double), double (*y)(double),
                           double a, double b, int n); 
extern void mglPlot2dCurve(double (*x)(double), double (*y)(double),
                           double a, double b, int n, double width); 
extern void mglPlot2dCurve(double (*f)(double t, double x),
                           double left_end, double right_end,
                           int ncells, double t);

extern void mglPlot2dCurve(void (*F)(double t, double &x, double &y), 
                           double a, double b, int n, double width);
extern void mglPlot2dCurve(void (*F)(double t, double &x, double &y), 
                           double a, double b, int n); 

extern void mglPlotGrid(const double lo[2], const double hi[2],
                        int m, int n); 
extern void mglPlotGrid(const double lo[2], const double hi[2], 
                        int m, int n, double);

void mglPlot3dLattice(const double low[3], 
                      const double high[3],
                      int I, int J, int K);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPlotSurface(const double low[2], const double high[2],
                           int gridSize, const double value[]); 

extern void mglPlotSurfaceLattice(const double low[2], const double high[2], 
                                  int cells_number, const double value[]);

extern void mglPlotIsoSurface(double (*f)(double, double, double), 
                              double iso_value, const double low[3],
                              const double high[3], int cells_number);

extern void mglPlotImplicitSurface(const double low[3], const double high[3], 
                                   int cells_number, const double value[]);
extern void mglPlotImplicitSurface(const double low[3], const double high[3],
                                   const double value[8]); 

extern void mglPlotIsoSurface(const double low[3], const double high[3], 
                              const double value[8]);

// plot iso-surface within a hexahedron.
extern void mglPlotIsoSurface8(const double coord[8][3], 
                               const double value[8]); 
extern void mglPlotIsoSurface(const double coord[24], 
                              const double value[8]);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPlotIsoContour(const double x[], const double y[],
                              const double z[], int n, 
                              double iso_value); 
extern void mglPlotIsoContour(const double x[], const double y[], 
                              const double z[], int n, 
                              const double iso_value[], 
                              int k); 

extern void mglPlotIsoContour(const double x[], int m,
                              const double y[], int n, 
                              const double z[], double iso_value);
extern void mglPlotIsoContour(const double x[], int m, 
                              const double y[], int n, const double z[], 
                              const double iso_value[], int k);
extern void mglPlotIsoContours(const double x[], int m, 
                               const double y[], int n, const double z[],
                               const double iso_value[], int k); 

extern void mglPlotNodeDataByGrid3D(const double x[], const double y[],
                                    double **vect, int m);

extern void mglPlotNodeDataByGrid3D(const double x[], const double y[],
                                    double **vect, int m, double min_v, 
                                     double max_v);

extern void mglPlotNodeDataByIsoLines(const double x[], const double y[],
                                      double **vect, int m, int ell);

extern void mglPlotNodeDataByIsoLines(const double x[], const double y[],
                                      double **vect, int m, double min_v, 
                                      double max_v, int ell); 

extern void mglPlotIsoLines(const double x[], const double y[], 
                            double **vect, int m, double min_v,
                            double max_v, int ell); 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPlotIsoContour2(const double a[2], const double b[2],
                               const double c[2], const double d[2], 
                               const double value[4]); 
extern void mglPlotIsoContour3(const double a[3], const double b[3],
                               const double c[3], const double d[3],
                               const double value[4]);

extern void mglPlotIsoContours2(const double a[2], const double b[2], 
                                const double c[2], const double d[2], 
                                const double value[4], double lo,
                                double hi, int n);
extern void mglPlotIsoContours3(const double a[3], const double b[3], 
                                const double c[3],  const double d[3],
                                const double value[4], double lo,
                                double hi, int n); 

extern void mglPlotIsoContours2(const double a[2], const double b[2], 
                                const double c[2], const double d[2],
                                const double value[4], double lo, double hi,
                                int n, double cold, double hot); 
extern void mglPlotIsoContours3(const double a[3], const double b[3], 
                                const double c[3], const double d[3],
                                const double value[4], double lo, double hi,
                                int n, double cold, double ho);

/**
 *****************************************************************************
 *                                                                           *
 *             (low[0], high[1])_____(high[0], high[1])                      *
 *                    |                      |                               *
 *                    |                      |                               *
 *             (low[0], low[1])_______(high[0], low[1])                      *
 *                                                                           *
 *****************************************************************************
 **/
extern void mglPlotIsoContour(const double low[2], 
                              const double high[2],
                              const double v[4]); 

/**
 *****************************************************************************
 *                                                                           *
 *             (low[0], high[1])_____(high[0], high[1])                      *
 *                    |                      |                               *
 *                    |                      |                               *
 *             (low[0], low[1])_______(high[0], low[1])                      *
 *                                                                           *
 *****************************************************************************
 **/
extern void mglPlotIsoContour(const double low[2],
                              const double high[2], 
                              const double v[4], 
                              double color); 

/**
 *****************************************************************************
 *                                                                           *
 *               (c[0], c[1])_______(d[0], d[1])                             *
 *                     |                  |                                  *
 *                     |                  |                                  *
 *               (a[0], a[1])_______(b[0], b[1])                             *
 *                                                                           *
 *****************************************************************************
 **/
extern void mglPlotIsoContour(const double a[2], const double b[2],
                              const double c[2], const double d[2],
                              const double v[4]); 

/**
 *****************************************************************************
 *                                                                           *
 *               (c[0], c[1])_______(d[0], d[1])                             *
 *                     |                  |                                  *
 *                     |                  |                                  *
 *               (a[0], a[1])_______(b[0], b[1])                             *
 *                                                                           *
 *****************************************************************************
 **/
extern void mglPlotIsoContour(const double a[2], const double b[2],
                              const double c[2], const double d[2], 
                              const double v[4], double color); 

/**
 *****************************************************************************
 *                                                                           *
 *               (c[0], c[1])_______(d[0], d[1])                             *
 *                     |                  |                                  *
 *                     |                  |                                  *
 *               (a[0], a[1])_______(b[0], b[1])                             *
 *                                                                           *
 *****************************************************************************
 **/
extern void mglPlotIsoContours(const double a[2], const double b[2], 
                               const double c[2], const double d[2],
                               const double v[4], double min_u, 
                               double max_u, int n); 

/**
 *****************************************************************************
 *                                                                           *
 *             (c[0], c[1], c[2])_______(d[0], d[1], d[2])                   *
 *                      |                        |                           *
 *                      |                        |                           *
 *             (a[0], a[1], a[2])_______(b[0], b[1], b[2])                   *
 *                                                                           *
 *****************************************************************************
 **/
extern void mglPlotIsoContour3(const double a[3], const double b[3],
                               const double c[3], const double d[3], 
                               const double v[4]); 

/**
 *****************************************************************************
 *                                                                           *
 *             (c[0], c[1], c[2])_______(d[0], d[1], d[2])                   *
 *                      |                        |                           *
 *                      |                        |                           *
 *             (a[0], a[1], a[2])_______(b[0], b[1], b[2])                   *
 *                                                                           *
 *****************************************************************************
 **/
extern void mglPlotIsoContour3(const double a[3], const double b[3], 
                               const double c[3], const double d[3],
                               const double v[4], double color);

/**
 *****************************************************************************
 *                                                                           *
 *             (c[0], c[1], c[2])_______(d[0], d[1], d[2])                   *
 *                      |                        |                           *
 *                      |                        |                           *
 *             (a[0], a[1], a[2])_______(b[0], b[1], b[2])                   *
 *                                                                           *
 *****************************************************************************
 **/
extern void mglPlotIsoContours3(const double a[3], const double b[3],
                                const double c[3], const double d[3], 
                                const double v[4], double min_u, 
                                double max_u, int n);

/**
 *****************************************************************************
 *****************************************************************************
 **/
extern void mglPlotIsolines(const double coord[3][2], const double u[3], 
                            double min_u, double max_u, const int n);

/**
 *****************************************************************************
 *****************************************************************************
 **/
extern void mglPlotIsolines(double ax, double ay, double bx, double by, 
                            double cx, double cy, double ua, double ub,
                            double uc, double min_u, double max_u,
                            const int isoline_num); 

/**
 *****************************************************************************
 *****************************************************************************
 **/
extern void mglPlotIsoContours(double ax, double ay, double bx, double by,
                               double cx, double cy, double ua, double ub, 
                               double uc, double min_u, double max_u, 
                               const int isoline_num);

/**
 *****************************************************************************
 * draw iso-contours on a two-dimensional triangle.                          *
 *                                                                           *
 * mglPlotIsolinesOnTriangle.                                                *
 *****************************************************************************
 **/
extern void mglPlotIsoContours2(const double coord0[2], const double coord1[2], 
                                const double coord2[2], const double u[3], 
                                double min_u, double max_u,
                                const int n);

/**
 *****************************************************************************
 * draw iso-contours on a three-dimensional triangle.                        *
 *****************************************************************************
 **/
extern void mglPlotIsoContours3(const double coord0[3], const double coord1[3],
                                const double coord2[3], const double u[3],
                                double min_u, double max_u, 
                                const int n); 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglPlotIsoSurface(const double coord[4][3], const double data[4]); 

extern void mglPlotIsoSurfaces(const double coord[4][3], const double u[4], 
                               double min_u, double max_u, int n);

extern void mglPlotIsoSurfaces(const double a[3], const double b[3], 
                               const double c[3], const double d[3],
                               const double data[4], double min_u,
                               double max_u, int n); 

extern void mglPlotIsoSurface(const double a[3], const double b[3], 
                              const double c[3], const double d[3], 
                              const double data[4]); 

extern void mglPlotIsolines(const double coord[4][3], const double data[4],
                            const double x[4],  double min_u,
                            double max_u, int n);

extern void mglPlotColormap(const double coord[4][3], const double data[4],
                            const double x[4]);

//*****************************************************************************

void mglPlotNodeDataByColormap(const double x[], const double y[], 
                               double **vect, int I, int J);

void mglPlotNodeDataByColormap(const double x[], const double y[], 
                               double **vect, int I, int J,
                               double min_v, double max_v); 

void mglPlotNodeDataByIsolines(const double x[], const double y[],
                               double **vect, int I, int J,
                               int isovalue_num);

void mglPlotNodeDataByIsolines(const double x[], const double y[],
                               double **vect, int I, int J,
                               double min_v, double max_v,
                               int isovalue_num); 

void mglPlotNodeDataByIsoLines(const double x[], const double y[], 
                               double **vect, int I, int J, 
                               int isovalue_num); 

void mglPlotNodeDataByIsoLines(const double x[], const double y[],
                               double **vect, int I, int J, 
                               double min_v, double max_v, 
                               int isovalue_num); 

void mglPlotNodeDataByIsocontours(const double x[], const double y[],
                                  double **vect, int I, int J, 
                                  int isovalue_num); 

void mglPlotNodeDataByIsocontours(const double x[], const double y[], 
                                  double **vect, int I, int J,
                                  double min_v, double max_v, 
                                  int isovalue_num);

void mglPlotNodeDataByIsoContours(const double x[], const double y[],
                                  double **vect, int I, int J,
                                  int isovalue_num);

void mglPlotNodeDataByIsoContours(const double x[], const double y[],
                                  double **vect, int I, int J, 
                                  double min_v, double max_v, 
                                  int isovalue_num);

void mglPlotInteriorNodeDataByColormap(const double *x, const double *y, 
                                       bool **interior, double **vect,
                                       int m, int n);

void mglPlotInteriorNodeDataByColormap(const double *x, const double *y,
                                       bool **interior, double **vect,
                                       int m, int n, double min_v, 
                                       double max_v); 

void mglPlotExteriorNodeDataByColormap(const double *x, const double *y, 
                                       bool **interior, double **vect,
                                       int m, int n); 

void mglPlotExteriorNodeDataByColormap(const double *x, const double *y,
                                       bool **interior, double **vect,
                                       int m, int n, double min_v,
                                       double max_v);

void mglPlotInteriorNodeDataByIsolines(const double *x, const double *y,
                                       bool **interior, double **vect,
                                       int m, int n, int isovalue_num);

void mglPlotExteriorNodeDataByIsolines(const double *x, const double *y,
                                       bool **interior, double **vect,
                                       int m, int n, int isovalue_num); 

void mglPlotInteriorNodeDataByIsoLines(const double *x, const double *y,
                                       bool **interior, double **vect,
                                       int m, int n, int isovalue_num);

void mglPlotExteriorNodeDataByIsoLines(const double *x, const double *y,
                                       bool **interior, double **vect, 
                                       int m, int n, int isovalue_num); 

void mglPlotInteriorNodeDataByIsocontours(const double *x, const double *y, 
                                           bool **interior, double **vect, 
                                          int m, int n, int isovalue_num);

void mglPlotExteriorNodeDataByIsocontours(const double *x, const double *y,
                                          bool **interior, double **vect, 
                                          int m, int n, int isovalue_num);

void mglPlotInteriorNodeDataByIsolines(const double *x, const double *y,
                                       bool **interior, double **vect,
                                       int m, int n, double min_v, 
                                       double max_v, int isovalue_num); 

void mglPlotExteriorNodeDataByIsolines(const double *x, const double *y, 
                                       bool **interior, double **vect, 
                                       int m, int n, double min_v,
                                       double max_v, int isovalue_num); 

void mglPlotInteriorNodeDataByIsoLines(const double *x, const double *y, 
                                       bool **interior, double **vect,
                                       int m, int n, double min_v, 
                                       double max_v, int isovalue_num); 

void mglPlotExteriorNodeDataByIsoLines(const double *x, const double *y, 
                                       bool **interior, double **vect, 
                                       int m, int n, double min_v, 
                                       double max_v, int isovalue_num);

void mglPlotInteriorNodeDataByIsocontours(const double *x, const double *y,
                                           bool **interior, double **vect, 
                                          int m, int n, double min_v,
                                          double max_v, int isovalue_num);

void mglPlotExteriorNodeDataByIsocontours(const double *x, const double *y,
                                          bool **interior, double **vect, 
                                          int m, int n, double min_v,
                                          double max_v, int isovalue_num); 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern int mglPlotIsoline(const double x[], const double y[], 
                          double **data, double v, int I, 
                          int J, int m, int n); 

extern int mglPlotLevelSet(const double x[], const double y[], 
                           double **data, double v, int I,
                            int J, int m, int n);

extern int mglPlotZeroLevelSet(const double x[], const double y[], 
                               double **data, int I, int J,
                               int m, int n);

extern void mglPlotZeroLevelSet(const double *x, const double *y, 
                                int I, int J, double **v); 

extern void mglPlotLevelSet(const double *x, const double *y,
                            int I, int J, double **v, double vth);

extern void mglPlotIsoline(const double *x, const double *y,
                           int I, int J, double **v, double vth); 

extern int mglPlotIsoline(double (*V)(double, double), double v,
                          const double low[2], const double high[2], 
                          int I, int J, double delta);

extern void mglPlotLevelSet(const double low[2], const double high[2],
                            double **data, int I, int J);

extern int mglPlotZeroLevelSet(double (*V)(double, double), 
                               const double low[2], const double high[2], 
                               int I, int J, double delta);

extern void mglPlotZeroLevelSet(const double low[2], const double high[2],
                                double **data, int I, int J); 

extern void mglPlotLineShadedDomain(const double low[2], const double high[2], 
                                    double (*F)(double x, double y), 
                                    int line_num = 20); 

extern void mglPlotLineShadedDomain(const double low[2], const double high[2],
                                    bool **interior, int I, int J, 
                                    int line_num = 20);

extern void mglPlotLineShadedDomain(const double low[2], const double high[2],
                                    double (*F)(double x, double y),
                                    int deg, int line_num); 

extern void mglPlotLineShadedDomain(const double low[2], const double high[2],
                                    bool **interior, int I, int J, 
                                    int deg, int line_num);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern void mglCreateLatexFile(const char *filename,
                               const char *caption,
                               int n = 40); 

extern void mglPlotVectorField(double x0, double y0,
                               double (*f)(double, double), 
                               double (*g)(double, double), 
                               const double low[3], 
                               const double high[3], 
                               int cells_number);

extern void mglPlotRingMesh(double cx, double cy,
                            double ri, double re,
                            int m, int n); 

extern void mglPlotTriangularMesh(double ax, double ay, 
                                  double bx, double by, 
                                  double cx, double cy, 
                                  int depth, double width);

extern void mglPlotTriangularMesh(double ax, double ay, 
                                  double bx, double by,
                                  double cx, double cy,
                                  int depth);

extern void mglPlotTriangle(const double center[2], double radius);

extern void mglPlotTriangle(const double center[2], double radius, 
                            double width);

extern void mglPlotTriangle(double ox, double oy, double radius);

extern void mglPlotTriangle(double ox, double oy, double radius, 
                             double width); 

extern void mglPlotTriangularMesh(const double center[2], double radius, 
                                  int depth);

extern void mglPlotTriangularMesh(const double center[2], double radius, 
                                  int depth, double width); 

extern void mglPlotTriangularMesh(double ox, double oy, double radius, 
                                  int depth); 

extern void mglPlotTriangularMesh(double ox, double oy, double radius,
                                  int depth, double width);

extern void mglPlotTriangularMesh(double ax, double ay,
                                  double bx, double by, 
                                  double cx, double cy, 
                                  int depth);

extern void mglPlotHexagonalMesh(double cx, double cy, 
                                 double radius, int depth,
                                 double theta = 0.0);
extern void mglPlotHexagonalMesh(const double center[2],
                                 double radius, int depth,
                                 double theta = 0.0);

extern void mglPlot2dHexagon(double cx, double cy,
                             double radius, double theta = 0.0); 
extern void mglPlot2dHexagon(const double center[2],
                             double radius, double theta = 0.0); 

extern void mglPlot2dDashLattice(double low0, double low1, 
                                 double high0, double high1,
                                 int m, int n, double width); 
extern void mglPlot2dDashLattice(const double low[2], 
                                 const double high[2],
                                 int m, int n, double width);

extern void mglPlot2dDashLattice(double low0, double low1, 
                                 double high0, double high1, 
                                 int m, int n); 
extern void mglPlot2dDashLattice(const double low[2],
                                 const double high[2],
                                 int m, int n); 

extern void mglPlot2dLattice(const double low[2],
                             const double high[2],
                             int m, int n); 
extern void mglPlot2dLattice(double low0, double low1,
                             double high0, double high1,
                             int m, int n); 
extern void mglPlot2dLattice(double low0, double low1,
                             double high0, double high1,
                             int m, int n, double width); 
extern void mglPlot2dLattice(const double low[2],
                             const double high[2],
                             int m, int n, double width); 

extern void mglPlot2dQuadLattice(const double a[2], const double b[2], 
                                 const double c[2], const double d[2], 
                                 int m, int n);

extern void mglPlot2dDashQuadLattice(const double a[2], const double b[2],
                                     const double c[2], const double d[2],
                                     int m, int n); 

extern void mglPlotMatrix(const double *mat, int n2);

// MathGL multiple windows output;

extern void mglPlot2d(const char x_name[], const char y_name[], 
                      const double x[], const double y[], int n);
extern void mglPlot2d(const char name[], const double x[], 
                      const double y[], int n); 

// MathGL Image Dumping;

extern int mglImport(void); 
extern int mglImageDump(void); 
extern int mglDumpImage(void); 
extern int mglImportImage(void); 

extern int mglDumpImage(int); 
extern int mglDumpImage(const char []);
extern int mglDumpImage(const char [], int);

extern void mglRenameEpsOutput(const char []);

extern int mglXPMImageDump(const char name[], int);
extern int mglDumpXPMImage(const char name[], int);
extern int mglDumpXPM(const char name[], int); 
extern int mglXPMImageDump(const char name[]);
extern int mglDumpXPMImage(const char name[]); 
extern int mglDumpXPM(const char name[]); 

extern int mglImportGIFImage(void);
extern int mglImportGifImage(void); 

extern int mglXPMImageDump(void);
extern int mglDumpXPMImage(void);
extern int mglDumpXPMImage(int);
extern int mglDumpXPM(void); 
extern int mglDumpXPM(int); 

extern int mglDumpXWD(void);
extern int mglXWDImageDump(void);
extern int mglDumpXWDImage(void);
extern int mglDumpXWD(const char img[]); 
extern int mglXWDImageDump(const char img[]);
extern int mglDumpXWDImage(const char img[]); 

extern int mglDumpGIF(void); 
extern int mglGIFImageDump(void);
extern int mglDumpGIFImage(void); 
extern int mglDumpGIF(const char img[]); 
extern int mglGIFImageDump(const char img[]);
extern int mglDumpGIFImage(const char img[]);

extern int mglDumpDrawingAreaWindow(void); 
extern int mglDumpMainGuiWindow(void); 

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#ifdef USING_NAMESPACES 
namespace EPS
{
#endif

extern void resetShadingParameters(double a, double b); 

extern void setFrontMaterialColor(double r, double g, double b); 
extern void setBackMaterialColor(double r, double g, double b); 

extern void resetLighteningDirection(double x, double y, double z);
extern void initLighteningDirection(double x, double y, double z); 

extern void resetLighteningDirection(double ldir[3]); 
extern void initLighteningDirection(double ldir[3]); 

extern void disableInfiniteLight(void); 
extern void enableInfiniteLight(void); 

extern void disablePolygonList(void);
extern void enablePolygonList(void); 

extern void disableLightening(void);
extern void enableLightening(void); 

#ifdef USING_NAMESPACES
}
#endif

#endif 
