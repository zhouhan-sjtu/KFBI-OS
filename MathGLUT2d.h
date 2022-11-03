
/**
 *************************************************************************** 
 * The OpenGL-based C/C++ codes were written for students by Wenjun Ying,  *
 * School of Mathematical Sciences and Institute of Natural Sciences,      *
 * Shanghai Jiao Tong University, Shanghai 200240.                         *
 ***************************************************************************
 **/

#ifndef __MathGLUT2d_h_IS_INCLUDED__
#define __MathGLUT2d_h_IS_INCLUDED__

//*****************************************************************************
namespace mgl { // start of the namespace "mgl"
//*****************************************************************************

extern void display(void); 
extern void initialize(int argc, char *argv[]); 

extern void createWindow(int width, int height,
                          const char title[]);

extern void setDrawingArea(double lx, double ly, 
                           double hx, double hy);

extern void enterMainLoop(void); 

extern void setSpecialFunc(int (*)(int)); 
extern void setMotionFunc(int (*)(int, int));
extern void setKeyboardFunc(int (*)(unsigned char));
extern void setMouseFunc(int (*)(int, int, int, int)); 

//*****************************************************************************
} // end of the namespace "mgl"
//*****************************************************************************

#endif
