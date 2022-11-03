
/**
 *************************************************************************** 
 * The OpenGL-based C/C++ codes were written for students by Wenjun Ying,  *
 * School of Mathematical Sciences and Institute of Natural Sciences,      *
 * Shanghai Jiao Tong University, Shanghai 200240.                         *
 ***************************************************************************
 **/

#include <cmath> 

#ifdef MAC_OS_X 
#include <GLUT/glut.h> 
#else 
#include <GL/glut.h>
#endif 

#include "NewMathGL.h" 

//*****************************************************************************
namespace mgl { // start of the namespace "mgl"
//*****************************************************************************

static bool mouse_in_motion = false; 
static bool left_button_down = false; 

static double low[2] = {-1.1, -1.1};
static double high[2] = {1.1, 1.1};

static double history_low[100][2];
static double history_high[100][2]; 

static int max_zooming_times = 100;
static int zooming_times = 0; 

static double temp[2][2]; 

static int win_width = 640; 
static int win_height = 640; 

static int win_x_offset = 0;
static int win_y_offset = 0; 

static int viewport_width = 500; 
static int viewport_height = 500; 

static double aspect_ratio = 1.0; 

static int (*mgl_special_func)(int key) = 0; 

static int (*mgl_motion_func)(int, int) = 0; 

static int (*mgl_key_func)(unsigned char ch) = 0;

static int (*mgl_mouse_func)(int, int, int, int) = 0; 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void display(void) 
{
  glClearColor(1.0, 1.0, 1.0, 1.0); 
  glClear(GL_COLOR_BUFFER_BIT);

  mglPlotAll();

  if (mouse_in_motion) {
    glColor3d(0.0, 0.0, 0.0); 
    glLineStipple(4, 0x9999); 
    glEnable(GL_LINE_STIPPLE); 
    glBegin(GL_LINE_LOOP);
      glVertex2d(temp[0][0], temp[0][1]);
      glVertex2d(temp[0][0], temp[1][1]);
      glVertex2d(temp[1][0], temp[1][1]); 
      glVertex2d(temp[1][0], temp[0][1]);
    glEnd(); 
    glDisable(GL_LINE_STIPPLE);
  }

  glFlush();
  glutSwapBuffers();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void resetDrawingArea(const double a[2], const double b[2])
{
  double lo[2], hi[2];

  if (a[0] < b[0]) {
    lo[0] = a[0];
    hi[0] = b[0];
  } else {
    lo[0] = b[0]; 
    hi[0] = a[0]; 
  }
  if (a[1] < b[1]) {
    lo[1] = a[1]; 
    hi[1] = b[1];
  } else {
    lo[1] = b[1];
    hi[1] = a[1]; 
  }

  double c[2], r[2]; 

  c[0] = 0.5 * (hi[0] + lo[0]);
  c[1] = 0.5 * (hi[1] + lo[1]);

  r[0] = 0.5 * (hi[0] - lo[0]); 
  r[1] = 0.5 * (hi[1] - lo[1]); 

  if ((r[1] / r[0]) < aspect_ratio) {
    r[1] = r[0] * aspect_ratio; 
    lo[1] = c[1] - r[1]; 
    hi[1] = c[1] + r[1]; 
  } else {
    r[0] = r[1] / aspect_ratio;
    lo[0] = c[0] - r[0]; 
    hi[0] = c[0] + r[0];
  }

  for (int d = 0; d < 2; d++) {
    low[d] = lo[d];
    high[d] = hi[d]; 
  }

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glOrtho(low[0], high[0], low[1], high[1], - 1.0, 1.0); 

  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline bool checkIfLargeEnough(const double lo[2], const double hi[2])
{
  double rx = fabs(hi[0] - lo[0]) / (high[0] - low[0]); 
  double ry = fabs(hi[1] - lo[1]) / (high[1] - low[1]); 
  if ((rx > 0.02) || (ry > 0.02)) {
    return true;
  } else {
    return false;
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mouse_func(int button, int state, int a, int b)
{
  if (mgl_mouse_func) {
    int value = mgl_mouse_func(button, state, a, b);
    if (value != 0) {
      return;
    }
  }

  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      double x = (a - win_x_offset) / (double)viewport_width; 
      double y = 1.0 - (b - win_y_offset) / (double)viewport_height;
      x = low[0] + (high[0] - low[0]) * x;
      y = low[1] + (high[1] - low[1]) * y; 
      temp[0][0] = x; 
      temp[0][1] = y; 
      left_button_down = true;
    } else if (state == GLUT_UP) {
      if (zooming_times < max_zooming_times) {
        if (mouse_in_motion && checkIfLargeEnough(temp[0], temp[1])) {
          for (int d = 0; d < 2; d++) {
            history_low[zooming_times][d] = low[d]; 
            history_high[zooming_times][d] = high[d]; 
          }
          zooming_times++; 
          resetDrawingArea(temp[0], temp[1]); 
        }
      }
      left_button_down = false;
      mouse_in_motion = false; 
      glutPostRedisplay();
    }
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void motion_func(int a, int b) 
{
  if (mgl_motion_func) {
    int value = mgl_motion_func(a, b);
    if (value != 0) {
      return; 
    }
  }

  if (left_button_down) {
    double x = (a - win_x_offset) / (double)viewport_width;
    double y = 1.0 - (b - win_y_offset) / (double)viewport_height; 
    x = low[0] + (high[0] - low[0]) * x; 
    y = low[1] + (high[1] - low[1]) * y; 
    temp[1][0] = x; 
    temp[1][1] = y; 
    mouse_in_motion = true; 
    glutPostRedisplay();
  }
}

////:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//void run(void) 
//{
//  double t = 0.0; 
//  double T = 1.0;
//  double dt = 0.01; 
//
//  t += dt;
//
//  mglNewpage();
//  //plotRandomObjects(); 
//  mglFlush();
//
//  glutIdleFunc(0);
//}
//
////:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//void start(int value) 
//{
//  mglNewpage(); 
//  //plotMultipleCircles(100); 
//  mglFlush(); 
//}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void key_func(unsigned char ch, int a, int b) 
{
  if (mgl_key_func) {
    if (ch == 8) {
      if (zooming_times > 0) {
        zooming_times--; 
        for (int d = 0; d < 2; d++) {
          low[d] = history_low[zooming_times][d]; 
          high[d] = history_high[zooming_times][d];
        }
        resetDrawingArea(low, high);
        glutPostRedisplay();
        return; 
      } else {
        int value = mgl_key_func(ch); 
        if (value != 0) {
          return; 
        }
      }
    } else {
      int value = mgl_key_func(ch);
      if (value != 0) {
        return; 
      }
    }
  }

  if ((ch == 'q') || (ch == 'Q') || (ch == 27)) {
    exit(0);
  }

  if (ch == ' ') {
    glutPostRedisplay();
  } else if (ch == 13) {
  } else if (ch == 'i') {
  } else if (ch == 'a') {
    //glutIdleFunc(run); 
  } else if (ch == 'p') {
  } else if (ch == 'l') {
  } else if (ch == 'c') {
  } else if (ch == 't') {
  } else if (ch == 'r') {
  } else if (ch == 's') {
    //glutTimerFunc(10, start, 1); 
  } else if (ch == 'g') {
  } else if (ch == 'R') {
  } else if (ch == 8) {
    if (zooming_times > 0) {
      zooming_times--; 
      for (int d = 0; d < 2; d++) {
        low[d] = history_low[zooming_times][d];
        high[d] = history_high[zooming_times][d];
      }
      resetDrawingArea(low, high);
      glutPostRedisplay();
      return;
    }
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void reshape(int w, int h)
{
  win_width = w; 
  win_height = h; 

  if (h > (w*aspect_ratio)) {
    viewport_width = w; 
    viewport_height = w * aspect_ratio; 
    int j = (h - w * aspect_ratio) / 2; 
    glViewport(0, j, viewport_width, viewport_height); 
  } else {
    viewport_height = h;
    viewport_width = h / aspect_ratio;
    int i = (w - h / aspect_ratio) / 2; 
    glViewport(i, 0, viewport_width, viewport_height); 
  }

  win_x_offset = (win_width - viewport_width) / 2;
  win_y_offset = (win_height - viewport_height) / 2; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool checkIfMatched(const double p[2], const double q[2], 
                    const double u[2], const double v[2])
{
	double e = fabs(p[0] - u[0]) + fabs(p[1] - u[1]);
	double f = fabs(q[0] - v[0]) + fabs(q[1] - v[1]);
  if ((e < 1.0E-8) && (f < 1.0E-8)) {
    return true;
  } else {
    return false; 
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void special_func(int key, int a, int b) 
{
  if (mgl_special_func) {
    int value = mgl_special_func(key);
    if (value != 0) {
      return; 
    }
  }

  if (key == GLUT_KEY_PAGE_DOWN) {

    if (zooming_times < max_zooming_times) {

      double cx = 0.5 * (high[0] + low[0]);
      double cy = 0.5 * (high[1] + low[1]);
      double rx0 = 0.5 * (high[0] - low[0]); 
      double ry0 = 0.5 * (high[1] - low[1]);

      double rx = 0.8 * rx0;
      double ry = 0.8 * ry0; 

      double lo[2], hi[2]; 

      lo[0] = cx - rx; 
      lo[1] = cy - ry;
      hi[0] = cx + rx; 
      hi[1] = cy + ry;

      if (zooming_times > 0) {
        int k = zooming_times - 1;
        if (checkIfMatched(lo, hi, history_low[k], history_high[k])) {
          zooming_times--;
        } else {
          for (int d = 0; d < 2; d++) {
            history_low[zooming_times][d] = low[d];
            history_high[zooming_times][d] = high[d]; 
          }
          zooming_times++;
        }
      } else {
        for (int d = 0; d < 2; d++) {
          history_low[zooming_times][d] = low[d]; 
          history_high[zooming_times][d] = high[d];
        }
        zooming_times++; 
      }

      resetDrawingArea(lo, hi); 
      glutPostRedisplay(); 
    }

  } else if (key == GLUT_KEY_PAGE_UP) {

    if (zooming_times < max_zooming_times) {

      double cx = 0.5 * (high[0] + low[0]); 
      double cy = 0.5 * (high[1] + low[1]);
      double rx0 = 0.5 * (high[0] - low[0]);
      double ry0 = 0.5 * (high[1] - low[1]);

      double rx = 1.25 * rx0;
      double ry = 1.25 * ry0; 

      double lo[2], hi[2];

      lo[0] = cx - rx; 
      lo[1] = cy - ry; 
      hi[0] = cx + rx; 
      hi[1] = cy + ry; 

      if (zooming_times > 0) {
        int k = zooming_times - 1;
        if (checkIfMatched(lo, hi, history_low[k], history_high[k])) {
          zooming_times--;
        } else {
          for (int d = 0; d < 2; d++) {
            history_low[zooming_times][d] = low[d];
            history_high[zooming_times][d] = high[d];
          }
          zooming_times++; 
        }
      } else {
        for (int d = 0; d < 2; d++) {
          history_low[zooming_times][d] = low[d]; 
          history_high[zooming_times][d] = high[d]; 
        }
        zooming_times++;
      }

      resetDrawingArea(lo, hi);
      glutPostRedisplay();
    }

  } else if (key == GLUT_KEY_F9) {

    if (zooming_times > 0) {
      zooming_times = 0;
      for (int d = 0; d < 2; d++) {
        low[d] = history_low[zooming_times][d]; 
        high[d] = history_high[zooming_times][d];
      }
      resetDrawingArea(low, high); 
      glutPostRedisplay();
    }

  } else {

  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void initialize(int argc, char *argv[]) 
{
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_DOUBLE); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void createWindow(int width, int height, const char title[])
{
	win_width = width;
	win_height = height;

  viewport_width = win_width; 
  viewport_height = win_height; 

  aspect_ratio = static_cast<double>(viewport_height) / viewport_width;

  glutInitWindowSize(win_width, win_height);

  glutCreateWindow(title); 

  glutReshapeFunc(reshape); 
  glutDisplayFunc(display);
  glutKeyboardFunc(key_func); 

  glutMouseFunc(mouse_func); 
  glutMotionFunc(motion_func); 
  glutSpecialFunc(special_func);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void setDrawingArea(double lx, double ly, double hx, double hy)
{
	low[0] = lx;	low[1] = hx;
	high[0] = ly;	high[1] = hy;

  glOrtho(lx, ly, hx, hy, - 1.0, 1.0);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void enterMainLoop(void)
{
  glutMainLoop();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void setSpecialFunc(int (*F)(int)) 
{
  mgl_special_func = F;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void setMotionFunc(int (*F)(int, int)) 
{
  mgl_motion_func = F;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void setKeyboardFunc(int (*F)(unsigned char)) 
{
  mgl_key_func = F;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void setMouseFunc(int (*F)(int, int, int, int))
{
  mgl_mouse_func = F;
}

//*****************************************************************************
} // end of the namespace "mgl"
//*****************************************************************************

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
