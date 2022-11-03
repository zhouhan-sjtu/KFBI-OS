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

#include <cmath> 
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream> 
#include <iostream> 

#include <time.h> 
#include <float.h>
#include <limits.h> 
#include <unistd.h> 

#include "Const.h" 
#include "NewMathGL.h"

#ifdef MAC_OS_X 
#include <GLUT/glut.h> 
#else 
#include <GL/glut.h>
#endif 

#ifdef USE_GL2PS
#include "gl2ps.h"
#endif

class PlottingElement
{
  public :
    PlottingElement(void); 
    virtual ~PlottingElement(void);
    PlottingElement(const PlottingElement &c); 
    PlottingElement(int type, char *buf, int size); 

    PlottingElement &operator = (const PlottingElement &); 

    void draw(void) const; 

    PlottingElement *getPrev(void); 
    PlottingElement *getNext(void); 

    void setPrev(PlottingElement *ptr);
    void setNext(PlottingElement *ptr);

    const PlottingElement *getPrev(void) const;
    const PlottingElement *getNext(void) const;

  private :
    int _type; 
    int _size;
    char *_buffer; 

    PlottingElement *_prev; 
    PlottingElement *_next;
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement::PlottingElement(const PlottingElement &c)
{
  _prev = 0;
  _next = 0; 

  _type = c._type;
  _size = c._size; 

  if (0 != c._buffer) {
    _buffer = new char[_size];
    for (int i = 0; i < _size; i++) {
      _buffer[i] = c._buffer[i]; 
    }
  } else {
    _buffer = 0;
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement &PlottingElement::operator = (const PlottingElement &c)
{
  if (this != &c) {
    _prev = 0; 
    _next = 0;

    if (0 != _buffer) {
      delete[] _buffer; 
      _buffer = 0;
      _size = 0; 
    }

    _size = c._size;
    _type = c._type; 

    if (0 != c._buffer) {
      _buffer = new char[_size];
      for (int i = 0; i < _size; i++) {
        _buffer[i] = c._buffer[i];
      }
    } else {
      _buffer = 0;
    }
  }
  return *this; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement::PlottingElement(int type, char *buf, int size)
{
  _type = type;
  _size = size;
  _buffer = buf; 

  _prev = 0; 
  _next = 0; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement::PlottingElement(void)
{
  _prev = 0; 
  _next = 0;

  _size = 0;
  _type = - 1;
  _buffer = 0; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement::~PlottingElement(void)
{
  _prev = 0;

  PlottingElement *ptr = _next;
  while (ptr != 0) {
    PlottingElement *t = ptr->_next; 
    ptr->_next = 0;
    delete ptr;
    ptr = t; 
  }
  ptr = 0;

  // //

  _size = 0;
  _type = 0; 

  if (0 != _buffer) {
    delete[] _buffer;
    _buffer = 0;
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement *PlottingElement::getPrev(void) 
{
  return _prev;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement *PlottingElement::getNext(void) 
{
  return _next; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void PlottingElement::setPrev(PlottingElement *ptr) 
{
  _prev = ptr; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void PlottingElement::setNext(PlottingElement *ptr) 
{
  _next = ptr;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const PlottingElement *PlottingElement::getPrev(void) const 
{
  return _prev;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const PlottingElement *PlottingElement::getNext(void) const
{
  return _next;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawIsolines(double coord[3][2], const double u[3], 
                  double min_u, double max_u,
                  const int isoval_num) 
{
  double v[3];
  int index[3]; 
  bool status = false;
  double ax, ay, bx, by;
  double iso_value = min_u; 
  const double delta = (max_u - min_u) / isoval_num; 
  const double denom = static_cast<double>(isoval_num);
  for (int k = 0; k <= isoval_num; k++) {
    for (int j = 0; j < 3; j++) {
      v[j] = u[j] - iso_value;
      if (v[j] < 0.0) {
        index[j] = 0;
      } else {
        index[j] = 1;
      }
    }
    if (index[0] == 0) {
      if (index[1] == 0) {
        if (index[2] == 1) {
          double alpha0 = v[2] / (v[2] - v[0]);
          double beta0 = 1.0 - alpha0;
          ax = alpha0 * coord[0][0] + beta0 * coord[2][0];
          ay = alpha0 * coord[0][1] + beta0 * coord[2][1]; 
          double alpha1 = v[2] / (v[2] - v[1]);
          double beta1 = 1.0 - alpha1;
          bx = alpha1 * coord[1][0] + beta1 * coord[2][0];
          by = alpha1 * coord[1][1] + beta1 * coord[2][1];
          double r = sqrt((bx - ax) * (bx - ax) + (by - ay) * (by - ay)); 
          if (r > EPSILON) {
            status = true;
          }
        }
      } else {
        if (index[2] == 0) {
          double alpha0 = v[1] / (v[1] - v[0]);
          double beta0 = 1.0 - alpha0;
          ax = alpha0 * coord[0][0] + beta0 * coord[1][0]; 
          ay = alpha0 * coord[0][1] + beta0 * coord[1][1]; 
          double alpha2 = v[1] / (v[1] - v[2]);
          double beta2 = 1.0 - alpha2; 
          bx = alpha2 * coord[2][0] + beta2 * coord[1][0]; 
          by = alpha2 * coord[2][1] + beta2 * coord[1][1]; 
          double r = sqrt((bx - ax) * (bx - ax) + (by - ay) * (by - ay));
          if (r > EPSILON) {
            status = true; 
          }
        } else {
          double alpha1 = v[0] / (v[0] - v[1]);
          double beta1 = 1.0 - alpha1; 
          ax = alpha1 * coord[1][0] + beta1 * coord[0][0]; 
          ay = alpha1 * coord[1][1] + beta1 * coord[0][1];
          double alpha2 = v[0] / (v[0] - v[2]); 
          double beta2 = 1.0 - alpha2; 
          bx = alpha2 * coord[2][0] + beta2 * coord[0][0];
          by = alpha2 * coord[2][1] + beta2 * coord[0][1];
          double r = sqrt((bx - ax) * (bx - ax) + (by - ay) * (by - ay)); 
          if (r > EPSILON) {
            status = true;
          }
        }
      }
    } else {
      if (index[1] == 0) {
        if (index[2] == 0) {
          double alpha1 = v[0] / (v[0] - v[1]);
          double beta1 = 1.0 - alpha1; 
          ax = alpha1 * coord[1][0] + beta1 * coord[0][0]; 
          ay = alpha1 * coord[1][1] + beta1 * coord[0][1]; 
          double alpha2 = v[0] / (v[0] - v[2]);
          double beta2 = 1.0 - alpha2;
          bx = alpha2 * coord[2][0] + beta2 * coord[0][0]; 
          by = alpha2 * coord[2][1] + beta2 * coord[0][1]; 
          double r = sqrt((bx - ax) * (bx - ax) + (by - ay) * (by - ay)); 
          if (r > EPSILON) {
            status = true; 
          }
        } else {
          double alpha0 = v[1] / (v[1] - v[0]); 
          double beta0 = 1.0 - alpha0;
          ax = alpha0 * coord[0][0] + beta0 * coord[1][0];
          ay = alpha0 * coord[0][1] + beta0 * coord[1][1];
          double alpha2 = v[1] / (v[1] - v[2]); 
          double beta2 = 1.0 - alpha2;
          bx = alpha2 * coord[2][0] + beta2 * coord[1][0];
          by = alpha2 * coord[2][1] + beta2 * coord[1][1]; 
          double r = sqrt((bx - ax) * (bx - ax) + (by - ay) * (by - ay));
          if (r > EPSILON) {
            status = true; 
          }
        }
      } else {
        if (index[2] == 0) {
          double alpha0 = v[2] / (v[2] - v[0]);
          double beta0 = 1.0 - alpha0; 
          ax = alpha0 * coord[0][0] + beta0 * coord[2][0];
          ay = alpha0 * coord[0][1] + beta0 * coord[2][1];
          double alpha1 = v[2] / (v[2] - v[1]);
          double beta1 = 1.0 - alpha1;
          bx = alpha1 * coord[1][0] + beta1 * coord[2][0]; 
          by = alpha1 * coord[1][1] + beta1 * coord[2][1];
          double r = sqrt((bx - ax) * (bx - ax) + (by - ay) * (by - ay));
          if (r > EPSILON) {
            status = true; 
          }
        }
      }
    }
    if (status == true) {
      status = false;
      extern void mglSetColor(double); 
      mglSetColor(k / denom); 
      //mglPlot2dLine(ax, ay, bx, by);
      glBegin(GL_LINES);
        glVertex2d(ax, ay); 
        glVertex2d(bx, by);
      glEnd(); 
    }
    iso_value += delta;
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void PlottingElement::draw(void) const 
{
  //std::cout << "type = " << _type << std::endl;

  if (_type == 0) { // set point size

    double *data = (double *)_buffer; 
    glPointSize(data[0]);

#ifdef USE_GL2PS
  	gl2psPointSize(data[0]);
#endif

  } else if (_type == 1) { // set line width

    double *data = (double *)_buffer; 
    glLineWidth(data[0]);

#ifdef USE_GL2PS
  	gl2psLineWidth(data[0]);
#endif

  } else if (_type == 2) { // set RGB color

    double *data = (double *)_buffer; 

    double r = data[0]; 
    double g = data[1];
    double b = data[2]; 

    glColor3d(r, g, b);

  } else if (_type == 3) { // point

    double *data = (double *)_buffer; 

    glBegin(GL_POINTS); 
      glVertex2d(data[0], data[1]);
    glEnd(); 

  } else if (_type == 4) { // line segment

    double *data = (double *)_buffer;

    glBegin(GL_LINES); 
      glVertex2d(data[0], data[1]); 
      glVertex2d(data[2], data[3]); 
    glEnd(); 

  } else if (_type == 5) { // dash line segment

    double *data = (double *)_buffer; 

    glLineStipple(4, 0x9999); 
    glEnable(GL_LINE_STIPPLE); 

#ifdef USE_GL2PS
  	gl2psEnable(GL2PS_LINE_STIPPLE);
#endif

    glBegin(GL_LINES);
      glVertex2d(data[0], data[1]);
      glVertex2d(data[2], data[3]);
    glEnd(); 
    glDisable(GL_LINE_STIPPLE); 

#ifdef USE_GL2PS
  	gl2psDisable(GL2PS_LINE_STIPPLE);
#endif

  } else if (_type == 6) { // circle

    double *data = (double *)_buffer;

    double cx = data[0];
    double cy = data[1];
    double r = data[2];

    int n = 400; 
    int n1 = n + 1;

    double delta = (M_PI + M_PI) / n; 
    double theta = 0.0;

    glBegin(GL_LINE_STRIP);
      for (int i = 0; i < n1; i++) {
        double x = cx + r * cos(theta);
        double y = cy + r * sin(theta);
        glVertex2d(x, y); 
        theta += delta;
      }
    glEnd();

  } else if (_type == 7) { // dash circle

    double *data = (double *)_buffer;

    double cx = data[0]; 
    double cy = data[1]; 
    double r = data[2]; 

    int n = 400;
    int n1 = n + 1;

    double delta = (M_PI + M_PI) / n;
    double theta = 0.0;

    glBegin(GL_LINES);
      for (int i = 0; i < n; i++) {
        double x = cx + r * cos(theta);
        double y = cy + r * sin(theta); 
        int l = i % 8;
        if ((l >= 0) && (l < 4)) {
          glVertex2d(x, y); 
        }
        theta += delta;
      }
    glEnd();

  } else if (_type == 8) { // triangle

    double *data = (double *)_buffer; 

    glBegin(GL_LINE_LOOP);
      glVertex2d(data[0], data[1]); 
      glVertex2d(data[2], data[3]);
      glVertex2d(data[4], data[5]);
    glEnd();

  } else if (_type == 9) { // rectangle

    double *data = (double *)_buffer; 

    double lx = data[0]; 
    double ly = data[1];
    double hx = data[2]; 
    double hy = data[3]; 

    glBegin(GL_LINE_LOOP);
      glVertex2d(lx, ly); 
      glVertex2d(hx, ly);
      glVertex2d(hx, hy); 
      glVertex2d(lx, hy);
    glEnd();

  } else if (_type == 10) { // quadrilateral

    double *data = (double *)_buffer;

    glBegin(GL_LINE_LOOP);
      for (int i = 0, j = 0; i < 4; i++, j += 2) {
        glVertex2d(data[j], data[j + 1]); 
      }
    glEnd(); 

  } else if (_type == 11) { // solid circle

    double *data = (double *)_buffer; 

    double cx = data[0];
    double cy = data[1];
    double r = data[2];

    int n = 400; 

    double delta = (M_PI + M_PI) / n; 
    double theta = 0.0;

    glBegin(GL_POLYGON); 
      for (int i = 0; i < n; i++) {
        double x = cx + r * cos(theta);
        double y = cy + r * sin(theta);
        glVertex2d(x, y); 
        theta += delta;
      }
    glEnd();

  } else if (_type == 12) { // solid triangle

    double *data = (double *)_buffer; 

    glBegin(GL_TRIANGLES);
      glVertex2d(data[0], data[1]);
      glVertex2d(data[2], data[3]); 
      glVertex2d(data[4], data[5]);
    glEnd(); 

  } else if (_type == 13) { // solid rectangle

    double *data = (double *)_buffer;

    double lx = data[0];
    double ly = data[1];
    double hx = data[2]; 
    double hy = data[3]; 

    //glRectd(lx, ly, hx, hy);

    glBegin(GL_QUADS);
      glVertex2d(lx, ly);
      glVertex2d(hx, ly);
      glVertex2d(hx, hy); 
      glVertex2d(lx, hy);
    glEnd(); 

  } else if (_type == 14) { // solid quadrilateral

    double *data = (double *)_buffer;

    glBegin(GL_QUADS); 
      glVertex2d(data[0], data[1]); 
      glVertex2d(data[2], data[3]);
      glVertex2d(data[4], data[5]);
      glVertex2d(data[6], data[7]); 
    glEnd(); 

  } else if (_type == 15) { // polygon

    int *num = (int *)_buffer;

    int n = num[0];

    int offset = sizeof(int);

    double *data = (double *)(_buffer + offset); 

    glBegin(GL_LINE_STRIP);
      for (int i = 0, j = 0; i < n; i++, j += 2) {
        glVertex2d(data[j], data[j + 1]);
      }
      glVertex2d(data[0], data[1]); 
    glEnd();

  } else if (_type == 16) { // solid polygon

    int *num = (int *)_buffer; 

    int n = num[0]; 

    int offset = sizeof(int); 

    double *data = (double *)(_buffer + offset);

    glBegin(GL_POLYGON); 
      for (int i = 0, j = 0; i < n; i++, j += 2) {
        glVertex2d(data[j], data[j + 1]);
      }
      glVertex2d(data[0], data[1]); 
    glEnd(); 

  } else if (_type == 17) { // multiple points

    int *num = (int *)_buffer; 

    int n = num[0];

    int offset = sizeof(int);

    double *data = (double *)(_buffer + offset); 

    glBegin(GL_POINTS);
      for (int i = 0, j = 0; i < n; i++, j += 2) {
        glVertex2d(data[j], data[j + 1]); 
      }
    glEnd();

  } else if (_type == 18) { // multiple lines

    int *num = (int *)_buffer;

    int n = num[0]; 

    int offset = sizeof(int);

    double *data = (double *)(_buffer + offset);

    glBegin(GL_LINE_STRIP); 
      for (int i = 0, j = 0; i < n; i++, j += 2) {
        glVertex2d(data[j], data[j + 1]);
      }
    glEnd();

  } else if (_type == 19) { // string

    int *num = (int *)_buffer;

    int n = num[0]; 

    int offset = sizeof(int);

    double *coord = (double *)(_buffer + offset);

    double x = coord[0];
    double y = coord[1]; 

    offset += 2 * sizeof(double);

    char *str = (char *)(_buffer + offset); 

#ifdef USE_GL2PS
  	const char *fonts[] =
  	  { "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic",
  	    "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique",
  	    "Courier", "Courier-Bold", "Courier-Oblique", "Courier-BoldOblique",
  	    "Symbol", "ZapfDingbats" };
#endif

    glRasterPos2d(x, y); 
#ifdef USE_GL2PS
		//double angle = 0.0;
  	//gl2psTextOpt(str, fonts[0], 24, GL2PS_TEXT_BL, angle);
  	gl2psText(str, fonts[0], 24);
#endif
    for (int i = 0; i < n; i++) {
      //glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, str[i]);
    }

  } else if (_type == 20) { // isolines

    int *num = (int *)_buffer;

    int isoval_n = num[0];

    int offset = sizeof(int); 

    double *dptr = (double *)(_buffer + offset); 

    double coord[3][2], u[3], min_u, max_u; 

    coord[0][0] = *dptr++;
    coord[0][1] = *dptr++; 
    coord[1][0] = *dptr++; 
    coord[1][1] = *dptr++; 
    coord[2][0] = *dptr++;
    coord[2][1] = *dptr++; 

    u[0] = *dptr++; 
    u[1] = *dptr++;
    u[2] = *dptr++; 

    min_u = *dptr++; 
    max_u = *dptr++;

    drawIsolines(coord, u, min_u, max_u, isoval_n); 

  } else if (_type == 21) { // solid triangle with multiple colors

    extern void mglSetColor(double); 

    double *data = (double *)_buffer; 

    double red, green, blue; 

    glBegin(GL_TRIANGLES);
      mglSetColor(data[2]);
      glVertex2d(data[0], data[1]); 
      mglSetColor(data[5]); 
      glVertex2d(data[3], data[4]);
      mglSetColor(data[8]);
      glVertex2d(data[6], data[7]);
    glEnd();

  }
}

//*****************************************************************************

class PlottingList 
{
  public :
    PlottingList(void); 
    virtual ~PlottingList(void);
    PlottingList(const PlottingList &);

    PlottingList &operator = (PlottingList &c); 

    void clear(void);
    void draw(void) const;

    int getLength(void) const;

    PlottingElement *getHead(void);
    PlottingElement *getTail(void); 

    void append(PlottingElement *ptr); 

    const PlottingElement *getHead(void) const; 
    const PlottingElement *getTail(void) const;

    void setColormap(const double map[2]);
    void getColormap(double map[2]) const; 

    void setColormap(double lo, double hi);
    void getColormap(double &lo, double &hi) const;

  private :
    int _length; 

    double _colormap[2]; 

    PlottingElement *_head;
    PlottingElement *_tail; 
}; 

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int PlottingList::getLength(void) const
{
  return _length; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingList::PlottingList(void) 
{
  _head = 0; 
  _tail = 0;

  _colormap[0] = 0.0; 
  _colormap[1] = 1.0; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingList::~PlottingList(void)
{
  _tail = 0; 
  if (0 != _head) {
    delete _head; 
    _head = 0; 
  }

  _length = 0; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingList::PlottingList(const PlottingList &c)
{
  _head = 0; 
  _tail = 0;

  _colormap[0] = c._colormap[0]; 
  _colormap[1] = c._colormap[1]; 

  const PlottingElement *ptr = c.getHead(); 
  if (0 != ptr) {
    _head = new PlottingElement(*ptr); 

    PlottingElement *q = _head;

    ptr = ptr->getNext(); 
    while (ptr != 0) {
      PlottingElement *t = new PlottingElement(*ptr); 
      ptr = ptr->getNext();
      q->setNext(t); 
      q = t;
    }
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingList &PlottingList::operator = (PlottingList &c) 
{
  if (this != &c) {

    clear();

    _colormap[0] = c._colormap[0];
    _colormap[1] = c._colormap[1]; 

    const PlottingElement *ptr = c.getHead(); 
    if (0 != ptr) {
      _head = new PlottingElement(*ptr);

      PlottingElement *q = _head; 

      ptr = ptr->getNext();
      while (ptr != 0) {
        PlottingElement *t = new PlottingElement(*ptr); 
        ptr = ptr->getNext(); 
        q->setNext(t);
        q = t; 
      }
    }
  }
  return *this;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement *PlottingList::getHead(void) 
{
  return _head;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline PlottingElement *PlottingList::getTail(void) 
{
  return _tail; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const PlottingElement *PlottingList::getHead(void) const
{
  return _head;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline const PlottingElement *PlottingList::getTail(void) const 
{
  return _tail;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void PlottingList::append(PlottingElement *ptr) 
{
  _length++; 

  if (_head == 0) {
    _head = ptr;
    _tail = ptr; 
  } else {
    _tail->setNext(ptr); 
    _tail = ptr; 
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void PlottingList::clear(void)
{
  if (0 != _head) {
    delete _head;
    _head = 0; 
  }
  _tail = 0;

  _head = 0; 

  _length = 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void PlottingList::setColormap(const double map[2])
{
  _colormap[0] = map[0];
  _colormap[1] = map[1]; 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void PlottingList::getColormap(double map[2]) const 
{
  map[0] = _colormap[0];
  map[1] = _colormap[1];
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void PlottingList::setColormap(double lo, double hi) 
{
  _colormap[0] = lo; 
  _colormap[1] = hi; 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void PlottingList::getColormap(double &lo, double &hi) const 
{
  lo = _colormap[0];
  hi = _colormap[1];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void PlottingList::draw(void) const
{
  PlottingElement *ptr = _head;
  while (0 != ptr) {
    ptr->draw(); 
    ptr = ptr->getNext(); 
  }
}

//*****************************************************************************

static PlottingList plotting_list; 

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotAll(void) 
{
  plotting_list.draw(); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglClear(void)
{
  plotting_list.clear(); 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglNewpage(void) 
{
  plotting_list.clear();
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglNewPage(void)
{
  plotting_list.clear();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglSetColormap(double lo, double hi)
{
  plotting_list.setColormap(lo, hi);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglResetColormap(double lo, double hi) 
{
  plotting_list.setColormap(lo, hi);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglSetColormap(const double colormap[2])
{
  plotting_list.setColormap(colormap);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglResetColormap(const double colormap[2])
{
  plotting_list.setColormap(colormap); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglSetPointSize(double sz) 
{
  int type = 0;

  int size = sizeof(double); 

  char *buf = new char[size];

  double *data = (double *)buf;

  data[0] = sz; 

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglSetLineWidth(double wd) 
{
  int type = 1;

  int size = sizeof(double); 

  char *buf = new char[size];

  double *data = (double *)buf;

  data[0] = wd; 

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglSetColor(double r, double g, double b) 
{
  int type = 2; 

  int size = 3*sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf; 

  data[0] = r; 
  data[1] = g; 
  data[2] = b; 

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotPoint(double x, double y)
{
  int type = 3; 

  int size = 2*sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf; 
  data[0] = x;  data[1] = y; 

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotLine(double x1, double y1, double x2, double y2) 
{
  int type = 4; 

  int size = 4*sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf;
  data[0] = x1;  data[1] = y1;  data[2] = x2;  data[3] = y2; 

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e); 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotDashLine(double x1, double y1, double x2, double y2)
{
  int type = 5; 

  int size = 4*sizeof(double); 

  char *buf = new char[size]; 

  double *data = (double *)buf; 
  data[0] = x1;  data[1] = y1;  data[2] = x2;  data[3] = y2; 

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotCircle(double x, double y, double r) 
{
  int type = 6;

  int size = 3*sizeof(double); 

  char *buf = new char[size]; 

  double *data = (double *)buf; 
  data[0] = x;  data[1] = y;  data[2] = r;

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e); 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotDashCircle(double x, double y, double r) 
{
  int type = 7; 

  int size = 3*sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf;
  data[0] = x;  data[1] = y;  data[2] = r; 

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotTriangle(double x1, double y1,
                     double x2, double y2, 
                     double x3, double y3)
{
  int type = 8;

  int size = 6*sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf; 

  data[0] = x1;  data[1] = y1; 
  data[2] = x2;  data[3] = y2; 
  data[4] = x3;  data[5] = y3; 

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e); 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotRectangle(double lx, double ly, double hx, double hy)
{
  int type = 9; 

  int size = 4*sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf;

  data[0] = lx;  data[1] = ly;
  data[2] = hx;  data[3] = hy; 

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotQuadrilateral(double x1, double y1,
                          double x2, double y2,
                          double x3, double y3, 
                          double x4, double y4)
{
  int type = 10;

  int size = 8*sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf; 

  data[0] = x1;  data[1] = y1; 
  data[2] = x2;  data[3] = y2; 
  data[4] = x3;  data[5] = y3; 
  data[6] = x4;  data[7] = y4;

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotSolidCircle(double x, double y, double r) 
{
  int type = 11; 

  int size = 3*sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf; 
  data[0] = x;  data[1] = y;  data[2] = r;

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e); 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotSolidTriangle(double x1, double y1, 
                          double x2, double y2, 
                          double x3, double y3) 
{
  int type = 12;

  int size = 6*sizeof(double);

  char *buf = new char[size]; 

  double *data = (double *)buf;

  data[0] = x1;  data[1] = y1;
  data[2] = x2;  data[3] = y2; 
  data[4] = x3;  data[5] = y3;

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e); 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotSolidRectangle(double lx, double ly, double hx, double hy) 
{
  int type = 13;

  int size = 4*sizeof(double); 

  char *buf = new char[size];

  double *data = (double *)buf; 

  data[0] = lx;  data[1] = ly;
  data[2] = hx;  data[3] = hy; 

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e); 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotSolidQuadrilateral(double x1, double y1, 
                               double x2, double y2, 
                               double x3, double y3,
                               double x4, double y4)
{
  int type = 14; 

  int size = 8*sizeof(double); 

  char *buf = new char[size]; 

  double *data = (double *)buf;

  data[0] = x1;  data[1] = y1; 
  data[2] = x2;  data[3] = y2;
  data[4] = x3;  data[5] = y3;
  data[6] = x4;  data[7] = y4;

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotPolygon(const double x[], const double y[], int n)
{
  int type = 15; // wire polygon

  int p = sizeof(int);

  int size = p + (n + n) * sizeof(double); 

  char *buf = new char[size];

  int *num = (int *) buf; 

  double *data = (double *) (buf + p);

  num[0] = n;

  for (int i = 0, j = 0; i < n; i++) {
    data[j++] = x[i];  data[j++] = y[i];
  }

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotSolidPolygon(const double x[], const double y[], int n)
{
  int type = 16; // solid polygon

  int p = sizeof(int); 

  int size = p + (n + n) * sizeof(double);

  char *buf = new char[size];

  int *num = (int *) buf; 

  double *data = (double *) (buf + p);

  num[0] = n; 

  for (int i = 0, j = 0; i < n; i++) {
    data[j++] = x[i];  data[j++] = y[i];
  }

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotPoints(const double x[], const double y[], int n)
{
  int type = 17;  // multiple points

  int p = sizeof(int); 

  int size = p + (n + n) * sizeof(double);

  char *buf = new char[size]; 

  int *num = (int *) buf;

  double *data = (double *) (buf + p);

  num[0] = n;

  for (int i = 0, j = 0; i < n; i++) {
    data[j++] = x[i];  data[j++] = y[i];
  }

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotLines(const double x[], const double y[], int n) 
{
  int type = 18; // multiple lines

  int p = sizeof(int);

  int size = p + (n + n) * sizeof(double); 

  char *buf = new char[size]; 

  int *num = (int *)buf; 

  double *data = (double *) (buf + p);

  num[0] = n;

  for (int i = 0, j = 0; i < n; i++) {
    data[j++] = x[i];  data[j++] = y[i];
  }

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotString(double x, double y, const char *str)
{
  int type = 19; // string

  int p1 = sizeof(int);
  int p2 = 2 *sizeof(double); 

  int len = strlen(str); 

  int size = p1 + p2 + len * sizeof(char); 

  char *buf = new char[size]; 

  int *num = (int *)buf; 

  double *coord = (double *) (buf + p1); 

  num[0] = len;

  coord[0] = x;
  coord[1] = y;

  char *data = (char *) (buf + p1 + p2); 

  for (int i = 0; i < len; i++) {
    data[i] = str[i]; 
  }

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotString(double x, double y, int idx)
{
  char str[32]; 
  sprintf(str, "%d", idx);
  mglPlotString(x, y, str); 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotString(double x, double y, double value)
{
  char str[32]; 
  sprintf(str, "%f", value); 
  mglPlotString(x, y, str); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglFlush(void)
{
  glutPostRedisplay();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotIsolines(const double coord[3][2], const double u[3],
                     double min_u, double max_u, int isoval_n)
{
  int type = 20; // isolines

  int p = sizeof(int); 

  int size = p + 11 * sizeof(double);

  char *buf = new char[size];

  int *num = (int *)buf;

  double *dptr = (double *) (buf + p);

  num[0] = isoval_n; 

  *dptr++ = coord[0][0]; 
  *dptr++ = coord[0][1]; 
  *dptr++ = coord[1][0]; 
  *dptr++ = coord[1][1]; 
  *dptr++ = coord[2][0];
  *dptr++ = coord[2][1];

  *dptr++ = u[0];
  *dptr++ = u[1];
  *dptr++ = u[2]; 

  *dptr++ = min_u; 
  *dptr++ = max_u; 

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotSolidTriangle(double x1, double y1, double v1, 
                          double x2, double y2, double v2,
                          double x3, double y3, double v3) 
{
  int type = 21; // a solid triangle with multiple colors

  int size = 9 * sizeof(double);

  char *buf = new char[size];

  double *data = (double *)buf; 

  double map[2]; 

  plotting_list.getColormap(map); 

  double c1 = (v1 - map[0]) / (map[1] - map[0]);
  double c2 = (v2 - map[0]) / (map[1] - map[0]); 
  double c3 = (v3 - map[0]) / (map[1] - map[0]); 

  data[0] = x1;  data[1] = y1;  data[2] = c1; 
  data[3] = x2;  data[4] = y2;  data[5] = c2;
  data[6] = x3;  data[7] = y3;  data[8] = c3;

  PlottingElement *e = new PlottingElement(type, buf, size); 
  plotting_list.append(e);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglPlotSolidTriangle(const double coord[3][2],
                          const double u[3]) 
{
  int type = 21; // a solid triangle with multiple colors

  int size = 9 * sizeof(double);

  char *buf = new char[size]; 

  double *data = (double *)buf; 

  double map[2];

  plotting_list.getColormap(map); 

  double c1 = (u[0] - map[0]) / (map[1] - map[0]);
  double c2 = (u[1] - map[0]) / (map[1] - map[0]);
  double c3 = (u[2] - map[0]) / (map[1] - map[0]);

  data[0] = coord[0][0];  data[1] = coord[0][1];  data[2] = c1; 
  data[3] = coord[1][0];  data[4] = coord[1][1];  data[5] = c2;
  data[6] = coord[2][0];  data[7] = coord[2][1];  data[8] = c3;

  PlottingElement *e = new PlottingElement(type, buf, size);
  plotting_list.append(e);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void allocateColor(double c, double &r, double &g, double &b)
{
  int k; 
  double d;
  if (c <= 0.0) {
    k = - 1; 
    d = - 4.0 * c; 

    if (d > 1.0) {
      d = 1.0;
    }
  } else if (c >= 1.0) {
    k = 4; 
    d = 4.0 * (c - 1.0);

    if (d > 1.0) {
      d = 1.0;
    }
  } else {
    double a = 4 * c; 
    k = (int) a;
    d = a - k; 
  }

  switch (k) {
    case (-1) : {
      r = 0.75 * d; 
      g = 0.0; 
      b = 1.0; 
      break;
    }
    case 0 : {
      r = 0.0; 
      g = d; 
      b = 1.0;
      break; 
    }
    case 1 : {
      r = 0.0;
      g = 1.0; 
      b = 1.0 - d; 
      break;
    }
    case 2 : {
      r = d; 
      g = 1.0; 
      b = 0.0; 
      break; 
    }
    case 3 : {
      r = 1.0;
      g = 1.0 - d;
      b = 0.0;
      break; 
    }
    default : {
      r = 1.0; 
      g = 0;
      b = 0; 
      break;
    }
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void mglSetColor(double c)
{
  double r, g, b;

  allocateColor(c, r, g, b); 

  glColor3d(r, g, b);
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
