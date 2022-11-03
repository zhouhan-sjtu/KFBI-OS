// Copyright (c) 2013, Wenjun Ying, Department of Mathematics and Institute of 
// Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240, P.R. China.

// The graphics package was developed for the course of numerical analysis and
// scientific computing that the author taught in the Zhiyuan College of 
// Shanghai Jiao Tong University in the Spring of 2013. 

// This file is free software; as a special exception the author gives unlimited
// permission to copy and/or distribute it, with or without modifications, as 
// long as this notice is preserved.

// This file is distributed in the hope that it will be useful, but WITHOUT ANY 
// WARRANTY, to the extent permitted by law; without even the implied warranty 
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#ifndef __MinMax_h_IS_INCLUDED__
#define __MinMax_h_IS_INCLUDED__ 

#include <cmath>

template <typename T>
T Max(T a, T b)
{
  return a > b ? a : b; 
}

template <typename T>
T Min(T a, T b) 
{
  return a > b ? b : a; 
}

template <typename T>
T max(T a, T b) 
{
  return a > b ? a : b; 
}

template <typename T>
T min(T a, T b) 
{
  return a > b ? b : a;
}

template <typename T>
T Max(T a, T b, T c) 
{
  if (a > b) {
    return a > c ? a : c; 
  } else {
    return b > c ? b : c;
  }
}

template <typename T>
T max(T a, T b, T c)
{
  if (a > b) {
    return a > c ? a : c; 
  } else {
    return b > c ? b : c;
  }
}

template <typename T>
T Min(T a, T b, T c)
{
  if (a < b) {
    return a < c ? a : c; 
  } else {
    return b < c ? b : c; 
  }
}

template <typename T>
T min(T a, T b, T c) 
{
  if (a < b) {
    return a < c ? a : c;
  } else {
    return b < c ? b : c; 
  }
}

template <typename T>
T Max(T a, T b, T c, T d)
{
  if (a > b) {
    return Max(a, c, d); 
  } else {
    return Max(b, c, d); 
  }
}

template <typename T>
T max(T a, T b, T c, T d) 
{
  if (a > b) {
    return Max(a, c, d);
  } else {
    return Max(b, c, d);
  }
}

template <typename T>
T Min(T a, T b, T c, T d) 
{
  if (a < b) {
    return Min(a, c, d); 
  } else {
    return Min(b, c, d); 
  }
}

template <typename T>
T min(T a, T b, T c, T d) 
{
  if (a < b) {
    return min(a, c, d);
  } else {
    return min(b, c, d); 
  }
}

template <typename T>
T min(const T data[], int n, int &idx)
{
  idx = 0; 
  T r = data[0];
  for (int i = 1; i < n; i++) {
    if (data[i] < r) {
      r = data[i]; 
      idx = i;
    }
  }
  return r;
}

template <typename T>
T min(const T data[], int n) 
{
  T r = data[0]; 
  for (int i = 1; i < n; i++) {
    if (data[i] < r) {
      r = data[i];
    }
  }
  return r; 
}

template <typename T>
T Min(const T data[], int n) 
{
  T r = data[0];
  for (int i = 1; i < n; i++) {
    if (data[i] < r) {
      r = data[i];
    }
  }
  return r;
}

template <typename T>
T max(const T data[], int n)
{
  T r = data[0]; 
  for (int i = 1; i < n; i++) {
    if (data[i] > r) {
      r = data[i]; 
    }
  }
  return r;
}

template <typename T>
T max(const T data[], int n, int &idx) 
{
  idx = 0;
  T r = data[0];
  for (int i = 1; i < n; i++) {
    if (data[i] > r) {
      r = data[i]; 
      idx = i; 
    }
  }
  return r; 
}

template <typename T>
T Max(const T data[], int n)
{
  T r = data[0]; 
  for (int i = 1; i < n; i++) {
    if (data[i] > r) {
      r = data[i];
    }
  }
  return r;
}

template <typename T>
T min(T **const u, int n)
{
  T r = u[0][0];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      T s = u[i][j];
      if (s < r) {
        r = s;
      }
    }
  }
  return r;
}

template <typename T>
T max(T **const u, int n) 
{
  T r = u[0][0]; 
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      T s = u[i][j]; 
      if (s > r) {
        r = s;
      }
    }
  }
  return r; 
}

template <typename T>
T min(T **const u, int m, int n) 
{
  T r = u[0][0];
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      T s = u[i][j];
      if (s < r) {
        r = s; 
      }
    }
  }
  return r; 
}

template <typename T>
T max(T **const u, int m, int n)
{
  T r = u[0][0]; 
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      T s = u[i][j]; 
      if (s > r) {
        r = s; 
      }
    }
  }
  return r; 
}

template <typename T>
T max(T **u, const int size[2])
{
  T max_u = u[0][0];
  for (int i = 0; i < size[0]; i++) {
    for (int j = 0; j < size[1]; j++) {
      T a = u[i][j]; 
      if (a > max_u) {
        max_u = a; 
      }
    }
  }
  return max_u; 
}

template <typename T>
T min(T **u, const int size[2]) 
{
  T min_u = u[0][0];
  for (int i = 0; i < size[0]; i++) {
    for (int j = 0; j < size[1]; j++) {
      T a = u[i][j]; 
      if (a < min_u) {
        min_u = a; 
      }
    }
  }
  return min_u; 
}

#endif 
