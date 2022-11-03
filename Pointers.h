// Copyright (c) Wenjun Ying, School of Mathematical Sciences and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is free software; as a special exception the author gives
// unlimited permission to copy and/or distribute it, with or without
// modifications, as long as this notice is preserved.

// This file is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#ifndef __Pointers_h_IS_INCLUDED__
#define __Pointers_h_IS_INCLUDED__ 

#include <complex>

inline double *allocateVector(double (*f)(double x),
                              double a, double b, int n)
{
  double *v = new double[n]; 

  double x = a;  v[0] = f(x);
  double h = (b - a) / (n - 1.0);
  for (int i = 1; i < n; i++) {
    x += h;  v[i] = f(x); 
  }

  return v;
}

inline double *allocateVector(double a, double b, int n)
{
  double *v = new double[n]; 

  v[0] = a; 
  double h = (b - a) / (n - 1.0); 
  for (int i = 1; i < n; i++) {
    v[i] = v[i - 1] + h;
  }

  return v;
}

inline double *allocateVector(int n) 
{
  double *v = new double[n]; 
  for (int i = 0; i < n; i++) {
    v[i] = 0; 
  }
  return v; 
}

inline int *allocateIntVector(int n) 
{
  int *v = new int[n]; 
  for (int i = 0; i < n; i++) {
    v[i] = 0; 
  }
  return v; 
}

inline bool *allocateBoolVector(int n) 
{
  bool *v = new bool[n];
  for (int i = 0; i < n; i++) {
    v[i] = false; 
  }
  return v;
}

inline double **allocateMatrix(int n)
{
  double **A = new double *[n]; 
  A[0] = new double[n*n]; 
  for (int i = 1; i < n; i++) {
    A[i] = A[i - 1] + n;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }
  return A; 
}

inline double **allocateMatrix2(int n)
{
  double **A = new double *[n];
  for (int i = 0; i < n; i++) {
    A[i] = new double[n]; 
    for (int j = 0; j < n; j++) {
      A[i][j] = 0.0; 
    }
  }
  return A; 
}

inline int **allocateIntMatrix(int n) 
{
  int **A = new int *[n];
  A[0] = new int[n*n]; 
  for (int i = 1; i < n; i++) {
    A[i] = A[i - 1] + n; 
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }
  return A; 
}

inline bool **allocateBoolMatrix(int n) 
{
  bool **A = new bool *[n]; 
  A[0] = new bool[n * n];
  for (int i = 1; i < n; i++) {
    A[i] = A[i - 1] + n; 
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = false; 
    }
  }
  return A;
}

inline double **allocateMatrix(int m, int n) 
{
  double **A = new double *[m];
  A[0] = new double[m*n]; 
  for (int i = 1; i < m; i++) {
    A[i] = A[i - 1] + n;
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }
  return A;
}

inline int **allocateIntMatrix(int m, int n)
{
  int **A = new int *[m]; 
  A[0] = new int[m*n]; 
  for (int i = 1; i < m; i++) {
    A[i] = A[i - 1] + n; 
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = 0; 
    }
  }
  return A;
}

inline bool **allocateBoolMatrix(int m, int n) 
{
  bool **A = new bool *[m];
  A[0] = new bool[m * n];
  for (int i = 1; i < m; i++) {
    A[i] = A[i - 1] + n;
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = false; 
    }
  }
  return A;
}

inline void freeVector(int *&v) 
{
  delete[] v;
  v = 0;
}

inline void freeVector(bool *&v) 
{
  delete[] v;
  v = 0;
}

inline void freeVector(double *&v) 
{
  delete[] v;
  v = 0;
}

inline void freeIntVector(int *&v)
{
  delete[] v; 
  v = 0;
}

inline void freeBoolVector(bool *&v)
{
  delete[] v;
  v = 0;
}

inline void freeMatrix(int **&A) 
{
  delete[] A[0];
  delete[] A;
  A = 0; 
}

inline void freeMatrix(bool **&A) 
{
  delete[] A[0];
  delete[] A; 
  A = 0; 
}

inline void freeMatrix(double **&A)
{
  delete[] A[0];
  delete[] A; 
  A = 0; 
}

inline void freeIntMatrix(int **&A)
{
  delete[] A[0];
  delete[] A;
  A = 0;
}

inline void freeBoolMatrix(bool **&A)
{
  delete[] A[0];
  delete[] A; 
  A = 0;
}

inline void freeMatrix2(double **&A, int n) 
{
  for (int i = 0; i < n; i++) {
    delete[] A[i];
    A[i] = 0;
  }
  delete[] A; 
  A = 0; 
}

inline std::complex<double> *allocateComplexVector(int n)
{
  std::complex<double> *v = new std::complex<double>[n]; 
  for (int i = 0; i < n; i++) {
    v[i] = 0;
  }
  return v;
}

inline std::complex<double> **allocateComplexMatrix(int n)
{
  std::complex<double> **A = new std::complex<double> *[n]; 
  A[0] = new std::complex<double>[n*n];
  for (int i = 1; i < n; i++) {
    A[i] = A[i - 1] + n; 
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = 0; 
    }
  }
  return A;
}

inline void freeComplexMatrix(std::complex<double> **&A) 
{
  delete[] A[0];
  delete[] A;
  A = 0; 
}

inline void freeMatrix(std::complex<double> **&A) 
{
  delete[] A[0]; 
  delete[] A;
  A = 0;
}

inline void freeVector(std::complex<double> *&v) 
{
  delete[] v; 
  v = 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
void allocateGridNodeData(T **&data, int I, int J)
{
  int I1 = I + 1;
  int J1 = J + 1; 

  data = new T*[I1]; 
  data[0] = new T[I1 * J1]; 
  for (int j = 0; j < J1; j++) {
    data[0][j] = 0; 
  }
  for (int i = 1; i < I1; i++) {
    data[i] = data[i - 1] + J1; 
    for (int j = 0; j < J1; j++) {
      data[i][j] = 0; 
    }
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
void freeGridNodeData(T **&data) 
{
  delete[] data[0];
  delete[] data;
  data = 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
void freeGridNodeData(T **&data, int I, int J) 
{
  delete[] data[0];
  delete[] data;
  data = 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
void allocateGridNodeData(T ***&data, int I, int J, int K)
{
  int I1 = I + 1; 
  int J1 = J + 1;
  int K1 = K + 1; 

  data = new T**[I1]; 
  for (int i = 0; i < I1; i++) {
    data[i] = new T*[J1]; 
    for (int j = 0; j < J1; j++) {
      data[i][j] = new T[K1]; 
      for (int k = 0; k < K1; k++) {
        data[i][j][k] = 0; 
      }
    }
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
void freeGridNodeData(T ***&data, int I, int J, int K) 
{
  if (0 != data) {
    int I1 = I + 1; 
    int J1 = J + 1;
    for (int i = 0; i < I1; i++) {
      for (int j = 0; j < J1; j++) {
        delete[] data[i][j]; 
      }
      delete[] data[i];
    }
    delete[] data;
    data = 0;
  }
}

#endif
