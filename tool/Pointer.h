/*=============================================================================
*   
*   Filename : Pointer.h
*   Creator : Han Zhou
*   Date : 06/15/21
*   Description : 
*
=============================================================================*/
   
#ifndef __Pointers_h_IS_INCLUDED__
#define __Pointers_h_IS_INCLUDED__ 

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T *allocateVector(size_t n)
{
  T *v = new T[n]; 
  for (size_t i = 0; i < n; i++) {
    v[i] = 0;
  }
  return v;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void freeVector(T *&v) 
{
  delete[] v;
  v = 0;
}

//*****************************************************************************
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T **allocateMatrix(size_t n) 
{
  T **A = new T*[n];
  for (size_t i = 0; i < n; i++) {
    A[i] = new T[n];
    for (size_t j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }
  return A; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T **allocateMatrix(size_t m, size_t n) 
{
  T **A = new T*[m];
  for (size_t i = 0; i < m; i++) {
    A[i] = new T[n];
    for (size_t j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }
  return A; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void freeMatrix(T **&A, size_t n)
{
  for (size_t i = 0; i < n; i++) {
    delete[] A[i]; 
  }
  delete[] A;
  A = 0; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void freeMatrix(T **&A, size_t m, size_t n)
{
  for (size_t i = 0; i < m; i++) {
    delete[] A[i]; 
  }
  delete[] A;
  A = 0; 
}

//*****************************************************************************
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T ***allocateTensor3d(size_t n) 
{
  T ***A = new T**[n];
  for (size_t i = 0; i < n; i++) {
    A[i] = new T*[n];
    for (size_t j = 0; j < n; j++) {
			A[i][j] = new T[n];
			for(size_t k = 0; k < n; k++){
      	A[i][j][k] = 0;
			}
    }
  }
  return A; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T ***allocateTensor3d(size_t m, size_t n, size_t l) 
{
  T ***A = new T**[m];
  for (size_t i = 0; i < m; i++) {
    A[i] = new T*[n];
    for (size_t j = 0; j < n; j++) {
			A[i][j] = new T[l];
			for(size_t k = 0; k < l; k++){
      	A[i][j][k] = 0;
			}
    }
  }
  return A; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void freeTensor3d(T ***&A, size_t n)
{
  for (size_t i = 0; i < n; i++) {
		for(size_t j = 0; j < n; j++){
    	delete[] A[i][j]; 
		}
    delete[] A[i]; 
  }
  delete[] A;
  A = 0; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void freeTensor3d(T ***&A, size_t m, size_t n, size_t l)
{
  for (size_t i = 0; i < m; i++) {
		for(size_t j = 0; j < n; j++){
			delete[] A[i][j];
		}
    delete[] A[i]; 
  }
  delete[] A;
  A = 0; 
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif 
