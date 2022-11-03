/*=============================================================================
*   
*   Filename : MyMatrix.H
*   Creator : Han Zhou
*   Date : 06/15/21
*   Description : 
*
=============================================================================*/
   
#ifndef _MYMATRIX_H
#define _MYMATRIX_H

#include <iostream>
#include <cfloat>

#include "MyVector.H"

template<typename T>
class Matrix
{
	public :
		Matrix(void);
		Matrix(size_t m, size_t n);
		Matrix(size_t m, size_t n, const T &);
		virtual ~Matrix(void);

	public :

		Matrix &operator = (const Matrix &m);

		size_t rows(void) const;
		size_t cols(void) const;

		void allocate(size_t m, size_t n);
		void reallocate(size_t m, size_t n);
		void clear(void);
		void fill(const T &a);

		void copyPointerFrom(const Matrix &m);

#ifdef DATA_DEBUG

		const Vector<T>& operator[](size_t i) const;
		Vector<T>& operator[](size_t i);

		void print(void) const;

#else 

		T* operator[](size_t i) const;
		T*& operator[](size_t i);

#endif

		T max(void) const;
		T min(void) const;

		T max_norm(void) const;
		T l1_norm(void) const;
		T l2_norm(void) const;
		T l1_norm_scaled(void) const;
		T l2_norm_scaled(void) const;

		//T max_error(const Matrix &u) const;
		//T l1_error_scaled(const Matrix &u) const;
		//T l2_error_scaled(const Matrix &u) const;

		//T inner_product(const Matrix &u) const;


	private :

#ifdef DATA_DEBUG
		Vector<T> *_data;
#else
		T **_data;
#endif

		size_t _row;
		size_t _col;

};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Matrix<T>::Matrix(void) : _row(0), _col(0), _data(0)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Matrix<T>::Matrix(size_t m, size_t n) : _row(m), _col(n)
{

#ifdef DATA_DEBUG

	_data = new Vector<T>[_row];
	for(size_t i = 0; i < _row; i++){
		_data[i].reallocate(_col);
	}

#else

	_data = new T*[m];
	for(size_t i = 0; i < m; i++){
		_data[i] = new T[n];
	}

#endif

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Matrix<T>::Matrix(size_t m, size_t n, const T &a) : _row(m), _col(n)
{
#ifdef DATA_DEBUG

	_data = new Vector<T>[_row];
	for(size_t i = 0; i < _row; i++){
		_data[i].reallocate(_col);
	}

#else

	_data = new T*[m];
	for(size_t i = 0; i < m; i++){
		_data[i] = new T[n];
	}

#endif

	for(size_t i = 0; i < m; i++){
		for(size_t j = 0; j < n; j++){
			_data[i][j] = a;
		}
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Matrix<T>::allocate(size_t m, size_t n)
{
	_row = m;
	_col = n;

#ifdef DATA_DEBUG

	_data = new Vector<T>[_row];
	for(size_t i = 0; i < _row; i++){
		_data[i].reallocate(_col);
	}

#else

	_data = new T*[m];
	for(size_t i = 0; i < m; i++){
		_data[i] = new T[n];
	}

#endif
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Matrix<T>::reallocate(size_t m, size_t n)
{
#ifdef DATA_DEBUG

	if(0 != _data) {
		delete[] _data;
	}

	_row = m;
	_col = n;

	_data = new Vector<T>[_row];
	for(size_t i = 0; i < _row; i++){
		_data[i].reallocate(_col);
	}

#else 

	if(0 != _data) {
		for(size_t i = 0; i < _row; i++){
			delete[] _data[i];
		}
		delete[] _data;
	}

	_row = m;
	_col = n;

	_data = new T*[m];
	for(size_t i = 0; i < m; i++){
		_data[i] = new T[n];
	}

#endif

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Matrix<T> &Matrix<T>::operator = (const Matrix<T> &m)
{
	if (this != &m) {
		reallocate(m.rows(), m.cols());
#ifdef USE_OPENMP
#pragma omp parallel for collapse(2)
#endif
		for(size_t i = 0; i < _row; i++){
			for(size_t j = 0; j < _col; j++){
				_data[i][j] = m._data[i][j];
			}
		}
	}
	return *this;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Matrix<T>::copyPointerFrom(const Matrix<T> &m)
{
	if (this != &m) {
		clear();
		_data = m._data;
		_row = m._row;
		_col = m._col;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Matrix<T>::~Matrix(void)
{
	clear();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline size_t Matrix<T>::rows(void) const
{
	return _row;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline size_t Matrix<T>::cols(void) const
{
	return _col;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Matrix<T>::fill(const T &a)
{
#ifdef USE_OPENMP
#pragma omp parallel for collapse(2)
#endif
	for(size_t i = 0; i < _row; i++){
		for(size_t j = 0; j < _col; j++){
			_data[i][j] = a;
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Matrix<T>::clear(void)
{

#ifdef DATA_DEBUG

	if(0 != _data) {
		delete[] _data;
		_data = 0;
		_row = _col = 0;
	}

#else

	if(0 != _data) {
		for(size_t i = 0; i < _row; i++){
			delete[] _data[i];
		}
		delete[] _data;
		_data = 0;
		_row = _col = 0;
	}

#endif

}

#ifdef DATA_DEBUG
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline const Vector<T>& Matrix<T>::operator[](size_t i) const
{
	if(i >= _row) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}

	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>& Matrix<T>::operator[](size_t i)
{
	if(i >= _row) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}

	return _data[i];
}

#else
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T* Matrix<T>::operator[](size_t i) const
{
	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T*& Matrix<T>::operator[](size_t i)
{
	return _data[i];
}
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Matrix<T>::max(void) const
{
	T max_u = _data[0][0];
	for(size_t i = 0; i < _row; i++){
		for(size_t j = 0; j < _col; j++){
			max_u = max_u > _data[i][j] ? max_u : _data[i][j];
		}
	}
	return max_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Matrix<T>::min(void) const
{
	T min_u = _data[0][0];
	for(size_t i = 0; i < _row; i++){
		for(size_t j = 0; j < _col; j++){
			min_u = min_u < _data[i][j] ? min_u : _data[i][j];
		}
	}
	return min_u;
}

#ifdef DATA_DEBUG
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Matrix<T>::print(void) const
{

	if (std::is_same<T, double>::value || 
			std::is_same<T, float>::value) {

  	for (size_t i = 0; i < _row; i++) {
  	  for (size_t j = 0; j < _col; j++) {
  	    T u = _data[i][j]; 
  	    if (fabs(u) < DBL_EPSILON) {
  	      //printf("          "); 
  	      printf("         0"); 
  	    } else {
  	      printf("%8.4f  ", u);
  	    }
  	  }
  	  printf("\n");
  	}
  	printf("\n"); 

	} else if (std::is_same<T, int>::value || 
						 std::is_same<T, unsigned int>::value || 
						 std::is_same<T, long unsigned int>::value) {

  	for (size_t i = 0; i < _row; i++) {
  	  for (size_t j = 0; j < _col; j++) {
  	    T u = _data[i][j]; 
  	    //if (u == 0) {
  	    //  printf("        "); 
  	    //} else {
  	    //  printf("%6d  ", u);
  	    //}
  	      printf("%6d  ", u);
  	  }
  	  printf("\n");
  	}
  	printf("\n"); 

	}

}
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Matrix<T>::max_norm(void) const
{
	T max_u = fabs(_data[0][0]);

	for(size_t i = 0; i < _row; i++){
		for(size_t j = 0; j < _col; j++){
			T u = fabs(_data[i][j]);
			max_u = max_u > u ? max_u : u;
		}
	}
	return max_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Matrix<T>::l1_norm(void) const
{
	T sum = 0.0;
	for(size_t i = 0; i < _row; i++){
		for(size_t j = 0; j < _col; j++){
			sum += fabs(_data[i][j]);
		}
	}
	return sum;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Matrix<T>::l2_norm(void) const
{
	T sum = 0.0;
	for(size_t i = 0; i < _row; i++){
		for(size_t j = 0; j < _col; j++){
			sum += _data[i][j] * _data[i][j];
		}
	}
	return sqrt(sum);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Matrix<T>::l1_norm_scaled(void) const
{
	T sum = 0.0;
	for(size_t i = 0; i < _row; i++){
		for(size_t j = 0; j < _col; j++){
			sum += fabs(_data[i][j]);
		}
	}
	return sum / (_row*_col);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Matrix<T>::l2_norm_scaled(void) const
{
	T sum = 0.0;
	for(size_t i = 0; i < _row; i++){
		for(size_t j = 0; j < _col; j++){
			sum += _data[i][j] * _data[i][j];
		}
	}
	return sqrt(sum / (_row*_col));
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

