/*=============================================================================
*   
*   Filename : MyTensor.H
*   Creator : Han Zhou
*   Date : 06/15/21
*   Description : 
*
=============================================================================*/
   
#ifndef _MYTENSOR_H
#define _MYTENSOR_H

#include <iostream>

#include "MyMatrix.H"

//=============================================================================
//*****************************************************************************

template<typename T>
class Tensor
{
	public :
		Tensor(void);
		Tensor(size_t n);
		Tensor(size_t I, size_t J, size_t K);
		virtual ~Tensor(void);

	public :

#ifdef DATA_DEBUG
		const Matrix<T>& operator[](size_t i) const;
		Matrix<T>& operator[](size_t i);
#else 
		T** operator[](size_t i) const;
		T**& operator[](size_t i);
#endif

		Tensor &operator = (const Tensor &m);

		void get_dim(size_t dim[3]) const;
		void get_dim(size_t&, size_t&, size_t&) const;

		void allocate(size_t I, size_t J, size_t K);
		void reallocate(size_t I, size_t J, size_t K);
		void clear(void);
		void fill(T a);

		T max(void) const;
		T min(void) const;

	private :

#ifdef DATA_DEBUG
		Matrix<T> *_data;
#else
		T ***_data;
#endif

		size_t _dim1;
		size_t _dim2;
		size_t _dim3;
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Tensor<T>::Tensor(void) : _dim1(0), _dim2(0), _dim3(0), _data(0)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Tensor<T>::Tensor(size_t n) : _dim1(n), _dim2(n), _dim3(n)
{
	allocate(n, n, n);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Tensor<T>::Tensor(size_t I, size_t J, size_t K)
 : _dim1(I), _dim2(J), _dim3(K)
{
	allocate(I, J, K);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Tensor<T>::allocate(size_t I, size_t J, size_t K)
{
	_dim1 = I;
	_dim2 = J;
	_dim3 = K;

#ifdef DATA_DEBUG

	_data = new Matrix<T>[_dim1];
	for(size_t i = 0; i < _dim1; i++){
		_data[i].allocate(_dim2, _dim3);
	}

#else

	_data = new T**[_dim1];
	for(size_t i = 0; i < _dim1; i++){
		_data[i] = new T*[_dim2];
		for(size_t j = 0; j < _dim2; j++){
			_data[i][j] = new T[_dim3];
		}
	}

#endif
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Tensor<T>::reallocate(size_t I, size_t J, size_t K)
{
#ifdef DATA_DEBUG

	if(0 != _data) {
		delete[] _data;
	}

	_dim1 = I;
	_dim2 = J;
	_dim3 = K;

	_data = new Matrix<T>[_dim1];
	for(size_t i = 0; i < _dim1; i++){
		_data[i].allocate(_dim2, _dim3);
	}

#else 

	if(0 != _data) {
		for(size_t i = 0; i < _dim1; i++){
			for(size_t j = 0; j < _dim2; j++){
				delete[] _data[i][j];
			}
			delete[] _data[i];
		}
		delete[] _data;
	}

	_dim1 = I;
	_dim2 = J;
	_dim3 = K;

	_data = new T**[I];
	for(size_t i = 0; i < I; i++){
		_data[i] = new T*[J];
		for(size_t j = 0; j < J; j++){
			_data[i][j] = new T[K];
		}
	}

#endif

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Tensor<T>::~Tensor(void)
{
	clear();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Tensor<T>::get_dim(size_t dim[3]) const
{
	dim[0] = _dim1;
	dim[1] = _dim2;
	dim[2] = _dim3;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Tensor<T>::get_dim(size_t &n1, size_t &n2, size_t &n3) const
{
	n1 = _dim1;
	n2 = _dim2;
	n3 = _dim3;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Tensor<T>::fill(T a)
{
#ifdef USE_OPENMP
#pragma omp parallel for collapse(3)
#endif
	for(size_t i = 0; i < _dim1; i++){
		for(size_t j = 0; j < _dim2; j++){
			for(size_t k = 0; k < _dim3; k++){
				_data[i][j][k] = a;
			}
		}
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Tensor<T>::clear(void)
{
#ifdef DATA_DEBUG

	if(0 != _data) {
		delete[] _data;
		_data = 0;
		_dim1 = _dim2 = _dim3 = 0;
	}


#else

	if(0 != _data) {
		for(size_t i = 0; i < _dim1; i++){
			for(size_t j = 0; j < _dim2; j++){
				delete[] _data[i][j];
			}
			delete[] _data[i];
		}
		delete[] _data;	
		_data = 0;
		
		_dim1 = _dim2 = _dim3 = 0;
	}

#endif

}

#ifdef DATA_DEBUG
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline const Matrix<T>& Tensor<T>::operator[](size_t i) const
{
	if(i >= _dim1) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}

	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Matrix<T>& Tensor<T>::operator[](size_t i)
{
	if(i >= _dim1) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}

	return _data[i];
}

#else
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T** Tensor<T>::operator[](size_t i) const
{
	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T**& Tensor<T>::operator[](size_t i)
{
	return _data[i];
}
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Tensor<T>::max(void) const
{
	T max_u = _data[0][0][0];

	for(size_t i = 0; i < _dim1; i++){
		for(size_t j = 0; j < _dim2; j++){
			for(size_t k = 0; k < _dim3; k++){
				max_u = max_u > _data[i][j][k] ? max_u : _data[i][j][k];
			}
		}
	}
	return max_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Tensor<T>::min(void) const
{
	T min_u = _data[0][0][0];

	for(size_t i = 0; i < _dim1; i++){
		for(size_t j = 0; j < _dim2; j++){
			for(size_t k = 0; k < _dim3; k++){
				min_u = min_u < _data[i][j][k] ? min_u : _data[i][j][k];
			}
		}
	}
	return min_u;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Tensor<T> &Tensor<T>::operator = (const Tensor<T> &m)
{
	if (this != &m) {

		if(_dim1 != m._dim1 || _dim2 != m._dim2 || _dim3 != m._dim3) {
			reallocate(m._dim1, m._dim2, m._dim3);
		}

#ifdef USE_OPENMP
#pragma omp parallel for collapse(3)
#endif
		for(size_t i = 0; i < _dim1; i++){
			for(size_t j = 0; j < _dim2; j++){
				for(size_t k = 0; k < _dim3; k++){
					_data[i][j][k] = m._data[i][j][k];
				}
			}
		}
	}
	return *this;
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

