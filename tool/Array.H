/*=============================================================================
*   
*   Filename : Array.H
*   Creator : Han Zhou
*   Date : 06/15/21
*   Description : 
*
=============================================================================*/
 
#ifndef _ARRAY_H
#define _ARRAY_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

template<typename T, size_t N>
class Array
{
	public :
		Array(void);
		virtual ~Array(void){};

		Array(const Array& a);

	public :

		T operator[](size_t) const;
		T& operator[](size_t);

		Array &operator = (const Array &v);

#ifdef DATA_DEBUG
		void fill(T a);
		size_t size(void) const;
#endif
		void print(void) const;
		void normalize(void);
		T inner_product(const Array &u) const;

	private :

		T _data[N];

#ifdef DATA_DEBUG
		size_t _size;
#endif
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline Array<T, N>::Array(void)
{
#ifdef DATA_DEBUG
	_size = N;
#endif
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline Array<T, N>::Array(const Array &a)
{
#ifdef DATA_DEBUG
	_size = N;
#endif
	for(size_t i = 0; i < N; i++){
		_data[i] = a._data[i];
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifdef DATA_DEBUG
template<typename T, size_t N>
inline void Array<T, N>::fill(T a)
{
	for(size_t i = 0; i < _size; i++){
		_data[i] = a;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline size_t Array<T, N>::size(void) const
{
	return _size;
}
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline T Array<T, N>::operator[](size_t i) const
{
#ifdef DATA_DEBUG
	if(i >= _size) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}
#endif

	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline T& Array<T, N>::operator[](size_t i)
{

#ifdef DATA_DEBUG
	if(i >= _size) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}
#endif

	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline void Array<T, N>::normalize(void)
{
	T sum = 0;
	for(size_t i = 0; i < N; i++){
		sum += _data[i] * _data[i];
	}

#ifdef DATA_DEBUG
	if(sum < 1.0E-15) {
		std::cout << "zero array." << std::endl;
		exit(1);
	}
#endif

	T r_n = 1.0 / sqrt(sum);
	for(size_t i = 0; i < N; i++){
		_data[i] *= r_n;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline void Array<T, N>::print(void) const
{
	std::cout << "(";
	for(size_t i = 0; i < N - 1; i++){
		std::cout << _data[i] << ", ";
	}
	std::cout << _data[N - 1] << ")" << std::endl;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline T Array<T, N>::inner_product(const Array<T, N> &u) const
{
#ifdef DATA_DEBUG
	if(u.size() != _size) {
		std::cout << "wrong size." << std::endl;
	}
#endif
	T sum = 0;
	for(size_t i = 0; i < N; i++){
		sum += _data[i] * u[i];
	}
	return sum;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, size_t N>
inline Array<T, N> &Array<T, N>::operator = (const Array<T, N> &v)
{
	if(this != &v) {

#ifdef DATA_DEBUG
	_size = v._size;
#endif

		for(size_t i = 0; i < N; i++){
			_data[i] = v._data[i];
		}

	}

	return *this;
}


//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

