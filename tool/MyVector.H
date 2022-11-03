/*=============================================================================
*   
*   Filename : MyVector.H
*   Creator : Han Zhou
*   Date : 06/15/21
*   Description : 
*
=============================================================================*/
 
#ifndef _MYVECTOR_H
#define _MYVECTOR_H

#include <iostream>

template<typename T>
class Vector
{
	public :
		Vector(void);
		Vector(size_t n);
		Vector(double low, double high, size_t n);
		Vector(size_t n, const T&);
		Vector(const Vector &v);
		virtual ~Vector(void);

	public :

		T operator[](size_t) const;
		T& operator[](size_t);
		Vector &operator = (const Vector &v);

		size_t size(void) const;

		void allocate(size_t n);
		void reallocate(size_t n);
		void setAxis(double low, double high, size_t n);
		void fill(const T &a);
		void clear(void);

#ifdef DATA_DEBUG
		void print(void) const;
#endif

		T max(void) const;
		T min(void) const;
		T max_norm(void) const;
		T l1_norm(void) const;
		T l2_norm(void) const;
		T l1_norm_scaled(void) const;
		T l2_norm_scaled(void) const;

		T max_error(const Vector &u) const;
		T l1_error_scaled(const Vector &u) const;
		T l2_error_scaled(const Vector &u) const;

		T inner_product(const Vector &u) const;

		void normalize(void);

	private :

		T *_data;
		size_t _size;
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(void):_size(0), _data(0)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(size_t n):_size(n)
{
	_data = new T[_size];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(double low, double high, size_t n)
{
#ifdef DATA_DEBUG
	if (low > high) {
		std::cout << "invalid interval." << std::endl;
		exit(1);
	}
#endif

	_size = n + 1;
	_data = new T[_size];
	double h = (high - low) / n;
	_data[0] = low;
	for(size_t i = 1; i <= n; i++){
		_data[i] = _data[i - 1] + h;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(size_t n, const T &a):_size(n)
{
	_data = new T[_size];

	for(size_t i = 0; i < n; i++){
		_data[i] = a;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(const Vector<T> &v):_size(v._size)
{
	allocate(v._size);

	for(size_t i = 0; i < v._size; i++){
		_data[i] = v._data[i];
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::allocate(size_t n)
{
	_size = n;
	_data = new T[_size];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::reallocate(size_t n)
{
	if(0 != _data) {
		delete[] _data;
	}

	_size = n;
	_data = new T[_size];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::setAxis(double low, double high, size_t n)
{
	if(0 != _data) {
		delete[] _data;
	}

	_size = n + 1;
	_data = new T[_size];

	double ds = (high - low) / n;

	_data[0] = low;
	for(size_t i = 1; i <= n; i++){
		_data[i] = _data[i-1] + ds;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::~Vector(void)
{
	if(0 != _data) {
		delete[] _data;
		_data = 0;
	}
	_size = 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T> &Vector<T>::operator = (const Vector<T> &v)
{
	if(this != &v) {

		reallocate(v._size);
		for(size_t i = 0; i < _size; i++){
			_data[i] = v._data[i];
		}

	}

	return *this;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::clear(void)
{
	if(0 != _data) {
		delete[] _data;
		_data = 0;
	}
	_size = 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::fill(const T &a)
{
	for(size_t i = 0; i < _size; i++){
		_data[i] = a;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline size_t Vector<T>::size(void) const
{
	return _size;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::operator[](size_t i) const
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
template<typename T>
inline T& Vector<T>::operator[](size_t i)
{
#ifdef DATA_DEBUG
	if(i >= _size) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}
#endif
	return _data[i];
}

#ifdef DATA_DEBUG
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::print(void) const
{

	if (std::is_same<T, double>::value || 
			std::is_same<T, float>::value || 
			std::is_same<T, long double>::value) {

  	for (size_t i = 0; i < _size; i++) {
  	  if ((i > 0) && (i%8 == 0)) {
  	    printf("\n"); 
  	  }
  	  printf("%9.6f ", _data[i]); 
  	}
  	printf("\n");

	} else if (std::is_same<T, int>::value || 
						 std::is_same<T, unsigned int>::value || 
						 std::is_same<T, long unsigned int>::value) {

  	for (size_t i = 0; i < _size; i++) {
  	  if ((i > 0) && (i%16 == 0)) {
  	    printf("\n"); 
  	  }
  	  printf("%4d ", _data[i]); 
  	}
  	printf("\n"); 

	}
}
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::max(void) const
{
	T max_u = _data[0];
	for(size_t i = 1; i < _size; i++){
		max_u = max_u > _data[i] ? max_u : _data[i];
	}
	return max_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::min(void) const
{
	T min_u = _data[0];
	for(size_t i = 1; i < _size; i++){
		min_u = min_u < _data[i] ? min_u : _data[i];
	}
	return min_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::max_norm(void) const
{
	T max_u = fabs(_data[0]);
	for(size_t i = 1; i < _size; i++){
		T u = fabs(_data[i]);
		max_u = max_u > u ? max_u : u;
	}
	return max_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l1_norm(void) const
{
	T sum = 0.0;
	for(size_t i = 0; i < _size; i++){
		sum += fabs(_data[i]);
	}
	return sum;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l2_norm(void) const
{
	T sum = 0.0;
	for(size_t i = 0; i < _size; i++){
		sum += _data[i] * _data[i];
	}
	return sqrt(sum);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l1_norm_scaled(void) const
{
	T sum = 0.0;
	for(size_t i = 0; i < _size; i++){
		sum += fabs(_data[i]);
	}
	return sum / _size;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l2_norm_scaled(void) const
{
	T sum = 0.0;
	for(size_t i = 0; i < _size; i++){
		sum += _data[i] * _data[i];
	}
	return sqrt(sum / _size);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::normalize(void)
{
	T sum = 0.0;
	for(size_t i = 0; i < _size; i++){
		sum += _data[i] * _data[i];
	}

#ifdef DATA_DEBUG
	if(sum < 1.0E-15) {
		std::cout << "zero vector." << std::endl;
		exit(1);
	}
#endif

	T r_norm = 1.0 / sqrt(sum);

	for(size_t i = 0; i < _size; i++){
		_data[i] *= r_norm;
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::max_error(const Vector<T> &u) const
{
#ifdef DATA_DEBUG
	if(_size != u.size()) {
		std::cout << "wrong size. " << std::endl;
		exit(1);
	}
#endif

	T max_e = fabs(_data[0] - u[0]);
	for(size_t i = 1; i < _size; i++){
		T e = fabs(_data[i] - u[i]);
		max_e = max_e > e ? max_e : e;
	}
	return max_e;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l1_error_scaled(const Vector<T> &u) const
{
#ifdef DATA_DEBUG
	if(_size != u.size()) {
		std::cout << "wrong size. " << std::endl;
		exit(1);
	}
#endif

	T sum = 0.0;
	for(size_t i = 0; i < _size; i++){
		sum += fabs(_data[i] - u[i]);
	}
	return sum / _size;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l2_error_scaled(const Vector<T> &u) const
{
#ifdef DATA_DEBUG
	if(_size != u.size()) {
		std::cout << "wrong size. " << std::endl;
		exit(1);
	}
#endif

	T sum = 0.0;
	for(size_t i = 0; i < _size; i++){
		T err = _data[i] - u[i];
		sum += err * err;
	}
	return sqrt(sum / _size);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template< typename T >
inline T Vector< T >::inner_product(const Vector<T> &u) const
{
	T sum = 0;
	for(size_t i = 0; i < _size; i++){
		sum += u[i] * _data[i];
	}
	return sum;
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif
