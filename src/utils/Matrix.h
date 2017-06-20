/*
 * Matrix.h
 *
 *  Created on: 2009-03-30
 *      Author: kamil
 */

#ifndef MATRIX_H_
#define MATRIX_H_
#include <ostream>

namespace Vinhac {

template<class T,unsigned N, unsigned M=N>
class Matrix {
private:
	T elements[N][M];

public:
	Matrix();
	Matrix(const Matrix<T,N,M> &source);
	virtual ~Matrix();


	T* operator[](unsigned i);

	Matrix<T,N,M> operator+(Matrix<T,N,M> &another);
	Matrix<T,N,M> operator-(Matrix<T,N,M> &another);

	Matrix<T,M,N> transpose();

	Matrix<T,N,M> conjugate();



	template<class R,unsigned V,unsigned W>
	friend std::ostream& operator<<(std::ostream& out, Matrix<R,V,W> &matrix);

	template<class A, unsigned B, unsigned C, unsigned D>
	friend Matrix<A,B,D> operator*(Matrix<A,B,C> &m1, Matrix<A,C,D> &m2);


	friend class Matrix<T,M,N> ;


};


template<class T, unsigned N, unsigned M, unsigned K>
Matrix<T,N,K> operator*(Matrix<T,N,M> &m1, Matrix<T,M,K> &m2);




}

#include "Matrix.hpp"

#endif /* MATRIX_H_ */
