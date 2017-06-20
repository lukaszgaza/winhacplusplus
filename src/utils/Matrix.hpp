/*
 * Matrix.hpp
 *
 *  Created on: 2009-03-30
 *      Author: kamil
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <ostream>
#include <complex>


namespace Vinhac {

template<class T,unsigned N,unsigned M>
Matrix<T,N,M>::Matrix(){
	if(N==0 || M==0) throw "0 - size of matrix forbidden";
	for(unsigned i=0;i<N;i++)
		for(unsigned j=0;j<M;j++) elements[i][j]=0.0;
}


template<class T,unsigned N, unsigned M>
Matrix<T,N,M>::~Matrix(){
}


template<class T,unsigned N,unsigned M>
Matrix<T,N,M>::Matrix(const Matrix<T,N,M> &source){

	for(unsigned i=0; i<N ;i++)
		for(unsigned j=0; j<M ; j++)
			this->elements[i][j] = source.elements[i][j];


}

template<class T,unsigned N,unsigned M>
Matrix<T,N,M> Matrix<T,N,M>::operator+(Matrix<T,N,M> &another){
	Matrix<T,N,M> result;

	for(unsigned i=0; i<N ;i++)
			for(unsigned j=0; j<M ; j++)
				result.elements[i][j] = this->elements[i][j] + another.elements[i][j];

	return result;
}




template<class T,unsigned N,unsigned M>
Matrix<T,N,M> Matrix<T,N,M>::operator-(Matrix<T,N,M> &another){
	Matrix<T,N,M> result;

	for(unsigned i=0; i<N ;i++)
			for(unsigned j=0; j<M ; j++)
				result.elements[i][j] = this->elements[i][j] - another.elements[i][j];

	return result;
}


template<class T,unsigned N,unsigned M>
Matrix<T,N,M> Matrix<T,N,M>::conjugate(){
	Matrix<T,N,M> result;

		for(unsigned i=0; i<N ;i++)
				for(unsigned j=0; j<M ; j++)
					result.elements[i][j] = std::conj(this->elements[i][j]);

		return result;
}


template<class T,unsigned N,unsigned M>
Matrix<T,M,N> Matrix<T,N,M>::transpose(){

	Matrix<T,M,N> result;

	for(unsigned i=0; i<N ;i++)
				for(unsigned j=0; j<M ; j++)
					result.elements[j][i] = this->elements[i][j];

	return result;
}



template<class T, unsigned N, unsigned M, unsigned K>
Matrix<T,N,K> operator*(Matrix<T,N,M> &m1, Matrix<T,M,K> &m2){
	Matrix<T,N,K> result;

	for(unsigned i=0;i<N;i++)
		for(unsigned j=0;j<K;j++){
			result.elements[i][j]=0;

			for(unsigned k=0;k<M;k++) result.elements[i][j] += m1.elements[i][k]*m2.elements[k][j];


		}
	return result;

}


template<class T,unsigned N, unsigned M>
T* Matrix<T,N,M>::operator[](unsigned i){
	return elements[i];
}


template<class T,unsigned N, unsigned M>
std::ostream& operator<<(std::ostream& out, Matrix<T,N,M>& matrix){
	for(unsigned i=0; i<N ; i++){
		for(unsigned j=0;j<M;j++) out<<"["<<matrix.elements[i][j]<<"] ";
		out<<std::endl;
	}
	return out;
}

}


#endif /* MATRIX_HPP_ */
