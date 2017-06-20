/*
 * SpecialMatrix.h
 *
 *  Created on: 2009-04-28
 *      Author: kamil
 */

#ifndef SPECIALMATRIX_H_
#define SPECIALMATRIX_H_

#include <vector>
#include <ostream>
#include <complex>
#include <iostream>
namespace VINHAC {

unsigned mypow(unsigned x,unsigned exp);


class SpecialMatrix {
private:
	std::complex<double> *data;
	unsigned *ranges;
	unsigned dim;
	unsigned total;

	unsigned multi2one(std::vector<unsigned> &multi) const;
	std::vector<unsigned> one2multi(unsigned one) const;


public:
	SpecialMatrix(std::vector<unsigned> &ranges);
	SpecialMatrix(const SpecialMatrix &another);
	std::complex<double> get(std::vector<unsigned> index) const;
	void set(std::vector<unsigned> index,std::complex<double> val);
	void set(unsigned index[],std::complex<double> val);
//	virtual ~SpecialMatrix();
	 ~SpecialMatrix();

	SpecialMatrix& operator=(const SpecialMatrix &another);

	inline void setRange(unsigned *newRange){ranges = newRange;}
	inline unsigned* getRange(){return ranges;}

	friend SpecialMatrix operator*(double lhs,SpecialMatrix &rhs);


	friend SpecialMatrix operator*(SpecialMatrix &lhs,double rhs);


	friend SpecialMatrix operator*(std::complex<double> lhs,SpecialMatrix &rhs);


	friend SpecialMatrix operator*(SpecialMatrix &lhs,std::complex<double> rhs);

	friend std::ostream& operator<<(std::ostream& out,SpecialMatrix &matrix);


	friend SpecialMatrix operator*(double lhs,const SpecialMatrix &rhs);


	friend SpecialMatrix operator*(const SpecialMatrix &lhs,double rhs);


	friend SpecialMatrix operator*(std::complex<double> lhs,const SpecialMatrix &rhs);


	friend SpecialMatrix operator*(const SpecialMatrix &lhs,std::complex<double> rhs);

	std::complex<double> operator*(SpecialMatrix &rhs);

	SpecialMatrix operator&(SpecialMatrix &rhs);


	SpecialMatrix operator+(SpecialMatrix &rhs);

	std::complex<double> operator*(const SpecialMatrix &rhs) const;

	SpecialMatrix operator&(const SpecialMatrix &rhs) const;


	SpecialMatrix operator+(const SpecialMatrix &rhs) const;

	std::complex<double> sumUp() const;
};



SpecialMatrix operator*(double lhs,SpecialMatrix &rhs);

SpecialMatrix operator*(SpecialMatrix &lhs,double rhs);

SpecialMatrix operator*(double lhs,const SpecialMatrix &rhs);

SpecialMatrix operator*(const SpecialMatrix &lhs,double rhs);


std::ostream& operator<<(std::ostream& out,SpecialMatrix &matrix);





}


#endif /* SPECIALMATRIX_H_ */
