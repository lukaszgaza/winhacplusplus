/*
 * SpecialMatrix.cpp
 *
 *  Created on: 2009-04-28
 *      Author: kamil
 */

#include "SpecialMatrix.h"

namespace VINHAC {

unsigned mypow(unsigned x, unsigned exp) {
	unsigned result = 1;
	for (unsigned i = 0; i < exp; i++)
		result *= x;
	return result;

}

SpecialMatrix::SpecialMatrix(const std::vector<unsigned> &ranges) {

	dim = ranges.size();
	this->ranges = new unsigned[dim];
	total = 1;
	for (unsigned i = 0; i < ranges.size(); ++i) {
		this->ranges[i] = ranges[i];
		total = total * this->ranges[i];
	}

	this->data = new std::complex<double>[total];

}

SpecialMatrix::SpecialMatrix(const SpecialMatrix &another) {
	dim = another.dim;
	ranges = new unsigned[dim];
	total = 1;

	for (unsigned i = 0; i < dim; ++i) {
		this->ranges[i] = another.ranges[i];
		total = total * this->ranges[i];
	}

	this->data = new std::complex<double>[total];

	for (unsigned i = 0; i < total; ++i) {
		data[i] = another.data[i];
	}

}

SpecialMatrix& SpecialMatrix::operator=(const SpecialMatrix &another) {

	dim = another.dim;
	delete[] ranges;
	ranges = new unsigned[dim];

	total = 1;

	for (unsigned i = 0; i < dim; ++i) {
		this->ranges[i] = another.ranges[i];
		total = total * this->ranges[i];
	}
	delete[] data;
	this->data = new std::complex<double>[total];

	for (unsigned i = 0; i < total; ++i) {
		data[i] = another.data[i];
	}

	return *this;

}

std::complex<double> SpecialMatrix::get(const std::vector<unsigned>& index) const {
	unsigned idx = 0;
	unsigned offset = 1;
	for (unsigned i = 0; i < dim && i < index.size(); ++i) {
		idx += offset * index[i];
		offset = offset * ranges[i];
	}
	return data[idx];
}

void SpecialMatrix::set(const std::vector<unsigned>& index,
		const std::complex<double>& val) {
	unsigned idx = 0;
	unsigned offset = 1;
	for (unsigned i = 0; i < dim && i < index.size(); ++i) {
		idx += offset * index[i];
		offset = offset * ranges[i];
	}
	data[idx] = val;
}

void SpecialMatrix::set(unsigned index[], const std::complex<double>& val) {
	unsigned idx = 0;
	unsigned tmp = 1;

	for (unsigned i = 0; i < dim; ++i) {
		idx += tmp * index[i];
		tmp = tmp * ranges[i];
	}
	//std::cout<<idx<<std::endl;
	//std::cout<<"["<<index[0]<<" "<<index[1]<<" "<<index[2]<<"] "<<idx<<std::endl;
	data[idx] = val;
}

SpecialMatrix::~SpecialMatrix() {
	delete[] data;
	delete[] ranges;
}

SpecialMatrix operator*(double lhs, SpecialMatrix &rhs) {
	SpecialMatrix result(rhs);

	for (unsigned i = 0; i < result.total; ++i) {
		result.data[i] = result.data[i] * lhs;
	}

	return result;
}

SpecialMatrix operator*(double lhs, const SpecialMatrix &rhs) {
	SpecialMatrix result(rhs);

	for (unsigned i = 0; i < result.total; ++i) {
		result.data[i] = result.data[i] * lhs;
	}

	return result;
}

SpecialMatrix operator*(SpecialMatrix &lhs, double rhs) {
	SpecialMatrix result(lhs);

	for (unsigned i = 0; i < result.total; ++i) {
		result.data[i] = result.data[i] * rhs;
	}

	return result;
}

SpecialMatrix operator*(const SpecialMatrix &lhs, double rhs) {
	SpecialMatrix result(lhs);

	for (unsigned i = 0; i < result.total; ++i) {
		result.data[i] = result.data[i] * rhs;
	}

	return result;
}

SpecialMatrix operator*(std::complex<double> lhs, SpecialMatrix &rhs) {
	SpecialMatrix result(rhs);

	for (unsigned i = 0; i < result.total; ++i) {
		result.data[i] = result.data[i] * lhs;
	}

	return result;
}

SpecialMatrix operator*(std::complex<double> lhs, const SpecialMatrix &rhs) {
	SpecialMatrix result(rhs);

	for (unsigned i = 0; i < result.total; ++i) {
		result.data[i] = result.data[i] * lhs;
	}

	return result;
}

SpecialMatrix operator*(SpecialMatrix &lhs, std::complex<double> rhs) {
	SpecialMatrix result(lhs);

	for (unsigned i = 0; i < result.total; ++i) {
		result.data[i] = result.data[i] * rhs;
	}

	return result;
}

SpecialMatrix operator*(const SpecialMatrix &lhs, std::complex<double> rhs) {
	SpecialMatrix result(lhs);

	for (unsigned i = 0; i < result.total; ++i) {
		result.data[i] = result.data[i] * rhs;
	}

	return result;
}

SpecialMatrix SpecialMatrix::operator+(SpecialMatrix &rhs) {
	if (dim != rhs.dim)
		throw "SpecialMatrix operator + : cannot add matrices whith not equal dimension";

	for (unsigned i = 0; i < dim; ++i) {
		if (ranges[i] != rhs.ranges[i])
			throw "SpecialMatrix operator + : not equal ranges in one of dimensions";
	}

	SpecialMatrix result(*this);

	for (unsigned i = 0; i < result.total; ++i)
		result.data[i] += rhs.data[i];

	return result;
}

SpecialMatrix SpecialMatrix::operator+(const SpecialMatrix &rhs) const {
	if (dim != rhs.dim)
		throw "SpecialMatrix operator + : cannot add matrices whith not equal dimension";

	for (unsigned i = 0; i < dim; ++i) {
		if (ranges[i] != rhs.ranges[i])
			throw "SpecialMatrix operator + : not equal ranges in one of dimensions";
	}

	SpecialMatrix result(*this);

	for (unsigned i = 0; i < result.total; ++i)
		result.data[i] += rhs.data[i];

	return result;
}

SpecialMatrix SpecialMatrix::operator&(SpecialMatrix &rhs) {

	std::vector<unsigned> rangestmp;
	rangestmp.reserve(dim + rhs.dim + 2);

	for (unsigned i = 0; i < 2 && i < dim - 1; ++i)
		rangestmp.push_back(ranges[i]);
	for (unsigned i = 1; i < 3 && i < rhs.dim; ++i)
		rangestmp.push_back(rhs.ranges[i]);
	for (unsigned i = 2; i < dim - 1; ++i)
		rangestmp.push_back(ranges[i]);
	for (unsigned i = 3; i < rhs.dim; ++i)
		rangestmp.push_back(rhs.ranges[i]);

	SpecialMatrix result(rangestmp);

	for (unsigned i = 0; i < result.total; ++i)
		result.data[i] = 0;

	for (unsigned i = 0; i < this->total; ++i)
		for (unsigned j = 0; j < rhs.total; ++j) {
			std::vector<unsigned> idx1 = this->one2multi(i);
			std::vector<unsigned> idx2 = rhs.one2multi(j);
			if (idx1[this->dim - 1] != idx2[0])
				continue;
			std::vector<unsigned> idx;
			idx.reserve(dim + rhs.dim + 2);
			for (unsigned k = 0; k < 2 && k < this->dim - 1; ++k)
				idx.push_back(idx1[k]);
			for (unsigned k = 1; k < 3 && k < rhs.dim; ++k)
				idx.push_back(idx2[k]);
			for (unsigned k = 2; k < this->dim - 1; ++k)
				idx.push_back(idx1[k]);
			for (unsigned k = 3; k < rhs.dim; ++k)
				idx.push_back(idx2[k]);

			result.data[result.multi2one(idx)] += this->data[i] * rhs.data[j];

		}

	return result;

}

SpecialMatrix SpecialMatrix::operator&(const SpecialMatrix &rhs) const {

	std::vector<unsigned> rangestmp;
	rangestmp.reserve(dim + rhs.dim + 2);

	for (unsigned i = 0; i < 2 && i < dim - 1; ++i)
		rangestmp.push_back(ranges[i]);
	for (unsigned i = 1; i < 3 && i < rhs.dim; ++i)
		rangestmp.push_back(rhs.ranges[i]);
	for (unsigned i = 2; i < dim - 1; ++i)
		rangestmp.push_back(ranges[i]);
	for (unsigned i = 3; i < rhs.dim; ++i)
		rangestmp.push_back(rhs.ranges[i]);

	SpecialMatrix result(rangestmp);

	for (unsigned i = 0; i < result.total; ++i)
		result.data[i] = 0;

	for (unsigned i = 0; i < this->total; ++i)
		for (unsigned j = 0; j < rhs.total; ++j) {
			std::vector<unsigned> idx1 = this->one2multi(i);
			std::vector<unsigned> idx2 = rhs.one2multi(j);
			if (idx1[this->dim - 1] != idx2[0])
				continue;
			std::vector<unsigned> idx;
			idx.reserve(dim + rhs.dim + 2);
			for (unsigned k = 0; k < 2 && k < this->dim - 1; ++k)
				idx.push_back(idx1[k]);
			for (unsigned k = 1; k < 3 && k < rhs.dim; ++k)
				idx.push_back(idx2[k]);
			for (unsigned k = 2; k < this->dim - 1; ++k)
				idx.push_back(idx1[k]);
			for (unsigned k = 3; k < rhs.dim; ++k)
				idx.push_back(idx2[k]);

			result.data[result.multi2one(idx)] += this->data[i] * rhs.data[j];

		}

	return result;

}

std::complex<double> SpecialMatrix::sumUp() const {
	std::complex<double> result(0);

	for (unsigned i = 0; i < this->total; ++i)
		result += data[i] * std::conj(data[i]);

	return result;
}

unsigned SpecialMatrix::multi2one(const std::vector<unsigned> &multi) const {
	unsigned idx = 0;
	unsigned offset = 1;
	for (unsigned i = 0; i < dim && i < multi.size(); ++i) {
		idx += offset * multi[i];
		offset = offset * ranges[i];
	}
	return idx;
}

std::vector<unsigned> SpecialMatrix::one2multi(unsigned one) const {
	std::vector<unsigned> idx(dim);

	unsigned tmp = 0;
	unsigned offset = 1;

	for (unsigned i = 0; i < dim; ++i) {
		offset = offset * ranges[i];
		tmp = one % offset;
		idx[i] = (tmp / (offset / ranges[i]));
		one = one - tmp;
	}

	return idx;
}

double SpecialMatrix::operator*(const SpecialMatrix &rhs) const {
	std::set<PolarizationName::name> polarization;
	polarization.insert(PolarizationName::unpolarized);
	return this->multiplyAndSquare(rhs, polarization)[PolarizationName::unpolarized];
}

std::map<PolarizationName::name, double> SpecialMatrix::multiplyAndSquare(
		const SpecialMatrix &rhs,
		const std::set<PolarizationName::name> &polarization) const {
	if (this->ranges[dim - 1] != rhs.ranges[0])
		throw "SpecialMatrix::multiplyAndSquare: ranges at lambda aren't equal";
	if (this->ranges[dim - 1] != 3)
		throw "SpecialMatrix::multiplyAndSquare: calculating W polarizations with range of lambda other than 1-3 does not make sence";
	if (polarization.size() == 0)
		throw "SpecialMatrix::multiplyAndSquare: passing zero set of polarizations does not make sense";

	std::map<PolarizationName::name, double> result;
	for (std::set<PolarizationName::name>::iterator it = polarization.begin(); it
			!= polarization.end(); ++it) {
		result[*it] = 0;
	}
	std::complex<double> tmp(0);
	std::complex<double> tmp_arr[3];

	unsigned offset1 = total / ranges[dim - 1];
	unsigned offset2 = rhs.total / ranges[dim - 1];

	for (unsigned k = 0; k < offset1; ++k)
		for (unsigned j = 0; j < offset2; ++j) {
			if (polarization.find(PolarizationName::unpolarized)
					!= polarization.end()) {
				tmp = 0;
				for (unsigned i = 0; i < ranges[dim - 1]; ++i) {
					tmp += data[i * offset1 + k] * rhs.data[j * ranges[dim - 1]
							+ i];
				}
				result[PolarizationName::unpolarized]
						+= (tmp * std::conj(tmp)).real();
			}
			if (polarization.find(PolarizationName::longitudinally)
					!= polarization.end()) {
				unsigned i = 2;
				tmp = data[i * offset1 + k] * rhs.data[j * ranges[dim - 1] + i];
				result[PolarizationName::longitudinally] += (tmp * std::conj(
						tmp)).real();
			}
			if (polarization.find(PolarizationName::transversely)
					!= polarization.end()) {
				tmp = 0;
				for (unsigned i = 0; i < 2; ++i) {
					tmp += data[i * offset1 + k] * rhs.data[j * ranges[dim - 1]
							+ i];
				}
				result[PolarizationName::transversely]
						+= (tmp * std::conj(tmp)).real();
			}
			if (polarization.find(PolarizationName::left) != polarization.end()) {
				tmp = (data[0 * offset1 + k] + std::complex<double>(0, 1)
						* data[1 * offset1 + k]) * (rhs.data[j
						* ranges[dim - 1] + 0] - std::complex<double>(0, 1)
						* rhs.data[j * ranges[dim - 1] + 1]) / 2.0;
				result[PolarizationName::left] += (tmp * std::conj(tmp)).real();
			}
			if (polarization.find(PolarizationName::right)
					!= polarization.end()) {
				tmp = (data[0 * offset1 + k] - std::complex<double>(0, 1)
						* data[1 * offset1 + k]) * (rhs.data[j
						* ranges[dim - 1] + 0] + std::complex<double>(0, 1)
						* rhs.data[j * ranges[dim - 1] + 1]) / 2.0;
				result[PolarizationName::right]
						+= (tmp * std::conj(tmp)).real();
			}
			if (polarization.find(
					PolarizationName::int_long_tran)
					!= polarization.end()) {
				for (unsigned i = 0; i < ranges[dim - 1]; ++i) {
					tmp_arr[i] = data[i * offset1 + k] * rhs.data[j
							* ranges[dim - 1] + i];
				}
				result[PolarizationName::int_long_tran]
						+= 2.0 * ((tmp_arr[0] + tmp_arr[1]) * std::conj(
								tmp_arr[2])).real();
			}
			if (polarization.find(
					PolarizationName::int_left_right)
					!= polarization.end()) {
				tmp_arr[0] = (data[0 * offset1 + k] + std::complex<double>(0, 1)
										* data[1 * offset1 + k]) * (rhs.data[j
										* ranges[dim - 1] + 0] - std::complex<double>(0, 1)
										* rhs.data[j * ranges[dim - 1] + 1]) / 2.0; // Left
				tmp_arr[1] = (data[0 * offset1 + k] - std::complex<double>(0, 1)
										* data[1 * offset1 + k]) * (rhs.data[j
										* ranges[dim - 1] + 0] + std::complex<double>(0, 1)
										* rhs.data[j * ranges[dim - 1] + 1]) / 2.0; // Right
				result[PolarizationName::int_left_right]
						+= 2.0 * (tmp_arr[0] * std::conj(tmp_arr[1])).real();
			}
		}

	return result;
}

std::ostream& operator<<(std::ostream& out, SpecialMatrix &matrix) {
	for (unsigned i = 0; i < matrix.total; ++i) {
		std::cout << "one dim idx: " << i << " , multidim idx:";
		std::vector<unsigned> idx = matrix.one2multi(i);
		for (unsigned j = 0; j < idx.size(); ++j)
			std::cout << " " << idx[j];
		std::cout << " , value: " << matrix.data[i];

		std::cout << std::endl;

	}

	return out;
}

}

