#include <string>
#include <vector>
#include <exception>
#include <iostream>
#include "PDF_LHA.h"

using namespace std;
//using namespace LHAPDF;


VINHAC::PDF_LHA::PDF_LHA() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::PDF_LHA()"<<std::endl;

#endif
	//throw "Not yet implemented";
}

VINHAC::PDF_LHA::~PDF_LHA() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::~PDF_LHA()"<<std::endl;

#endif

}
// Initialization Single-set
void VINHAC::PDF_LHA::initPDFByName(const string &name, LHAPDF::SetType type,
		int memset) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::initPDFByName() 1"<<std::endl;

#endif
	LHAPDF::initPDFByName(name, type, memset);
}

void VINHAC::PDF_LHA::initPDFByName(const string &filename, int memset) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::initPDFByName() 2"<<std::endl;

#endif
	LHAPDF::initPDFByName(filename, memset);
}
void VINHAC::PDF_LHA::initPDFSet(const string &path) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::initPDFSet()"<<std::endl;

#endif
	LHAPDF::initPDFSet(path);
}

void VINHAC::PDF_LHA::initPDFSetByName(const string &name, LHAPDF::SetType type) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::initPDFSetByName() 1"<<std::endl;

#endif
	LHAPDF::initPDFSetByName(name, type);
}

void VINHAC::PDF_LHA::initPDFSetByName(const string &filename) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::initPDFSetByName() 2"<<std::endl;

#endif
	LHAPDF::initPDFSetByName(filename);
}
//! Initialization - Multiple
void VINHAC::PDF_LHA::initPDFByNameM(int nset, const string &name,
		LHAPDF::SetType type, int memset) {
	LHAPDF::initPDFByNameM(nset, name, type, memset);
}
void VINHAC::PDF_LHA::initPDFByNameM(int nset, const string &filename,
		int memset) {
	LHAPDF::initPDFByNameM(nset, filename, memset);
}

void VINHAC::PDF_LHA::initPDFSetM(int nset, const string &path) {
	LHAPDF::initPDFSetM(nset, path);
}

void VINHAC::PDF_LHA::initPDFSetByNameM(int nset, const string &filename) {
	LHAPDF::initPDFSetByNameM(nset, filename);
}
void VINHAC::PDF_LHA::setPDFPath(const string &path) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::setPDFPath()"<<std::endl;

#endif
	LHAPDF::setPDFPath(path);
}
//Single-set functions
void VINHAC::PDF_LHA::initPDF(int memset) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::initPDF()"<<std::endl;

#endif
	LHAPDF::initPDF(memset);
}
void VINHAC::PDF_LHA::getDescription() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: PDF_LHA::getDescription()"<<std::endl;

#endif
	LHAPDF::getDescription();
}

vector<double> VINHAC::PDF_LHA::xfx(double x, double Q) {
	return LHAPDF::xfx(x, Q);
}
void VINHAC::PDF_LHA::xfx(double x, double Q, double *results) {
	return LHAPDF::xfx(x, Q, results);
}
double VINHAC::PDF_LHA::xfx(double x, double Q, int fl) {
	return LHAPDF::xfx(x, Q, fl);
}
// Nuclear: Extra a param for atomic mass number.
vector<double> VINHAC::PDF_LHA::xfxa(double x, double Q, double a) {
	return LHAPDF::xfxa(x, Q, a);
}
void VINHAC::PDF_LHA::xfxa(double x, double Q, double a, double *results) {
	return LHAPDF::xfxa(x, Q, a, results);
}
double VINHAC::PDF_LHA::xfxa(double x, double Q, double a, int fl) {
	return LHAPDF::xfxa(x, Q, a, fl);
}

//Other
int VINHAC::PDF_LHA::numberPDF() {
	return LHAPDF::numberPDF();
}
double VINHAC::PDF_LHA::alphasPDF(double Q) {
	return LHAPDF::alphasPDF(Q);
}
int VINHAC::PDF_LHA::getOrderPDF(int nset) {
	return LHAPDF::getOrderPDFM(nset);
}

int VINHAC::PDF_LHA::getOrderAlphaS(int nset) {
	return LHAPDF::getOrderAlphaSM(nset);
}

double VINHAC::PDF_LHA::getQMass(int f) {
	return LHAPDF::getQMass(f);
}
double VINHAC::PDF_LHA::getThreshold(int f) {
	return LHAPDF::getThreshold(f);
}
double VINHAC::PDF_LHA::getLam4(int m) {
	return LHAPDF::getLam4(m);
}
double VINHAC::PDF_LHA::getLam5(int m) {
	return LHAPDF::getLam5(m);
}

double VINHAC::PDF_LHA::getXmin(int m) {
	LHAPDF::getXmin(m);
}

double VINHAC::PDF_LHA::getXmax(int m) {
	return LHAPDF::getXmax(m);
}

double VINHAC::PDF_LHA::getQ2min(int m) {
	return LHAPDF::getQ2min(m);
}

double VINHAC::PDF_LHA::getQ2max(int m) {
	return LHAPDF::getQ2max(m);
}

