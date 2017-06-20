#ifndef __PDF_LHA_h__
#define __PDF_LHA_h__

#include <string>
#include <vector>
#include <exception>
#include "LHAPDF/LHAPDF.h"
#include "PDF.h"
using namespace std;

namespace VINHAC {
class PDF;
class PDF_LHA;
}

namespace VINHAC {
/**
 * \brief class handling Parton Distribution Functions from package LHAPDF
 */
class PDF_LHA: public PDF {
public:
	PDF_LHA();
	~PDF_LHA();

public:
	//! Initialization Single-set
	void initPDF();

	LHAPDF::PDFSetInfo initPDFByName(const string &name, LHAPDF::SetType type, int memset);

	void initPDFByName(const string &filename, int memset);

	void initPDFSet(const string &path);

	void initPDFSetByName(const string &name, LHAPDF::SetType type);

	void initPDFSetByName(const string &filename);
	//! Initialization Mulitple

	void initPDFByNameM(int nset, const string &name, LHAPDF::SetType type,
			int memset);

	void initPDFByNameM(int nset, const string &filename, int memset);

	void initPDFSetM(int nset, const string &path);

	void initPDFSetByNameM(int nset, const string &name, LHAPDF::SetType type);

	void initPDFSetByNameM(int nset, const string &filename);

	void setPDFPath(const string &path);

	//! Single-set functions

	void initPDF(int memset);

	void getDescription();

	vector<double> xfx(double x, double Q);

	void xfx(double x, double Q, double *results);

	double xfx(double x, double Q, int fl);
	// Nuclear: Extra a param for atomic mass number.

	vector<double> xfxa(double x, double Q, double a);

	void xfxa(double x, double Q, double a, double *results);

	double xfxa(double x, double Q, double a, int fl);
	// Other

	int numberPDF();

	double alphasPDF(double Q);

	int getOrderPDF(int nset);

	int getOrderAlphaS(int nset);

	double getQMass(int f);

	double getThreshold(int f);

	double getLam4(int m);

	double getLam5(int m);

	double getXmin(int m);

	double getXmax(int m);

	double getQ2min(int m);

	double getQ2max(int m);
};
}

#endif
