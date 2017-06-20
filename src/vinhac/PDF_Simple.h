/*
 * PDF_Simple.h
 *
 *  Created on: 15 Jul 2009
 *      Author: siudmy
 */

#ifndef PDF_SIMPLE_H_
#define PDF_SIMPLE_H_
#include "PDF.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#endif /* PDF_SIMPLE_H_ */

namespace VINHAC {
class PDF;
class PDF_Simple;
}

namespace VINHAC {
/**
 * \brief class handling simple unrealistic parametrization of PDF
 *
 * This is usefull for testing.
 */
class PDF_Simple: public PDF {

private:
	double m_Xmin;
	double m_Xmax;
	double m_Q2min;
	double m_Q2max;

public:
	PDF_Simple() {
		cout << "PDF simple constructor" << endl;
	}
	;
	~PDF_Simple() {
		cout << "SIMPLE PDF destructor" << endl;
	}
	;
	//! Initialization Single-set

	vector<double> xfx(double x, double q);

	double xfx(double x, double q, int fl);
	/*
	 public: virtual void initPDF() = 0;
	 */
	double getXmin(int m = 0) {
		return m_Xmin;
	}
	;
	double getXmax(int m = 0) {
		return m_Xmax;
	}
	;
	double getQ2min(int m = 0) {
		return m_Q2min;
	}
	;
	double getQ2max(int m = 0) {
		return m_Q2max;
	}
	;

	void setXmin(double Xmin) {
		m_Xmin = Xmin;
	}
	;
	void setXmax(double Xmax) {
		m_Xmax = Xmax;
	}
	;
	void setQ2min(double Q2min) {
		m_Q2min = Q2min;
	}
	void setQ2max(double Q2max) {
		m_Q2max = Q2max;
	}

};
}
