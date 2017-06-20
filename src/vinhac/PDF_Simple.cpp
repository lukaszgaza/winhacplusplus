/*
 * PDF_Simple.cpp
 *
 *  Created on: 15 Jul 2009
 *      Author: siudmy
 */

#include "PDF_Simple.h"
#include "VinhacException.h"

std::vector<double> VINHAC::PDF_Simple::xfx(double x, double q) {
	std::vector<double> result;
	for (int i = 0; i <= 12; ++i) {
		result.push_back(xfx(x, q, i - 6));
	}

	return result;
}

double VINHAC::PDF_Simple::xfx(double x, double q, int fl) {
	if (abs(fl) < 7) {

		if (fl != 1 && fl != 2 && abs(fl) != 6) {
			return 1.0 / x / 1.0e4;
		} else if (fl == 1)
			return 1.0 / x / 1.0e4;
		else if (fl == 2)
			return 2.0 / x / 1.0e4;
		else
			return 0;
	} else {
		throw VinhacException("PDF SIMPLE wrong flavour");
	}

}
