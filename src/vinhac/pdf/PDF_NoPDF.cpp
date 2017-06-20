/*
 * PDF_noPDF.cpp
 *
 *  Created on: 15 Jul 2009
 *      Author: siudmy
 */

#include "PDF_NoPDF.h"
#include "../core/VinhacException.h"

double VINHAC::PDF_NoPDF::xfx(double x, double q, int fl) {
	if (abs(fl) < 7) {
		return x;
	} else {
		throw VinhacException("PDF no PDF wrong flavour");
	}
	return x;
}
