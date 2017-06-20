/*
 * SpinorialFunction.h
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#ifndef SPINORIALFUNCTION_H_
#define SPINORIALFUNCTION_H_

#include <complex>
#include "../../../utils/CLHEP/Vector/Vector/LorentzVector.h"
using namespace std;

namespace VINHAC {

class SpinorialFunction {
public:
	SpinorialFunction();
	virtual ~SpinorialFunction();

	static complex<double> spinorialFunction(complex<double> hip1[2],
			const CLHEP::HepLorentzVector& a, complex<double> hip2[2],
			const int& ialfa);
	static complex<double> spinorialFunction(complex<double> hip1[2],
			const CLHEP::HepLorentzVector& a, const CLHEP::HepLorentzVector& b,
			const CLHEP::HepLorentzVector& c, complex<double> hip2[2],
			const int& ialfa);
};

}

#endif /* SPINORIALFUNCTION_H_ */
