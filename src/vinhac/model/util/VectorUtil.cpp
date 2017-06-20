/*
 * VectorUtil.cpp
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#include "VectorUtil.h"
#include "../../core/VinhacException.h"

namespace VINHAC {

VectorUtil::VectorUtil() {

}

VectorUtil::~VectorUtil() {

}

void VINHAC::VectorUtil::VecPol(CLHEP::HepLorentzVector& epsilonW,
		const CLHEP::HepLorentzVector& p, const int& lambda) {
	///////////////////////////////////////////////////////////////////////////////////////
	//  Calculation of polarization lambda state vectors of a vector boson in the
	//  rectangular basis, see K. Hagiwara and D. Zeppenfeld,
	//  Nucl. Phys. B274 (1986) 1, eq. (3.47).
	//  INPUT:
	//     lambda   - polarization component i.e. 1,2 or 3
	//  OUTPUT:
	//     eps -  lambda polarization state i.e. 4-vector
	//-------------------------------------------------------------------------------------
	//  Written by: Andrzej Siodmok                                       date: 25.04.2005
	//  Last update: 25.04.2005                                           by: A.S.
	///////////////////////////////////////////////////////////////////////////////////////

	epsilonW.set(0);
	double qt2, qt, aq2, aq, am, ws3;
	qt2 = p.perp2(); // x^2 + y^2
	qt = sqrt(qt2); // sqrt(x^2 + y^2)
	aq2 = p.z();
	aq2 *= aq2;
	aq2 += qt2;
	aq = sqrt(aq2); // x^2 + y^2 + z^2
	am = sqrt(fabs(p.e() * p.e() - aq2));// (*this).restMass();

	// lambda polarization state = 1
	if (lambda == 1) {
		if (aq < 1e-8) { // <--- small Rho() approximation
			epsilonW.setX(1.0);
			// cout << "test_male aq" << endl;
		} else if (qt < 1e-8) { // <---  small Perp() approximation
			epsilonW.setX(p.z() / aq);
			//cout << "test male qt" << endl;
		} else { // eq. (3.47a)
			epsilonW.setX(p.x() / qt * p.z() / aq);
			epsilonW.setY(p.y() / qt * p.z() / aq);
			epsilonW.setZ(-qt / aq);
			//epsilonW.SetT(0.); nie trzeba bo na poczatku wyzerowany
		}
	}
	// lambda polarization state = 2
	else if (lambda == 2) {
		if (aq < 1e-8) { // <--- small Rho() approximation
			epsilonW.setY(1e0);
		} else if (qt < 1e-8) { // <--- small Perp() approximation
			epsilonW.setY(1.0);
		} else { // eq. (3.47b)
			epsilonW.setX(-p.y() / qt);
			epsilonW.setY(p.x() / qt);
			epsilonW.setZ(0.);
			epsilonW.setT(0.);
		}
	}
	// lambda polarization state = 3
	else if (lambda == 3) {
		if (aq < 1e-8) {
			if (am > 1e-5) {
				epsilonW.setZ(1.);
				//cout << "test1" << endl;
			}
		} else if (qt < 1e-8) {
			if (am > 1e-5) {
				epsilonW.setZ(p.t() / am / aq * p.z());
				epsilonW.setT(aq / am);
				//cout << "test2" << endl;
			}
		} else {
			if (am > 1e-5) {
				ws3 = p.t();
				ws3 /= am;
				ws3 /= aq;
				epsilonW.setX(ws3 * p.x());
				epsilonW.setY(ws3 * p.y());
				epsilonW.setZ(ws3 * p.z());
				epsilonW.setT(ws3 * aq2 / p.t());
				//cout << "test3" << endl;
			}
		}
	} else {
		throw VinhacException("Zla wartosc lambda");

	}
}

// Additional heg vectors stuff
void VINHAC::VectorUtil::HelEig(const CLHEP::HepLorentzVector& p,
		const int& hel, complex<double> chip[2]) {
	///////////////////////////////////////////////////////////////////////////////////////
	//  This function provides a helicity eigenstate (Pauli spinor) for
	//  a fermion of helicity ihel and 4-momentum p.
	//  INPUT:
	//     hel, p - helicity of a fermion
	//  OUTPUT:
	//     chip - 2-component (complex) Pauli spinor
	//-------------------------------------------------------------------------------------
	//  Written by: Andrzej Siodmok                                       date: 25.04.2005
	//  Last update: 25.04.2005                                           by: A.S.
	///////////////////////////////////////////////////////////////////////////////////////

	double pmp3, fnor, temp1, temp2;
	pmp3 = p.rho();// <--- Rho() = sqrt ( x^2 + y^2 + z^2)
	pmp3 += p.z();
	// helicity of fermion = -1
	if (hel == -1) {
		if (pmp3 > 1e-10) {
			fnor = 1. / sqrt(2 * p.rho() * pmp3); // norm of Pauli spinor from eq.(3.22b)
			temp1 = -fnor * p.x();
			temp2 = fnor * p.y();
			chip[0] = complex<double> (temp1, temp2);
			temp1 = fnor * pmp3;
			temp2 = 0;
			chip[1] = complex<double> (temp1, temp2);
		} else { // when p_z = -|\vec{p}|
			chip[0] = complex<double> (-1., 0.);
			chip[1] = complex<double> (0., 0.);
		}
	} else if (hel == 1) {
		// helicity of fermion = 1
		if (pmp3 > 1e-10) {
			fnor = 1. / sqrt(2 * p.rho() * pmp3);
			temp1 = fnor * pmp3;
			temp2 = 0;
			chip[0] = complex<double> (temp1, temp2);
			temp1 = fnor * p.x();
			temp2 = fnor * p.y();
			chip[1] = complex<double> (temp1, temp2); //eq. (3.22b)
		} else { // when p_z = -|\vec{p}|
			chip[0] = complex<double> (0., 0.);
			chip[1] = complex<double> (1., 0.);
		}
	} else {
		throw VinhacException(" heleig: Wrong fermion helicity");
	}
}



void VINHAC::VectorUtil::MxV2dC(complex<double> matrix[2][2],
		complex<double> vec[2], complex<double> res[2]) {

	///////////////////////////////////////////////////////////////////////////////////////
	// Multiply 2-dim. complex matrix by complex vector.
	// INPUT:
	//    matrix, vec - 2x2 matrix and 2-component vector, resp.
	// OUTPUT:
	//    res - resulting 2-component complex vector
	//-------------------------------------------------------------------------------------
	//  Written by: Andrzej Siodmok                                       date: 25.04.2005
	//  Last update: 25.04.2005                                           by: A.S.
	///////////////////////////////////////////////////////////////////////////////////////
	complex<double> temp[2];
	temp[0] = vec[0];
	temp[1] = vec[1];
	res[0] = matrix[0][0] * temp[0] + matrix[0][1] * temp[1];
	res[1] = matrix[1][0] * temp[0] + matrix[1][1] * temp[1];
}

void VINHAC::VectorUtil::MvtoWm(const CLHEP::HepLorentzVector& p, complex<
		double> as[2][2], const int& alpha) {

	///////////////////////////////////////////////////////////////////////////////////////
	//  Translate Minkowski 4-vector into 2x2 complex Matrix in Weyl basis.
	//  INPUT:
	//     alpha - index alpha=-/+ of a projetion operator in Weyl basis
	//  OUTPUT:
	//     as - 2x2 complex matrix (a-slash)_alpha eq.3.11
	//-------------------------------------------------------------------------------------
	//  Written by: Andrzej Siodmok                                       date: 25.04.2005
	//  Last update: 25.04.2005                                           by: A.S.
	///////////////////////////////////////////////////////////////////////////////////////

	as[0][0] = complex<double> ((p.t() - alpha * p.z()), 0.);
	as[0][1] = complex<double> ((-alpha * p.x()), (alpha * p.y()));
	as[1][0] = complex<double> ((-alpha * p.x()), (-alpha * p.y()));
	as[1][1] = complex<double> ((p.t() + alpha * p.z()), 0.);
}


}
