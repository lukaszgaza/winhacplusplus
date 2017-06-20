/*
 * SpinorialFunction.cpp
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#include "SpinorialFunction.h"

#include "VectorUtil.h"

namespace VINHAC {

SpinorialFunction::SpinorialFunction() {

}

SpinorialFunction::~SpinorialFunction() {

}

complex<double> VINHAC::SpinorialFunction::spinorialFunction(
		complex<double> hip1[2], const CLHEP::HepLorentzVector& a, complex<
				double> hip2[2], const int& ialfa) {

	complex<double> as[2][2];
	complex<double> vx[2];

	VectorUtil::MvtoWm(a, as, ialfa);

	VectorUtil::MxV2dC(as, hip2, vx);

	return conj(hip1[0]) * vx[0] + conj(hip1[1]) * vx[1];
}

complex<double> VINHAC::SpinorialFunction::spinorialFunction(
		complex<double> hip1[2], const CLHEP::HepLorentzVector& a,
		const CLHEP::HepLorentzVector& b, const CLHEP::HepLorentzVector& c,
		complex<double> hip2[2], const int& ialfa) {
	//!----------------------------------------------------------------------!
	//! This function provides a value of the spinorial function             !
	//!           S(p1,a1,,a2,a3,p2)_{lambda1,lambda2}^alpha                 !
	//! given in Ref: K. Hagiwara & D. Zeppenfeld, Nucl.Phys. B274 (1986) 1. !
	//! INPUT: chip1,chip2 - Pauli spinors of fermions: f1(p1,lambda1),      !
	//!                                                 f2(p2,lambda2)       !
	//!        a1,a2,a3 - 4-vectors in Minkowski space                       !
	//!        ialf - index alpha=-/+ of a projetion operator in Weyl basis  !
	//!----------------------------------------------------------------------!
	//! Written by: Wieslaw Placzek            date: 19.06.1997              !
	//! Last update: 30.08.2002                by: W.P.                      !
	//!----------------------------------------------------------------------!
	complex<double> as[2][2];
	complex<double> bs[2][2];
	complex<double> cs[2][2];
	complex<double> vx[2];
	//! Translate Minkowki 4-vectors into 2x2 complex matrices in Weyl basis
	VectorUtil::MvtoWm(a, as, ialfa);
	VectorUtil::MvtoWm(b, bs, -ialfa);
	VectorUtil::MvtoWm(c, cs, ialfa);
	//! Multiply Weyl matrix by Pauli spinor, etc.
	VectorUtil::MxV2dC(cs, hip2, vx);
	VectorUtil::MxV2dC(bs, vx, vx);
	VectorUtil::MxV2dC(as, vx, vx);
	//! S(p1,a1,a2,a3,p2)_{lambda1,lambda2}^alpha
	return conj(hip1[0]) * vx[0] + conj(hip1[1]) * vx[1];

}

}
