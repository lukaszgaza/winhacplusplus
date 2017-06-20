/*
 * Transformation.cpp
 *
 *  Created on: 24-08-2011
 *      Author: sobol
 */

#include "Transformation.h"
#include "../core/VinhacException.h"
#include "../input/DataSource.h"

namespace VINHAC {

Transformation::Transformation() {

}

Transformation::~Transformation() {

}

CLHEP::HepLorentzVector VINHAC::Transformation::BOSTDQ(const int& MODE,
		const CLHEP::HepLorentzVector& Q, const CLHEP::HepLorentzVector& P) {
	CLHEP::HepLorentzVector result;
	BOSTDQ(MODE, Q, P, result);
	return result;
}

void VINHAC::Transformation::BOSTDQ(const int& MODE,
		const CLHEP::HepLorentzVector& Q, const CLHEP::HepLorentzVector& P,
		CLHEP::HepLorentzVector& R) {

#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::BOSTDQ()" << std::endl;

#endif
	/*
	 C     *******************************
	 C BOOST ALONG ARBITRARY AXIS (BY RONALD KLEISS).
	 C P BOOSTED INTO R  FROM ACTUAL FRAME TO REST FRAME OF Q
	 C FORTH (MODE = 1) OR BACK (MODE = -1).
	 C Q MUST BE A TIMELIKE, P MAY BE ARBITRARY.
	 */
	R.set(0);
	double FAC = 0;

	double AMQ = Q.m();

	if (MODE == -1) {
		R.setE(
				(P.x() * Q.x() + P.y() * Q.y() + P.z() * Q.z() + P.e() * Q.e())
						/ AMQ);
		FAC = (R.e() + P.e()) / (Q.e() + AMQ);
	} else if (MODE == 1) {
		R.setE(
				(-P.x() * Q.x() - P.y() * Q.y() - P.z() * Q.z() + P.e() * Q.e())
						/ AMQ);
		FAC = -(R.e() + P.e()) / (Q.e() + AMQ);
	} else {
		throw VinhacException(" ++++++++ WRONG MODE IN BOOST3 ");
	}

	R.setX(P.x() + FAC * Q.x());
	R.setY(P.y() + FAC * Q.y());
	R.setZ(P.z() + FAC * Q.z());
}

void VINHAC::Transformation::Rot2le(const double& the, const double& phi,
		CLHEP::HepLorentzVector& pvec) {
	//     ****************************************
	//Euler rotation and boosts for leptonic W-decays

	RXTOD2(the, pvec); //! rotation z-x
	RXTOD3(phi, pvec); //! rotation x-y
}

void VINHAC::Transformation::RXTOD2(const double& PHI,
		CLHEP::HepLorentzVector& PVEC) {

	double CS = cos(PHI);
	double SN = sin(PHI);
	double x, y, z, e;

	x = (CS * PVEC.x() + SN * PVEC.z());
	y = (PVEC.y());
	z = (-SN * PVEC.x() + CS * PVEC.z());
	e = (PVEC.e());
	PVEC.set(x, y, z, e);
}

void VINHAC::Transformation::RXTOD3(const double& PHI,
		CLHEP::HepLorentzVector& PVEC) {

	double CS = cos(PHI);
	double SN = sin(PHI);
	double x, y, z, e;

	x = (CS * PVEC.x() - SN * PVEC.y());
	y = (SN * PVEC.x() + CS * PVEC.y());
	z = (PVEC.z());
	e = (PVEC.e());
	PVEC.set(x, y, z, e);
}

double VINHAC::Transformation::angfix(double x, double y) {
	/*    *******************
	 * CALCULATES ANGLE IN (0,2*PI) RANGE OUT OF X-Y
	 *     ***********************/

	double theta;
	if (fabs(y) < fabs(x)) {
		theta = atan(fabs(y / x));
		if (x <= 0.0)
			theta = DataSource::get().pi - theta;
	} else {
		theta = acos(x / sqrt(x * x + y * y));
	}
	if (y < 0.0)
		theta = 2.0 * DataSource::get().pi - theta;
	return theta;
}

void VINHAC::Transformation::transformationToCMSofWNucleonSystem(
		CLHEP::HepLorentzVector& P, CLHEP::HepLorentzVector& quark,
		CLHEP::HepLorentzVector& antiQuark, CLHEP::HepLorentzVector& boson,
		CLHEP::HepLorentzVector& lepton, CLHEP::HepLorentzVector& neutrino) {
	/*************************************************************************
	 * Lorentz tranformations of quarks, W and leptons 4-momenta to CMS of   *
	 * W-Nucleon system with +z axis pointing along W direction.             *
	 * Input/Output: P - nucleon 4-momentum                                  *
	 *               pqa, paq - 4-momenta of quark and antiquark, resp.      *
	 *               Q,qle,qnu - 4-momenta of W, charged lepton and neutrino *
	 *-----------------------------------------------------------------------*
	 * Written by Wieslaw Placzek,                          Paris, July 2005 *
	 * Last update: 21.07.2005         by: WP                                *
	 *************************************************************************/

	//--- Boost vector
	CLHEP::HepLorentzVector QN = P + boson;
	// Boost to the CMS of W-Nucleon
	P = BOSTDQ(1, QN, P);
	quark = BOSTDQ(1, QN, quark);
	antiQuark = BOSTDQ(1, QN, antiQuark);
	boson = BOSTDQ(1, QN, boson);
	lepton = BOSTDQ(1, QN, lepton);
	neutrino = BOSTDQ(1, QN, neutrino);
	// Rotate W onto xz-plane
	if ((boson.x() * boson.x() + boson.y() * boson.y()) > 1e-10) {
		double Phixy = angfix(boson.x(), boson.y());
		if (Phixy > 1e-8) {
			RXTOD3(-Phixy, P);
			RXTOD3(-Phixy, quark);
			RXTOD3(-Phixy, antiQuark);
			RXTOD3(-Phixy, boson);
			RXTOD3(-Phixy, lepton);
			RXTOD3(-Phixy, neutrino);
		}
	}
	// Rotate W in the xz-plane onto z-axis
	if ((boson.x() * boson.x() + boson.z() * boson.z()) > 1e-10) {
		double Phizx = angfix(boson.z(), boson.x());
		if (Phizx > 1e-8) {
			RXTOD2(-Phizx, P);
			RXTOD2(-Phizx, quark);
			RXTOD2(-Phizx, antiQuark);
			RXTOD2(-Phizx, boson);
			RXTOD2(-Phizx, lepton);
			RXTOD2(-Phizx, neutrino);
		}
	}
}

void VINHAC::Transformation::transformationToBosonRestFrameZAxisAsLABXAxisAccordingToBeamWithSmallerPolarAngle(
		CLHEP::HepLorentzVector& beamA, CLHEP::HepLorentzVector& beamB,
		CLHEP::HepLorentzVector& quark, CLHEP::HepLorentzVector& antiQuark,
		CLHEP::HepLorentzVector& boson, CLHEP::HepLorentzVector& lepton,
		CLHEP::HepLorentzVector& neutrino) {
	/*************************************************************************
	 * Lorentz tranformations of quarks, W and leptons 4-momenta to W rest   *
	 * frame with +z axis pointing in W direction in LAB frame and +x axis   *
	 * oriented according to beam with smaller polar angle in this frame.    *
	 * INPUT/OUTPUT: Q,qle,qnu - 4-momenta of W, charged lepton and neutrino *
	 *-----------------------------------------------------------------------*
	 * Written by Wieslaw Placzek,                        Paris, August 2011 *
	 * Last update: 22.08.2011         by: WP                                *
	 *************************************************************************/

	double Phixy = 0.0;
	double Phizx = 0.0;

	// W 4-momentum in LAB
	CLHEP::HepLorentzVector bosonLab = boson;
	// Boost to the W rest frame
	beamA = BOSTDQ(1, boson, beamA);
	beamB = BOSTDQ(1, boson, beamB);
	quark = BOSTDQ(1, boson, quark);
	antiQuark = BOSTDQ(1, boson, antiQuark);
	lepton = BOSTDQ(1, boson, lepton);
	neutrino = BOSTDQ(1, boson, neutrino);
	boson = BOSTDQ(1, boson, boson);
	// Rotate to make +z axis pointing into W direction in LAB frame
	// 1. Rotation in around z-axis to put xz plane on W direction in LAB
	if (bosonLab.x() * bosonLab.x() + bosonLab.y() * bosonLab.y() > 1e-10) {
		Phixy = angfix(bosonLab.x(), bosonLab.y());
		if (Phixy > 1e-8) {
			RXTOD3(-Phixy, beamA);
			RXTOD3(-Phixy, beamB);
			RXTOD3(-Phixy, quark);
			RXTOD3(-Phixy, antiQuark);
			RXTOD3(-Phixy, bosonLab);
			RXTOD3(-Phixy, lepton);
			RXTOD3(-Phixy, neutrino);
		}
	}
	// 2. Rotation in xz-plane to put z-axis onto W direction in LAB
	if (bosonLab.x() * bosonLab.x() + bosonLab.z() * bosonLab.z() > 1e-10) {
		Phizx = angfix(bosonLab.z(), bosonLab.x());
		if (Phizx > 1e-8) {
			RXTOD2(-Phizx, beamA);
			RXTOD2(-Phizx, beamB);
			RXTOD2(-Phizx, quark);
			RXTOD2(-Phizx, antiQuark);
			RXTOD2(-Phizx, lepton);
			RXTOD2(-Phizx, neutrino);
		}
	}
	// Calculate cosine polar angles of beams momenta
	double costh1 = beamA.z() / sqrt(
			beamA.x() * beamA.x() + beamA.y() * beamA.y() + beamA.z()
					* beamA.z());
	double costh2 = beamB.z() / sqrt(
			beamB.x() * beamB.x() + beamB.y() * beamB.y() + beamB.z()
					* beamB.z());
	// Caculate azimuthal angle for a beam with a smaller polar angle
	if (costh1 > costh2) {
		Phixy = angfix(beamA.x(), beamA.y());
	} else {
		Phixy = angfix(beamB.x(), beamB.y());
	}
	// Rotate all momenta by -Phixy around +z axis
	if (Phixy > 1e-8) {
		RXTOD3(-Phixy, quark);
		RXTOD3(-Phixy, antiQuark);
		RXTOD3(-Phixy, lepton);
		RXTOD3(-Phixy, neutrino);
	}
}

void VINHAC::Transformation::generateFixedPT(
		CLHEP::HepLorentzVector& beamA, CLHEP::HepLorentzVector& beamB,
		CLHEP::HepLorentzVector& quark, CLHEP::HepLorentzVector& antiQuark,
		CLHEP::HepLorentzVector& boson, CLHEP::HepLorentzVector& lepton,
		CLHEP::HepLorentzVector& neutrino) {

	// W 4-momentum in LAB
	CLHEP::HepLorentzVector bosonLab = boson;
	// Boost to the W rest frame
	beamA = BOSTDQ(1, boson, beamA);
	beamB = BOSTDQ(1, boson, beamB);
	quark = BOSTDQ(1, boson, quark);
	antiQuark = BOSTDQ(1, boson, antiQuark);
	lepton = BOSTDQ(1, boson, lepton);
	neutrino = BOSTDQ(1, boson, neutrino);
	boson = BOSTDQ(1, boson, boson);

	double px = 10.0;
	double py = 5.0;
	double pz = (bosonLab.pz()>=0?1:-1)*sqrt(bosonLab.e()*bosonLab.e()-bosonLab.m2()-px*px-py*py);
	bosonLab.setX(px);
	bosonLab.setY(py);
	bosonLab.setZ(pz);


	beamA = BOSTDQ(-1, bosonLab, beamA);
	beamB = BOSTDQ(-1, bosonLab, beamB);
	quark = BOSTDQ(-1, bosonLab, quark);
	antiQuark = BOSTDQ(-1, bosonLab, antiQuark);
	lepton = BOSTDQ(-1, bosonLab, lepton);
	neutrino = BOSTDQ(-1, bosonLab, neutrino);
	boson = BOSTDQ(-1, bosonLab, boson);

}

void VINHAC::Transformation::generateRandomPT(
		CLHEP::HepLorentzVector& beamA, CLHEP::HepLorentzVector& beamB,
		CLHEP::HepLorentzVector& quark, CLHEP::HepLorentzVector& antiQuark,
		CLHEP::HepLorentzVector& boson, CLHEP::HepLorentzVector& lepton,
		CLHEP::HepLorentzVector& neutrino,
		TRND *randomEngine) {

	// W 4-momentum in LAB
	CLHEP::HepLorentzVector bosonLab = boson;
	// Boost to the W rest frame
	beamA = BOSTDQ(1, boson, beamA);
	beamB = BOSTDQ(1, boson, beamB);
	quark = BOSTDQ(1, boson, quark);
	antiQuark = BOSTDQ(1, boson, antiQuark);
	lepton = BOSTDQ(1, boson, lepton);
	neutrino = BOSTDQ(1, boson, neutrino);
	boson = BOSTDQ(1, boson, boson);

	double phi = 2*M_PI*randomEngine->Flat();
	double pT = -10*log(1-randomEngine->Flat());
	double px = pT*cos(phi);
	double py = pT*sin(phi);
	double pz = (bosonLab.pz()>=0?1:-1)*sqrt(bosonLab.e()*bosonLab.e()-bosonLab.m2()-px*px-py*py);
	bosonLab.setX(px);
	bosonLab.setY(py);
	bosonLab.setZ(pz);


	beamA = BOSTDQ(-1, bosonLab, beamA);
	beamB = BOSTDQ(-1, bosonLab, beamB);
	quark = BOSTDQ(-1, bosonLab, quark);
	antiQuark = BOSTDQ(-1, bosonLab, antiQuark);
	lepton = BOSTDQ(-1, bosonLab, lepton);
	neutrino = BOSTDQ(-1, bosonLab, neutrino);
	boson = BOSTDQ(-1, bosonLab, boson);

}


}
