/*
 * Transformation.h
 *
 *  Created on: 24-08-2011
 *      Author: sobol
 */

#ifndef TRANSFORMATION_H_
#define TRANSFORMATION_H_

#include "../../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include "../../utils/FOAM/TRND.h"

namespace VINHAC {

class Transformation {
public:
	Transformation();
	virtual ~Transformation();

	static void BOSTDQ(const int& MODE, const CLHEP::HepLorentzVector& Q,
			const CLHEP::HepLorentzVector& P, CLHEP::HepLorentzVector& R);
	static CLHEP::HepLorentzVector BOSTDQ(const int& MODE, const CLHEP::HepLorentzVector& Q,
				const CLHEP::HepLorentzVector& P);

	static void Rot2le(const double& the, const double& phi,
			CLHEP::HepLorentzVector& pvec);
	static void RXTOD2(const double& PHI, CLHEP::HepLorentzVector& PVEC);
	static void RXTOD3(const double& PHI, CLHEP::HepLorentzVector& PVEC);
	static double angfix(double x,double y);

	static void transformationToCMSofWNucleonSystem(
			CLHEP::HepLorentzVector& P, CLHEP::HepLorentzVector& quark,
			CLHEP::HepLorentzVector& antiQuark, CLHEP::HepLorentzVector& boson,
			CLHEP::HepLorentzVector& lepton, CLHEP::HepLorentzVector& neutrino);

	static void transformationToBosonRestFrameZAxisAsLABXAxisAccordingToBeamWithSmallerPolarAngle(
			CLHEP::HepLorentzVector& beamA,
			CLHEP::HepLorentzVector& beamB,
			CLHEP::HepLorentzVector& quark,
			CLHEP::HepLorentzVector& antiQuark, CLHEP::HepLorentzVector& boson,
			CLHEP::HepLorentzVector& lepton, CLHEP::HepLorentzVector& neutrino);

	static void generateFixedPT(
			CLHEP::HepLorentzVector& beamA,
			CLHEP::HepLorentzVector& beamB,
			CLHEP::HepLorentzVector& quark,
			CLHEP::HepLorentzVector& antiQuark, CLHEP::HepLorentzVector& boson,
			CLHEP::HepLorentzVector& lepton, CLHEP::HepLorentzVector& neutrino);

	static void generateRandomPT(
				CLHEP::HepLorentzVector& beamA,
				CLHEP::HepLorentzVector& beamB,
				CLHEP::HepLorentzVector& quark,
				CLHEP::HepLorentzVector& antiQuark, CLHEP::HepLorentzVector& boson,
				CLHEP::HepLorentzVector& lepton, CLHEP::HepLorentzVector& neutrino
				, TRND *randomEngine);
};

}

#endif /* TRANSFORMATION_H_ */
