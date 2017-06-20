/*
 * MatrixElement.h
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#ifndef MATRIXELEMENT_H_
#define MATRIXELEMENT_H_

#include "../core/Event.h"
#include "../input/DataSource.h"
#include <map>
#include "../core/enum/PolarizationDefinition.h"

namespace VINHAC {

class MatrixElement {
public:
	MatrixElement();
	virtual ~MatrixElement();

	static double matrixElement0(const CLHEP::HepLorentzVector& quark,
			const CLHEP::HepLorentzVector& antiquark,
			const CLHEP::HepLorentzVector& boson,
			const CLHEP::HepLorentzVector& lepton,
			const CLHEP::HepLorentzVector& antilepton, double ckmElement);
	static std::map<PolarizationName::name,double> matrixElement0(
			const CLHEP::HepLorentzVector& quark,
			const CLHEP::HepLorentzVector& antiquark,
			const CLHEP::HepLorentzVector& boson,
			const CLHEP::HepLorentzVector& lepton,
			const CLHEP::HepLorentzVector& antilepton, double ckmElement , std::set<PolarizationName::name> polarizations);
	static double matrixElementFSR1(const CLHEP::HepLorentzVector& quark,
			const CLHEP::HepLorentzVector& antiquark,
			const CLHEP::HepLorentzVector& boson,
			const CLHEP::HepLorentzVector& lepton,
			const CLHEP::HepLorentzVector& antilepton, const CLHEP::HepLorentzVector& photon,
			int chargeOfBoson,
			double ckmElement);
};

}

#endif /* MATRIXELEMENT_H_ */
