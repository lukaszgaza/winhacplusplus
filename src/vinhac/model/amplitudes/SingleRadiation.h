/*
 * SingleRadiation.h
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#ifndef SINGLERADIATION_H_
#define SINGLERADIATION_H_

#include "../../../utils/SpecialMatrix.h"
#include "../../core/Event.h"
#include "../../input/DataSource.h"

namespace VINHAC {

class SingleRadiation {
public:
	SingleRadiation();
	virtual ~SingleRadiation();

	static SpecialMatrix radiativeWdecay(const CLHEP::HepLorentzVector& boson,
			const CLHEP::HepLorentzVector& lepton,
			const CLHEP::HepLorentzVector& antilepton,
			const CLHEP::HepLorentzVector& photon,
			int chargeOfBoson);
};

}

#endif /* SINGLERADIATION_H_ */
