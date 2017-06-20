/*
 * Born.h
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#ifndef BORN_H_
#define BORN_H_

#include "../../../utils/SpecialMatrix.h"
#include "../../core/Event.h"
#include "../../input/DataSource.h"

namespace VINHAC {

class Born {
public:
	Born();
	virtual ~Born();

	static SpecialMatrix bornWproduction(const CLHEP::HepLorentzVector& quark,
			const CLHEP::HepLorentzVector& antiquark,
			const CLHEP::HepLorentzVector& boson, double ckmElement);
	static SpecialMatrix bornWdecay(const CLHEP::HepLorentzVector& boson,
			const CLHEP::HepLorentzVector& lepton,
			const CLHEP::HepLorentzVector& antilepton);
};

}

#endif /* BORN_H_ */
