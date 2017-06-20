/*
 * SingleRadiation.cpp
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#include "SingleRadiation.h"
#include "../util/VectorUtil.h"
#include "../util/SpinorialFunction.h"

namespace VINHAC {

SingleRadiation::SingleRadiation() {

}

SingleRadiation::~SingleRadiation() {
}

VINHAC::SpecialMatrix VINHAC::SingleRadiation::radiativeWdecay(
		const CLHEP::HepLorentzVector& boson,
		const CLHEP::HepLorentzVector& lepton,
		const CLHEP::HepLorentzVector& antilepton,
		const CLHEP::HepLorentzVector& photon, int chargeOfBoson) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::radiativeWdecay()" << std::endl;

#endif
	double e(sqrt(4 * DataSource::get().pi / DataSource::get().invAlphaQED));
	double colourFactor(1);
	complex<double> factorRadiativeWdecay(0,
			-e * e * colourFactor / sqrt(2.0 * DataSource::get().sinThetaW2));

	vector<unsigned> ranges(4);
	ranges[0] = 3;
	ranges[1] = 2;
	ranges[2] = 2;
	ranges[3] = 2;

	VINHAC::SpecialMatrix result(ranges);

	//! Electric charge of fermion f2 in units of e
	int Qf1 = (chargeOfBoson - 1) / 2;
	int Qf2 = Qf1 - chargeOfBoson;

	//! zero/non-zero fermion charges
	bool LQf1 = Qf1 != 0;
	bool LQf2 = Qf2 != 0;

	CLHEP::HepLorentzVector epsilonW;
	CLHEP::HepLorentzVector epsilonK;
	//! Calculation of the radiative W-decay amplitudes
	for (int bos = 0; bos < 3; ++bos) { //              !<-- loop over W polarizations
		VectorUtil::VectorUtil::VecPol(epsilonW, boson, bos + 1);
		for (int l = 0; l < 2; ++l) {//         !<-- loop over f1 helicities

			double omegalept = VectorUtil::VectorUtil::omega(-(2 * l - 1),
					lepton.e(), lepton.getV().mag());
			//! Pauli spinor for f1
			complex<double> hipl[2];

			VectorUtil::VectorUtil::HelEig(lepton, l * 2 - 1, hipl);

			for (int al = 0; al < 2; ++al) {//      !<-- loop over f2~ helicities

				double omegaantilept = VectorUtil::omega(al * 2 - 1,
						antilepton.e(), antilepton.getV().mag());
				//! Pauli spinor for f2~
				complex<double> hipal[2];

				VectorUtil::HelEig(antilepton, -al * 2 + 1, hipal);
				complex<double> Sf1 = SpinorialFunction::spinorialFunction(
						hipl, epsilonW, hipal, -1);
				for (int kap = 0; kap < 2; ++kap) { //DO kap = 1,2 !<-- loop over photon polarizations
					VectorUtil::VecPol(epsilonK, photon, kap + 1);

					double radf = static_cast<double> (-chargeOfBoson);
					radf *= (boson * epsilonK);
					radf /= (boson * photon);

					if (LQf1)
						radf += static_cast<double> (Qf1) * (lepton * epsilonK)
								/ (lepton * photon);
					if (LQf2)
						radf -= static_cast<double> (Qf2) * (antilepton
								* epsilonK) / (antilepton * photon);
					complex<double> amr = radf;
					amr *= Sf1;

					amr -= static_cast<double> (chargeOfBoson) * (photon
							* epsilonW) / (boson * photon)
							* SpinorialFunction::spinorialFunction(hipl,
									epsilonK, hipal, -1) + static_cast<double> (chargeOfBoson) * (epsilonK
							* epsilonW) / (boson * photon)
							* SpinorialFunction::spinorialFunction(hipl,
									photon, hipal, -1);

					if (LQf1)
						amr += static_cast<double> (Qf1) / 2.0 / (lepton
								* photon)
								* SpinorialFunction::spinorialFunction(hipl,
										epsilonK, photon, epsilonW, hipal, -1);
					if (LQf2)
						amr -= static_cast<double> (Qf2) / 2.0 / (antilepton
								* photon)
								* SpinorialFunction::spinorialFunction(hipl,
										epsilonW, photon, epsilonK, hipal, -1);

					ranges[0] = bos;
					ranges[1] = l;
					ranges[2] = al;
					ranges[3] = kap;

					complex<double> tmp = factorRadiativeWdecay;
					tmp *= static_cast<double> (al * 2 - 1);
					tmp *= omegalept;
					tmp *= omegaantilept;
					tmp *= amr;

					result.set(ranges, tmp);
				}// ! photon
			}// ! f2~
		}// ! f1
	}//! W

#ifdef DEBUG_FLAG_AMPL
	fcomplex amp[3][2][2][2];
	double p1[4];
	double p2[4];
	double q[4];
	double pk[4];
	double one = 1.0;
	for (int i = 0; i < 4; ++i) {
		p1[i] = lepton[i];
		p2[i] = antilepton[i];
		q[i] = boson[i];
		pk[1] = photon[i];
	}
	wdecrad_(DataSource::get().sinThetaW2, &one, &one, &QW, &Qf1, q, p1, p2, pk, amp);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				for (int l = 0; l < 2; ++l) {
					std::cout << "DEBUG: ModelHandler::radWdecay() fortran ["
					<< i << "," << j << "," << k << "," << l << "]("
					<< amp[i][j][k][l].r << "," << amp[i][j][k][l].i
					<< ")" << std::endl;
					ranges[0] = i;
					ranges[1] = j;
					ranges[2] = k;
					ranges[3] = l;
					std::cout << "DEBUG: ModelHandler::radWdecay() cpp [" << i
					<< "," << j << "," << k << "," << l << "]("
					<< result.get(ranges) << ")" << std::endl;
				}
			}
		}
	}

#endif
	return result;
}

}
