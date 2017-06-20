/*
 * Born.cpp
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#include "Born.h"
#include "../util/SpinorialFunction.h"
#include "../util/VectorUtil.h"

#ifdef DEBUG_FLAG
#include "../../../utils/FortranFunction.h"
#endif

namespace VINHAC {

Born::Born() {
}

Born::~Born() {
}

VINHAC::SpecialMatrix VINHAC::Born::bornWproduction(
		const CLHEP::HepLorentzVector& quark,
		const CLHEP::HepLorentzVector& antiquark,
		const CLHEP::HepLorentzVector& boson, double ckmElement) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::bornWproduction()" << std::endl;

#endif
	double e(sqrt(4 * DataSource::get().pi / DataSource::get().invAlphaQED));
	complex<double> factor(0, - ckmElement * e / sqrt(2.0 * DataSource::get().sinThetaW2));

	vector<unsigned> ranges(3);
	ranges[0] = 2;
	ranges[1] = 2;
	ranges[2] = 3;

	VINHAC::SpecialMatrix result(ranges);

	CLHEP::HepLorentzVector epsilonW;
	for (int bos = 0; bos < 3; ++bos) {

		VectorUtil::VecPol(epsilonW, boson, bos + 1);

		for (int q = 0; q < 2; ++q) {
			double omegaq = VectorUtil::omega(-(2 * q - 1), quark.e(),
					quark.getV().mag());

			complex<double> hipq[2];

			VectorUtil::HelEig(quark, q * 2 - 1, hipq);

			for (int aq = 0; aq < 2; ++aq) {
				double omegaaq = VectorUtil::omega(aq * 2 - 1, antiquark.e(),
						antiquark.getV().mag());

				complex<double> hipaq[2];

				VectorUtil::HelEig(antiquark, -aq * 2 + 1, hipaq);

				ranges[0] = q;
				ranges[1] = aq;
				ranges[2] = bos;

				complex<double> tmp = factor;
				tmp *= static_cast<double> (aq * 2 - 1);
				tmp *= omegaq;
				tmp *= omegaaq;
				tmp *= SpinorialFunction::spinorialFunction(hipaq, epsilonW,
						hipq, -1);
				result.set(ranges, tmp);

			}

		}

	}

#ifdef DEBUG_FLAG_AMPL
	//fcomplex amp[2][2][3];
	double p1[4];
	double p2[4];
	double q[4];
	for (int i = 0; i < 4; ++i) {
		p1[i] = quark[i];
		p2[i] = antiquark[i];
		q[i] = boson[i];
	}
	//wprobor_(&DataSource::get().sinThetaW2, &ckmElement, p1, p2, q, amp);

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 3; ++k) {
			//	std::cout << "DEBUG: ModelHandler::bornWproduction() fortran ["
			//	<< i << "," << j << "," << k << "](" << amp[i][j][k].r
			//	<< "," << amp[i][j][k].i << ")" << std::endl;
				ranges[0] = i;
				ranges[1] = j;
				ranges[2] = k;
				std::cout << "DEBUG: ModelHandler::bornWproduction() cpp ["
				<< i << "," << j << "," << k << "](" << result.get(
						ranges) << ")" << std::endl;
			}
		}
	}

#endif

	return result;
}

VINHAC::SpecialMatrix VINHAC::Born::bornWdecay(
		const CLHEP::HepLorentzVector& boson,
		const CLHEP::HepLorentzVector& lepton,
		const CLHEP::HepLorentzVector& antilepton) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::bornWdecay()" << std::endl;

#endif

	double colourFactor(1);
	double e(sqrt(4 * DataSource::get().pi / DataSource::get().invAlphaQED));
	complex<double> factorBornWdecay(0,
			-e * colourFactor / sqrt(2.0 * DataSource::get().sinThetaW2));


	vector<unsigned> ranges(2);
	ranges[0] = 3;
	ranges[1] = 2;
	ranges[2] = 2;

	VINHAC::SpecialMatrix result(ranges);

	CLHEP::HepLorentzVector epsilonW;
	for (int bos = 0; bos < 3; ++bos) {

		VectorUtil::VecPol(epsilonW, boson, bos + 1);

		for (int l = 0; l < 2; ++l) {
			double omegalept = VectorUtil::omega(-(2 * l - 1), lepton.e(),
					lepton.getV().mag());

			complex<double> hipl[2];

			VectorUtil::HelEig(lepton, l * 2 - 1, hipl);

			for (int al = 0; al < 2; ++al) {
				double omegaantilept = VectorUtil::omega(al * 2 - 1,
						antilepton.e(), antilepton.getV().mag());

				complex<double> hipal[2];

				VectorUtil::HelEig(antilepton, -al * 2 + 1, hipal);

				ranges[1] = l;
				ranges[2] = al;
				ranges[0] = bos;

				complex<double> tmp = factorBornWdecay;
				tmp *= static_cast<double> (al * 2 - 1);
				tmp *= omegalept;
				tmp *= omegaantilept;
				tmp *= SpinorialFunction::spinorialFunction(hipl, epsilonW,
						hipal, -1);
				result.set(ranges, tmp);

			}

		}

	}
#ifdef DEBUG_FLAG_AMPL
	//fcomplex amp[3][2][2];
	double p1[4];
	double p2[4];
	double q[4];
	double one = 1.0;
	for (int i = 0; i < 4; ++i) {
		p1[i] = lepton[i];
		p2[i] = antilepton[i];
		q[i] = boson[i];
	}
	//wdecbor_(&DataSource::get().sinThetaW2, &one, &one, q, p1, p2, amp);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
			//	std::cout << "DEBUG: ModelHandler::bornWdecay() fortran [" << i
			//	<< "," << j << "," << k << "](" << amp[i][j][k].r
			//	<< "," << amp[i][j][k].i << ")" << std::endl;
				ranges[0] = i;
				ranges[1] = j;
				ranges[2] = k;
				std::cout << "DEBUG: ModelHandler::bornWdecay() cpp [" << i
				<< "," << j << "," << k << "](" << result.get(ranges)
				<< ")" << std::endl;
			}
		}
	}

#endif

	return result;

}

}
