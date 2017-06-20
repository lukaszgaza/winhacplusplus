/*
 * ModelHandler.cpp
 *
 *  Created on: 26 Aug 2009
 *      Author: siudmy
 */
#include "ModelHandler.h"
#include "../utils/SANC/wh_ewcdec.h"

#ifdef DEBUG_FLAG
#include "../utils/FortranFunction.h"
#endif

//#define DEBUG_FLAG_AMPL

VINHAC::ModelHandler::ModelHandler(DataSource *ds,
		VINHAC::HardProcessHandler* p_Hard) :
	ds(ds), p_HardMatrix(p_Hard), colourFactor(1)

{
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::ModelHandler()" << std::endl;

#endif

	e = sqrt(4 * ds->pi / ds->invAlphaQED);
	normalizationFactor0 = 1.0 / 8.0 / (2.0 * 2.0 * ds->pi * ds->pi);
	normalizationFactor1 = 1.0 / 16.0 / pow(2.0 * ds->pi, 5);

	if (ds->keyRad == 2) {
		double mw = ds->MW;
		int id = 11;
		int KeyGmu = ds->keyGmu;
		double aMW = ds->MW;
		double aMZ = ds->MZ;
		double aMH = ds->MH;
		double Xmfe[20];
		double Xckm[3][3];
		for(unsigned i = 0 ; i<20 ;i++){
			Xmfe[i]=ds->masses[i+1];
		}

		for(unsigned i = 0 ;i<ds->ckm.size();++i){
			for(unsigned j = 0 ; j<ds->ckm[i].size(); ++j){
				Xckm[j][i]=ds->ckm[i][j];
			}
		}

		deltaWeak[11] = wh_delewdec_(&mw,&id,&KeyGmu,
				&aMW,&aMZ, &aMH,Xmfe,
				Xckm);
		id=13;
		deltaWeak[13] = wh_delewdec_(&mw,&id,&KeyGmu,
				&aMW,&aMZ, &aMH,Xmfe,
				Xckm);
		id=15;
		deltaWeak[15] = wh_delewdec_(&mw,&id,&KeyGmu,
				&aMW,&aMZ, &aMH,Xmfe,
				Xckm);
	} else {
		deltaWeak[11] = 0.0;
		deltaWeak[13] = 0.0;
		deltaWeak[15] = 0.0;
	}

}


void VINHAC::ModelHandler::prepareParticles(VINHAC::Event& event) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::prepareParticles()" << std::endl;

#endif

	CLHEP::HepLorentzVector pqu, paq, QQ, qlp, qlm, q1, q2, qk;
	// get initial quark and anti quark

	int iu;
	int id;
	if (event.getBoson().getPDGid() == -24) {
		iu = abs(event.getAntiQuark().getPDGid()) / 2 - 1;
		id = event.getQuark().getPDGid() / 2 + 1 - 1;
	} else {
		id = abs(event.getAntiQuark().getPDGid()) / 2 + 1 - 1;
		iu = event.getQuark().getPDGid() / 2 - 1;
	}

	ckmElement = ds->ckm[iu][id];
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::prepareParticles() # ckm[" << iu << "]["
			<< id << "] " << ds->ckm[iu][id] << std::endl;

#endif
}

double VINHAC::ModelHandler::matrixElement0(Event& event) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::matrixElement0()" << std::endl;

#endif
	double result;
	SpecialMatrix prod = bornWproduction(event);
	SpecialMatrix dec = bornWdecay(event);

	result = (prod * dec).real();

	double initialSpinFactor = 0.25;
	double colourFac = 1.0 / 3.0;

	double s = event.getBoson().getFourMomentumBosonRestFrame().mag2();
	result = result / inverseWPropagator(s, ds);
	result = result * initialSpinFactor * colourFac;

#ifdef DEBUG_FLAG
	int keyWid = 0;
	double amW = ds->MW;
	double GaW = ds->GW;
	double sW2 = ds->sinThetaW2;
	double Uij = ckmElement;
	double p1[4];
	double p2[4];
	double Q[4];
	double q1[4];
	double q2[4];
	for (int i = 0; i < 4; ++i) {
		p1[i] = event.getQuark().getFourMomentumBosonRestFrame()[i];
		p2[i] = event.getAntiQuark().getFourMomentumBosonRestFrame()[i];
		Q[i] = event.getBoson().getFourMomentumBosonRestFrame()[i];
		if (event.getLepton().getPDGid() > 0) {
			q1[i] = event.getLepton().getFourMomentumBosonRestFrame()[i];
			q2[i] = event.getNeutrino().getFourMomentumBosonRestFrame()[i];
		} else {
			q2[i] = event.getLepton().getFourMomentumBosonRestFrame()[i];
			q1[i] = event.getNeutrino().getFourMomentumBosonRestFrame()[i];
		}
	}
	double fortranxmat0 = wh_mat0_(&keyWid, &amW, &GaW, &sW2, &Uij, p1, p2, Q,
			q1, q2);

	std::cout << "DEBUG: ModelHandler::matrixElement0() # fortran: "
			<< fortranxmat0 << std::endl;
	std::cout << "DEBUG: ModelHandler::matrixElement0() # cpp: " << result
			<< std::endl;
	std::cout << "DEBUG: ModelHandler::matrixElement0() # cpp/fortran: "
			<< result / fortranxmat0 << std::endl;
#endif

	return result;
}

VINHAC::SpecialMatrix VINHAC::ModelHandler::bornWproduction(Event& event) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::bornWproduction()" << std::endl;

#endif
	vector<unsigned> ranges;
	ranges.push_back(2);
	ranges.push_back(2);
	ranges.push_back(3);

	VINHAC::SpecialMatrix result(ranges);

	double e = sqrt(4 * ds->pi / ds->invAlphaQED);

	complex<double> factor(0, -e * ckmElement / sqrt(2.0 * ds->sinThetaW2));

	for (int bos = 0; bos < 3; ++bos) {

		CLHEP::HepLorentzVector epsilonW = VecPol(
				event.getBoson().getFourMomentumBosonRestFrame(), bos + 1);

		for (int q = 0; q < 2; ++q) {
			double
					omegaq =
							omega(
									-(2 * q - 1),
									event.getQuark().getFourMomentumBosonRestFrame().e(),
									event.getQuark().getFourMomentumBosonRestFrame().getV().mag());

			complex<double> hipq[2];

			HelEig(event.getQuark().getFourMomentumBosonRestFrame(), q * 2 - 1,
					hipq);

			for (int aq = 0; aq < 2; ++aq) {
				double
						omegaaq =
								omega(
										aq * 2 - 1,
										event.getAntiQuark().getFourMomentumBosonRestFrame().e(),
										event.getAntiQuark().getFourMomentumBosonRestFrame().getV().mag());

				complex<double> hipaq[2];

				HelEig(event.getAntiQuark().getFourMomentumBosonRestFrame(),
						-aq * 2 + 1, hipaq);

				ranges[0] = q;
				ranges[1] = aq;
				ranges[2] = bos;

				complex<double> tmp = factor * static_cast<double> (aq * 2 - 1)
						* omegaq * omegaaq * spinorialFunction(hipaq, epsilonW,
						hipq, -1);
				result.set(ranges, tmp);

			}

		}

	}

#ifdef DEBUG_FLAG_AMPL
	fcomplex amp[2][2][3];
	double p1[4];
	double p2[4];
	double q[4];
	for (int i = 0; i < 4; ++i) {
		p1[i] = event.getQuark().getFourMomentumBosonRestFrame()[i];
		p2[i] = event.getAntiQuark().getFourMomentumBosonRestFrame()[i];
		q[i] = event.getBoson().getFourMomentumBosonRestFrame()[i];
	}
	wprobor_(&ds->sinThetaW2, &ckmElement, p1, p2, q, amp);

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 3; ++k) {
				std::cout << "DEBUG: ModelHandler::bornWproduction() fortran ["
				<< i << "," << j << "," << k << "](" << amp[i][j][k].r
				<< "," << amp[i][j][k].i << ")" << std::endl;
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

VINHAC::SpecialMatrix VINHAC::ModelHandler::bornWdecay(Event& event) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::bornWdecay()" << std::endl;

#endif
	Particle* lept;
	Particle* antilept;

	if (event.getLepton().getPDGid() > 0) {
		lept = &event.getLepton();
		antilept = &event.getNeutrino();
	} else {
		antilept = &event.getLepton();
		lept = &event.getNeutrino();
	}

	vector<unsigned> ranges;
	ranges.push_back(3);
	ranges.push_back(2);
	ranges.push_back(2);

	VINHAC::SpecialMatrix result(ranges);

	double e = sqrt(4 * ds->pi / ds->invAlphaQED);

	complex<double> factor(0, -e * colourFactor / sqrt(2.0 * ds->sinThetaW2));

	for (int bos = 0; bos < 3; ++bos) {

		CLHEP::HepLorentzVector epsilonW = VecPol(
				event.getBoson().getFourMomentumBosonRestFrame(), bos + 1);

		for (int l = 0; l < 2; ++l) {
			double omegalept = omega(-(2 * l - 1),
					lept->getFourMomentumBosonRestFrame().e(),
					lept->getFourMomentumBosonRestFrame().getV().mag());

			complex<double> hipl[2];

			HelEig(lept->getFourMomentumBosonRestFrame(), l * 2 - 1, hipl);

			for (int al = 0; al < 2; ++al) {
				double omegaantilept = omega(al * 2 - 1,
						antilept->getFourMomentumBosonRestFrame().e(),
						antilept->getFourMomentumBosonRestFrame().getV().mag());

				complex<double> hipal[2];

				HelEig(antilept->getFourMomentumBosonRestFrame(), -al * 2 + 1,
						hipal);

				ranges[1] = l;
				ranges[2] = al;
				ranges[0] = bos;

				complex<double> tmp = factor * static_cast<double> (al * 2 - 1)
						* omegalept * omegaantilept * spinorialFunction(hipl,
						epsilonW, hipal, -1);
				result.set(ranges, tmp);

			}

		}

	}
#ifdef DEBUG_FLAG_AMPL
	fcomplex amp[3][2][2];
	double p1[4];
	double p2[4];
	double q[4];
	double one = 1.0;
	for (int i = 0; i < 4; ++i) {
		p1[i] = lept->getFourMomentumBosonRestFrame()[i];
		p2[i] = antilept->getFourMomentumBosonRestFrame()[i];
		q[i] = event.getBoson().getFourMomentumBosonRestFrame()[i];
	}
	wdecbor_(&ds->sinThetaW2, &one, &one, q, p1, p2, amp);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				std::cout << "DEBUG: ModelHandler::bornWdecay() fortran [" << i
				<< "," << j << "," << k << "](" << amp[i][j][k].r
				<< "," << amp[i][j][k].i << ")" << std::endl;
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

VINHAC::SpecialMatrix VINHAC::ModelHandler::radiativeWdecay(Event& event,
		Particle& photon) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::radiativeWdecay()" << std::endl;

#endif
	Particle* lept;
	Particle* antilept;

	if (event.getLepton().getPDGid() > 0) {
		lept = &event.getLepton();
		antilept = &event.getNeutrino();
	} else {
		antilept = &event.getLepton();
		lept = &event.getNeutrino();
	}

	vector<unsigned> ranges;
	ranges.push_back(3);
	ranges.push_back(2);
	ranges.push_back(2);
	ranges.push_back(2);

	VINHAC::SpecialMatrix result(ranges);

	double e = sqrt(4.0 * ds->pi / ds->invAlphaQED);

	complex<double> factor(0, -e * e * colourFactor
			/ sqrt(2.0 * ds->sinThetaW2));

	double QW = 1.0;
	if (event.getBoson().getPDGid() < 0) {
		QW = -1.0;
	}

	//! Electric charge of fermion f2 in units of e
	double Qf1 = (QW - 1.0) / 2.0;
	double Qf2 = Qf1 - QW;

	//! zero/non-zero fermion charges
	bool LQf1 = abs(Qf1) > 1.0E-10;
	bool LQf2 = abs(Qf2) > 1.0E-10;

	//! Calculation of the radiative W-decay amplitudes
	for (int bos = 0; bos < 3; ++bos) { //              !<-- loop over W polarizations
		CLHEP::HepLorentzVector epsilonW = VecPol(
				event.getBoson().getFourMomentumBosonRestFrame(), bos + 1);
		for (int l = 0; l < 2; ++l) {//         !<-- loop over f1 helicities

			double omegalept = omega(-(2 * l - 1),
					lept->getFourMomentumBosonRestFrame().e(),
					lept->getFourMomentumBosonRestFrame().getV().mag());
			//! Pauli spinor for f1
			complex<double> hipl[2];

			HelEig(lept->getFourMomentumBosonRestFrame(), l * 2 - 1, hipl);

			for (int al = 0; al < 2; ++al) {//      !<-- loop over f2~ helicities

				double omegaantilept = omega(al * 2 - 1,
						antilept->getFourMomentumBosonRestFrame().e(),
						antilept->getFourMomentumBosonRestFrame().getV().mag());
				//! Pauli spinor for f2~
				complex<double> hipal[2];

				HelEig(antilept->getFourMomentumBosonRestFrame(), -al * 2 + 1,
						hipal);
				complex<double> Sf1 = spinorialFunction(hipl, epsilonW, hipal,
						-1);
				for (int kap = 0; kap < 2; ++kap) { //DO kap = 1,2 !<-- loop over photon polarizations
					CLHEP::HepLorentzVector epsilonK = VecPol(
							photon.getFourMomentumBosonRestFrame(), kap + 1);

					double radf = -QW
							* (event.getBoson().getFourMomentumBosonRestFrame()
									* epsilonK)
							/ (event.getBoson().getFourMomentumBosonRestFrame()
									* photon.getFourMomentumBosonRestFrame());

					if (LQf1)
						radf
								= radf
										+ Qf1
												* (lept->getFourMomentumBosonRestFrame()
														* epsilonK)
												/ (lept->getFourMomentumBosonRestFrame()
														* photon.getFourMomentumBosonRestFrame());
					if (LQf2)
						radf
								= radf
										- Qf2
												* (antilept->getFourMomentumBosonRestFrame()
														* epsilonK)
												/ (antilept->getFourMomentumBosonRestFrame()
														* photon.getFourMomentumBosonRestFrame());
					complex<double> amr = radf * Sf1;

					amr = amr - QW * (photon.getFourMomentumBosonRestFrame()
							* epsilonW)
							/ (event.getBoson().getFourMomentumBosonRestFrame()
									* photon.getFourMomentumBosonRestFrame())
							* spinorialFunction(hipl, epsilonK, hipal, -1) + QW
							* (epsilonK * epsilonW)
							/ (event.getBoson().getFourMomentumBosonRestFrame()
									* photon.getFourMomentumBosonRestFrame())
							* spinorialFunction(hipl,
									photon.getFourMomentumBosonRestFrame(),
									hipal, -1);

					if (LQf1)
						amr
								= amr
										+ Qf1 / 2.0
												/ (lept->getFourMomentumBosonRestFrame()
														* photon.getFourMomentumBosonRestFrame())
												* spinorialFunction(
														hipl,
														epsilonK,
														photon.getFourMomentumBosonRestFrame(),
														epsilonW, hipal, -1);
					if (LQf2)
						amr
								= amr
										- Qf2 / 2.0
												/ (antilept->getFourMomentumBosonRestFrame()
														* photon.getFourMomentumBosonRestFrame())
												* spinorialFunction(
														hipl,
														epsilonW,
														photon.getFourMomentumBosonRestFrame(),
														epsilonK, hipal, -1);

					ranges[0] = bos;
					ranges[1] = l;
					ranges[2] = al;
					ranges[3] = kap;

					complex<double> tmp = factor * static_cast<double> (al * 2
							- 1) * omegalept * omegaantilept * amr;

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
		p1[i] = lept->getFourMomentumBosonRestFrame()[i];
		p2[i] = antilept->getFourMomentumBosonRestFrame()[i];
		q[i] = event.getBoson().getFourMomentumBosonRestFrame()[i];
		pk[1] = photon.getFourMomentumBosonRestFrame()[i];
	}
	wdecrad_(&ds->sinThetaW2, &one, &one, &QW, &Qf1, q, p1, p2, pk, amp);

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

double VINHAC::ModelHandler::matrixElementFSR1(Event& event, Particle& photon) {

	VINHAC::SpecialMatrix prod = bornWproduction(event);
	VINHAC::SpecialMatrix dec = radiativeWdecay(event, photon);

	double result = (prod * dec).real();

	double initialSpinFactor = 0.25;
	double colourFac = 1.0 / 3.0;

	double s = event.getBoson().getFourMomentumBosonRestFrame().mag2();
	result = result / inverseWPropagator(s, ds);
	result = result * initialSpinFactor * colourFac;

	//*! Correction for very collinear radiation (big cancellations between
	//*! various contributions to spin amplitudes lead to lost of numerical
	//*! precision in these regions!)
	double xmbo = matrixElement0(event);
	double QW = 1.0;
	if (event.getBoson().getPDGid() < 0) {
		QW = -1.0;
	}

	//! Electric charge of fermion f2 in units of e
	double Qf1 = (QW - 1.0) / 2.0;

	Particle* lept;
	Particle* antilept;

	if (event.getLepton().getPDGid() > 0) {
		lept = &event.getLepton();
		antilept = &event.getNeutrino();
	} else {
		antilept = &event.getLepton();
		lept = &event.getNeutrino();
	}
	CLHEP::HepLorentzVector ql;
	if (abs(Qf1) > 1.0E-10) {
		ql = lept->getFourMomentumBosonRestFrame();
	} else {
		ql = antilept->getFourMomentumBosonRestFrame();
	}
	double aml2 = ql * ql;
	double amQ2 = event.getBoson().getFourMomentumBosonRestFrame()
			* event.getBoson().getFourMomentumBosonRestFrame();
	double qlk = ql * photon.getFourMomentumBosonRestFrame();
	double QWk = event.getBoson().getFourMomentumBosonRestFrame()
			* photon.getFourMomentumBosonRestFrame();
	double plQ = ql * event.getBoson().getFourMomentumBosonRestFrame();
	//! Soft photon factor from 4-momenta
	double Sfac = (2.0 * plQ - QWk / qlk * aml2 - qlk / QWk * amQ2) / qlk / QWk;

	//! Soft photon factor from spin amplitudes
	double xsof = 0.0;
	for (int kap = 0; kap < 2; ++kap) {//              !<-- loop over photon polarizations
		//! Polarization vectors of photon in rectangular basis
		CLHEP::HepLorentzVector epsilonK = VecPol(
				photon.getFourMomentumBosonRestFrame(), kap + 1);
		xsof = xsof + pow(((ql * epsilonK) / qlk
				- (event.getBoson().getFourMomentumBosonRestFrame() * epsilonK)
						/ QWk), 2);
	}
	double xmsof = 4 * ds->pi * ds->alphaQED * xmbo * xsof;
	double xmbos = 4 * ds->pi * ds->alphaQED * xmbo * Sfac;
	//! Difference of soft-photon factors from 4-momenta and spin amplitudes.
	double xmsdif = xmbos - xmsof;
	//*!......................................................................
	//*! Corrected matrix element
	result = result + xmsdif;
#ifdef DEBUG_FLAG
	int keyWid = 0;
	double amW = ds->MW;
	double GaW = ds->GW;
	double sW2 = ds->sinThetaW2;
	double Uij = ckmElement;
	double p1[4];
	double p2[4];
	double Q[4];
	double q1[4];
	double q2[4];
	double pk[4];
	for (int i = 0; i < 4; ++i) {
		p1[i] = event.getQuark().getFourMomentumBosonRestFrame()[i];
		p2[i] = event.getAntiQuark().getFourMomentumBosonRestFrame()[i];
		Q[i] = event.getBoson().getFourMomentumBosonRestFrame()[i];
		pk[i] = photon.getFourMomentumBosonRestFrame()[i];
		if (event.getLepton().getPDGid() > 0) {
			q1[i] = event.getLepton().getFourMomentumBosonRestFrame()[i];
			q2[i] = event.getNeutrino().getFourMomentumBosonRestFrame()[i];
		} else {
			q2[i] = event.getLepton().getFourMomentumBosonRestFrame()[i];
			q1[i] = event.getNeutrino().getFourMomentumBosonRestFrame()[i];
		}
	}
	double fortranxmat1 = wh_matfsr1_(&keyWid, &amW, &GaW, &sW2, &Uij, &QW, p1,
			p2, Q, q1, q2, pk);

	std::cout << "DEBUG: ModelHandler::matrixElementFSR1() # fortran: "
			<< fortranxmat1 << std::endl;
	std::cout << "DEBUG: ModelHandler::matrixElementFSR1() # cpp: " << result
			<< std::endl;
	std::cout << "DEBUG: ModelHandler::matrixElementFSR1() # cpp/fortran: "
			<< result / fortranxmat1 << std::endl;
#endif
	return result;
}

bool VINHAC::ModelHandler::eventEvolve(VINHAC::Event& event) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::eventEvolve()" << std::endl;

#endif
	///////////////////////////////////////////////////////////////////////
	// Model weights for the process:                                    //
	//                   q + q~ --> Z --> l + l~ + n*gamma               //
	// OUTPUT: WtMdl - the best model weight for Mode=0, dummy otherwise //
	//-------------------------------------------------------------------//
	///////////////////////////////////////////////////////////////////////


	double weight = 1.0;

	// Check if prievious weights were nonzero

#ifdef DEBUG_FLAG

	std::cout
			<< "DEBUG: ModelHandler::eventEvolve() event.currentCommonWtMultiplication()"
			<< event.currentCommonWtMultiplication() << std::endl;

#endif

	double WtCrud = event.currentCommonWtMultiplication();
	if (WtCrud == 0.0) {
		//event.addModelWeight("ModelBest",WtSet);
		map<string, double> modelWeights = event.getModelWeightsMap();
		for (map<string, double>::iterator it = modelWeights.begin(); it
				!= modelWeights.end(); ++it) {
			event.addModelWeight((*it).first, 0.0);
		}
		return false;
	}

	prepareParticles(event);

	double xmat0 = matrixElement0(event);

	double sqq = event.getBoson().getFourMomentumBosonRestFrame().e()
			* event.getBoson().getFourMomentumBosonRestFrame().e();
	double dsigma0 = normalizationFactor0 / sqq * xmat0;

	double dsigmaCrud = 2.0 * differentialXsectionCrude(sqq, event);
	event.setDifferentialXsectionCrude(dsigmaCrud);

	if (ds->keyRad == 0) {

		weight *= dsigma0 / dsigmaCrud;
		event.addModelWeight("ModelBorn", weight);
		event.addModelWeight("ModelBest", weight);
#ifdef DEBUG_FLAG
		std::cout << "DEBUG: ModelHandler::eventEvolve() # dsig0 dsigcru: "
				<< dsigma0 << " " << dsigmaCrud << endl;
		std::cout << "DEBUG: ModelHandler::eventEvolve() # ModelBorn weight: "
				<< weight << std::endl;

#endif
		return weight != 0.0;

	} else {
		double deltaQED = 0.0;
		double deltaEW = 0.0;

		//---------------------------------------------------------------------
		//------------------- O(alpha) EW CORRECTIONS--------------------------
		//---------------------------------------------------------------------
		// Virtual+soft QED-like corrections for leptonic W decays
		if (ds->keyRad <= 2) {
			deltaQED = ModelHandler::deltaQED(
					event.getBoson().getFourMomentumBosonRestFrame().e(),
					ds->masses[abs(event.getLepton().getPDGid())]);

			deltaEW = deltaQED + deltaWeak[abs(event.getLepton().getPDGid())];
		} else {
			throw VinhacException("keyRad > 2 not implemented yet");
		}

		// LL approximation
		double alfpi = ds->alphaQED / ds->pi;
		double gamLL = 2.0 * alfpi * (2 * log(
				event.getBoson().getFourMomentumBosonRestFrame().e()
						/ ds->masses[abs(event.getLepton().getPDGid())]) - 1.0);
		double deltaLL = 0.25 * gamLL;

		//======================= BETA0 BETA0 BETA0 =============================
		// beta0:
		double beta00 = dsigma0;
		double beta01_QED = beta00 * (1 + deltaQED);
		double beta01_EW = beta00 * (1 + deltaEW);
		double beta01_LL = beta00 * (1 + deltaLL);

		//======================= BETA1 BETA1 BETA1 =============================
		// beta1:
		double beta11 = 0.0;
		double beta11_LL = 0.0;

		double dsig1 = 0.0;
		double Stil = 0.0;
		for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
			// Single hard-photon matrix element
			// 1. FSR only
			double xmat1 = matrixElementFSR1(event, event.getPhotons()[i]);
			if (ds->keyRad >= 3) {
				// 2. FSR + iterference
				throw VinhacException("keyRad >= 3 not implemented yet");
			}

			//--- Normalized hard photon cross section
			dsig1 = normalizationFactor1 / sqq * xmat1;

			// Soft factor S-tilde:

			if (ds->keyRad <= 2) {
				// 1. FSR only
				Stil = StPair(sqq, pow(ds->masses[abs(
						event.getLepton().getPDGid())], 2),
						event.getBoson().getFourMomentumBosonRestFrame(),
						event.getLepton().getFourMomentumBosonRestFrame(),
						event.getPhotons()[i].getFourMomentumBosonRestFrame());
			} else {
				// 2. FSR + iterference
				throw VinhacException("Keyrad >2 not implemented yet");
			}

			// beta1 O(alpha1)
			double beta1i = dsig1 / Stil - beta00;
			beta11 = beta11 + beta1i;
			// LL approximation
			double z = 2.0
					* event.getPhotons()[i].getFourMomentumBosonRestFrame().e()
					/ event.getBoson().getFourMomentumBosonRestFrame().e();
			double RadLL = pow(z, 2) / 2.0 / (1.0 - z);
			double beta1i_LL = RadLL * beta00;
			beta11_LL = beta11_LL + beta1i_LL;

		} //  ! <-- photon loop


		// Check on betas - if not huge
		double Huge = 1.0E3;
		double FacWcr = WtCrud / dsigmaCrud;

		if (abs(beta00 * FacWcr) > Huge) {
			//Nbt0 = Nbt0 + 1
			// Print-outs if LevPrn > 0

			beta00 = 0.0;
			//TODO Printout;


		}
		double beta11tmp = beta11;

		if (abs((beta01_EW + beta11) * FacWcr) > Huge) {
			//Nbt1_EW = Nbt1_EW + 1
			// Print-outs if LevPrn > 1

			beta01_EW = 0.0;
			beta11 = 0.0;
			dsig1 = 0.0;

		}

		if (abs((beta01_QED + beta11tmp) * FacWcr) > Huge) {
			//Nbt1_QED = Nbt1_QED + 1
			// Print-outs if LevPrn > 1

			beta01_QED = 0.0;
			beta11 = 0.0;
			dsig1 = 0.0;
		}

		if (abs((beta01_LL + beta11_LL) * FacWcr) > Huge) {
			//Nbt1_LL = Nbt1_LL + 1
			// Print-outs if LevPrn > 1

			beta01_LL = 0.0;
			beta11_LL = 0.0;
		}

		//=======================================================================
		//=============== FIXED ORDER, NO EXPONENTIATION ========================
		//=======================================================================
		double xs00 = 0.0;
		double xs0gQED = 0.0;
		double xs0gEW = 0.0;
		double xs1g = 0.0;
		if (event.getPhotons().size() == 0) {
			double BtiReb = 0.0;
			double delvsEW = 0.0;
			double delvsQED = 0.0;
			if (ds->keyRad <= 2) {
				//--- YFS FormFactor and QED/EW corrections
				BtiReb = YWldec(
						event.getBoson().getFourMomentumBosonRestFrame(),
						event.getLepton().getFourMomentumBosonRestFrame(),
						ds->EgMin);
				delvsQED = VirSofQED(
						event.getBoson().getFourMomentumBosonRestFrame().e(),
						ds->masses[abs(event.getLepton().getPDGid())],
						ds->EgMin);
				delvsEW = delvsQED
						+ deltaWeak[abs(event.getLepton().getPDGid())];
			} else {
				throw VinhacException("keyrad > 2 not implemented yet");
			}
			double YFSfmf = exp(BtiReb);
			// Born
			xs00 = beta00 / YFSfmf;
			// Born + virtual + soft
			xs0gQED = xs00 * (1 + delvsQED);// ! pure QED
			xs0gEW = xs00 * (1 + delvsEW);// ! EW corrections

		} else if (event.getPhotons().size() == 1) {
			// 1 real hard photon
			//--- YFS FormFactor
			double BtiReb = 0.0;
			if (ds->keyRad <= 2) {
				BtiReb = YWldec(
						event.getBoson().getFourMomentumBosonRestFrame(),
						event.getLepton().getFourMomentumBosonRestFrame(),
						ds->EgMin);
			} else {
				throw VinhacException("keyrad > 2 not implemented yet");
			}
			double YFSfmf = exp(BtiReb);
			xs1g = dsig1 / Stil / YFSfmf;

		}

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++ WEIGHTS, WEIGHTS, WEIGHTS ++++++++++++++++++++++++
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// Principal weights:
		//----------------------------------------------------------------------
		//--- O(alpha^0)_exp
		event.addModelWeight("beta00", beta00 / dsigmaCrud);
		//--- O(alpha^1)_exp: QED only
		event.addModelWeight("beta01_QED + beta11", (beta01_QED + beta11)
				/ dsigmaCrud);
		//--- O(alpha^1)_exp: EW corrections
		event.addModelWeight("beta01_EW + beta11", (beta01_EW + beta11)
				/ dsigmaCrud);
		event.addModelWeight("ModelBest", (beta01_EW + beta11) / dsigmaCrud);
		weight = (beta01_EW + beta11) / dsigmaCrud;
		//--- O(alpha^1)_exp: LL QED only
		event.addModelWeight("beta01_LL + beta11_LL", (beta01_LL + beta11_LL)
				/ dsigmaCrud);
		//----------------------------------------------------------------------
		// Individual betas...
		//--- O(alpha^0) beta00
		event.addModelWeight("beta00", beta00 / dsigmaCrud);
		//--- O(alpha^1) beta01_QED
		event.addModelWeight("beta01_QED", beta01_QED / dsigmaCrud);
		//--- O(alpha^1) beta11
		event.addModelWeight("beta11", beta11 / dsigmaCrud);
		//--- O(alpha^1) beta01_EW
		event.addModelWeight("beta01_EW", beta01_EW / dsigmaCrud);
		//--- O(alpha^1) beta01_LL
		event.addModelWeight("beta01_LL", beta01_LL / dsigmaCrud);
		//--- O(alpha^1) beta11_LL
		event.addModelWeight("beta11_LL", beta11_LL / dsigmaCrud);
		//========================================================================
		//================ FIXED ORDER, NO EXPONENTIATION ========================
		//========================================================================
		//--- Born
		event.addModelWeight("xs00", xs00 / dsigmaCrud);
		//--- O(alpha^1)_QED
		event.addModelWeight("xs0gQED + xs1g", (xs0gQED + xs1g) / dsigmaCrud);
		//--- O(alpha^1)_EW
		event.addModelWeight("xs0gEW + xs1g", (xs0gEW + xs1g) / dsigmaCrud);
		//--- Born + virtual_QED + soft photon corr.
		event.addModelWeight("xs0gQED", xs0gQED / dsigmaCrud);
		//--- Born + virtual_EW + soft photon corr.
		event.addModelWeight("xs0gEW", xs0gEW / dsigmaCrud);
		//--- One real hard photon
		event.addModelWeight("xs1g", xs1g / dsigmaCrud);

		//=======================================================================
		//***********************************************************************

		// Return best weight
		return weight != 0.0;
		//***********************************************************************
	}

	return weight != 0.0;
}

CLHEP::HepLorentzVector VINHAC::ModelHandler::VecPol(CLHEP::HepLorentzVector p,
		int lambda) {
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

	CLHEP::HepLorentzVector eps; // costructor Sets (0., 0., 0., 0.)
	double qt2, qt, aq2, aq, am, ws3;
	qt2 = p.perp2(); // x^2 + y^2
	qt = sqrt(qt2); // sqrt(x^2 + y^2)
	aq2 = qt2 + p.z() * p.z();
	aq = sqrt(aq2); // x^2 + y^2 + z^2
	am = sqrt(fabs(p.e() * p.e() - aq2));// (*this).restMass();

	// lambda polarization state = 1
	if (lambda == 1) {
		if (aq < 1e-8) { // <--- small Rho() approximation
			eps.setX(1.0);
			// cout << "test_male aq" << endl;
		} else if (qt < 1e-8) { // <---  small Perp() approximation
			eps.setX(p.z() / aq);
			//cout << "test male qt" << endl;
		} else { // eq. (3.47a)
			eps.setX(p.x() / qt * p.z() / aq);
			eps.setY(p.y() / qt * p.z() / aq);
			eps.setZ(-qt / aq);
			//eps.SetT(0.); nie trzeba bo na poczatku wyzerowany
		}
	}
	// lambda polarization state = 2
	else if (lambda == 2) {
		if (aq < 1e-8) { // <--- small Rho() approximation
			eps.setY(1e0);
		} else if (qt < 1e-8) { // <--- small Perp() approximation
			eps.setY(1.0);
		} else { // eq. (3.47b)
			eps.setX(-p.y() / qt);
			eps.setY(p.x() / qt);
			eps.setZ(0.);
			eps.setT(0.);
		}
	}
	// lambda polarization state = 3
	else if (lambda == 3) {
		if (aq < 1e-8) {
			if (am > 1e-5) {
				eps.setZ(1.);
				//cout << "test1" << endl;
			}
		} else if (qt < 1e-8) {
			if (am > 1e-5) {
				eps.setZ(p.t() / am / aq * p.z());
				eps.setT(aq / am);
				//cout << "test2" << endl;
			}
		} else {
			if (am > 1e-5) {
				ws3 = p.t() / am / aq;
				eps.setX(ws3 * p.x());
				eps.setY(ws3 * p.y());
				eps.setZ(ws3 * p.z());
				eps.setT(ws3 * aq2 / p.t());
				//cout << "test3" << endl;
			}
		}
	} else {
		throw VinhacException("Zla wartosc lambda");

	}
	return eps;
}

// Additional heg vectors stuff
void VINHAC::ModelHandler::HelEig(CLHEP::HepLorentzVector p, int hel, complex<
		double> chip[2]) {
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
	chip[0] = complex<double> (0., 0.);
	chip[1] = complex<double> (0., 0.);
	pmp3 = p.rho() + p.z(); // <--- Rho() = sqrt ( x^2 + y^2 + z^2)
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
		}
		else { // when p_z = -|\vec{p}|
			chip[0] = complex<double> (0., 0.);
			chip[1] = complex<double> (1., 0.);
		}
	} else {
		throw VinhacException(" heleig: Wrong fermion helicity");
	}
}

complex<double> VINHAC::ModelHandler::spinorialFunction(
		complex<double> hip1[2], CLHEP::HepLorentzVector a,
		complex<double> hip2[2], int ialfa) {

	complex<double> as[2][2];
	complex<double> vx[2];

	MvtoWm(a, as, ialfa);

	MxV2dC(as, hip2, vx);

	return conj(hip1[0]) * vx[0] + conj(hip1[1]) * vx[1];
}

complex<double> VINHAC::ModelHandler::spinorialFunction(
		complex<double> hip1[2], CLHEP::HepLorentzVector a,
		CLHEP::HepLorentzVector b, CLHEP::HepLorentzVector c,
		complex<double> hip2[2], int ialfa) {
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
	MvtoWm(a, as, ialfa);
	MvtoWm(b, bs, -ialfa);
	MvtoWm(c, cs, ialfa);
	//! Multiply Weyl matrix by Pauli spinor, etc.
	MxV2dC(cs, hip2, vx);
	MxV2dC(bs, vx, vx);
	MxV2dC(as, vx, vx);
	//! S(p1,a1,a2,a3,p2)_{lambda1,lambda2}^alpha
	return conj(hip1[0]) * vx[0] + conj(hip1[1]) * vx[1];

}

void VINHAC::ModelHandler::MxV2dC(complex<double> matrix[2][2],
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

void VINHAC::ModelHandler::MvtoWm(CLHEP::HepLorentzVector& p,
		complex<double> as[2][2], int alpha) {

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

double VINHAC::ModelHandler::inverseWPropagator(double s, DataSource *ds) {
	double result = 0.0;
	if (ds->keyWid == 0) {
		result = (s - (ds->MW) * (ds->MW)) * (s - (ds->MW) * (ds->MW))
				+ ((ds->MW) * (ds->GW)) * ((ds->MW) * (ds->GW));
	} else {
		result = (s - (ds->MW) * (ds->MW)) * (s - (ds->MW) * (ds->MW)) + (s * s
				* (ds->GW) * (ds->GW)) / ((ds->MW) * (ds->MW));
	}

	return result;
}

double VINHAC::ModelHandler::differentialXsectionCrude(double s, Event& event) {

	double invProp = inverseWPropagator(s, ds);
	double dsig;
	double QW = 1.0;
	double costheta = event.getCosTheta();
	if (event.getBoson().getPDGid() < 0) {
		QW = -1.0;
	}
	dsig = (ckmElement / 8.0 / ds->sinThetaW2 / ds->invAlphaQED) * (ckmElement
			/ 8.0 / ds->sinThetaW2 / ds->invAlphaQED) * s / invProp * (1.0 - QW
			* costheta);
	if (ds->keyPol == 0) {
		dsig *= (1.0 - QW * costheta);
	}

	return dsig / 3.0;
}

double VINHAC::ModelHandler::YWldec(CLHEP::HepLorentzVector Q,
		CLHEP::HepLorentzVector ql, double aKmax) {

	// Dummy photon mass
	double amg = 1.0E-20;
	//!
	double aM = Q.mag();
	double aml = ql.mag();
	double t = (Q - ql).mag2();
	// YFS IR factor Y, i.e. sum of virtual and real photon IR functions
	double Bvirt = ReBtug(t, aM, aml, amg);
	double Breal = Btilde(Q, ql, aM, aml, aKmax, amg);

	return Bvirt + Breal;
}

double VINHAC::ModelHandler::ReBtug(double t, double xm1, double xm2,
		double amg) {


	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / alfinv / pi;
	double epst = 1.E-6;
	double am1, am2;

	if (xm1 > xm2) {
		am1 = xm1;
		am2 = xm2;
	} else {
		am1 = xm2;
		am2 = xm1;
	}
	double am1s = pow(am1, 2);
	double am2s = pow(am2, 2);
	double am12 = am1 * am2;
	double amds = (am1 - am2) * (am1 + am2);
	// IR regulator LOG
	double RegLog = -2 * log(sqrt(am12) / amg);
	double ReB2pi;
	//--------------------------------------------------------------------
	// Non-zero t case:
	if (abs(t) > epst) {
		double xnu = (am1s + am2s - t) / 2;
		double xla = sqrt((xnu + am12) * (xnu - am12));
		double xlanu = xla + xnu;
		// Function A(p1,p2)
		double A = log(xlanu / am12) / xla;
		// Function A1(p1,p2)
		double A1 = amds / t * log(am1 / am2) - 2 * pow(xla, 2) / t * A - 2;
		// Function A3(p1,p2)
		double tix = t / (2 * xla * (xlanu - am2s));
		double eta = am2s * tix;
		double zet = xlanu * tix;
		double xlg1 = log(xlanu / am1s);
		double xlg2 = log(xlanu / am2s);
		double xar1 = abs((xla - xnu + am1s) / t);
		double xar2 = abs((xlanu - am2s) / am2s);
		double xlet = log(abs(eta)) * log(1 + eta);
		double xlze = log(abs(zet)) * log(1 + zet);
		double A3 = log(2 * xla / am12) * A + (0.25 * (xlg1 + 2 * log(xar1))
				* xlg1 + 0.25 * (xlg2 - 2 * log(xar2)) * xlg2 + 0.5 * (xlet
				- xlze) + DILOGY(-eta) - DILOGY(-zet)) / xla;
		//	Function ReB(s)
		ReB2pi = (xnu * A - 1) * RegLog + 0.5 * A1 - xnu * A3;
		//--------------------------------------------------------------------
		// Zero t case:
	} else {
		ReB2pi = ((am1s + am2s) / amds * log(am1 / am2) - 1) * (RegLog + 0.5);
	}
	return alfpi * ReB2pi;
}

double VINHAC::ModelHandler::DILOGY(double X) {
	//C-------------------------------------------- REMARKS ---------------
	//C DILOGARITHM FUNCTION: DILOG(X)=INT( -LN(1-Z)/Z ) , 0 < Z < X .
	//C THIS IS THE CERNLIB VERSION.
	//C--------------------------------------------------------------------
	//IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	double Z = -1.644934066848226;
	double T, S;

	if (X < -1.0) {
		T = 1.0 / X;
		S = -0.5;
		Z = Z - 0.5 * pow(log(abs(X)), 2);
	} else if (X <= 0.5) {
		T = X;
		S = 0.5;
		Z = 0.0;
	} else if (X == 1.0) {
		return 1.644934066848226;
	} else if (X <= 2.0) {
		T = 1.0 - X;
		S = -0.5;
		Z = 1.644934066848226 - log(X) * log(abs(T));
	} else {
		Z = 3.289868133696453;
		T = 1.0 / X;
		S = -0.5;
		Z = Z - 0.5 * pow(log(abs(X)), 2);
	}
	double Y = 2.666666666666667 * T + 0.666666666666667;
	double B = 0.000000000000001;
	double A = Y * B + 0.000000000000004;
	B = Y * A - B + 0.000000000000011;
	A = Y * B - A + 0.000000000000037;
	B = Y * A - B + 0.000000000000121;
	A = Y * B - A + 0.000000000000398;
	B = Y * A - B + 0.000000000001312;
	A = Y * B - A + 0.000000000004342;
	B = Y * A - B + 0.000000000014437;
	A = Y * B - A + 0.000000000048274;
	B = Y * A - B + 0.000000000162421;
	A = Y * B - A + 0.000000000550291;
	B = Y * A - B + 0.000000001879117;
	A = Y * B - A + 0.000000006474338;
	B = Y * A - B + 0.000000022536705;
	A = Y * B - A + 0.000000079387055;
	B = Y * A - B + 0.000000283575385;
	A = Y * B - A + 0.000001029904264;
	B = Y * A - B + 0.000003816329463;
	A = Y * B - A + 0.000014496300557;
	B = Y * A - B + 0.000056817822718;
	A = Y * B - A + 0.000232002196094;
	B = Y * A - B + 0.001001627496164;
	A = Y * B - A + 0.004686361959447;
	B = Y * A - B + 0.024879322924228;
	A = Y * B - A + 0.166073032927855;
	A = Y * A - B + 1.935064300869969;
	return S * T * (A - B) + Z;
}

double VINHAC::ModelHandler::Btilde(CLHEP::HepLorentzVector p1,
		CLHEP::HepLorentzVector p2, double am1, double am2, double aKmax,
		double amgam) {


	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / alfinv / pi;
	double epsb = 1.0E-10;

	double E1 = p1.e();
	double E2 = p2.e();
	double am12 = am1 * am2;
	double p1p2 = p1 * p2;
	if (p1p2 - am12 < 1.0E-10)
		return 0;
	double xlam = sqrt((p1p2 - am12) * (p1p2 + am12));
	// Function A(p1,p2)
	double A = 1.0 / xlam * log((p1p2 + xlam) / am12);
	// betas and big logs - limit beta->0 has to be treated carefully!
	double bt1s = abs((1.0 - am1 / E1) * (1.0 + am1 / E1));
	double bet1 = 0, b1ln = 0;
	if (bt1s > epsb) {
		bet1 = sqrt(bt1s);
		b1ln = log((1.0 + bet1) * E1 / am1) / bet1;
	} else {
		b1ln = 1.0;
	}
	double bt2s = abs((1.0 - am2 / E2) * (1.0 + am2 / E2));
	double bet2 = 0, b2ln = 0;
	if (bt2s > epsb) {
		bet2 = sqrt(bt2s);
		b2ln = log((1 + bet2) * E2 / am2) / bet2;
	} else {
		b2ln = 1.0;
	}
	// Function A4(p1,p2)
	double A4 = ModelHandler::A4anal(p1p2, E1, E2, am1, am2);

	// B-tilde(p1,p2;aKmax,amgam)
	double Btian = (p1p2 * A - 1) * 2 * log(2 * aKmax / amgam) + b1ln + b2ln
			+ p1p2 * A4;
	return alfpi * Btian;
}

double VINHAC::ModelHandler::A4anal(double p1p2, double E1, double E2,
		double am1, double am2) {


	bool RF1, RF2;
	double A4=0.0, A4p, A4z;
	//
	double Q2 = 2 * p1p2 - pow(am1, 2) - pow(am2, 2);
	// Check if momenta in rest frame of p1 or p2
	RF1 = (1 - am1 / E1) > 1.0E-10;
	RF2 = (1 - am2 / E2) > 1.0E-10;
	// Q2 <> 0
	if (abs(Q2) > 1.0E-10) {
		//     ==========================
		//... NOT in p1/p2 rest frame - usual case
		if (RF1 && RF2) {
			// Q2 > 0 (usual case)
			if (Q2 > 0.0) {
				A4 = ModelHandler::A4Q2pos(p1p2, E1, E2, am1, am2);
				// Q2 < 0 (e.g. radiative leptonic W-deays):
			} else if (Q2 < 0.0) {
				A4 = ModelHandler::A4Q2neg(p1p2, E1, E2, am1, am2);
				A4p = ModelHandler::A4Q2pos(p1p2, E1, E2, am1, am2);
				A4z = ModelHandler::A4Q2zer(E1, E2, am1, am2);
			}
			//... In p1 rest frame
		} else if (RF2) {
			A4 = ModelHandler::A4Q2fRF(E2, am2, am1);
			// In p2 rest frame
		} else if (RF1) {
			A4 = ModelHandler::A4Q2fRF(E1, am1, am2);
		} else {
			A4 = 0.0;
		}
		// Q2=0 case (e.g. leptonic W-decays):
	} else {
		//     ====
		if (RF1 && RF2) {
			A4 = ModelHandler::A4Q2zer(E1, E2, am1, am2);
			//... In p1 rest frame
		} else if (RF2) {
			A4 = ModelHandler::A4Q2zRF(am2, am1);
			//... In p2 rest frame
		} else if (RF1) {
			A4 = ModelHandler::A4Q2zRF(am1, am2);
		} else {
			A4 = 0.0;
		}
	}
	//     =====
	return A4;
}

double VINHAC::ModelHandler::A4Q2pos(double p1p2, double E1, double E2,
		double am1, double am2) {



	double p1s = pow(E1, 2) - pow(am1, 2);
	double p2s = pow(E2, 2) - pow(am2, 2);
	if (p1s < p2s) {
		double tmp = am1;
		am1 = am2;
		am2 = tmp;
		tmp = E1;
		E1 = E2;
		E2 = tmp;
	}
	double Ep = E1 + E2;
	double Em = E1 - E2;
	double sm = am1 + am2;
	double dm = am1 - am2;
	double Q2 = 2 * p1p2 - pow(am1, 2) - pow(am2, 2);
	double xl = sqrt((Q2 + pow(sm, 2)) * (Q2 + pow(dm, 2)));
	double xq = sqrt(Q2 + pow(Em, 2));
	double qp = xq + Em;
	double qm = xq - Em;
	double et0 = sqrt(pow(E2, 2) - pow(am2, 2));
	if (p1p2 > E1 * E2)
		et0 = -et0;
	double et1 = sqrt(pow(E1, 2) - pow(am1, 2)) + xq;
	double y1 = 0.5 * ((xq - Ep) + (sm * dm + xl) / qp);
	double y2 = y1 - xl / qp;
	double y3 = 0.5 * ((xq + Ep) + (sm * dm + xl) / qm);
	double y4 = y3 - xl / qm;
	// Some auxiliary functions
	double Eln;
	if (abs(Em) > 1.0E-10) {
		Eln = log(abs(qm / qp)) * (log(abs((et1 - y1) * (et1 - y4) / (et1 - y2)
				/ (et1 - y3))) - log(abs((et0 - y1) * (et0 - y4) / (et0 - y2)
				/ (et0 - y3))));
	} else {
		Eln = 0.0;
	}
	double Vet0 = ModelHandler::Yijeta(y1, y4, et0) + ModelHandler::Yijeta(y2,
			y1, et0) + ModelHandler::Yijeta(y3, y2, et0)
			- ModelHandler::Yijeta(y3, y4, et0) + 0.5 * etaln(y1, y2, y3, y4,
			et0) * etaln(y2, y3, y1, y4, et0);
	double Vet1 = ModelHandler::Yijeta(y1, y4, et1) + ModelHandler::Yijeta(y2,
			y1, et1) + ModelHandler::Yijeta(y3, y2, et1)
			- ModelHandler::Yijeta(y3, y4, et1) + 0.5 * etaln(y1, y2, y3, y4,
			et1) * etaln(y2, y3, y1, y4, et1);
	// Function A4(p1,p2)
	return (Eln + Vet1 - Vet0) / xl;
}

double VINHAC::ModelHandler::Yijeta(double yi, double yj, double eta) {

	return 2 * DILOGY((yj - yi) / (eta - yi)) + 0.5 * pow(log(abs((eta - yi)
			/ (eta - yj))), 2);
}

double VINHAC::ModelHandler::A4Q2neg(double p1p2, double E1, double E2,
		double am1, double am2) {


	// Some auxiliary variables
	if (am1 < am2) {
		double tmp = am1;
		am1 = am2;
		am2 = tmp;
		tmp = E1;
		E1 = E2;
		E2 = tmp;
	}
	double p1m = sqrt((E1 - am1) * (E1 + am1));
	double p2m = sqrt((E2 - am2) * (E2 + am2));

	double Ep = E1 + E2;
	double Em = E1 - E2;
	double sm = am1 + am2;
	double dm = am1 - am2;
	double Q2 = 2 * p1p2 - pow(am1, 2) - pow(am2, 2);
	double xl = sqrt((Q2 + pow(sm, 2)) * (Q2 + pow(dm, 2)));
	double xq = sqrt(Q2 + pow(Em, 2));
	double qp = xq + Em;
	double et0 = p2m / xq;
	double et1 = p1m / xq - 1;
	double y1 = -(xq + Ep + (sm * dm + xl) / Q2 * qp) / 2 / xq;
	double y2 = y1 + xl / Q2 * qp / xq;
	double y3 = -(xq - Ep + (sm * dm + xl) / qp) / 2 / xq;
	double y4 = y3 + xl / qp / xq;
	// Some auxiliary functions
	double Eln;
	if (abs(Em) > 1.0E-10) {
		Eln = log(abs(pow(qp, 2) / Q2)) * (etaln(y1, y4, y2, y3, et1) - etaln(
				y1, y4, y2, y3, et0));
	} else {
		Eln = 0.0;
	}
	double Vet0 = Yijeta(y1, y4, et0) + Yijeta(y2, y1, et0) + Yijeta(y3, y2,
			et0) - Yijeta(y3, y4, et0) + 0.5 * etaln(y1, y2, y3, y4, et0)
			* etaln(y2, y3, y1, y4, et0);
	double Vet1 = Yijeta(y1, y4, et1) + Yijeta(y2, y1, et1) + Yijeta(y3, y2,
			et1) - Yijeta(y3, y4, et1) + 0.5 * etaln(y1, y2, y3, y4, et1)
			* etaln(y2, y3, y1, y4, et1);
	// Function A4(p1,p2)
	double A4 = (Eln + Vet1 - Vet0) / xl;
	if (A4 < 1.0E6) {

	} else {
		//TODO    print*,'E1,am1=',E1,am1
		//    print*,'E2,am2=',E2,am2
		//   print*,'p1p2,E1E2=',p1p2,E1*E2
		//   print*,'et0,et1=',et0,et1
		//   print*,'y1,y2  =',y1,y2
		//   print*,'y3,y4  =',y3,y4
		//   print*,'Eln,Vet0,Vet1=',Eln,Vet0,Vet1
		//   print*,'A4=',A4
	}
	return A4;
}

double VINHAC::ModelHandler::A4Q2zer(double E1, double E2, double am1,
		double am2) {


	// Some auxiliary variables

	if (am1 < am2) {
		double tmp = am1;
		am1 = am2;
		am2 = tmp;
		tmp = E1;
		E1 = E2;
		E2 = tmp;
	}

	double am2s = pow(am2, 2);
	double Em = E1 - E2;
	double dm = (am1 + am2) * (am1 - am2);
	double xem = dm / 2 / pow(Em, 2);
	double p1m = sqrt((E1 - am1) * (E1 + am1));
	double p2m = sqrt((E2 - am2) * (E2 + am2));
	double et0 = p2m / Em;
	double et1 = p1m / Em - 1;
	double y1 = E2 / Em;
	double y2 = y1 - xem;
	double y3 = -y1 + 2 * am2s / dm;
	// eta_i - y_j, y_j - y_i
	double e0y1 = -am2s / Em / (E2 + p2m);
	double e0y2 = e0y1 + xem;
	double e0y3 = et0 - y3;
	double e1y1 = et1 - y1;
	double e1y2 = et1 - y2;
	double e1y3 = et1 - y3;
	double y2y1 = -xem;
	double y3y2 = y3 - y2;
	// Some auxiliary functions
	double Eln = -log(xem) * log(abs(e1y1 / e0y1 * e0y2 / e1y2 * e0y3 / e1y3));
	double Vet0 = 0.5 * pow(log(abs(e0y1 * e0y2 / e0y3)), 2) + log(abs(e0y1))
			* log(abs(e0y1) / pow(e0y2, 2)) + 2 * DILOGY(y2y1 / e0y1) + 2
			* DILOGY(y3y2 / e0y2);
	double Vet1 = 0.5 * pow(log(abs(e1y1 * e1y2 / e1y3)), 2) + log(abs(e1y1))
			* log(abs(e1y1) / pow(e1y2, 2)) + 2 * DILOGY(y2y1 / e1y1) + 2
			* DILOGY(y3y2 / e1y2);
	// Function A4(p1,p2)
	double A4 = (Eln + Vet1 - Vet0) / dm;
	if (A4 < 1.0E5) {

	} else {

		//TODO
		/*c        A4 = 0
		 print*,'am1,am2=',am1,am2
		 print*,'E1,E2=',E1,E2
		 print*,'xem,et0,et1=',xem,et0,et1
		 print*,'y1,y2,y3=',y1,y2,y3
		 print*,'e0y1,e0y2,e0y3=',e0y1,e0y2,e0y3
		 print*,'e1y1,e1y2,e1y3=',e1y1,e1y2,e1y3
		 print*,'y2y1,y3y1=',y2y1,y3y1
		 print*,'Eln,Vet0,Vet1=',Eln,Vet0,Vet1
		 print*,'A4=',A4
		 A4qzero = A4*/
	}
	return A4;
}

double VINHAC::ModelHandler::A4Q2fRF(double E1, double am1, double am2) {


	// Some auxiliary variables
	double am1s = pow(am1, 2);
	double p1 = sqrt((E1 - am1) * (E1 + am1));
	double Epp = E1 + p1;
	double Emp = am1s / Epp;
	double xmm = am2 - Emp;
	double xmp = (am2 * (am2 - 2 * E1) + am1s) / xmm;
	double twp = 2 * p1;
	double tmp = twp * am2;
	double ymp = am2 * Epp - am1s;
	double ymm = am2 * Emp - am1s;
	// Function A4(p1,p2)

	double A4 = (log(abs(xmm / xmp)) * log(Epp / am2) - 2 * log(twp * xmm / am2
			/ am1) * log(Epp / am1) + 2 * DILOGY(Emp / am2) - 2 * DILOGY(Epp
			/ am2) + DILOGY(-xmp / twp) - DILOGY(xmm / twp) + DILOGY(ymp / tmp)
			- DILOGY(-ymm / tmp)) / tmp;
	return A4;
}

double VINHAC::ModelHandler::A4Q2zRF(double am1, double am2) {


	//
	// Some auxiliary variables
	double xm12 = (am1 - am2) * (am1 + am2);
	// Function A4(p1,p2)
	double A4 = -2 * (pow(log(am1 / am2), 2) + DILOGY(xm12 / pow(am1, 2)))
			/ xm12;
	return A4;
}

double VINHAC::ModelHandler::deltaQED(double amW, double aml) {

	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / pi / alfinv;
	//
	//	 Based on: Bardin & Passarino, "The standard model ...", Oxford, 1999.
	//	*WP      del = alfpi * ( LOG(amW/aml) - 1.5 )
	//	 Based on: Marciano & Sirlin, Phys. Rev. D8 (1973) 3612.
	double del = alfpi * (log(amW / aml) + 0.5);
	return del;
}

double VINHAC::ModelHandler::StPair(double am1s, double am2s,
		CLHEP::HepLorentzVector p1, CLHEP::HepLorentzVector p2,
		CLHEP::HepLorentzVector pk) {


	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / pi / alfinv;

	double p1p2 = p1 * p2;
	double p1k = p1 * pk;
	double p2k = p2 * pk;
	double SoFa = (p1p2 / p2k - am1s / p1k) / p1k + (p1p2 / p1k - am2s / p2k)
			/ p2k;

	return alfpi / pi / 4.0 * SoFa;
}

double VINHAC::ModelHandler::VirSofQED(double amW, double aml, double aKmax) {


	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / pi / alfinv;

	double Blog = log(amW / aml);
	double viso = 2.0 * (Blog - 1.0) * log(2.0 * aKmax / amW) + 1.5 * Blog
			- pow(pi, 2.0) / 6.0 + 1.0;
	return alfpi * viso;
}
