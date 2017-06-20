/*
 * MatrixElement.cpp
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#include "MatrixElement.h"
#include "amplitudes/Born.h"
#include "amplitudes/SingleRadiation.h"
#include "util/CommonUtil.h"
#include "util/VectorUtil.h"

namespace VINHAC {

MatrixElement::MatrixElement() {
}

MatrixElement::~MatrixElement() {

}

double VINHAC::MatrixElement::matrixElement0(
		const CLHEP::HepLorentzVector& quark,
		const CLHEP::HepLorentzVector& antiquark,
		const CLHEP::HepLorentzVector& boson,
		const CLHEP::HepLorentzVector& lepton,
		const CLHEP::HepLorentzVector& antilepton, double ckmElement) {
	std::set<PolarizationName::name> polarization;
	polarization.insert(PolarizationName::unpolarized);
	return matrixElement0(quark, antiquark, boson, lepton, antilepton, ckmElement, polarization)[PolarizationName::unpolarized];
}

std::map<PolarizationName::name,double> VINHAC::MatrixElement::matrixElement0(
		const CLHEP::HepLorentzVector& quark,
		const CLHEP::HepLorentzVector& antiquark,
		const CLHEP::HepLorentzVector& boson,
		const CLHEP::HepLorentzVector& lepton,
		const CLHEP::HepLorentzVector& antilepton, double ckmElement , std::set<PolarizationName::name> polarizations) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::matrixElement0()" << std::endl;

#endif
	std::map<PolarizationName::name,double>  result;
	SpecialMatrix prod = Born::bornWproduction(quark, antiquark, boson, ckmElement);
	SpecialMatrix dec = Born::bornWdecay(boson,lepton,antilepton);

	result = (prod.multiplyAndSquare(dec,polarizations));

	double initialSpinFactor = 0.25;
	double colourFac = 1.0 / 3.0;

	double s = boson.mag2();
	double invPropagator = CommonUtil::inverseWPropagator(s);
	for(std::map<PolarizationName::name,double>::iterator it = result.begin(); it!= result.end(); ++it){
		result[it->first] /= invPropagator;
		result[it->first] *= initialSpinFactor;
		result[it->first] *= colourFac;
	}

#ifdef DEBUG_FLAG
	int keyWid = 0;
	double amW = DataSource::get().MW;
	double GaW = DataSource::get().GW;
	double sW2 = DataSource::get().sinThetaW2;
	double Uij = ckmElement;
	double p1[4];
	double p2[4];
	double Q[4];
	double q1[4];
	double q2[4];
	for (int i = 0; i < 4; ++i) {
		p1[i] = event.getQuark().getFourMomentumBosonRestFrame()[i];
		p2[i] = event.getAntiQuark().getFourMomentumBosonRestFrame()[i];
		Q[i] = boson[i];
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

double VINHAC::MatrixElement::matrixElementFSR1(
		const CLHEP::HepLorentzVector& quark,
		const CLHEP::HepLorentzVector& antiquark,
		const CLHEP::HepLorentzVector& boson,
		const CLHEP::HepLorentzVector& lepton,
		const CLHEP::HepLorentzVector& antilepton, const CLHEP::HepLorentzVector& photon,
		int chargeOfBoson,
		double ckmElement) {

	VINHAC::SpecialMatrix prod = Born::bornWproduction(quark,antiquark,boson, ckmElement);
	VINHAC::SpecialMatrix dec = SingleRadiation::radiativeWdecay(boson, lepton, antilepton, photon,
			chargeOfBoson);

	double result = (prod * dec);

	const double initialSpinFactor = 0.25;
	const double colourFac = 1.0 / 3.0;

	double s = boson.mag2();
	result /= CommonUtil::inverseWPropagator(s);
	result *= initialSpinFactor;
	result *= colourFac;

	//*! Correction for very collinear radiation (big cancellations between
	//*! various contributions to spin amplitudes lead to lost of numerical
	//*! precision in these regions!)
	double xmbo = matrixElement0(quark, antiquark, boson, lepton , antilepton, ckmElement);

	//! Electric charge of fermion f2 in units of e
	int Qf1 = (chargeOfBoson - 1) / 2;

	CLHEP::HepLorentzVector ql = antilepton;
	if (Qf1 != 0) {
		ql = lepton;
	}

	double aml2(ql * ql);
	double amQ2(
			boson
					* boson);
	double qlk(ql * photon);
	double QWk(
			boson
					* photon);
	double plQ(ql * boson);
	//! Soft photon factor from 4-momenta
	double Sfac((2.0 * plQ - QWk / qlk * aml2 - qlk / QWk * amQ2) / qlk / QWk);

	//! Soft photon factor from spin amplitudes
	double xsof(0.0);
	CLHEP::HepLorentzVector epsilonK;
	for (int kap = 0; kap < 2; ++kap) {//              !<-- loop over photon polarizations
		//! Polarization vectors of photon in rectangular basis
		VectorUtil::VecPol(epsilonK, photon,
				kap + 1);
		xsof = xsof + pow(
				((ql * epsilonK) / qlk
						- (boson
								* epsilonK) / QWk), 2);
	}
	double xmsof(4);
	xmsof *= DataSource::get().pi;
	xmsof *= DataSource::get().alphaQED;
	xmsof *= xmbo;
	xmsof *= xsof;
	double xmbos(4);
	xmbos *= DataSource::get().pi;
	xmbos *= DataSource::get().alphaQED;
	xmbos *= xmbo;
	xmbos *= Sfac;
	//! Difference of soft-photon factors from 4-momenta and spin amplitudes.
	double xmsdif(xmbos);
	xmsdif -= xmsof;
	//*!......................................................................
	//*! Corrected matrix element
	result += xmsdif;
#ifdef DEBUG_FLAG
	int keyWid = 0;
	double amW = DataSource::get().MW;
	double GaW = DataSource::get().GW;
	double sW2 = DataSource::get().sinThetaW2;
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
		Q[i] = boson[i];
		pk[i] = photon[i];
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

}
