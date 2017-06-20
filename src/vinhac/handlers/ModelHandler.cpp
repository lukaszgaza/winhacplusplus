/*
 * ModelHandler.cpp
 *
 *  Created on: 26 Aug 2009
 *      Author: siudmy
 */
#include "ModelHandler.h"
#include "../../utils/SANC/wh_ewcdec.h"
#include "../input/DataSource.h"
#include "../core/Event.h"
#include "../core/Particle.h"
#include "../core/VinhacException.h"
#include "../model/MatrixElement.h"
#include "../model/util/YFSUtil.h"
#include "../model/util/CommonUtil.h"

#ifdef DEBUG_FLAG
#include "../utils/FortranFunction.h"
#endif

//#define DEBUG_FLAG_AMPL

VINHAC::ModelHandler::ModelHandler(VINHAC::HardProcessHandler* p_Hard) :
			p_HardMatrix(p_Hard),
			normalizationFactor0(
					1.0 / 8.0 / (2.0 * 2.0 * DataSource::get().pi
							* DataSource::get().pi)),
			normalizationFactor1(
					1.0 / 16.0 / pow(2.0 * DataSource::get().pi, 5)),
			alfpi(DataSource::get().alphaQED / DataSource::get().pi),
			overWeightedEvents(0)

{
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::ModelHandler()" << std::endl;

#endif

	if (DataSource::get().keyRad == 2) {
		double mw = DataSource::get().MW;
		int id = 11;
		int KeyGmu = DataSource::get().keyGmu;
		double aMW = DataSource::get().MW;
		double aMZ = DataSource::get().MZ;
		double aMH = DataSource::get().MH;
		double Xmfe[20];
		double Xckm[3][3];
		for (unsigned i = 0; i < 20; i++) {
			Xmfe[i] = DataSource::get().masses[i + 1];
		}

		for (unsigned i = 0; i < DataSource::get().ckm.size(); ++i) {
			for (unsigned j = 0; j < DataSource::get().ckm[i].size(); ++j) {
				Xckm[j][i] = DataSource::get().ckm[i][j];
			}
		}

		deltaWeak[11] = wh_delewdec_(&mw, &id, &KeyGmu, &aMW, &aMZ, &aMH, Xmfe,
				Xckm);
		id = 13;
		deltaWeak[13] = wh_delewdec_(&mw, &id, &KeyGmu, &aMW, &aMZ, &aMH, Xmfe,
				Xckm);
		id = 15;
		deltaWeak[15] = wh_delewdec_(&mw, &id, &KeyGmu, &aMW, &aMZ, &aMH, Xmfe,
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
	// get initial quark and anti quark

	int iu;
	int id;
	if (event.getBoson().getPDGid() == -24) {
		iu = abs(event.getAntiQuark().getPDGid()) / 2 - 1;
		id = event.getQuark().getPDGid() / 2;
	} else {
		id = abs(event.getAntiQuark().getPDGid()) / 2;
		iu = event.getQuark().getPDGid() / 2 - 1;
	}

	ckmElement = DataSource::get().ckm[iu][id];

	quark = event.getQuark().getFourMomentum(DataSource::get().modelReferenceFrame);
	antiquark = event.getAntiQuark().getFourMomentum(DataSource::get().modelReferenceFrame);
	boson = event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame);
	if (event.getLepton().getPDGid() > 0) {
		lepton = event.getLepton().getFourMomentum(DataSource::get().modelReferenceFrame);
		antilepton = event.getNeutrino().getFourMomentum(DataSource::get().modelReferenceFrame);
	} else {
		antilepton = event.getLepton().getFourMomentum(DataSource::get().modelReferenceFrame);
		lepton = event.getNeutrino().getFourMomentum(DataSource::get().modelReferenceFrame);
	}
	chargeOfBoson = 1;
	if (event.getBoson().getPDGid() < 0) {
		chargeOfBoson = -1;
	}

#ifdef DEBUG_FLAG

	std::cout << "DEBUG: ModelHandler::prepareParticles() # ckm[" << iu << "]["
	<< id << "] " << DataSource::get().ckm[iu][id] << std::endl;

#endif
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


	std::map<PolarizationName::name, double> weight;

	for (set<PolarizationName::name>::const_iterator polarization =
			DataSource::get().polarizations.begin(); polarization
			!= DataSource::get().polarizations.end(); ++polarization) {
		weight[*polarization] = 1.0;
	}

	// Check if prievious weights were nonzero

#ifdef DEBUG_FLAG

	std::cout
	<< "DEBUG: ModelHandler::eventEvolve() event.currentCommonWtMultiplication()"
	<< event.currentCommonWtMultiplication() << std::endl;

#endif

	double WtCrud(event.currentCommonWtMultiplication());
	if (WtCrud == 0.0) {
		//event.addModelWeight("ModelBest",WtSet);
		const map<pair<WeightName::name, PolarizationName::name> , double>
				& modelWeights = event.getModelWeightsMap();
		for (map<pair<WeightName::name, PolarizationName::name> , double>::const_iterator
				it = modelWeights.begin(); it != modelWeights.end(); ++it) {
			event.addModelWeight((*it).first, 0.0);
		}
		event.addModelWeight(DataSource::get().mainWeight, 0.0);
		return false;
	}

	prepareParticles(event);

	const std::map<PolarizationName::name, double>& xmat0 =
			MatrixElement::matrixElement0(quark, antiquark, boson, lepton,
					antilepton, ckmElement, DataSource::get().polarizations);

	double sqq(
			event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame)
					* event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame));
	std::map<PolarizationName::name, double> dsigma0;
	for (std::map<PolarizationName::name, double>::const_iterator it =
			xmat0.begin(); it != xmat0.end(); ++it) {
		dsigma0[it->first] = normalizationFactor0 * it->second / sqq;
	}

	double dsigmaCrud(2.0);
	dsigmaCrud *= differentialXsectionCrude(sqq, event);
	event.setDifferentialXsectionCrude(dsigmaCrud);

	if (DataSource::get().keyRad == 0) {
		for (set<PolarizationName::name>::const_iterator polarization =
				DataSource::get().polarizations.begin(); polarization
				!= DataSource::get().polarizations.end(); ++polarization) {
			weight[*polarization] *= dsigma0[*polarization] / dsigmaCrud;
			event.addModelWeight(
					make_pair(WeightName::xs00, *polarization),
					weight[*polarization]);
		}

#ifdef DEBUG_FLAG
		std::cout << "DEBUG: ModelHandler::eventEvolve() # dsig0 dsigcru: "
		<< dsigma0 << " " << dsigmaCrud << endl;
		std::cout << "DEBUG: ModelHandler::eventEvolve() # ModelBorn weight: "
		<< weight << std::endl;

#endif
		return weight[DataSource::get().mainWeight.second] != 0.0;

	} else {
		double deltaQED(0.0);
		double deltaEW(0.0);

		//---------------------------------------------------------------------
		//------------------- O(alpha) EW CORRECTIONS--------------------------
		//---------------------------------------------------------------------
		// Virtual+soft QED-like corrections for leptonic W decays
		if (DataSource::get().keyRad <= 2) {
			deltaQED
					= YFSUtil::deltaQED(
							event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame).m(),
							DataSource::get().masses[abs(
									event.getLepton().getPDGid())]);

			deltaEW = deltaQED;
			deltaEW += deltaWeak[abs(event.getLepton().getPDGid())];
		} else {
			throw VinhacException("keyRad > 2 not implemented yet");
		}

		// LL approximation
		double gamLL(2.0);
		gamLL *= alfpi;
		gamLL *= (2 * log(
				event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame).m()
						/ DataSource::get().masses[abs(
								event.getLepton().getPDGid())]) - 1.0);
		double deltaLL(0.25);
		deltaLL *= gamLL;

		//======================= BETA0 BETA0 BETA0 =============================
		// beta0:
		double beta00(dsigma0[PolarizationName::unpolarized]);
		double beta01_QED(beta00);
		beta01_QED *= (1 + deltaQED);
		double beta01_EW(beta00);
		beta01_EW *= (1 + deltaEW);
		double beta01_LL(beta00);
		beta01_LL *= (1 + deltaLL);

		//======================= BETA1 BETA1 BETA1 =============================
		// beta1:
		double beta11(0.0);
		double beta11_LL(0.0);

		double dsig1(0.0);
		double Stil(0.0);
		for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
			// Single hard-photon matrix element
			// 1. FSR only
			double
					xmat1(
							MatrixElement::matrixElementFSR1(
									quark,
									antiquark,
									boson,
									lepton,
									antilepton,
									event.getPhotons()[i].getFourMomentum(DataSource::get().modelReferenceFrame),
									chargeOfBoson, ckmElement));
			if (DataSource::get().keyRad >= 3) {
				// 2. FSR + iterference
				throw VinhacException("keyRad >= 3 not implemented yet");
			}

			//--- Normalized hard photon cross section
			dsig1 = normalizationFactor1;
			dsig1 /= sqq;
			dsig1 *= xmat1;

			// Soft factor S-tilde:

			if (DataSource::get().keyRad <= 2) {
				// 1. FSR only
				Stil = YFSUtil::StPair(
						sqq,
						pow(
								DataSource::get().masses[abs(
										event.getLepton().getPDGid())], 2),
						event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame),
						event.getLepton().getFourMomentum(DataSource::get().modelReferenceFrame),
						event.getPhotons()[i].getFourMomentum(DataSource::get().modelReferenceFrame));
			} else {
				// 2. FSR + iterference
				throw VinhacException("Keyrad >2 not implemented yet");
			}

			// beta1 O(alpha1)
			double beta1i(dsig1);
			beta1i /= Stil;
			beta1i -= beta00;
			beta11 += beta1i;
			// LL approximation
			double z(2.0);
			z *= event.getPhotons()[i].getFourMomentum(DataSource::get().modelReferenceFrame).m();
			z /= event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame).m();
			double RadLL = pow(z, 2) / 2.0 / (1.0 - z);
			double beta1i_LL(RadLL);
			beta1i_LL *= beta00;
			beta11_LL += beta1i_LL;

		} //  ! <-- photon loop


		// Check on betas - if not huge
		double Huge(1.0E3);
		double FacWcr(WtCrud);
		FacWcr /= dsigmaCrud;

		if (abs(beta00 * FacWcr) > Huge) {
			overWeightedEvents++;
			beta00 = 0.0;

			if(DataSource::get().printOutLevel == PrintOutLevel::Detailed){
				printHighWeightWaring(event, "beta00",
					beta00,dsigmaCrud, beta01_EW, beta01_QED,
					beta11, beta01_LL, WtCrud, FacWcr);
			}
		}
		double beta11tmp = beta11;

		if (abs((beta01_EW + beta11) * FacWcr) > Huge) {
			overWeightedEvents++;
			beta01_EW = 0.0;
			beta11 = 0.0;
			dsig1 = 0.0;

			if(DataSource::get().printOutLevel == PrintOutLevel::Detailed){
				printHighWeightWaring(event, "beta01_EW",
					beta00,dsigmaCrud, beta01_EW, beta01_QED,
					beta11, beta01_LL, WtCrud, FacWcr);
			}
		}

		if (abs((beta01_QED + beta11tmp) * FacWcr) > Huge) {
			overWeightedEvents++;

			beta01_QED = 0.0;
			beta11 = 0.0;
			dsig1 = 0.0;

			if(DataSource::get().printOutLevel == PrintOutLevel::Detailed){
				printHighWeightWaring(event, "beta01_QED",
					beta00,dsigmaCrud, beta01_EW, beta01_QED,
					beta11, beta01_LL, WtCrud, FacWcr);
			}
		}

		if (abs((beta01_LL + beta11_LL) * FacWcr) > Huge) {
			overWeightedEvents++;

			beta01_LL = 0.0;
			beta11_LL = 0.0;

			if(DataSource::get().printOutLevel == PrintOutLevel::Detailed){
				printHighWeightWaring(event, "beta01_LL",
					beta00,dsigmaCrud, beta01_EW, beta01_QED,
					beta11, beta01_LL, WtCrud, FacWcr);
			}
		}

		//=======================================================================
		//=============== FIXED ORDER, NO EXPONENTIATION ========================
		//=======================================================================
		double xs00(0.0);
		double xs0gQED(0.0);
		double xs0gEW(0.0);
		double xs1g(0.0);
		if (event.getPhotons().size() == 0) {
			double BtiReb(0.0);
			double delvsEW(0.0);
			double delvsQED(0.0);
			if (DataSource::get().keyRad <= 2) {
				//--- YFS FormFactor and QED/EW corrections
				BtiReb = YFSUtil::YWldec(
						event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame),
						event.getLepton().getFourMomentum(DataSource::get().modelReferenceFrame),
						DataSource::get().EgMin);
				delvsQED = YFSUtil::VirSofQED(
						event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame).m(),
						DataSource::get().masses[abs(
								event.getLepton().getPDGid())],
						DataSource::get().EgMin);
				delvsEW = delvsQED
						+ deltaWeak[abs(event.getLepton().getPDGid())];
			} else {
				throw VinhacException("keyrad > 2 not implemented yet");
			}
			double YFSfmf(exp(BtiReb));
			// Born
			xs00 = beta00;
			xs00 /= YFSfmf;
			// Born + virtual + soft
			xs0gQED = xs00;
			xs0gQED *= (1 + delvsQED);// ! pure QED
			xs0gEW = xs00;
			xs0gEW *= (1 + delvsEW);// ! EW corrections

		} else if (event.getPhotons().size() == 1) {
			// 1 real hard photon
			//--- YFS FormFactor
			double BtiReb(0.0);
			if (DataSource::get().keyRad <= 2) {
				BtiReb = YFSUtil::YWldec(
						event.getBoson().getFourMomentum(DataSource::get().modelReferenceFrame),
						event.getLepton().getFourMomentum(DataSource::get().modelReferenceFrame),
						DataSource::get().EgMin);
			} else {
				throw VinhacException("keyrad > 2 not implemented yet");
			}
			double YFSfmf(exp(BtiReb));
			xs1g = dsig1;
			xs1g /= Stil;
			xs1g /= YFSfmf;

		}

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++ WEIGHTS, WEIGHTS, WEIGHTS ++++++++++++++++++++++++
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// Principal weights:
		//----------------------------------------------------------------------
		//--- O(alpha^0)_exp
		event.addModelWeight(WeightName::beta00, beta00 / dsigmaCrud);
		//--- O(alpha^1)_exp: QED only
		event.addModelWeight(WeightName::beta01_QED_plus_beta11,
				(beta01_QED + beta11) / dsigmaCrud);
		//--- O(alpha^1)_exp: EW corrections
		event.addModelWeight(WeightName::beta01_EW_plus_beta11,
				(beta01_EW + beta11) / dsigmaCrud);
		weight[PolarizationName::unpolarized] = (beta01_EW + beta11) / dsigmaCrud;
		//--- O(alpha^1)_exp: LL QED only
		event.addModelWeight(WeightName::beta01_LL_plus_beta11_LL,
				(beta01_LL + beta11_LL) / dsigmaCrud);
		//----------------------------------------------------------------------
		// Individual betas...
		//--- O(alpha^0) beta00
		event.addModelWeight(WeightName::beta00, beta00 / dsigmaCrud);
		//--- O(alpha^1) beta01_QED
		event.addModelWeight(WeightName::beta01_QED, beta01_QED / dsigmaCrud);
		//--- O(alpha^1) beta11
		event.addModelWeight(WeightName::beta11, beta11 / dsigmaCrud);
		//--- O(alpha^1) beta01_EW
		event.addModelWeight(WeightName::beta01_EW, beta01_EW / dsigmaCrud);
		//--- O(alpha^1) beta01_LL
		event.addModelWeight(WeightName::beta01_LL, beta01_LL / dsigmaCrud);
		//--- O(alpha^1) beta11_LL
		event.addModelWeight(WeightName::beta11_LL, beta11_LL / dsigmaCrud);
		//========================================================================
		//================ FIXED ORDER, NO EXPONENTIATION ========================
		//========================================================================
		//--- Born
		event.addModelWeight(WeightName::xs00, xs00 / dsigmaCrud);
		//--- O(alpha^1)_QED
		event.addModelWeight(WeightName::xs0gQED_plus_xs1g,
				(xs0gQED + xs1g) / dsigmaCrud);
		//--- O(alpha^1)_EW
		event.addModelWeight(WeightName::xs0gEW_plus_xs1g,
				(xs0gEW + xs1g) / dsigmaCrud);
		//--- Born + virtual_QED + soft photon corr.
		event.addModelWeight(WeightName::xs0gQED, xs0gQED / dsigmaCrud);
		//--- Born + virtual_EW + soft photon corr.
		event.addModelWeight(WeightName::xs0gEW, xs0gEW / dsigmaCrud);
		//--- One real hard photon
		event.addModelWeight(WeightName::xs1g, xs1g / dsigmaCrud);

		//=======================================================================
		//***********************************************************************

		// Return best weight
		return weight[PolarizationName::unpolarized] != 0.0;
		//***********************************************************************
	}

	return weight[PolarizationName::unpolarized] != 0.0;
}

double VINHAC::ModelHandler::differentialXsectionCrude(const double& s,
		Event& event) {

	double invProp = CommonUtil::inverseWPropagator(s);
	double dsig;
	double QW = 1.0;
	double costheta = event.getCosTheta();
	if (event.getBoson().getPDGid() < 0) {
		QW = -1.0;
	}
	dsig = (ckmElement / 8.0 / DataSource::get().sinThetaW2
			/ DataSource::get().invAlphaQED) * (ckmElement / 8.0
			/ DataSource::get().sinThetaW2 / DataSource::get().invAlphaQED) * s
			/ invProp * (1.0 - QW * costheta) * (1.0 - QW * costheta);

	return dsig / 3.0;
}

void VINHAC::ModelHandler::printHighWeightWaring(Event& event, string name,
		double beta00, double DsiCru, double beta01_EW, double beta01_QED,
		double beta11, double beta01_LL, double WtCrud, double FacWcr){
	cout<<"*********************************************"<<endl
	<<" >>>>>> Huge "<<name<<" !!! - reset to 0."<<endl
	<<" No. of overweighted events = "<<overWeightedEvents<<endl
	<<" beta00,DsiCru = "<<beta00<<" "<<DsiCru<<endl
	<<" beta01_EW,beta11 = "<<beta01_EW<<" "<<beta11<<endl
	<<" beta01_QED,beta01_LL = "<<beta01_QED<<" "<<beta01_LL<<endl
	<<" WtCrud,beta00*FacWcr = "<<WtCrud<<" "<<beta00*FacWcr<<endl
	<<event.printInFrame(DataSource::get().modelReferenceFrame)<<endl
	<<"*********************************************"<<endl;
}
