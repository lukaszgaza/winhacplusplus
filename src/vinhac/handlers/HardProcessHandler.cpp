#include "HardProcessHandler.h"
#include "../input/DataSource.h"
#include "../core/Event.h"
#include "../core/VinhacException.h"
#include "../model/util/YFSUtil.h"
#include "../util/Transformation.h"
using namespace std;

VINHAC::HardProcessHandler::~HardProcessHandler() {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::~HardProcessHandler()"
	<< std::endl;

#endif

}

VINHAC::HardProcessHandler::HardProcessHandler()
	 {

#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::HardProcessHandler() 3"
	<< std::endl;

#endif

	KeyWid = 0;

	//set name
	setHandlerName("HardProcessHandler");

	//calculate gW coupling constant
	gW = 1.0 / (sqrt((1.0 - (DataSource::get().MW / DataSource::get().MZ) * (DataSource::get().MW / DataSource::get().MZ)) * (1 - (1.0
			- (DataSource::get().MW / DataSource::get().MZ) * (DataSource::get().MW / DataSource::get().MZ)))));

	Qe = sqrt(4 * DataSource::get().pi * DataSource::get().alphaQED);
	this->Ecm = DataSource::get().Ecm;

	epsk = 2.0 * DataSource::get().EgMin / DataSource::get().Ecm;

#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::HardProcessHandler() Ecm = "
	<< Ecm << std::endl;

#endif
}

bool VINHAC::HardProcessHandler::eventEvolve(VINHAC::Event& event) {

	if (event.currentCommonWtMultiplication() == 0.0) {
		event.addCommonWeight(WeightName::HardProcessWeight, 0.0);
		return false;
	}

	double weight = prepareParticles(event);

	weight *= generateLepton(event);

	generateLeptonAngles(event);

	weight *= generateLeptonMomenta(event);
	transformationToLAB(event);

	event.addCommonWeight(WeightName::HardProcessWeight, weight);
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::eventEvolve() : weight " << weight
	<< std::endl;

#endif
	return weight != 0.0;
}

void VINHAC::HardProcessHandler::generateLeptonAngles(Event& event) {

	phi = 2 * DataSource::get().pi * p_randomGenerator->Flat();
	double costheta = 0.0;

	costheta = 1.0 - 2.0 * pow(p_randomGenerator->Flat(), 1.0 / 3.0);


	if (event.getBoson().getPDGid() < 0) {
		costheta = -costheta;
	}

	event.setPhi(phi);
	event.setCosTheta(costheta);
	theta = acos(costheta);

#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::generateLeptonAngles() # phi: "
	<< phi * 180.0 / DataSource::get().pi << " # theta : " << theta * 180.0 / DataSource::get().pi
	<< std::endl;

#endif

}

double VINHAC::HardProcessHandler::generateLepton(Event& event) {
	int charge = 0;

	if (event.getBoson().getPDGid() == W) {

		charge = +1;

	} else if (event.getBoson().getPDGid() == -W) {

		charge = -1;

	} else {
		return 0.0;
	}

	int leptonId = 0;
	int neutrinoId = 0;
	if (DataSource::get().finalParticles.size() == 1) {
		leptonId = -charge * DataSource::get().finalParticles[0];
	} else {
		int
				idx =
						static_cast<int> (static_cast<double> (DataSource::get().finalParticles.size())
								* p_randomGenerator->Flat());
		leptonId = -charge * DataSource::get().finalParticles[idx];
	}

	if (leptonId > 0) {
		neutrinoId = -(leptonId + 1);
	} else {
		neutrinoId = -(leptonId - 1);
	}

	event.getLepton().setPDGid(leptonId);
	event.getNeutrino().setPDGid(neutrinoId);

	event.getLepton().addParent(&event.getBoson());
	event.getNeutrino().addParent(&event.getBoson());

#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::generateLepton() # leptonId="
	<< leptonId << " # neutrino Id=" << neutrinoId << std::endl;

#endif

	return 1.0;
}

double VINHAC::HardProcessHandler::prepareParticles(VINHAC::Event& event) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::prepareParticles()" << std::endl;

#endif

	event.getBoson().setFourMomentum(event.getQuark().getFourMomentum());
	event.getBoson().getFourMomentum()
			+= event.getAntiQuark().getFourMomentum();
	double sqq = event.getBoson().getFourMomentum().mag2();
	double amqq = sqrt(sqq);
	event.getBoson().getFourMomentumBosonRestFrame().set(0.0, 0.0, 0.0, amqq);

	event.getBoson().addParent(&event.getQuark());
	event.getBoson().addParent(&event.getAntiQuark());

	int quarkPDGid = event.getQuark().getPDGid();
	int antiQuarkPDGid = event.getAntiQuark().getPDGid();

	int charge = DataSource::get().quarksCharges[abs(quarkPDGid)] - DataSource::get().quarksCharges[abs(
			antiQuarkPDGid)];

	if (charge == 3) {
		event.getBoson().setPDGid(24);
	} else if (charge == -3) {
		event.getBoson().setPDGid(-24);
	} else {
		event.getBoson().setPDGid(24);
		return 0.0;
	}

	bool notFound = true;
	for (vector<int>::iterator it = DataSource::get().intermediateParticles.begin(); it
			!= DataSource::get().intermediateParticles.end(); ++it) {
		if (event.getBoson().getPDGid() == (*it)) {
			notFound = false;
			break;
		}
	}

	if (notFound)
		return 0.0;

#ifdef DEBUG_FLAG

	std::cout << "DEBUG: HardProcessHandler::prepareParticles() # boson id: "
	<< event.getBoson().getPDGid() << std::endl;
	std::cout << "DEBUG: HardProcessHandler::prepareParticles() # boson 4mom: "
	<< event.getBoson().getFourMomentum() << std::endl;
	std::cout
	<< "DEBUG: HardProcessHandler::prepareParticles() # boson 4mom_wrf: "
	<< event.getBoson().getFourMomentumBosonRestFrame() << std::endl;

#endif

	return 1.0;
}

double VINHAC::HardProcessHandler::generateLeptonMomenta(Event& event) {

	if (DataSource::get().keyRad == 0) {

		const CLHEP::HepLorentzVector& quarkMomenta =
				event.getQuark().getFourMomentum();
		const CLHEP::HepLorentzVector& antiquarkMomenta =
				event.getAntiQuark().getFourMomentum();

		double M2 = (quarkMomenta + antiquarkMomenta).mag2();
		double leptonMass = DataSource::get().masses[abs(event.getLepton().getPDGid())];
		double neutrinoMass = DataSource::get().masses[abs(event.getNeutrino().getPDGid())];

		if (M2 <= pow((leptonMass + neutrinoMass), 2)) {
			return 0.0;
		}

		double M = sqrt(M2);
		double leptonEnergy = (M2 + leptonMass * leptonMass - neutrinoMass
				* neutrinoMass) / 2 / M;
		double leptonMomentum = sqrt((M2 - (leptonMass + neutrinoMass)
				* (leptonMass + neutrinoMass)) * ((M2 - (leptonMass
				- neutrinoMass) * (leptonMass - neutrinoMass)))) / 2 / M;

		CLHEP::HepLorentzVector& leptonMomenta =
				event.getLepton().getFourMomentumBosonRestFrame();
		leptonMomenta.set(leptonMomentum * sin(theta) * cos(phi),
				leptonMomentum * sin(theta) * sin(phi), leptonMomentum * cos(
						theta), leptonEnergy);

		CLHEP::HepLorentzVector& neutrinoMomenta =
				event.getNeutrino().getFourMomentumBosonRestFrame();
		neutrinoMomenta.set(0, 0, 0, M - leptonEnergy);

		for (int i = 0; i < 3; ++i) {
			neutrinoMomenta[i] = -leptonMomenta[i];
		}

#ifdef DEBUG_FLAG

		std::cout
		<< "DEBUG: HardProcessHandler::generateLeptonMomenta() # ql_wrf :"
		<< leptonMomenta << std::endl;
		std::cout
		<< "DEBUG: HardProcessHandler::generateLeptonMomenta() # qn_wrf :"
		<< neutrinoMomenta << std::endl;

#endif

		return 1.0;
	} else {

		double weight = 0.0;

		event.getPhotons().clear();

		double Ebeam = 0.5 * DataSource::get().Ecm;
		double x1 = event.getQuark().getFourMomentum().e() / Ebeam;
		double x2 = event.getAntiQuark().getFourMomentum().e() / Ebeam;
		double xx = max(x1, x2);
		delta = epsk / xx;

		double radiationWeight = generateMultiphotonRadiation(event);

#ifdef DEBUG_FLAG
		std::cout
		<< "DEBUG: HardProcessHandler::generateLeptonMomenta() # radiationW :"
		<< radiationWeight << std::endl;
#endif
		if (radiationWeight == 0.0)
			return 0.0;

		double weightYff = YWlRem(event);

		//! Weight compensating for approximations in S-factor
		double WtSfa = 1.0;
		if (DataSource::get().keyRad > 2 && event.getPhotons().size() > 0) {

			throw VinhacException("keyrad > 2 in S-factor");
		}

		weight = radiationWeight * weightYff * WtSfa;

		// Check on Wt - if not huge
		if (abs(weight) > 1.0E6) {
			weight = 0.0;

			if(DataSource::get().printOutLevel == PrintOutLevel::Detailed){
				cout<<"*********************************************"<<endl
				<<" >>>>>> HardProcessHandler: Huge weight !!! - reset to 0."<<endl
				<<" WtYFS = "<<weight<<endl
				<<" WtRad,WtYff,WtSfa = "<<radiationWeight<<" "<<weightYff<<" "<<WtSfa<<endl
				<<"*********************************************"<<endl;
			}
		}

		return weight;
	}

	return 1.0;
}

double VINHAC::HardProcessHandler::YWlRem(Event& event) {
	//!     ***************************************************************
	//!----------------------------------------------------------------------!
	//! A Weight compensating for simplified YFS formfactor                  !
	//! used at the generation level for radiative W-decays.                 !
	//!                                                                      !
	//! INPUT: KeyRad -- radiative corrections switch                        !
	//!        amd,amu -- masses of initial quarks d (s, b) and u (c, t)     !
	//!        amW,aml  -- masses of W (invariant) and final-state lepton    !
	//!        qd,qu   -- masses of initial quarks d (s, b) and u (c, t)     !
	//!        Q,ql    -- 4-momenta of W and lepton, resp.                   !
	//!        Emin    -- soft photon energy cut-off [in GeV]                !
	//! OUPUT: WtYff   -- compensation weight                                !
	//!----------------------------------------------------------------------!
	//! Written by: Wieslaw Placzek                        CERN,  July 2002  !
	//! Last update: 30.07.2007           by: W.P.                           !
	//!----------------------------------------------------------------------!
	//IMPLICIT REAL*8 (A-H,O-Z)
	//PARAMETER( pi = 3.1415926535897932d0, alfinv=137.03599911d0)
	//PARAMETER( alfpi=  1/pi/alfinv, alfa= 1d0/alfinv)
	//DIMENSION qd(4),qu(4),Q(4),ql(4)
	//!
	//! ... YFS formfactor weight
	double aM = event.getBoson().getFourMomentumBosonRestFrame().mag();
	double Emin = delta * aM / 2;
	double epsg = 2 * Emin / aM;
	double Yscut = yfsCut(aM, DataSource::get().masses[abs(event.getLepton().getPDGid())],
			event.getBoson().getFourMomentumBosonRestFrame(),
			event.getLepton().getFourMomentumBosonRestFrame(), epsg);

	double Yfull = 0.0;
	if (DataSource::get().keyRad <= 2) {
		//! ... FSR only
		Yfull = VINHAC::YFSUtil::YWldec(
				event.getBoson().getFourMomentumBosonRestFrame(),
				event.getLepton().getFourMomentumBosonRestFrame(), Emin);
	} else {
		//! ... FSR + Interference
		//Yfull = YWDecInt(amd, amu, amW, aml, qd, qu, Q, ql, Emin);
		throw VinhacException("keyRad > 2 not implemented yet");
	}
	double Yrem = Yfull - Yscut;
	double WtY = exp(Yrem);
	//! Total weight
	return WtY;
}

double VINHAC::HardProcessHandler::generateMultiphotonRadiation(Event& event) {

	double svar = event.getBoson().getFourMomentumBosonRestFrame().mag2();
	double aM = sqrt(svar);
	double aml = DataSource::get().masses[abs(event.getLepton().getPDGid())];
	double aml2 = aml * aml;
	double Enec = (svar + aml2) / 2.0 / aM; // ! overvalued lepton energy for crude
	double amc2 = (aml / Enec) * (aml / Enec);
	double betc = (svar - aml2) / (svar + aml2);

	double alf1 = 1.0 / DataSource::get().pi * DataSource::get().alphaQED;
	double gamf2 = 2.0 * alf1 / betc * log(aM / aml);
	double Averg = gamf2 * log(1.0 / delta);
	vector<double> rr;

	int nPhotons = PoissGenerator(Averg, 100, rr);

	vector<double> dis0;
	dis0.reserve(nPhotons);
	vector<double> cgx;
	cgx.reserve(nPhotons);
	vector<double> sgx;
	sgx.reserve(nPhotons);
	vector<double> xk;
	xk.reserve(nPhotons);
	double sprim;
	sumOfPhotonsMomentum.setX(0);
	sumOfPhotonsMomentum.setY(0);
	sumOfPhotonsMomentum.setZ(0);
	sumOfPhotonsMomentum.setE(0);

	if (nPhotons == 0) {
		sprim = svar;
		//-----------------------------------------------------------------------
		//     Reject events outside phase space
		//-----------------------------------------------------------------------
		double smini = pow((DataSource::get().masses[abs(event.getLepton().getPDGid())]
				+ DataSource::get().masses[abs(event.getNeutrino().getPDGid())]), 2);
		if (sprim < smini)
			return 0;
	} else {
		// Photon energy
		double xsum = 0.0;

		for (int i = 0; i < nPhotons; ++i) {
			double tmp = pow(delta, rr[i]);
			xk.push_back(pow(delta, rr[i]));
			xsum = xsum + tmp;
		}
		if (xsum >= 2.0) {
			// Event outside phase space (too much radiation)
			return 0.0;
		}

		for (int i = 0; i < nPhotons; ++i) {
			event.getPhotons().resize(event.getPhotons().size() + 1);
			Particle& photon = event.getPhotons().back();
			photon.setPDGid(22);
			photon.setTag("Photon");

			double cg, sg, dis;
			generatePhotonAngularDistribution(aM, aml, amc2, betc, cg, sg, dis);
			dis0.push_back(dis);

			CLHEP::HepLorentzVector& photonMom =
					photon.getFourMomentumBosonRestFrame();
			//     define photon momenta (in units of aM/2)
			double phig = 2 * DataSource::get().pi * p_randomGenerator->Flat();
			photonMom.setX(xk[i] * sg * cos(phig));
			photonMom.setY(xk[i] * sg * sin(phig));
			photonMom.setZ(xk[i] * cg);
			photonMom.setE(xk[i]);

			sumOfPhotonsMomentum += photonMom;
			cgx.push_back(cg);
			sgx.push_back(sg);
		}
		//-----------------------------------------------------------------------
		//     Reject events with too much radiation (outside phase space)
		//-----------------------------------------------------------------------
		double xmk2 = sumOfPhotonsMomentum.mag2();
		double yy = 1.0 - sumOfPhotonsMomentum.e() + xmk2 / 4.0;
		sprim = svar * yy;
		double smini = pow((DataSource::get().masses[abs(event.getLepton().getPDGid())]
				+ DataSource::get().masses[abs(event.getNeutrino().getPDGid())]), 2);
		if (sprim < smini) {
			return 0;
		}
		//-----------------------------------------------------------------------
		//    Rescale properly all photon momenta
		//-----------------------------------------------------------------------
		double ener = aM / 2.0;

		sumOfPhotonsMomentum = sumOfPhotonsMomentum * ener;
		for (unsigned i = 0; i < event.getPhotons().size(); ++i) {

			event.getPhotons()[i].getFourMomentumBosonRestFrame() *= ener;
			//Insufficient
			//event.getPhotons()[i].setFourMomentumBosonRestFrame(
			//		event.getPhotons()[i].getFourMomentumBosonRestFrame()
			//				* ener);
		}
	}

	double WtKin = setLeptonMomentaInRadiativeDecay(event, sprim, aM,
			sumOfPhotonsMomentum);

	if (WtKin == 0.0)
		return 0.0;
	// Velocity of charged lepton
	double bet1 = event.getLepton().getFourMomentumBosonRestFrame().z()
			/ event.getLepton().getFourMomentumBosonRestFrame().e();

	// Weight correponding to phase-space integration over energy-momentum
	// conservation delta-function
	double
			WtPS =
					2 * event.getLepton().getFourMomentumBosonRestFrame().z()
							/ (aM - sumOfPhotonsMomentum.e()
									+ event.getLepton().getFourMomentumBosonRestFrame().e()
											/ event.getLepton().getFourMomentumBosonRestFrame().z()
											* sumOfPhotonsMomentum.z());

	//-----------------------------------------------------------------------
	//     Mass weight for theta distribution
	//-----------------------------------------------------------------------
	//     Mass weight compensates for s'->s and droping terms -m**2/(k.q)**2
	//     Care is taken of machine rounding errors.
	//     delt RECALCULATED out of angles sgx(i),cgx(i)
	//     with TRUE lepton energy, using EXACT formulas
	double amd1 = pow((DataSource::get().masses[abs(event.getLepton().getPDGid())]
			/ event.getLepton().getFourMomentumBosonRestFrame().e()), 2);

	double Wt2 = 1.0;
	ygr.clear();
	zet.clear();
	for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
		double delt = 0.0;
		if (cgx[i] > 0.0) {
			delt = amd1 / (1 + bet1) + bet1 * pow(sgx[i], 2) / (1 + cgx[i]);
		} else {
			delt = 1.0 - bet1 * cgx[i];
		}
		double dist1 = 2 / delt * (1.0 - amd1 / 2 / delt - delt / 2);
		double wtm = dist1 / dis0[i];
		if (wtm < 1.0E-90)
			wtm = 0.0;
		Wt2 = Wt2 * wtm;
		//     Finaly define Sudakov variables
		ygr.push_back(xk[i] * delt
				* event.getLepton().getFourMomentumBosonRestFrame().e() / aM);
		zet.push_back(xk[i]);

	}

	// The total phase space integral for crude x-section,
	// Note that delta is a lower limit on x-variables in crude generation
	double Fphs = 2 * alf1 / betc * log(aM / aml) * log(1.0 / delta);
	// YFS formfactor cut-off dependend part in W rest frame
	// Note: The same cut-off delta as for "crude formfactor" is used here!
	double Fyfs = yfsCut(aM, aml,
			event.getBoson().getFourMomentumBosonRestFrame(),
			event.getLepton().getFourMomentumBosonRestFrame(), delta);
	double delb = Fyfs + Fphs;
	Wt2 = Wt2 * exp(delb);

	//-----------------------------------------------------------------------
	// Rotations in WRF to get z-axis along the beams, with +z pointing in
	// the initial quark direction.
	//-----------------------------------------------------------------------
	Transformation::Rot2le(theta, phi, event.getLepton().getFourMomentumBosonRestFrame());
	Transformation::Rot2le(theta, phi, event.getNeutrino().getFourMomentumBosonRestFrame());
	Transformation::Rot2le(theta, phi, sumOfPhotonsMomentum);
	for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
		Transformation::Rot2le(theta, phi,
				event.getPhotons()[i].getFourMomentumBosonRestFrame());
	}

#ifdef DEBUG_FLAG

	std::cout
	<< "DEBUG: HardProcessHandler::generateMultiphotonRadiation() # nPhotons :"
	<< nPhotons << std::endl;

#endif

	return Wt2 * WtPS;
}

double VINHAC::HardProcessHandler::yfsCut(const double& am1, const double& am2,
		const CLHEP::HepLorentzVector& p1, const CLHEP::HepLorentzVector& p2,
		const double& epsg) {
	//*     ***********************************
	//!----------------------------------------------------------------------!
	//! This function provides a value of the soft photon cut-off dependent  !
	//! part of the YFS infrared function Y = Btilde + ReB for any pair of   !
	//! charged particles.                                                   !
	//! INPUT: am1,am2 - particles masses                                    !
	//!        q1,q2   - particles 4-momenta                                 !
	//!        epsg    - dimensionless soft photon cut-off (defined in       !
	//!                  a given frame)                                      !
	//!----------------------------------------------------------------------!
	//! Written by:  Wieslaw Placzek                 Knoxville, Aug. 1997    !
	//! Last update:  4.08.1997                by: W.P.                      !
	//!----------------------------------------------------------------------!

	double alfpi = DataSource::get().alphaQED / DataSource::get().pi;

	double am12 = am1 * am2;
	double p1p2 = p1 * p2;
	double xlam = sqrt((p1p2 - am12) * (p1p2 + am12));
	double A = log((p1p2 + xlam) / am12) / xlam;
	return -2.0 * alfpi * (p1p2 * A - 1.0) * log(1.0 / epsg);
}

double VINHAC::HardProcessHandler::setLeptonMomentaInRadiativeDecay(
		Event& event, const double& sprim, const double& aM,
		const CLHEP::HepLorentzVector& phSum) {
	double WtKin = 0.0;
	double aml1s = pow(DataSource::get().masses[abs(event.getLepton().getPDGid())], 2);
	double aml2s = pow(DataSource::get().masses[abs(event.getNeutrino().getPDGid())], 2);
	double E12 = aM - sumOfPhotonsMomentum.e();
	double phs3 = sumOfPhotonsMomentum.z();
	double eph3 = (E12 - phs3) * (E12 + phs3);
	double rsp = sprim + aml1s - aml2s;
	double dlt = pow(rsp, 2) - 4 * aml1s * eph3;
	// Outside phase space
	if (dlt <= 0.0)
		return 0;
	double q1m = (E12 * sqrt(dlt) - rsp * phs3) / 2.0 / eph3;
	// Outside phase space
	if (q1m <= 0.0)
		return 0;
	double E1 = sqrt(pow(q1m, 2) + aml1s);

	// Charged lepton 4-momentum
	CLHEP::HepLorentzVector& q1 =
			event.getLepton().getFourMomentumBosonRestFrame();
	q1.setZ(q1m);
	q1.setE(E1);
	// Neutrino 4-momentum
	CLHEP::HepLorentzVector& q2 =
			event.getNeutrino().getFourMomentumBosonRestFrame();
	q2.setX(-sumOfPhotonsMomentum.x());
	q2.setY(-sumOfPhotonsMomentum.y());
	q2.setZ(-sumOfPhotonsMomentum.z() - q1m);
	q2.setE(E12 - E1);
	// Kinematics weight
	WtKin = 1.0;
	return WtKin;

}

void VINHAC::HardProcessHandler::generatePhotonAngularDistribution(
		const double& aM, const double& aml, const double& amx2,
		const double& beta, double& costhg, double& sinthg, double& dist) {

	// Generation of a random number
	double rn = p_randomGenerator->Flat();
	// Generation of cos(theta_gamma)
	double R2 = 2.0 * rn;
	double del = (1.0 + beta) / pow((aM / aml), R2);//             ! 1 - beta*costhg
	costhg = (1.0 - del) / beta;
	sinthg = sqrt(del * (2.0 - del) - amx2 * pow(costhg, 2));
	dist = 2.0 / del;
}

void VINHAC::HardProcessHandler::transformationToLAB(Event& event) {

	double epsb = 1.0E-16;

	const CLHEP::HepLorentzVector& bosonMom = event.getBoson().getFourMomentum();
	double beta2 = (bosonMom[0] * bosonMom[0] + bosonMom[1] * bosonMom[1]
			+ bosonMom[2] * bosonMom[2]) / (bosonMom[3] * bosonMom[3]);

	if (beta2 > epsb) {
		Transformation::BOSTDQ(-1, event.getBoson().getFourMomentum(),
				event.getLepton().getFourMomentumBosonRestFrame(),
				event.getLepton().getFourMomentum());
		Transformation::BOSTDQ(-1, event.getBoson().getFourMomentum(),
				event.getNeutrino().getFourMomentumBosonRestFrame(),
				event.getNeutrino().getFourMomentum());

		if (event.getPhotons().size() > 0) {

			for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
				Transformation::BOSTDQ(-1, event.getBoson().getFourMomentum(),
						event.getPhotons()[i].getFourMomentumBosonRestFrame(),
						event.getPhotons()[i].getFourMomentum());
			}
		}

	} else {
		event.getLepton().setFourMomentum(
				event.getLepton().getFourMomentumBosonRestFrame());
		event.getNeutrino().setFourMomentum(
				event.getNeutrino().getFourMomentumBosonRestFrame());
		if (event.getPhotons().size() > 0) {

			for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
				event.getPhotons()[i].setFourMomentum(
						event.getPhotons()[i].getFourMomentumBosonRestFrame());
			}
		}
	}

	if (event.getQuark().getParents().size()>0){
		if (event.getQuark().getParents()[0]->getTag() == "BeamB") {

	#ifdef DEBUG_FLAG
			std::cout
			<< "DEBUG: HardProcessHandler::transformationToLAB() # rotation ! :"
			<< std::endl;
	#endif
			{
				CLHEP::HepLorentzVector& tmp = event.getQuark().getFourMomentum();
				tmp[1] = -tmp[1];
				tmp[2] = -tmp[2];
			}

			{
				CLHEP::HepLorentzVector& tmp =
						event.getAntiQuark().getFourMomentum();
				tmp[1] = -tmp[1];
				tmp[2] = -tmp[2];
			}

			{
				CLHEP::HepLorentzVector& tmp =
						event.getNeutrino().getFourMomentum();
				tmp[1] = -tmp[1];
				tmp[2] = -tmp[2];
			}

			{
				CLHEP::HepLorentzVector& tmp = event.getLepton().getFourMomentum();
				tmp[1] = -tmp[1];
				tmp[2] = -tmp[2];
			}

			{
				CLHEP::HepLorentzVector&tmp = event.getBoson().getFourMomentum();
				tmp[1] = -tmp[1];
				tmp[2] = -tmp[2];
			}

			if (event.getPhotons().size() > 0) {

				for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
					{
						CLHEP::HepLorentzVector&tmp =
								event.getPhotons()[i].getFourMomentum();
						tmp[1] = -tmp[1];
						tmp[2] = -tmp[2];
					}
				}

			}
		}

		CLHEP::HepLorentzVector PP(
				event.getQuark().getParents()[0]->getFourMomentum());
		PP += event.getAntiQuark().getParents()[0]->getFourMomentum();

		beta2 = (PP[0] * PP[0] + PP[1] * PP[1] + PP[2] * PP[2]) / (PP[3] * PP[3]);

		if (beta2 > epsb) {
			Transformation::BOSTDQ(-1, PP, event.getQuark().getFourMomentum(),
					event.getQuark().getFourMomentum());
			Transformation::BOSTDQ(-1, PP, event.getAntiQuark().getFourMomentum(),
					event.getAntiQuark().getFourMomentum());
			Transformation::BOSTDQ(-1, PP, event.getLepton().getFourMomentum(),
					event.getLepton().getFourMomentum());
			Transformation::BOSTDQ(-1, PP, event.getNeutrino().getFourMomentum(),
					event.getNeutrino().getFourMomentum());
			Transformation::BOSTDQ(-1, PP, event.getBoson().getFourMomentum(),
					event.getBoson().getFourMomentum());
			if (event.getPhotons().size() > 0) {

				for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
					Transformation::BOSTDQ(-1, PP, event.getPhotons()[i].getFourMomentum(),
							event.getPhotons()[i].getFourMomentum());
				}
			}
		}

		if (event.getQuark().getParents()[0]->getTag() == "BeamA") {
			event.getBeamRestA().setFourMomentum(
					event.getBeamParticleA().getFourMomentum());
			event.getBeamRestA().getFourMomentum()
					-= event.getQuark().getFourMomentum();
			event.getBeamRestB().setFourMomentum(
					event.getBeamParticleB().getFourMomentum());
			event.getBeamRestB().getFourMomentum()
					-= event.getAntiQuark().getFourMomentum();
		} else if (event.getQuark().getParents()[0]->getTag() == "BeamB") {
			event.getBeamRestA().setFourMomentum(
					event.getBeamParticleA().getFourMomentum());
			event.getBeamRestA().getFourMomentum()
					-= event.getAntiQuark().getFourMomentum();
			event.getBeamRestB().setFourMomentum(
					event.getBeamParticleB().getFourMomentum());
			event.getBeamRestB().getFourMomentum()
					-= event.getQuark().getFourMomentum();
		}
	}
#ifdef DEBUG_FLAG
	std::cout << "DEBUG: HardProcessHandler::transformationToLAB() # beta2 :"
	<< beta2 << std::endl;
	std::cout << "DEBUG: HardProcessHandler::transformationToLAB() # ql :"
	<< event.getLepton().getFourMomentum() << " " << event.getLepton().getFourMomentum().m()
	<< std::endl;
	std::cout << "DEBUG: HardProcessHandler::transformationToLAB() # qn :"
	<< event.getNeutrino().getFourMomentum() << " "
	<< event.getNeutrino().getFourMomentum().m() << std::endl;
	std::cout << "DEBUG: HardProcessHandler::transformationToLAB() # pq :"
	<< event.getQuark().getParents()[0]->getFourMomentum() << " "
	<< event.getQuark().getParents()[0]->getFourMomentum().m() << std::endl;
	std::cout << "DEBUG: HardProcessHandler::transformationToLAB() # paq :"
	<< event.getAntiQuark().getParents()[0]->getFourMomentum() << " "
	<< event.getAntiQuark().getParents()[0]->getFourMomentum().m() << std::endl;
	std::cout << "DEBUG: HardProcessHandler::transformationToLAB() # q :"
	<< event.getQuark().getFourMomentum() << " " << event.getQuark().getFourMomentum().m()
	<< std::endl;
	std::cout << "DEBUG: HardProcessHandler::transformationToLAB() # aq :"
	<< event.getAntiQuark().getFourMomentum() << " "
	<< event.getAntiQuark().getFourMomentum().m() << std::endl;
	std::cout << "DEBUG: HardProcessHandler::transformationToLAB() # W :"
	<< event.getBoson().getFourMomentum() << " " << event.getBoson().getFourMomentum().m()
	<< std::endl;
	for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
		std::cout
		<< "DEBUG: HardProcessHandler::transformationToLAB() # photon :"
		<< event.getPhotons()[i].getFourMomentum() << " "
		<< event.getPhotons()[i].getFourMomentum().m() << std::endl;
	}
	for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
		std::cout
		<< "DEBUG: HardProcessHandler::transformationToLAB() # photon wrf:"
		<< event.getPhotons()[i].getFourMomentumBosonRestFrame() << " "
		<< event.getPhotons()[i].getFourMomentumBosonRestFrame().m() << std::endl;
	}
#endif

}

int VINHAC::HardProcessHandler::PoissGenerator(double average, int nmax,
		vector<double> &rr) {
#ifdef DEBUG_FLAG

	std::cout << "DEBUG: YFShandler::PoissGenerator()" << std::endl;

#endif
	const double Nlim = 100;
	double sum, y;
	int nn, Nfail = 0;
	int mult = 0;
	//-------------------------------------------------------------------------------------------
	// check if input parameters are correct
	if (average <= 0) {
		throw VinhacException("STOP in HardProcessHandler, PoissGenerator");
	}
	//make a loop over Nlim tries of generation
	for (Nfail = 0; Nfail < Nlim; Nfail++) {
		nn = 0;
		sum = 0.0;
		rr.clear();

		for (int i = 1; i <= nmax; i++) {
			double rand = p_randomGenerator->Flat();
			y = log(rand);
			sum = sum + y;
			nn++;
			rr.push_back(sum / (-average));
			//	cout << rr[i-1] << endl;
			if (sum < -average) {
				mult = nn - 1;
				return mult;
			}
		}
	}
	throw VinhacException("HardProcessHandler, PoissGenerator: TOO SMALL NMax");
}

