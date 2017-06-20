using namespace std;

#include "BeamHandler.h"
#include "../pdf/PDF_NoPDF.h"
#include "../core/VinhacException.h"
#include "../../utils/FortranFunction.h"
#include <algorithm>
#include "../../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include "../input/DataSource.h"
#include "ppWintegrand.h"
#include "../pdf/PDF_LHA.h"
#include "LHAPDF/LHAPDF.h"
#include "../pdf/PDF_Simple.h"
#include "../core/BeamParticle.h"
#include "../../utils/FOAM/TFOAM.h"
#include "ModelHandler.h"
#include "../core/Event.h"
#include "../../utils/FOAM/TRND.h"
#include "../model/util/CommonUtil.h"

//! constructor

VINHAC::BeamHandler::BeamHandler(TRND *p_randomGenerator) :
	pdfSubsetA(0), pdfSubsetB(0), pdflibNGroupA(-1), pdflibNGroupB(-1),
			pdflibNSetA(-1), pdflibNSetB(-1), pdfIdA(0), pdfIdB(0),
			p_randomGenerator(p_randomGenerator), wtFoam(0.0), XsCrud(0),
			XsCrudError(0), XsCrudWt(1.0), integrand(0), pdfObject(0), scalePDF(-1) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::BeamHandler() 1"<<std::endl;

#endif

	for (vector<int>::iterator it = DataSource::get().initialParticles.begin(); it
			!= DataSource::get().initialParticles.end(); ++it) {
		if ((*it) > 0) {
			if ((*it) % 2 == 1) {
				dQuarks.push_back((*it));
			} else {
				uQuarks.push_back((*it));
			}

		} else {
			if ((-(*it)) % 2 == 1) {
				dAntiQuarks.push_back((*it));
			} else {
				uAntiQuarks.push_back((*it));
			}

		}
	}

	sort(dQuarks.begin(), dQuarks.end(), comparator);
	sort(dAntiQuarks.begin(), dAntiQuarks.end(), comparator);
	sort(uQuarks.begin(), uQuarks.end(), comparator);
	sort(uAntiQuarks.begin(), uAntiQuarks.end(), comparator);

	isWplus = false;
	isWminus = false;
	for (vector<int>::iterator it = DataSource::get().intermediateParticles.begin(); it
			!= DataSource::get().intermediateParticles.end(); ++it) {
		if ((*it) == 24)
			isWplus = true;
		if ((*it) == -24)
			isWminus = true;

	}

	p_FoamX = 0;

	integrand = new ppWintegrand(this);

	//set handler name
	setHandlerName("BeamHandler");

	/////////////////////////BEAMS INITIALIZATION////////////////////////////////////
	/////////// PDFS
	if(DataSource::get().isPDFOn){
		if (DataSource::get().pdfInterface == PDFInterface::LHAPDF) {
			if(DataSource::get().printOutLevel>0){
				cout
				<< "========================LHAPDF=============================="
				<< endl;
			}
			if(DataSource::get().printOutLevel==PrintOutLevel::None){
				LHAPDF::setVerbosity(LHAPDF::SILENT);
			} else if(DataSource::get().printOutLevel==PrintOutLevel::Basic){
				LHAPDF::setVerbosity(LHAPDF::LOWKEY);
			} else if(DataSource::get().printOutLevel==PrintOutLevel::Detailed){
				LHAPDF::setVerbosity(LHAPDF::DEFAULT);
			}
			string pdfName = DataSource::get().pdfName;
			pdfNameA = pdfName;
			pdfNameB = pdfName;
			VINHAC::PDF_LHA* p_pdf = new VINHAC::PDF_LHA();
			p_pdf->alphasPDF(DataSource::get().MZ);
			LHAPDF::PDFSetInfo pdfinfo = p_pdf->initPDFByName(pdfName,
					LHAPDF::LHGRID, DataSource::get().pdfSubset);
			pdflibNGroupA = pdflibNGroupB = pdfinfo.pdflibNGroup;
			pdflibNSetA = pdflibNSetB = pdfinfo.pdflibNSet;
			pdfIdA = pdfIdB = pdfinfo.id;
			pdfSubsetA = pdfSubsetB = DataSource::get().pdfSubset;
			pdfObject = p_pdf;
		} else if(DataSource::get().pdfInterface ==  PDFInterface::Simple) {
			VINHAC::PDF_Simple * p_pdf = new PDF_Simple();
			p_pdf->setXmin(1e-6);
			p_pdf->setQ2min(1.69);
			p_pdf->setXmax(1.0);
			pdfObject = p_pdf;
		} else {
			throw VinhacException("Unsupported PDFInterface");
		}
	} else {
		VINHAC::PDF_NoPDF * p_pdf = new PDF_NoPDF();
		p_pdf->setXmin(1e-6);
		p_pdf->setQ2min(1.69);
		p_pdf->setXmax(1.0);
		pdfObject = p_pdf;
	}

	/////////// initialization of beamA
	// PDGid and Energy of beam particle
	string PDGidA = DataSource::get().beamA.find("PDGid")->second;
	string energyA = DataSource::get().beamA.find("energy")->second;

	double massA = DataSource::get().masses[abs(atoi(PDGidA.c_str()))];
	double energyAD = atof(energyA.c_str());
	// 4-momentum of the particle
	CLHEP::HepLorentzVector momentumA(0.0, 0.0, sqrt(energyAD * energyAD
			- massA * massA), energyAD);

	p_beamParticleA = new VINHAC::BeamParticle(atoi(PDGidA.c_str()),
					momentumA, pdfObject);

	if (!DataSource::get().isPDFOn && fabs(atof(PDGidA.c_str())) >= 7) {
			cout << "=== Bad beamA particle PDGid when no PDF used ====="
					<< endl;
			throw VinhacException();
	}

	/////////// initialization of beamB
	// PDGid and Energy of beam particle
	string PDGidB = DataSource::get().beamB.find("PDGid")->second;
	string energyB = DataSource::get().beamB.find("energy")->second;

	double massB = DataSource::get().masses[abs(atoi(PDGidB.c_str()))];
	double energyBD = atof(energyB.c_str());
	// 4-momentum of the particle
	CLHEP::HepLorentzVector momentumB(0.0, 0.0, -sqrt(energyBD * energyBD
			- massB * massB), energyBD);


	p_beamParticleB = new VINHAC::BeamParticle(atoi(PDGidB.c_str()),
			momentumB, pdfObject);

	if (!DataSource::get().isPDFOn && fabs(atof(PDGidB.c_str())) >= 7) {
		cout << "=== Bad beamB particle PDGid when no PDF used ====="
				<< endl;
		throw VinhacException();
	}

	//// Set Ecm and s
	Ecm = p_beamParticleA->getFourMomentum().e()
			+ p_beamParticleB->getFourMomentum().e();

	s = (p_beamParticleA->getFourMomentum()
			+ p_beamParticleB->getFourMomentum()).mag2();

	////Beam particle
	p_beamParticleA->setTag("BeamA"); //!< Beam particle A
	p_beamParticleB->setTag("BeamB"); //!< Beam particle B

	if (DataSource::get().isPDFOn) {
		initializeIntegrand();
		initializeFOAM();
	}

	DataSource::get().Ecm = this->Ecm;

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::BeamHandler() Ecm = "<<Ecm<<std::endl;

#endif
	setNormalizationFactor();
	setXsCrud();
}
bool VINHAC::BeamHandler::initializeFOAM() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::initializeFOAM()"<<std::endl;

#endif

	//=========================================================
	int nDim = 0; // SIMPLICAL   subspace
	int kDim = 2; // HYP-CUBICAL subspace
	int nCells = 5000;//5000; // Number of Cells
	int nSampl = 1000;//1000; // Number of MC events per cell in build-up
	int OptDrive = 2; // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	int OptEdge = 0; // (D=0) Option, vertices are included in the sampling (1) or not (0)
	//double MaxWtRej = 2.0; // Maximum wt for rejection, for OptRej=1
	//=========================================================
	p_FoamX = new TFOAM("p_FoamX"); // Create Simulator
	//=========================================================
	p_FoamX->SetChat(DataSource::get().printOutLevel);
	p_FoamX->SetnDim(nDim);
	p_FoamX->SetkDim(kDim);
	p_FoamX->SetnSampl(nSampl);
	p_FoamX->SetnCells(nCells);
	p_FoamX->SetOptDrive(OptDrive);
	p_FoamX->SetOptEdge(OptEdge);
	//p_FoamX->SetMaxWtRej(MaxWtRej);
	//p_FoamX->SetRho(Density);
	//===============================
	//p_randomGenerator = new CLHEP::RanluxEngine();

	p_FoamX->Initialize(p_randomGenerator, integrand);

	return true;
}

bool VINHAC::BeamHandler::initializeIntegrand() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::initializeIntegrand()"<<std::endl;

#endif

	integrand->setBeamParticleA(p_beamParticleA);
	integrand->setBeamParticleB(p_beamParticleB);
	integrand->setS(s);

	return true;
}

void VINHAC::BeamHandler::setNormalizationFactor() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::setNormalizationFactor()"<<std::endl;

#endif

	///// Overall normalization factor depending on IPS
	// Coupling 1/alfinv is used internally in matrix element calculations
	if (DataSource::get().keyGmu == 0) {
		//--- alpha-scheme
		//	cout << "GeV2nb = " << GeV2nb << endl;
		normalizationFactor = DataSource::get().invGeV2toNb;
		//cout << "GeV2nb" << DataSource::get().invGeV2toNb << endl;
	} else if (DataSource::get().keyGmu == 1) {
		//--- G_mu-scheme: alpha_Gmu - effective coupling constant at M_W scale
		double alpha_Gmu = sqrt(2.0) * DataSource::get().fermiConst * DataSource::get().MW * DataSource::get().MW
				* DataSource::get().sinThetaW2 / DataSource::get().pi;
		normalizationFactor = DataSource::get().invGeV2toNb * (alpha_Gmu * DataSource::get().invAlphaQED)
				* (alpha_Gmu * DataSource::get().invAlphaQED);
		//cout << "GeV2nb 2" << DataSource::get().invGeV2toNb << endl;
	} else {
		throw VinhacException("Wrong choise of alpha scheme!!!");
	}
	// Multiply FacNor by nuber of families in the final state
	normalizationFactor = normalizationFactor
			* static_cast<double> (DataSource::get().finalParticles.size());

}

void VINHAC::BeamHandler::setXsCrud() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::setXsCrud()"<<std::endl;

#endif

	if (!DataSource::get().isPDFOn) {

		double ckmElement = 0.0;
		int iu;
		int id;

		BeamParticle *quark;
		BeamParticle *antiQuark;
		if (p_beamParticleA->getPDGid() > 0 && p_beamParticleB->getPDGid() < 0) {
			quark = p_beamParticleA;
			antiQuark = p_beamParticleB;
		} else if (p_beamParticleA->getPDGid() < 0
				&& p_beamParticleB->getPDGid() > 0) {
			quark = p_beamParticleB;
			antiQuark = p_beamParticleA;
		} else {
			throw VinhacException("Bad choise of quarks when noPDF");
		}
		if (DataSource::get().quarksCharges[abs(antiQuark->getPDGid())] == 2) {
			iu = abs(antiQuark->getPDGid()) / 2 - 1;
			id = quark->getPDGid() / 2 + 1 - 1;
		} else {
			id = abs(antiQuark->getPDGid()) / 2 + 1 - 1;
			iu = quark->getPDGid() / 2 - 1;
		}
		ckmElement = DataSource::get().ckm[iu][id];

		double sqq = Ecm * Ecm;

		double sigCru = totalXsCrude(ckmElement, sqq);
		XsCrud = sigCru * normalizationFactor;

	} else {
		double SigCru = p_FoamX->GetPrimary();
#ifdef DEBUG_FLAG

		std::cout<<"DEBUG: BeamHandler::setXsCrud() # z foam : "<<SigCru<<std::endl;

#endif
		XsCrud = normalizationFactor * SigCru;

		DataSource::get().maxWeightRej = DataSource::get().maxWeightRej*XsCrud;
	}

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::setXsCrud() # XsCrud : "<<XsCrud<<std::endl;

#endif

	XsCrudWt = XsCrud;
	XsCrud = 1.0;
}

double VINHAC::BeamHandler::totalXsCrude(const double& Uij, const double& s) {
	double result = DataSource::get().pi / 12.0 * (DataSource::get().alphaQED * Uij / DataSource::get().sinThetaW2)
			* (DataSource::get().alphaQED * Uij / DataSource::get().sinThetaW2) * s
			/ CommonUtil::inverseWPropagator(s);

	return result / 3.0;
}

bool VINHAC::BeamHandler::eventEvolve(VINHAC::Event& event) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::eventEvolve()"<<std::endl;

#endif

	// Clean Beam Handler before generation of next event
	cleanBeamHandlerData();

	event.getBeamParticleA().setFourMomentum(p_beamParticleA->getFourMomentum());
	event.getBeamParticleB().setFourMomentum(p_beamParticleB->getFourMomentum());
	event.getBeamParticleA().setPDGid(p_beamParticleA->getPDGid());
	event.getBeamParticleB().setPDGid(p_beamParticleB->getPDGid());

	/// generate quarks particle
	double weight = generateQuarks(event);

	/// set parton kinematics
	weight *= setPartonsKinematics(event);

	//weight *= XsCrudWt;
	event.addCommonWeight(WeightName::XsCrudWt, XsCrudWt);
	event.addCommonWeight(WeightName::PartonsGeneration, weight);

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::eventEvolve() # weight: "<<weight<<std::endl;

#endif
	return weight != 0.0;
}

double VINHAC::BeamHandler::generateQuarks(Event& event) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::generateQuarks()"<<std::endl;

#endif
	if (!DataSource::get().isPDFOn) {
		//--- Parton level:
		int Idq1;
		int Idq2;
		Idq1 = p_beamParticleA->getPDGid();
		Idq2 = p_beamParticleB->getPDGid();
		if (Idq1 > 0 && Idq2 < 0) {
			event.getQuark().setPDGid(Idq1);
			event.getAntiQuark().setPDGid(Idq2);
		} else if (Idq1 < 0 && Idq2 > 0) {
			event.getAntiQuark().setPDGid(Idq1);
			event.getQuark().setPDGid(Idq2);
		} else {
			throw VinhacException("Bad partons in beams");
		}

		x1_parton = 1.0;
		x2_parton = 1.0;

		return 1.0;
	} else {
		event.pdfUsed = true;
		p_FoamX->MakeEvent();
		double mcVector[2];
		p_FoamX->GetMCvect(mcVector);
		double weightFoam = p_FoamX->GetMCwt();

		double eta = mcVector[0];
		double zet = mcVector[1];

		double x = integrand->getX_min() * exp(integrand->getXminimal() * zet);
		double alf_min = atan((integrand->getX_min() * x - integrand->getA())
				/ integrand->getB());
		double alf_max = atan((integrand->getX_max() * x - integrand->getA())
				/ integrand->getB());
		double delalf = alf_max - alf_min;
		double alf = alf_min + delalf * eta;
		double tau = integrand->getA() + integrand->getB() * tan(alf);

		x1_parton = min(1.0, x);
		x2_parton = min(1.0, tau / x1_parton);

		generateFlavours(event);
		event.pdfScale = scalePDF;
		event.quarkXpdf = xPDFa[event.getQuark().getPDGid() + 6];
		event.antiquarkXpdf = xPDFb[event.getAntiQuark().getPDGid() + 6];
		event.quarkPDFid = pdfIdA;
		event.antiquarkPDFid = pdfIdB;

#ifdef DEBUG_FLAG

		std::cout<<"DEBUG: BeamHandler::generateFlavours() # q aq : "<<event.getQuark().getPDGid()<<" "<<event.getAntiQuark().getPDGid()<<std::endl;

#endif

		return weightFoam;
	}

	return 0.0;
}

double VINHAC::BeamHandler::setPartonsKinematics(Event& event) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::setPartonsKinematics()"<<std::endl;

#endif

	////////////////////////////////////////////////////////////////////////
	// Fill quarks and Z-boson data records:                              //
	// At the parton level: +z axis pointing in the quark direction;      //
	// it is adjusted to the LAB frame.								      //
	////////////////////////////////////////////////////////////////////////


	double xqu = x1_parton;
	double xaq = x2_parton;
	event.quarkX = x1_parton;
	event.antiquarkX = x2_parton;

	event.getQuark().getFourMomentum().set(0);
	event.getAntiQuark().getFourMomentum().set(0);
	event.getQuark().getFourMomentumBosonRestFrame().set(0);
	event.getAntiQuark().getFourMomentumBosonRestFrame().set(0);

	int Idqu;
	int Idaq;
	// Set quark ID's, x's

	Idqu = event.getQuark().getPDGid();
	Idaq = event.getAntiQuark().getPDGid();
	double amqu = DataSource::get().masses[abs(Idqu)];
	double amaq = DataSource::get().masses[abs(Idaq)];

	if (abs(p_beamParticleA->getPDGid()) >= 10 && abs(
			p_beamParticleB->getPDGid()) >= 10) {
		// Set quark 4-momenta in beams CMS
		double Ebcm = 0.5 * Ecm;

		double ampro = DataSource::get().masses[abs(p_beamParticleA->getPDGid())];

		//--- 3. Bjorken x as a fraction of proton light-cone momentum in beams CMS
		double pcmp = Ebcm + sqrt(Ebcm * Ebcm - ampro * ampro);
		double pqup = xqu * pcmp;

		if (pqup <= amqu) {
			wtEvent = 0.0;
			return 0.0; // outside phase space
		}
		double pqum = (pqup - amqu) * (pqup + amqu) / 2.0 / pqup;
		double Equ = pqup - pqum;
		double paqp = xaq * pcmp;

		if (paqp <= amaq) {
			wtEvent = 0.0;
			return 0.0; // outside phase space
		}
		double paqm = (paqp - amaq) * (paqp + amaq) / 2.0 / paqp;
		double Eaq = paqp - paqm;

		// Quarks 4-momenta in CMS frame of beams
		// (x,y,z,t)
		event.getQuark().getFourMomentum().set(0.0, 0.0, pqum, Equ);
		event.getAntiQuark().getFourMomentum().set(0.0, 0.0, -paqm, Eaq);
	} else {
		if (p_beamParticleA->getPDGid() > 0 && p_beamParticleB->getPDGid() < 0) {
			event.getQuark().setFourMomentum(p_beamParticleA->getFourMomentum());
			event.getAntiQuark().setFourMomentum(
					p_beamParticleB->getFourMomentum());
		} else if (p_beamParticleA->getPDGid() < 0
				&& p_beamParticleB->getPDGid() > 0) {
			event.getQuark().setFourMomentum(p_beamParticleB->getFourMomentum());
			event.getAntiQuark().setFourMomentum(
					p_beamParticleA->getFourMomentum());
		} else {
			throw VinhacException("Bad initial data of beams");
		}
	}
	// Set quark, antiquark ID's and 4-momenta in their CMS (Z rest frame),
	// +z axis pointing in the quark direction.

	double sqq = (event.getQuark().getFourMomentum()
			+ event.getAntiQuark().getFourMomentum()).mag2();

	double amqq = sqrt(sqq);
	double Equ_wrf = (sqq + amqu * amqu - amaq * amaq) / amqq / 2;
	double pzqu_wrf = sqrt((Equ_wrf - amqu) * (Equ_wrf + amqu));
	// (x,y,z,t)
	event.getQuark().getFourMomentumBosonRestFrame().set(0.0, 0.0, pzqu_wrf,
			Equ_wrf);
	event.getAntiQuark().getFourMomentumBosonRestFrame().set(0.0, 0.0,
			-pzqu_wrf, amqq - Equ_wrf);

	wtEvent = 1.0;

	return 1.0;
}

VINHAC::BeamHandler::~BeamHandler() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::~BeamHandler()"<<std::endl;

#endif

	// Foam printout
	if (DataSource::get().printOutLevel == PrintOutLevel::Detailed && p_FoamX != 0) {
		double MCresult, MCerror;
		double eps = 0.0005;
		double Effic, WtMax, AveWt, Sigma;
		double IntNorm, Errel;
		double nCalls;

		p_FoamX->Finalize(IntNorm, Errel); // final printout
		p_FoamX->GetIntegMC(MCresult, MCerror); // get MC integral
		p_FoamX->GetWtParams(eps, AveWt, WtMax, Sigma); // get MC wt parameters
		nCalls = p_FoamX->GetnCalls();
		Effic = 0;
		if (WtMax > 0)
			Effic = AveWt / WtMax;
		cout
				<< "================================================================"
				<< endl;
		cout << " MCresult = " << MCresult << " +- " << MCerror << " RelErr= "
				<< MCerror / MCresult << endl;
		cout << " Dispersion/<wt>= " << Sigma << " AveWt" << AveWt << endl;
		cout << "      <wt>/WtMax= " << Effic << ",    for epsilon = " << eps
				<< endl;
		cout << " nCalls (initialization only) =   " << nCalls << endl;
		cout
				<< "================================================================"
				<< endl;
	}

	if (p_FoamX) {
		delete p_FoamX;
		p_FoamX = 0;
	}

	if (integrand) {
		delete integrand;
		integrand = 0;
	}

	if(pdfObject) {
		delete pdfObject;
		pdfObject = 0;
	}

	if (p_beamParticleA) {
		delete p_beamParticleA;
		p_beamParticleA = 0;
	}
	if (p_beamParticleB) {
		delete p_beamParticleB;
		p_beamParticleB = 0;
	}

}
void VINHAC::BeamHandler::cleanBeamHandlerData() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamHandler::cleanBeamHandlerData()"<<std::endl;

#endif

	wtFoam = 0;
	wtEvent = 0;
	//// Clean particles data

	p_beamParticleA->clearChildren();
	p_beamParticleB->clearChildren();

}

void VINHAC::BeamHandler::generateFlavours(Event& event) {

	if (p_randomGenerator->Flat() <= 0.5) {
		event.getQuark().addParent(p_beamParticleA);
		event.getAntiQuark().addParent(p_beamParticleB);
	} else {
		event.getQuark().addParent(p_beamParticleB);
		event.getAntiQuark().addParent(p_beamParticleA);
	}

	if (DataSource::get().keyCol == 0) {
		double sumWp = 0.0;
		double sumWm = 0.0;
		double cumWp[3][3];
		double cumWm[3][3];
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j) {
				cumWp[i][j] = 0.0;
				cumWm[i][j] = 0.0;
			}

		bool choosenWp = true;
		if (isWminus) {
			for (vector<int>::iterator it_d = dQuarks.begin(); it_d
					!= dQuarks.end(); ++it_d) {
				for (vector<int>::iterator it_u = uAntiQuarks.begin(); it_u
						!= uAntiQuarks.end(); ++it_u) {
					sumWm += xPDFa[(*it_d) + 6] * xPDFb[(*it_u) + 6]
							* DataSource::get().ckm2[-(*it_u) / 2 - 1][(*it_d) / 2];
					cumWm[-(*it_u) / 2 - 1][(*it_d) / 2] = sumWm;
				}
			}
		}
		if (isWplus) {
			for (vector<int>::iterator it_u = uQuarks.begin(); it_u
					!= uQuarks.end(); ++it_u) {
				for (vector<int>::iterator it_d = dAntiQuarks.begin(); it_d
						!= dAntiQuarks.end(); ++it_d) {
					sumWp += xPDFa[(*it_u) + 6] * xPDFb[(*it_d) + 6]
							* DataSource::get().ckm2[(*it_u) / 2 - 1][-(*it_d) / 2];
					cumWp[(*it_u) / 2 - 1][-(*it_d) / 2] = sumWp;
				}
			}
		}

		if (isWminus && isWplus) {
			double ranWmp = (sumWm + sumWp) * p_randomGenerator->Flat();
			if (ranWmp <= sumWm) {
				choosenWp = false;
			} else {
				choosenWp = true;
			}
		} else if (isWplus) {
			choosenWp = true;
		} else if (isWminus) {
			choosenWp = false;
		} else {
			throw VinhacException("Did not select either W+ or W-");
		}

		double random = p_randomGenerator->Flat();
		if (choosenWp) {
			for (int k = 0; k < 3; ++k)
				for (int l = 0; l < 3; ++l) {
					if (random * sumWp <= cumWp[k][l]) {
						event.getQuark().setPDGid(2 * (k + 1));
						event.getAntiQuark().setPDGid(-2 * (l + 1) + 1);
						return;
					}
				}
		} else {
			for (int k = 0; k < 3; ++k)
				for (int l = 0; l < 3; ++l) {
					if (random * sumWm <= cumWm[l][k]) {
						event.getQuark().setPDGid(2 * (k + 1) - 1);
						event.getAntiQuark().setPDGid(-2 * (l + 1));
						return;
					}
				}
		}
	}

}

