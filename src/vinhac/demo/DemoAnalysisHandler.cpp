#include <string>
#include <cstdlib>
using namespace std;

#include "DemoAnalysisHandler.h"
#include "../core/Event.h"
#include "../core/Particle.h"
#include "../input/DataSource.h"

// Creates histograms objects
VINHAC::DemoAnalysisHandler::DemoAnalysisHandler()
#ifdef ROOT_HIST
:
	numberOfBins(100),
	hst_W_rapidity("W: rapidity y", "W: rapidity y",
		numberOfBins, -6, 6),
	hst_cosTheta("lepton: cos(theta)",
		"lepton: cos(theta)", numberOfBins, -1, 1),
	hst_phi("lepton: phi",
			"lepton: phi", numberOfBins, -4, 4),
	hst_cosThetaWL("lepton: cos(thetaWL)",
			"lepton: cos(thetaWL)", numberOfBins, -1, 1),
	hst_phiWL("lepton: phiWL",
			"lepton: phiWL", numberOfBins, 0, 4),
	hst_leptonPT("lepton: p_T", "lepton: p_T", numberOfBins, 0, 50),
	hst_leptonPseudorapidity("lepton: pseudorapidity eta",
		"lepton: pseudorapidity eta", numberOfBins, -6, 6),
	hst_photon_multiplicity("photon: multiplicity",
		"photon: multiplicity", 10, 0, 10),
	hst_WTransverseMass(
		"W: transverse mass", "W: transverse mass", numberOfBins,
		0, 150),
	hst_WMass("W: mass", "W: mass",
		numberOfBins, 0, 160),
	hst_WPT(
		"W: transverse momentum", "W: transverse momentum", numberOfBins,
		0, 50),
	hst_photonsEnergy("photon: energy sum",
		"photon: energy sum", numberOfBins, -1, 50),
	hst_hardPhotonEnergy("photon: hard photon energy",
		"photon: hard photon energy", numberOfBins, -1, 50),
	hst_hardPhotonPt("photon: hard photon pt",
		"photon: hard photon pt", numberOfBins, -1, 50),
	hst_hardPhotonEta("photon: hard photon eta",
		"photon: hard photon eta", numberOfBins, -6, 6)
#endif
{

	// Forces histograms to store sum of weights squared
#ifdef ROOT_HIST
	hst_W_rapidity.Sumw2();
	hst_cosTheta.Sumw2();
	hst_phi.Sumw2();
	hst_cosThetaWL.Sumw2();
	hst_phiWL.Sumw2();
	hst_leptonPT.Sumw2();
	hst_leptonPseudorapidity.Sumw2();
	hst_photon_multiplicity.Sumw2();
	hst_WTransverseMass.Sumw2();
	hst_WMass.Sumw2();
	hst_WPT.Sumw2();
	hst_photonsEnergy.Sumw2();
	hst_hardPhotonEnergy.Sumw2();
	hst_hardPhotonPt.Sumw2();
	hst_hardPhotonEta.Sumw2();
#endif

	setHandlerName("DemoAnalysisHandler");

}
bool VINHAC::DemoAnalysisHandler::eventEvolve(VINHAC::Event& event) {

	// Extracting data from Event object
	Particle& boson = event.getBoson();
	Particle& lepton = event.getLepton();
	Particle& neutrino = event.getNeutrino();

	double yw = -50.0;
	double cosTheta = -2.0;
	double phi = -10.0;
	double cosThetaWL = -2.0;
	double phiWL = -10.0;
	double etal = -50.0;
	double pTl = -10.0;
	double pTn = -10.0;
	double wTransverseMass = -10.0;
	double wMass = -10.0;
	double wPT = -10.0;
	double photonsEnergy = -10.0;
	double hardPhotonE = -10.0;
	double hardPhotonPt = -10.0;
	double hardPhotonEta = -10.0;

	// Extracting main weight from event
	double mainWt = event.mainWt();

	// Construction of observables
	if (mainWt != 0.0) {

	double EnW = boson.getFourMomentum(DataSource::get().modelReferenceFrame).e();
	double pzW = boson.getFourMomentum(DataSource::get().modelReferenceFrame).pz();
	yw = 0.5 * log((EnW + pzW) / (EnW - pzW));
	const CLHEP::HepLorentzVector& pl = lepton.getFourMomentum(DataSource::get().modelReferenceFrame);
	const CLHEP::HepLorentzVector& pn = neutrino.getFourMomentum(DataSource::get().modelReferenceFrame);
	const CLHEP::HepLorentzVector& pW = boson.getFourMomentum(DataSource::get().modelReferenceFrame);
	cosTheta = pl.z() / pl.getV().mag();
	cosThetaWL = pl.getV().dot(pW.getV())/pW.getV().r()/pl.getV().r();
	phi = pl.getV().phi();
	phiWL = acos((pl.x()*pW.x()+pl.y()*pW.y())/pl.getV().rho()/pW.getV().rho());
	double theta = acos(cosTheta);
	etal = -log(tan(theta / 2));
	pTl = lepton.getFourMomentum(DataSource::get().modelReferenceFrame).getV().rho();
	pTn = neutrino.getFourMomentum(DataSource::get().modelReferenceFrame).getV().rho();

	wMass = boson.getFourMomentum(DataSource::get().modelReferenceFrame).m();
	wTransverseMass = sqrt(2 * (pTl * pTn - pl.x() * pn.x() - pl.y()
		* pn.y()));

	wPT= sqrt((lepton.getFourMomentum(DataSource::get().modelReferenceFrame).px()+neutrino.getFourMomentum(DataSource::get().modelReferenceFrame).px())*(lepton.getFourMomentum(DataSource::get().modelReferenceFrame).px()+neutrino.getFourMomentum(DataSource::get().modelReferenceFrame).px())+
				(lepton.getFourMomentum(DataSource::get().modelReferenceFrame).py()+neutrino.getFourMomentum(DataSource::get().modelReferenceFrame).py())*(lepton.getFourMomentum(DataSource::get().modelReferenceFrame).py()+neutrino.getFourMomentum(DataSource::get().modelReferenceFrame).py()));

	photonsEnergy = 0.0;
	for (unsigned i = 0; i < event.getPhotons().size(); ++i) {
		photonsEnergy += event.getPhotons()[i].getFourMomentum(DataSource::get().modelReferenceFrame).e();
	}

	if (event.getPhotons().size() > 0) {
		hardPhotonE =
			event.getPhotons()[0].getFourMomentum(DataSource::get().modelReferenceFrame).e();
		hardPhotonPt =
			event.getPhotons()[0].getFourMomentum(DataSource::get().modelReferenceFrame).getV().rho();

		double cosThetaPhoton =
			event.getPhotons()[0].getFourMomentum(DataSource::get().modelReferenceFrame).z()
			/ event.getPhotons()[0].getFourMomentum(DataSource::get().modelReferenceFrame).getV().mag();
		double thetaPhoton = acos(cosThetaPhoton);
		hardPhotonEta = -log(tan(thetaPhoton / 2));

	}

	}

	// Filling histograms
#ifdef ROOT_HIST
	hst_W_rapidity.Fill(yw, mainWt);
	hst_cosTheta.Fill(cosTheta, mainWt);
	hst_phi.Fill(phi, mainWt);
	hst_cosThetaWL.Fill(cosThetaWL, mainWt);
	hst_phiWL.Fill(phiWL, mainWt);
	hst_leptonPT.Fill(pTl, mainWt);
	hst_leptonPseudorapidity.Fill(etal, mainWt);
	hst_WTransverseMass.Fill(wTransverseMass, mainWt);
	hst_WMass.Fill(wMass, mainWt);
	hst_WPT.Fill(wPT, mainWt);
	hst_photon_multiplicity.Fill(event.getPhotons().size(), mainWt);
	hst_photonsEnergy.Fill(photonsEnergy, mainWt);
	hst_hardPhotonEnergy.Fill(hardPhotonE, mainWt);
	hst_hardPhotonPt.Fill(hardPhotonPt, mainWt);
	hst_hardPhotonEta.Fill(hardPhotonEta, mainWt);

#endif

	return true;
}

bool VINHAC::DemoAnalysisHandler::makeSnapshot(string snapshotPath) {

#ifdef ROOT_HIST
	string filePath = snapshotPath + "/histograms.root";
	string filePathBak = snapshotPath + "/histogramsBak.root";
	string mv_file_comand = "mv " + filePath + " " + filePathBak;
	system(mv_file_comand.c_str());

	TFile rootSanpshot(filePath.c_str(), "RECREATE", "Snapshot histograms");

	hst_W_rapidity.Write();
	hst_cosTheta.Write();
	hst_phi.Write();
	hst_cosThetaWL.Write();
	hst_phiWL.Write();
	hst_leptonPT.Write();
	hst_leptonPseudorapidity.Write();
	hst_photon_multiplicity.Write();
	hst_WTransverseMass.Write();
	hst_WMass.Write();
	hst_WPT.Write();
	hst_photonsEnergy.Write();
	hst_hardPhotonEnergy.Write();
	hst_hardPhotonPt.Write();
	hst_hardPhotonEta.Write();

	rootSanpshot.Write();
	rootSanpshot.Close();
#endif

	return true;
}
