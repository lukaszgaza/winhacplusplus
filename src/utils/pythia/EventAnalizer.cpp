/*
 * EventAnalizer.cpp
 *
 *  Created on: 2011-04-04
 *      Author: kamil
 */

#include "EventAnalizer.h"

EventAnalizer::EventAnalizer()
#ifdef ROOT_HIST
:
			numberOfBins(100),
			hst_W_rapidity("W: rapidity y", "W: rapidity y", numberOfBins, -6,
					6),
			hst_cosTheta("lepton: cos(theta)", "lepton: cos(theta)",
					numberOfBins, -1, 1),
			hst_leptonPT("lepton: p_T", "lepton: p_T", numberOfBins, 0, 50),
			hst_leptonPseudorapidity("lepton: pseudorapidity eta",
					"lepton: pseudorapidity eta", numberOfBins, -6, 6),
			hst_photon_multiplicity("photon: multiplicity",
					"photon: multiplicity", 10, 0, 10),
			hst_WTransverseMass("W: transverse mass", "W: transverse mass",
					numberOfBins, 0, 150),
			hst_WMass("W: mass", "W: mass",
					numberOfBins, 0, 160),
			hst_WPT("W: transverse momentum", "W: transverse momentum",
					numberOfBins, 0, 50),
			hst_photonsEnergy("photon: energy sum", "photon: energy sum",
					numberOfBins, -1, 50),
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
}

EventAnalizer::~EventAnalizer() {
}

void EventAnalizer::analize(Pythia8::Event &event, double weight) {

	// Extracting data from Event object
	Pythia8::Particle boson;
	Pythia8::Particle lepton;
	Pythia8::Particle neutrino;
	std::vector < Pythia8::Particle > photons;

	for (int i = 0; i < event.size(); ++i) {
		int absid = abs(event[i].id());
		if (absid == 24)
			boson = event[i];
		if (absid == 11 || absid == 13 || absid == 15)
			lepton = event[i];
		if (absid == 12 || absid == 14 || absid == 16)
			neutrino = event[i];
		if (absid == 22)
			photons.push_back(event[i]);
	}

	double yw = -50.0;
	double cosTheta = -2.0;
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
	double mainWt = weight;

	// Construction of observables
	if (mainWt != 0.0) {

		double EnW = boson.e();
		double pzW = boson.pz();
		yw = 0.5 * log((EnW + pzW) / (EnW - pzW));
		Pythia8::Vec4 pl = lepton.p();
		Pythia8::Vec4 pn = neutrino.p();
		cosTheta = pl.pz() / pl.pAbs();
		double theta = acos(cosTheta);
		etal = -log(tan(theta / 2));
		pTl = lepton.p().pT();
		pTn = neutrino.p().pT();

		wMass = boson.p().mCalc();
		wTransverseMass = sqrt(
				2 * (pTl * pTn - pl.px() * pn.px() - pl.py() * pn.py()));
		//wPT = boson.pT();
		wPT = sqrt(
				(lepton.p().px() + neutrino.p().px()) * (lepton.p().px()
						+ neutrino.p().px()) + (lepton.p().py()
						+ neutrino.p().py()) * (lepton.p().py()
						+ neutrino.p().py()));

		photonsEnergy = 0.0;
		int hardest = -1;
		double tmp = 0.0;
		for (unsigned i = 0; i < photons.size(); ++i) {
			photonsEnergy += photons[i].e();
			if (photons[i].e() > tmp) {
				hardest = i;
				tmp = photons[i].e();
			}
		}

		if (hardest >= 0) {
			hardPhotonE = photons[hardest].e();
			hardPhotonPt = photons[hardest].p().pT();

			double cosThetaPhoton = photons[hardest].pz()
					/ photons[hardest].p().pAbs();
			double thetaPhoton = acos(cosThetaPhoton);
			hardPhotonEta = -log(tan(thetaPhoton / 2));

		}

	}

	// Filling histograms
#ifdef ROOT_HIST
	hst_W_rapidity.Fill(yw, mainWt);
	hst_cosTheta.Fill(cosTheta, mainWt);
	hst_leptonPT.Fill(pTl, mainWt);
	hst_leptonPseudorapidity.Fill(etal, mainWt);
	hst_WTransverseMass.Fill(wTransverseMass, mainWt);
	hst_WMass.Fill(wMass, mainWt);
	hst_WPT.Fill(wPT, mainWt);
	hst_photon_multiplicity.Fill(photons.size(), mainWt);
	hst_photonsEnergy.Fill(photonsEnergy, mainWt);
	hst_hardPhotonEnergy.Fill(hardPhotonE, mainWt);
	hst_hardPhotonPt.Fill(hardPhotonPt, mainWt);
	hst_hardPhotonEta.Fill(hardPhotonEta, mainWt);

#endif
}

void EventAnalizer::scaleHistograms(double xsection) {
#ifdef ROOT_HIST
	hst_W_rapidity.Scale(xsection);
	hst_cosTheta.Scale(xsection);
	hst_leptonPT.Scale(xsection);
	hst_leptonPseudorapidity.Scale(xsection);
	hst_WTransverseMass.Scale(xsection);
	hst_WMass.Scale(xsection);
	hst_WPT.Scale(xsection);
	hst_photon_multiplicity.Scale(xsection);
	hst_photonsEnergy.Scale(xsection);
	hst_hardPhotonEnergy.Scale(xsection);
	hst_hardPhotonPt.Scale(xsection);
	hst_hardPhotonEta.Scale(xsection);

#endif
}

bool EventAnalizer::makeSnapshot(string snapshotPath) {

#ifdef ROOT_HIST
	string filePath = snapshotPath + "/histograms.root";
	string filePathBak = snapshotPath + "/histogramsBak.root";
	string mv_file_comand = "mv " + filePath + " " + filePathBak;
	system(mv_file_comand.c_str());

	TFile rootSanpshot(filePath.c_str(), "RECREATE", "Snapshot histograms");

	hst_W_rapidity.Write();
	hst_cosTheta.Write();
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
