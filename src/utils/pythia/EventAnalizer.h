/*
 * EventAnalizer.h
 *
 *  Created on: 2011-04-04
 *      Author: kamil
 */

#ifndef EVENTANALIZER_H_
#define EVENTANALIZER_H_

#include "Pythia.h"

#ifdef ROOT_HIST
#include "TH1D.h"
#include "TFile.h"
#include "TObject.h"
#endif

class EventAnalizer {
private:
#ifdef ROOT_HIST
	////ROOT histograms:
	TFile hfile;

	int numberOfBins;
	TH1D hst_W_rapidity;
	TH1D hst_cosTheta;
	TH1D hst_leptonPT;
	TH1D hst_leptonPseudorapidity;
	TH1D hst_photon_multiplicity;
	TH1D hst_WTransverseMass;
	TH1D hst_WMass;
	TH1D hst_WPT;

	TH1D hst_photonsEnergy;
	TH1D hst_hardPhotonEnergy;
	TH1D hst_hardPhotonPt;
	TH1D hst_hardPhotonEta;
#endif

public:
	EventAnalizer();
	virtual ~EventAnalizer();

	void analize(Pythia8::Event &event,double weight);
	bool makeSnapshot(string snapshotPath);
	void scaleHistograms(double xsection);


#ifdef ROOT_HIST
	inline TH1D& getWrapidityHist() {
		return hst_W_rapidity;
	}
	inline TH1D& getCosThetaHist() {
		return hst_cosTheta;
	}
	inline TH1D& getLeptonPTHist() {
		return hst_leptonPT;
	}
	inline TH1D& getLeptonPseudorapidityHist() {
		return hst_leptonPseudorapidity;
	}
	inline TH1D& getWTransverseMass() {
		return hst_WTransverseMass;
	}
	inline TH1D& getWMass() {
		return hst_WMass;
	}
	inline TH1D& getWPT() {
		return hst_WPT;
	}
	inline TH1D& getPhotonMultiplicity() {
		return hst_photon_multiplicity;
	}
	inline TH1D& getPhotonsEnergy() {
		return hst_photonsEnergy;
	}
	inline TH1D& getHardPhotonEnergy() {
		return hst_hardPhotonEnergy;
	}
	inline TH1D& getHardPhotonPt() {
		return hst_hardPhotonPt;
	}
	inline TH1D& getHardPhotonEta() {
		return hst_hardPhotonEta;
	}
#endif
};

#endif /* EVENTANALIZER_H_ */
