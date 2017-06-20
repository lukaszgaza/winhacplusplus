#include <string>
using namespace std;

#ifndef __DemoAnalysisHandler_h__
#define __DemoAnalysisHandler_h__

#include "Manager.h"
#include "Event.h"
#include "ProcessHandler.h"
#include "Particle.h"
#include <fstream>
#include "../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include <cstring>

#ifdef ROOT_HIST
#include "TH1D.h"
#include "TFile.h"
#include "TObject.h"
#endif

namespace VINHAC {
class Manager;
class Event;
class Particle;
class AnalisysHandler;
}

namespace VINHAC {

//!  demonstration class showing how to write your own handler.
/*!
 * It is showed here how to write your own handler by extending ProcessHandler
 * and you can pass it to generator using Manager (i.e. WINHAC) methods.
 * This example just cumulates some observables using ROOT histograms.
 * @see ProcessHandler
 * @see Manager
 * @see WINHAC
*/
class DemoAnalysisHandler: public VINHAC::ProcessHandler {
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

	TH1D hst_photonsEnergy;
	TH1D hst_hardPhotonEnergy;
	TH1D hst_hardPhotonPt;
	TH1D hst_hardPhotonEta;
#endif

	VINHAC::Manager* p_Manager;


public:
	/**
	* \brief evolves the event
	*
	* It does whole part of generation process which is to do in this handler
	* @see Event
	* @param event reference to current Event object
	* @return if weight was zero
	*/
	bool eventEvolve(VINHAC::Event& event);


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
	/**
	* \brief saves histograms to specified path
	*
	* @param snapshotPath path to folder where histograms have to be saved
	*/
	bool makeSnapshot(string snapshotPath);
	DemoAnalysisHandler();
	~DemoAnalysisHandler(){};

}; // class DemoAnalysisHandler
} // namespace VINHAC

#endif
