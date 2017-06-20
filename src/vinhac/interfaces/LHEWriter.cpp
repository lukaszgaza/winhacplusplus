/*
 * LHEWriter.cpp
 *
 *  Created on: 2010-11-04
 *      Author: kamil
 */

#include "LHEWriter.h"
#include "../core/Manager.h"
#include "../core/Event.h"
#include "../core/Particle.h"
#include "../handlers/BeamHandler.h"

namespace VINHAC {

LHEWriter::LHEWriter(const Manager& manager, std::ostream& os) :
	manager(manager), os(os) {
}
;

void LHEWriter::writeOpening() {
	// Write header.
	os << "<LesHouchesEvents version=\"1.0\">\n" << "<!--\n"
			<< "  File written by WINHAC++ \n" << "-->" << std::endl;
	os.flush();
}

void LHEWriter::writeInit() {
	// Write information on beams.
	os << "<init>\n"
	<< "  " << DataSource::get().beamA.find("PDGid")->second //ID of beam particle A according to the Particle Data Group
	<< "  " << DataSource::get().beamB.find("PDGid")->second //ID of beam particle B according to the Particle Data Group
	<< "  " << DataSource::get().beamA.find("energy")->second //energy in GeV of beam particle A
	<< "  " << DataSource::get().beamB.find("energy")->second //energy in GeV of beam particle B
	<< "  " << manager.p_BeamHandler->pdflibNGroupA //the author group for beam A, according to the Cernlib PDFlibspecification
	<< "  " << manager.p_BeamHandler->pdflibNGroupB //the author group for beam B, according to the Cernlib PDFlibspecification
	<< "  " << manager.p_BeamHandler->pdflibNSetA //the PDF set ID for beam A, according to the Cernlib PDFlib specification
	<< "  " << manager.p_BeamHandler->pdflibNSetB //the PDF set ID for beam A, according to the Cernlib PDFlib specification
	<< "  " << -4 //master switch dictating how the event weights are interpreted
	<< "  " << 1 // number of processes passed through interface
	<< "\n";

	// Write information on all the subprocesses.
	os << " " << 0 //subprocess xsection in pb this is mandatory for weight switch +/-2
	<< " " << 0 //subprocess xsection error in pb this is mandatory for weight switch +/-2
	<< " " << 0 //maximum weight , no meaning for weigth switch +/- 3,4 */
	<< " " << 222 //ID of the process
	<< "\n";


	os << "</init>" << std::endl;

	os.flush();
}

void LHEWriter::writeClosing() {
	// Write an end to the file.
	os << "</LesHouchesEvents>" << std::endl;
	os.flush();
}

void LHEWriter::writeEvent(Event & event) {
	// Write information on process as such.
	os << "<event>\n" << std::scientific << std::setprecision(6)
	   <<getNumberOfParticles(event)
	   << " " << 222 //ID of the process for this event
	   << " " << getEventWeight(event)
	   << " " << (event.isPDFUsed()?event.getPDFScale():-1) //scale of the event in GeV, as used for calculation of PDFs,
													        //If the scale has not been defined, this should be denoted by setting the scale to –1.
	   << " " << DataSource::get().alphaQED //the QED coupling
	   << " " << DataSource::get().alphaQCD //the QCD coupling
	<< "\n";

	Particle *partonA = 0;
	Particle *partonB = 0;
	double partonAx = 0.0;
	double partonAxfx = 0.0;
	double partonBx = 0.0;
	double partonBxfx = 0.0;

	if(event.getQuark().getFourMomentum().z()>0){
		partonA = &event.getQuark();
		partonAx = event.getQuarkX();
		partonAxfx = event.getQuarkXpdf();
		partonB = &event.getAntiQuark();
		partonBx = event.getAntiquarkX();
		partonBxfx = event.getAntiquarkXpdf();
	} else {
		partonB = &event.getQuark();
		partonBx = event.getQuarkX();
		partonBxfx = event.getQuarkXpdf();
		partonA = &event.getAntiQuark();
		partonAx = event.getAntiquarkX();
		partonAxfx = event.getAntiquarkXpdf();
	}

	// Write information on the particles, excluding zeroth.
	if(partonA->getPDGid()<0){
		printParticle(*partonA,-1,0,0,0,101);
	} else {
		printParticle(*partonA,-1,0,0,101);
	}
	if(partonB->getPDGid()<0){
		printParticle(*partonB,-1,0,0,0,101);
	} else {
		printParticle(*partonB,-1,0,0,101);
	}
	printParticle(event.getBoson(),2,1,2);
	printParticle(event.getLepton(),1,3,3);
	printParticle(event.getNeutrino(),1,3,3);
	for (unsigned i = 0; i < event.getPhotons().size(); i++) {
		printParticle(event.getPhotons()[i],1,3,3); //Parent of photon = boson & final lepton
	}

	if(event.isPDFUsed()){
		os << "#pdf"
		<< " " << partonA->getPDGid() // id1
		<< " " << partonB->getPDGid() // id2
		<< " " << partonAx // x1
		<< " " << partonBx // x2
		<< " " << event.getPDFScale() // scale PDF
		<< " " << partonAxfx // xpdf1
		<< " " << partonBxfx // xpdf2
		<< "\n";
	}
	os << "</event>" << std::endl;

	os.flush();
}

int LHEWriter::getNumberOfParticles(Event& event) {
	int result = event.getPhotons().size();

	result += 5; // easiest parton lvl process, q q_bar , W , lepton, nu. Hadron lvl is inluded in init information, it is the same for every event.

	return result;
}

double LHEWriter::getEventWeight(Event& event) {
	//Total best weight multiplied by crude
	//In LHE units of xsection are pb
	return event.mainWt()*manager.getCrudeXsection()*1000.0;
}

void LHEWriter::printParticle(Particle& particle, int statusPart,
		int mother1Part, int mother2Part,int colorTag,int antiColorTag) {

	os << " " << particle.getPDGid() // PDGid of Particle
	<< " " << statusPart // status according to LHE spec
	<< " " << mother1Part // index of mother particle 1
	<< " " << mother2Part // index of mother particle 2
		           << " " << colorTag // color flow
		           << " " << antiColorTag << std::setprecision(10) // anticolor flow
		           << " " << particle.getFourMomentum().x()
		           << " " << particle.getFourMomentum().y()
		           << " " << particle.getFourMomentum().z()
		           << " " << particle.getFourMomentum().e()
		           << " " << particle.getFourMomentum().mag() << std::setprecision(6) // mass
		    << " " << 0 // invariant lifetime ctau (distance from production to decay) in mm
		    << " " << 9; //cosine of the angle between the spin-vector of particle I and the 3-
						 //momentum of the decaying particle, specified in the lab frame
						 //Unknown or unpolarized particles should be given 9
	os << "\n";
}

}

