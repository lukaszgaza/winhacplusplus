#include <string>
#include <vector>
#include <exception>
#include <sstream>
using namespace std;

#include "Event.h"

VINHAC::Event::Event() {
	beamParticleA.setTag("BeamA");
	beamParticleB.setTag("BeamB");
	beamRestA.setTag("BeamRestA");
	beamRestB.setTag("BeamRestB");
	quark.setTag("Quark");
	antiQuark.setTag("Antiquark");
	boson.setTag("Boson");
	lepton.setTag("Lepton");
	neutrino.setTag("Neutrino");

	beamRestA.setPDGid(0);
	beamRestB.setPDGid(0);

}

void VINHAC::Event::printParticle(stringstream& ss, Particle& p) {
	int width = 20;
	//int precision = 13;

	ss.flags(ios::left);
	ss.width(12);
	ss << p.getTag();

	ss.flags(ios::right);
	ss.width(6);
	ss << p.getPDGid();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum().x();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum().y();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum().z();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum().e();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum().mag();

	ss << endl;
}

void VINHAC::Event::printParticleWrf(stringstream& ss, Particle& p) {
	int width = 20;
	//int precision = 13;

	ss.flags(ios::left);
	ss.width(12);
	ss << p.getTag();

	ss.flags(ios::right);
	ss.width(6);
	ss << p.getPDGid();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentumBosonRestFrame().x();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentumBosonRestFrame().y();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentumBosonRestFrame().z();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentumBosonRestFrame().e();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentumBosonRestFrame().mag();

	ss << endl;
}

string VINHAC::Event::printInLabFrame() {
	stringstream ss;

	int width = 18;

	ss
			<< "===================================================  (LAB frame)  ===================================================="
			<< endl;
	ss.flags(ios::left);
	ss.width(12);
	ss << "Particle";
	ss.flags(ios::right);
	ss.width(6);
	ss << "Code";
	ss.flags(ios::right);
	ss.width(width);
	ss << "px";
	ss.flags(ios::right);
	ss.width(width);
	ss << "py";
	ss.flags(ios::right);
	ss.width(width);
	ss << "pz";
	ss.flags(ios::right);
	ss.width(width);
	ss << "E";
	ss.flags(ios::right);
	ss.width(width);
	ss << "Mass";

	ss << endl;
	ss
			<< "----------------------------------------------------------------------------------------------------------------------"
			<< endl;
	printParticle(ss, beamParticleA);
	printParticle(ss, beamParticleB);
	printParticle(ss, beamRestA);
	printParticle(ss, beamRestB);
	printParticle(ss, quark);
	printParticle(ss, antiQuark);
	printParticle(ss, boson);
	printParticle(ss, lepton);
	printParticle(ss, neutrino);
	for (unsigned i = 0; i < photons.size(); ++i) {
		printParticle(ss, photons[i]);
	}
	ss
			<< "----------------------------------------------------------------------------------------------------------------------"
			<< endl;

	Particle lnu;
	lnu.setTag("Lep+Neu");
	lnu.setPDGid(0);
	lnu.setFourMomentum(lepton.getFourMomentum() + neutrino.getFourMomentum());

	if (photons.size() > 0) {
		lnu.setTag("Lep+Neu+Pho");
		for (unsigned i = 0; i < photons.size(); ++i) {
			lnu.setFourMomentum(lnu.getFourMomentum()
					+ photons[i].getFourMomentum());
		}
	}

	printParticle(ss, lnu);

	Particle shouldBe;
	shouldBe.setTag("Beams-Rests");
	shouldBe.setPDGid(0);
	shouldBe.setFourMomentum((beamParticleA.getFourMomentum()
			- beamRestA.getFourMomentum()) + (beamParticleB.getFourMomentum()
			- beamRestB.getFourMomentum()));

	printParticle(ss, shouldBe);
	ss
			<< "======================================================================================================================";

	return ss.str();
}

string VINHAC::Event::printInQuarksCMSFrame() {
	stringstream ss;

	int width = 18;

	ss
			<< "===================================================  (Quarks CMS)  ==================================================="
			<< endl;
	ss.flags(ios::left);
	ss.width(12);
	ss << "Particle";
	ss.flags(ios::right);
	ss.width(6);
	ss << "Code";
	ss.flags(ios::right);
	ss.width(width);
	ss << "px";
	ss.flags(ios::right);
	ss.width(width);
	ss << "py";
	ss.flags(ios::right);
	ss.width(width);
	ss << "pz";
	ss.flags(ios::right);
	ss.width(width);
	ss << "E";
	ss.flags(ios::right);
	ss.width(width);
	ss << "Mass";

	ss << endl;
	ss
			<< "----------------------------------------------------------------------------------------------------------------------"
			<< endl;
	printParticleWrf(ss, quark);
	printParticleWrf(ss, antiQuark);
	printParticleWrf(ss, boson);
	printParticleWrf(ss, lepton);
	printParticleWrf(ss, neutrino);
	for (unsigned i = 0; i < photons.size(); ++i) {
		printParticleWrf(ss, photons[i]);
	}
	ss
			<< "----------------------------------------------------------------------------------------------------------------------"
			<< endl;

	Particle lnu;
	lnu.setTag("Lep+Neu");
	lnu.setPDGid(0);
	lnu.setFourMomentumBosonRestFrame(lepton.getFourMomentumBosonRestFrame()
			+ neutrino.getFourMomentumBosonRestFrame());

	if (photons.size() > 0) {
		lnu.setTag("Lep+Neu+Pho");
		for (unsigned i = 0; i < photons.size(); ++i) {
			lnu.setFourMomentumBosonRestFrame(
					lnu.getFourMomentumBosonRestFrame()
							+ photons[i].getFourMomentumBosonRestFrame());
		}
	}

	printParticleWrf(ss, lnu);
	ss
			<< "======================================================================================================================";

	return ss.str();
}

bool VINHAC::Event::clean() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Event::clean()"<<std::endl;

#endif

	m_commonWeightsMap.clear();
	m_modelWeightsMap.clear();
	return true;
}

