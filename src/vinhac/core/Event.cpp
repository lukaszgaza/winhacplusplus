#include <string>
#include <vector>
#include <exception>
#include <sstream>
using namespace std;

#include "Event.h"
#include "../input/DataSource.h"
#include "../util/Commons.h"

VINHAC::Event::Event() :
	pdfUsed(false), quarkX(1), antiquarkX(1), quarkXpdf(1), antiquarkXpdf(1),
			pdfScale(-1), quarkPDFid(0), antiquarkPDFid(0) {
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
	photons.reserve(3);
}

void VINHAC::Event::printParticle(stringstream& ss, Particle& p, FrameName::name name) {
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
	ss << p.getFourMomentum(name).x();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum(name).y();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum(name).z();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum(name).e();

	ss.flags(ios::right | ios::showpoint | ios::fixed);
	ss.width(width);
	ss.precision(13);
	ss << p.getFourMomentum(name).mag();

	ss << endl;
}

void VINHAC::Event::printParticle(stringstream& ss, Particle& p) {
	printParticle(ss,p,FrameName::LAB);
}

void VINHAC::Event::printParticleWrf(stringstream& ss, Particle& p) {
	printParticle(ss,p,FrameName::bosonRestFrame);
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
	lnu.getFourMomentum() += lepton.getFourMomentum();
	lnu.getFourMomentum() += neutrino.getFourMomentum();
	//Insufficient
	//lnu.setFourMomentum(lepton.getFourMomentum() + neutrino.getFourMomentum());

	if (photons.size() > 0) {
		lnu.setTag("Lep+Neu+Pho");
		for (unsigned i = 0; i < photons.size(); ++i) {
			lnu.getFourMomentum() += photons[i].getFourMomentum();
			//Insufficient
			//lnu.setFourMomentum(lnu.getFourMomentum()
			//		+ photons[i].getFourMomentum());
		}
	}

	printParticle(ss, lnu);

	Particle shouldBe;
	shouldBe.setTag("Beams-Rests");
	shouldBe.setPDGid(0);

	shouldBe.getFourMomentum() += beamParticleA.getFourMomentum();
	shouldBe.getFourMomentum() -= beamRestA.getFourMomentum();
	shouldBe.getFourMomentum() += beamParticleB.getFourMomentum();
	shouldBe.getFourMomentum() -= beamRestB.getFourMomentum();
	//Insufficient
	//shouldBe.setFourMomentum((beamParticleA.getFourMomentum()
	//		- beamRestA.getFourMomentum()) + (beamParticleB.getFourMomentum()
	//		- beamRestB.getFourMomentum()));

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

	lnu.getFourMomentumBosonRestFrame()
			+= lepton.getFourMomentumBosonRestFrame();
	lnu.getFourMomentumBosonRestFrame()
			+= neutrino.getFourMomentumBosonRestFrame();
	//Insufficient
	//lnu.setFourMomentumBosonRestFrame(lepton.getFourMomentumBosonRestFrame()
	//		+ neutrino.getFourMomentumBosonRestFrame());

	if (photons.size() > 0) {
		lnu.setTag("Lep+Neu+Pho");
		for (unsigned i = 0; i < photons.size(); ++i) {

			lnu.getFourMomentumBosonRestFrame()
					+= photons[i].getFourMomentumBosonRestFrame();
			//Insufficient
			//lnu.setFourMomentumBosonRestFrame(lnu.getFourMomentumBosonRestFrame()
			//		+ photons[i].getFourMomentumBosonRestFrame());
		}
	}

	printParticleWrf(ss, lnu);
	ss
			<< "======================================================================================================================";

	return ss.str();
}

string VINHAC::Event::printInFrame(FrameName::name name) {
	stringstream ss;

	int width = 18;

	ss
			<< "===================================================  (Custom frame)  ================================================="
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
	printParticle(ss, beamParticleA, name);
	printParticle(ss, beamParticleB, name);
	printParticle(ss, beamRestA, name);
	printParticle(ss, beamRestB, name);
	printParticle(ss, quark, name);
	printParticle(ss, antiQuark, name);
	printParticle(ss, boson, name);
	printParticle(ss, lepton, name);
	printParticle(ss, neutrino, name);
	for (unsigned i = 0; i < photons.size(); ++i) {
		printParticle(ss, photons[i], name);
	}
	ss
			<< "----------------------------------------------------------------------------------------------------------------------"
			<< endl;

	Particle lnu;
	lnu.setTag("Lep+Neu");
	lnu.setPDGid(0);
	lnu.getFourMomentum(name) += lepton.getFourMomentum(name);
	lnu.getFourMomentum(name) += neutrino.getFourMomentum(name);
	//Insufficient
	//lnu.setFourMomentum(lepton.getFourMomentum() + neutrino.getFourMomentum());

	if (photons.size() > 0) {
		lnu.setTag("Lep+Neu+Pho");
		for (unsigned i = 0; i < photons.size(); ++i) {
			lnu.getFourMomentum(name) += photons[i].getFourMomentum(name);
			//Insufficient
			//lnu.setFourMomentum(lnu.getFourMomentum()
			//		+ photons[i].getFourMomentum());
		}
	}

	printParticle(ss, lnu, name);

	Particle shouldBe;
	shouldBe.setTag("Beams-Rests");
	shouldBe.setPDGid(0);

	shouldBe.getFourMomentum(name) += beamParticleA.getFourMomentum(name);
	shouldBe.getFourMomentum(name) -= beamRestA.getFourMomentum(name);
	shouldBe.getFourMomentum(name) += beamParticleB.getFourMomentum(name);
	shouldBe.getFourMomentum(name) -= beamRestB.getFourMomentum(name);
	//Insufficient
	//shouldBe.setFourMomentum((beamParticleA.getFourMomentum()
	//		- beamRestA.getFourMomentum()) + (beamParticleB.getFourMomentum()
	//		- beamRestB.getFourMomentum()));

	printParticle(ss, shouldBe, name);
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

void VINHAC::Event::setMainWt(double mainWt) {
	m_commonWeightsMap.clear();
	m_modelWeightsMap.clear();
	m_modelWeightsMap[DataSource::get().mainWeight] = mainWt;
}

double VINHAC::Event::mainWt() const {
	double result = 1.0;
	for (map<WeightName::name, double>::const_iterator it =
			m_commonWeightsMap.begin(); it != m_commonWeightsMap.end(); ++it) {
		result *= (*it).second;
	}

	map<pair<WeightName::name, PolarizationName::name> , double>::const_iterator
			mainWeight = m_modelWeightsMap.find(DataSource::get().mainWeight);
	if (mainWeight != m_modelWeightsMap.end()) {
		result *= mainWeight->second;
	}

	return result;
}

