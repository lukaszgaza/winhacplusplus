#include <string>
#include <vector>
#include <exception>
using namespace std;

#include "Particle.h"

using namespace VINHAC;

Particle::Particle() {
}

Particle::Particle(const Particle& target) :
	_PDGid(target._PDGid), p(target.p), p_brf2(target.p_brf2), m_tag(
			target.m_tag) {
	for (unsigned i = 0; i < target._parents.size(); ++i) {
		_parents.push_back(target._parents.at(i));
	}
	for (unsigned i = 0; i < target._children.size(); ++i) {
		_children.push_back(target._children.at(i));
	}

}

Particle::Particle(int PDGid) {
	setPDGid(PDGid);
	setFourMomentum(CLHEP::HepLorentzVector(0, 0, 0, 0));
	setFourMomentumBosonRestFrame(CLHEP::HepLorentzVector(0, 0, 0, 0));

}
Particle::Particle(int PDGid, CLHEP::HepLorentzVector momentum) {
	setPDGid(PDGid);
	setFourMomentum(momentum);
	setFourMomentumBosonRestFrame(CLHEP::HepLorentzVector(0, 0, 0, 0));

}
Particle::Particle(int PDGid, vector<VINHAC::Particle*> parents) {
	setPDGid(PDGid);
	_parents = parents;
	setFourMomentum(CLHEP::HepLorentzVector(0, 0, 0, 0));
	setFourMomentumBosonRestFrame(CLHEP::HepLorentzVector(0, 0, 0, 0));

}
Particle::Particle(int PDGid, vector<VINHAC::Particle*> parents,
		CLHEP::HepLorentzVector momentum) {
	setPDGid(PDGid);
	_parents = parents;
	setFourMomentum(momentum);
	setFourMomentumBosonRestFrame(CLHEP::HepLorentzVector(0, 0, 0, 0));
}

void Particle::clear() {
	setFourMomentum(CLHEP::HepLorentzVector(0, 0, 0, 0));
	setFourMomentumBosonRestFrame(CLHEP::HepLorentzVector(0, 0, 0, 0));
	setPDGid(0);
	_parents.clear();
	_children.clear();
}
