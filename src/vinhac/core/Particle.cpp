#include <string>
#include <vector>
#include <exception>
using namespace std;

#include "Particle.h"

using namespace VINHAC;

Particle::Particle() {
}

Particle::Particle(const Particle& target) :
	_PDGid(target._PDGid), m_tag(target.m_tag) {
	if (target.fourMomenta.find(FrameName::LAB) != target.fourMomenta.end()) {
		fourMomenta[FrameName::LAB]
				= target.fourMomenta.find(FrameName::LAB)->second;
	}
	if (target.fourMomenta.find(FrameName::bosonRestFrame)
			!= target.fourMomenta.end()) {
		fourMomenta[FrameName::bosonRestFrame] = target.fourMomenta.find(
				FrameName::bosonRestFrame)->second;
	}
	for (unsigned i = 0; i < target._parents.size(); ++i) {
		_parents.push_back(target._parents.at(i));
	}
	for (unsigned i = 0; i < target._children.size(); ++i) {
		_children.push_back(target._children.at(i));
	}

}

Particle::Particle(const int& PDGid) :
	_PDGid(PDGid) {
}

Particle::Particle(const int& PDGid, const CLHEP::HepLorentzVector& momentum) :
	_PDGid(PDGid) {
	fourMomenta[FrameName::LAB] = momentum;
}

Particle::Particle(const int& PDGid, const vector<VINHAC::Particle*>& parents) :
	_PDGid(PDGid) {
	_parents = parents;

}
Particle::Particle(const int& PDGid, const vector<VINHAC::Particle*>& parents,
		const CLHEP::HepLorentzVector& momentum) :
	_PDGid(PDGid) {
	fourMomenta[FrameName::LAB] = momentum;
	_parents = parents;
}

void Particle::clear() {
	getFourMomentum().set(0);
	getFourMomentumBosonRestFrame().set(0);
	_PDGid = 0;
	_parents.clear();
	_children.clear();
}
