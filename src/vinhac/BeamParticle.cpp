/*
 * \file BeamParticle.cpp
 * \brief
 * \date Created on: 2009-06-04
 * \author     Author: siudmy
 */

#include "BeamParticle.h"

using namespace VINHAC;

BeamParticle::BeamParticle() :
	pdf(0) {

}

BeamParticle::BeamParticle(int PDGid, PDF *pdf) :
	Particle(PDGid) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamParticle::BeamParticle() 1"<<std::endl;

#endif

	setPDF(pdf);
}
BeamParticle::BeamParticle(int PDGid, CLHEP::HepLorentzVector momentum,
		PDF *pdf) :
	Particle(PDGid, momentum) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamParticle::BeamParticle() 2"<<std::endl;

#endif

	setPDF(pdf);
}

BeamParticle::~BeamParticle() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: BeamParticle::~BeamParticle()"<<std::endl;

#endif

	if (pdf) {
		delete pdf;
		pdf = 0;
	}
}

