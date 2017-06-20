/*
 * FrameTransformationHandler.cpp
 *
 *  Created on: 26-08-2011
 *      Author: sobol
 */

#include "FrameTransformationHandler.h"
#include "../input/DataSource.h"
#include "../util/Transformation.h"

namespace VINHAC {

FrameTransformationHandler::FrameTransformationHandler(TRND *randomEngine):
	randomEngine(randomEngine){

}

FrameTransformationHandler::~FrameTransformationHandler() {
}

bool FrameTransformationHandler::eventEvolve(VINHAC::Event& event) {

	switch (DataSource::get().modelReferenceFrame) {
	case FrameName::randomPT:
	case FrameName::fixedPT:
	case FrameName::bosonRestFrameZAxisAsLABXAxisAccordingToBeamWithSmallerPolarAngle:
		copyLabMomentum(event.getBeamParticleA(),
				DataSource::get().modelReferenceFrame);
	case FrameName::bosonAndSecondBeamNucleonWithZAxisAlongW:
		copyLabMomentum(event.getBeamParticleB(),
				DataSource::get().modelReferenceFrame);
		copyLabMomentum(event.getQuark(), DataSource::get().modelReferenceFrame);
		copyLabMomentum(event.getAntiQuark(),
				DataSource::get().modelReferenceFrame);
		copyLabMomentum(event.getBoson(), DataSource::get().modelReferenceFrame);
		copyLabMomentum(event.getLepton(),
				DataSource::get().modelReferenceFrame);
		copyLabMomentum(event.getNeutrino(),
				DataSource::get().modelReferenceFrame);
		break;
	default:
		break;
	}

	switch (DataSource::get().modelReferenceFrame) {
	case FrameName::bosonAndSecondBeamNucleonWithZAxisAlongW:
		Transformation::transformationToCMSofWNucleonSystem(
				event.getBeamParticleB().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getQuark().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getAntiQuark().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getBoson().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getLepton().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getNeutrino().getFourMomentum(
						DataSource::get().modelReferenceFrame));
		break;
	case FrameName::bosonRestFrameZAxisAsLABXAxisAccordingToBeamWithSmallerPolarAngle:
		Transformation::transformationToBosonRestFrameZAxisAsLABXAxisAccordingToBeamWithSmallerPolarAngle(
				event.getBeamParticleA().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getBeamParticleB().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getQuark().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getAntiQuark().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getBoson().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getLepton().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getNeutrino().getFourMomentum(
						DataSource::get().modelReferenceFrame));
		break;
	case FrameName::fixedPT:
		Transformation::generateFixedPT(
				event.getBeamParticleA().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getBeamParticleB().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getQuark().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getAntiQuark().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getBoson().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getLepton().getFourMomentum(
						DataSource::get().modelReferenceFrame),
				event.getNeutrino().getFourMomentum(
						DataSource::get().modelReferenceFrame));
		break;
	case FrameName::randomPT:
			Transformation::generateRandomPT(
					event.getBeamParticleA().getFourMomentum(
							DataSource::get().modelReferenceFrame),
					event.getBeamParticleB().getFourMomentum(
							DataSource::get().modelReferenceFrame),
					event.getQuark().getFourMomentum(
							DataSource::get().modelReferenceFrame),
					event.getAntiQuark().getFourMomentum(
							DataSource::get().modelReferenceFrame),
					event.getBoson().getFourMomentum(
							DataSource::get().modelReferenceFrame),
					event.getLepton().getFourMomentum(
							DataSource::get().modelReferenceFrame),
					event.getNeutrino().getFourMomentum(
							DataSource::get().modelReferenceFrame),
							randomEngine);
			break;
	default:
		break;
	}

//	cout<<endl;
//	cout<<event.printInLabFrame()<<endl;
//	cout<<event.printInFrame(FrameName::randomPT)<<endl;

	return true;
}

void FrameTransformationHandler::copyLabMomentum(Particle& particle,
		FrameName::name name) {
	CLHEP::HepLorentzVector& target = particle.getFourMomentum(name);
	CLHEP::HepLorentzVector& source = particle.getFourMomentum();
	target = source;
}

}
