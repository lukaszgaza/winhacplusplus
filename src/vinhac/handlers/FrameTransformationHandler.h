/*
 * FrameTransformationHandler.h
 *
 *  Created on: 26-08-2011
 *      Author: sobol
 */

#ifndef FRAMETRANSFORMATIONHANDLER_H_
#define FRAMETRANSFORMATIONHANDLER_H_

#include "ProcessHandler.h"
#include "../../utils/FOAM/TRND.h"

namespace VINHAC {

class FrameTransformationHandler: public VINHAC::ProcessHandler {
private:
	TRND *randomEngine;
	void copyLabMomentum(Particle& particle,FrameName::name name);
public:
	FrameTransformationHandler(TRND *randomEngine);
	virtual ~FrameTransformationHandler();

	/**
	 * \brief evolves the event
	 *
	 * It does whole part of generation process which is to do in this handler
	 * @see Event
	 * @param event reference to current Event object
	 * @return if weight was zero
	 */
	bool eventEvolve(VINHAC::Event& event);
};

}

#endif /* FRAMETRANSFORMATIONHANDLER_H_ */
