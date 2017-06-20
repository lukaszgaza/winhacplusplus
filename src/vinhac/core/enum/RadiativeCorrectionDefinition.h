/*
 * FrameDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef RADIATIVECORRECTIONDEFINITION_H_
#define RADIATIVECORRECTIONDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct RadiativeCorrection {
	enum name {
		Born = 0,
		AboveQEDradiativeInWDecay = 1,
		AboveOAlphaEWCorrectionsInWDecay = 2
	};

	static name mapStringToEnum(const std::string& target){
		if(target=="Born"){
			return Born;
		} else if(target=="Above + QED radiative correction in W-decay"){
			return AboveQEDradiativeInWDecay;
		} else if(target=="Above + O(alpha) EW corrections in W-decay"){
			return AboveOAlphaEWCorrectionsInWDecay;
		}
		throw VinhacException("Unknown Radiative Correction name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){
			case Born:								return "Born";
			case AboveQEDradiativeInWDecay:			return "Above + QED radiative correction in W-decay";
			case AboveOAlphaEWCorrectionsInWDecay:	return "Above + O(alpha) EW corrections in W-decay";
		}
		throw VinhacException("Unknown Radiative Correction name");
	}
};

}

#endif /* RADIATIVECORRECTIONDEFINITION_H_ */
