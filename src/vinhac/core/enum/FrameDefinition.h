/*
 * FrameDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef FRAMEDEFINITION_H_
#define FRAMEDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct FrameName {
	enum name {
		bosonRestFrame,
		LAB,
		bosonAndSecondBeamNucleonWithZAxisAlongW,
		bosonRestFrameZAxisAsLABXAxisAccordingToBeamWithSmallerPolarAngle,
		fixedPT,
		randomPT
	};

	static name mapStringToEnum(const std::string& target){
		if(target=="Boson Rest Frame"){
			return bosonRestFrame;
		} else if(target=="LAB"){
			return LAB;
		} else if(target=="W and 2nd beam nucleon CMS with +z axis along W"){
			return bosonAndSecondBeamNucleonWithZAxisAlongW;
		} else if(target=="W rest frame with +z axis along W direction in LAB and +x axis oriented according to beam with smaller polar angle in this frame"){
			return bosonRestFrameZAxisAsLABXAxisAccordingToBeamWithSmallerPolarAngle;
		} else if(target=="Fixed PT"){
			return fixedPT;
		} else if(target=="Random PT"){
			return randomPT;
		}
		throw VinhacException("Unknown FrameName name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){
			case bosonRestFrame:			return "Boson Rest Frame";
			case LAB:						return "LAB";
			case bosonAndSecondBeamNucleonWithZAxisAlongW:
											return "W and 2nd beam nucleon CMS with +z axis along W";
			case bosonRestFrameZAxisAsLABXAxisAccordingToBeamWithSmallerPolarAngle:
											return "W rest frame with +z axis along W direction in LAB and +x axis oriented according to beam with smaller polar angle in this frame";
			case fixedPT:					return "Fixed PT";
			case randomPT:					return "Random PT";
		}
		throw VinhacException("Unknown FrameName name");
	}
};

}

#endif /* FRAMEDEFINITION_H_ */
