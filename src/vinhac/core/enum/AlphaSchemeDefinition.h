/*
 * FrameDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef ALPHASCHEMEDEFINITION_H_
#define ALPHASCHEMEDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct AlphaScheme {
	enum name {
		AlphaQED = 0,
		AlphaGmu = 1
	};

	static name mapStringToEnum(const std::string& target){
		if(target=="Alpha QED"){
			return AlphaQED;
		} else if(target=="Alpha Gmu"){
			return AlphaGmu;
		}
		throw VinhacException("Unknown alpha scheme name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){
			case AlphaQED:						return "Alpha QED";
			case AlphaGmu:						return "Alpha Gmu";
		}
		throw VinhacException("Unknown alpha scheme name");
	}
};

}

#endif /* ALPHASCHEMEDEFINITION_H_ */
