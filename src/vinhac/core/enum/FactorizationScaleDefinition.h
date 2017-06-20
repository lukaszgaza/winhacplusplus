/*
 * FrameDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef FACTORIZATIONSCALEDEFINITION_H_
#define FACTORIZATIONSCALEDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct FactorizationScale {
	enum name {
		invariantMass = 0,
		fixedWmass = 1
	};

	static name mapStringToEnum(const std::string& target){
		if(target=="Invariant mass"){
			return invariantMass;
		} else if(target=="Fixed W mass"){
			return fixedWmass;
		}
		throw VinhacException("Unknown factorization scale name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){
			case invariantMass:							return "Invariant mass";
			case fixedWmass:							return "Fixed W mass";
		}
		throw VinhacException("Unknown factorization scale name");
	}
};

}

#endif /* FACTORIZATIONSCALEDEFINITION_H_ */
