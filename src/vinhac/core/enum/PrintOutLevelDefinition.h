/*
 * FrameDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef PRINTOUTLEVELDEFINITION_H_
#define PRINTOUTLEVELDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct PrintOutLevel {
	enum name {
		None = 0,
		Basic = 1,
		Detailed = 2
	};

	static name mapStringToEnum(const std::string& target){
		if(target=="None"){
			return None;
		} else if(target=="Basic"){
			return Basic;
		} else if(target=="Detailed"){
			return Detailed;
		}
		throw VinhacException("Unknown PrintOutLevel name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){
			case None:							return "None";
			case Basic:							return "Basic";
			case Detailed:						return "Detailed";
		}
		throw VinhacException("Unknown PrintOutLevel name");
	}
};

}

#endif /* PRINTOUTLEVELDEFINITION_H_ */
