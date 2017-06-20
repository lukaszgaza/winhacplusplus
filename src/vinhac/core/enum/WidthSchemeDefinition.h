/*
 * FrameDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef WIDTHSCHEMEDEFINITION_H_
#define WIDTHSCHEMEDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct WidthScheme {
	enum name {
		fixed = 0,
		running = 1
	};

	static name mapStringToEnum(const std::string& target){
		if(target=="fixed"){
			return fixed;
		} else if(target=="running"){
			return running;
		}
		throw VinhacException("Unknown width scheme name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){
			case fixed:							return "fixed";
			case running:						return "running";
		}
		throw VinhacException("Unknown width scheme name");
	}
};

}

#endif /* WIDTHSCHEMEDEFINITION_H_ */
