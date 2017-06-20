/*
 * PolarizationDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef POLARIZATIONDEFINITION_H_
#define POLARIZATIONDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct PolarizationName {
	enum name {
		unpolarized,
		longitudinally,
		transversely,
		left,
		right,
		int_long_tran,
		int_left_right
	};

	static name mapStringToEnum(const std::string& target){
		if(target=="unpolarized"){
			return unpolarized;
		} else if(target=="longitudinally"){
			return longitudinally;
		} else if(target=="transversely"){
			return transversely;
		} else if(target=="left"){
			return left;
		} else if(target=="right"){
			return right;
		} else if(target=="int_left_right"){
			return int_left_right;
		} else if(target=="int_long_tran"){
			return int_long_tran;
		}
		throw VinhacException("Unknown polarization name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){
			case unpolarized:						return "unpolarized";
			case longitudinally:					return "longitudinally";
			case transversely:						return "transversely";
			case left:								return "left";
			case right:								return "right";
			case int_left_right:					return "int_left_right";
			case int_long_tran:						return "int_long_tran";

		}
		throw VinhacException("Unknown polarization name");
	}
};

}

#endif /* POLARIZATIONDEFINITION_H_ */
