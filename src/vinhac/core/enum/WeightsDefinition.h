/*
 * WeightsDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef WEIGHTSDEFINITION_H_
#define WEIGHTSDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct WeightName {
	enum name {
		XsCrudWt = 1,
		PartonsGeneration = 2,
		HardProcessWeight = 3,


		beta00 = 101,
		beta01_QED= 11,
		beta11 = 12,
		beta01_EW = 13,
		beta01_LL = 14,
		beta11_LL = 15,
		beta01_QED_plus_beta11 = 102,
		beta01_EW_plus_beta11 = 103,
		beta01_LL_plus_beta11_LL = 104,

		xs00 = 201,
		xs0gQED = 202,
		xs0gEW = 203,
		xs1g = 206,
		xs0gQED_plus_xs1g = 204,
		xs0gEW_plus_xs1g = 205
	};


	static name mapStringToEnum(const std::string& target){
		if(target=="XsCrudWt"){
			return XsCrudWt;
		} else if(target=="PartonsGeneration"){
			return PartonsGeneration;
		} else if(target=="HardProcessWeight"){
			return HardProcessWeight;
		} else if(target=="YFS O(alf0)"){
			return beta00;
		} else if(target=="beta01_QED"){
			return beta01_QED;
		} else if(target=="beta11"){
			return beta11;
		} else if(target=="beta01_EW"){
			return beta01_EW;
		} else if(target=="beta01_LL"){
			return beta01_LL;
		} else if(target=="beta11_LL"){
			return beta11_LL;
		} else if(target=="YFS O(alf1)_QED"){
			return beta01_QED_plus_beta11;
		} else if(target=="YFS O(alf1)_EW"){
			return beta01_EW_plus_beta11;
		} else if(target=="YFS O(alf1)_LL"){
			return beta01_LL_plus_beta11_LL;
		} else if(target=="Born"){
			return xs00;
		} else if(target=="O(alpha)_QED"){
			return xs0gQED_plus_xs1g;
		} else if(target=="O(alpha)_EW"){
			return xs0gEW_plus_xs1g;
		} else if(target=="Born + virtual_QED + soft"){
			return xs0gQED;
		} else if(target=="Born + virtual_EW + soft"){
			return xs0gEW;
		} else if(target=="One real hard photon"){
			return xs1g;
		}

		throw VinhacException("Unknown Weight name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){

		case XsCrudWt:							return "XsCrudWt";
		case PartonsGeneration:					return "PartonsGeneration";
		case HardProcessWeight:					return "HardProcessWeight";

		case beta00: 							return "YFS O(alf0)";
		case beta01_QED:						return "beta01_QED";
		case beta11:							return "beta11";
		case beta01_EW:							return "beta01_EW";
		case beta01_LL:							return "beta01_LL";
		case beta11_LL:							return "beta11_LL";
		case beta01_QED_plus_beta11:			return "YFS O(alf1)_QED";
		case beta01_EW_plus_beta11:				return "YFS O(alf1)_EW";
		case beta01_LL_plus_beta11_LL:  		return "YFS O(alf1)_LL";


		case xs00:								return "Born";
		case xs0gQED_plus_xs1g:					return "O(alpha)_QED";
		case xs0gEW_plus_xs1g:					return "O(alpha)_EW";
		case xs0gQED:							return "Born + virtual_QED + soft";
		case xs0gEW:							return "Born + virtual_EW + soft";
		case xs1g:								return "One real hard photon";
		}
		throw VinhacException("Unknown Weight name");
	}
};

}

#endif /* WEIGHTSDEFINITION_H_ */
