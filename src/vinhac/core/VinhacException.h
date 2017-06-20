/*
 * VinhacException.h
 *
 *  Created on: 2009-11-11
 *      Author: kamil
 */

#ifndef VINHACEXCEPTION_H_
#define VINHACEXCEPTION_H_

#include <string>

namespace VINHAC {

/**
 * \brief Simple class for exception in VINHAC type generator.
 * It has a message inside.
 */
class VinhacException {
private:
	std::string message;
public:
	VinhacException();

	/**
	 * \brief creates a message containing some text provided by parameter.
	 * This parameter should be convertible to string.
	 *
	 * @param message text provided as message
	 */
	template<typename T>
	VinhacException(T message) :
		message(message) {
	}
	;
	virtual ~VinhacException();

	/**
	 * @returns message
	 */
	std::string getMessage();
};

}

#endif /* VINHACEXCEPTION_H_ */
