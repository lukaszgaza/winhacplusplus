#include <string>
#include <vector>
#include <exception>
#include <iostream>
using namespace std;

#ifndef __ProcessHandler_h__
#define __ProcessHandler_h__
#include "Event.h"

namespace VINHAC {

/**
 * \brief Base class for object doing something with event.
 *
 * Every extending class has to implement eventEvolve method which does anything with Event object
 */
class ProcessHandler {

private:
	string handlerName;

public:

	/**
	* \brief evolves the event
	*
	* It does whole part of generation process which is to do in this handler
	* @see Event
	* @param event reference to current Event object
	* @return if weight was zero
	*/
	virtual bool eventEvolve(VINHAC::Event& event) = 0;

	virtual ~ProcessHandler() {
	}
	;
	inline string getHandlerName() {
		return handlerName;
	}
	inline void setHandlerName(string name) {
		handlerName = name;
	}
};
}

#endif
