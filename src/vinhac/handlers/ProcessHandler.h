using namespace std;
#include <string>
#include "../core/Event.h"

#ifndef __ProcessHandler_h__
#define __ProcessHandler_h__

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
	inline const string& getHandlerName() const {
		return handlerName;
	}
	inline void setHandlerName(const string& name) {
		handlerName = name;
	}
};
}

#endif
