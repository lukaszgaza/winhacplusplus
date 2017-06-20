#include <string>
#include <vector>
#include <exception>
using namespace std;

#ifndef __WINHAC_h__
#define __WINHAC_h__

#include "Manager.h"

namespace VINHAC {
class Manager;
class WINHAC;
}

namespace VINHAC {
/**
 * \brief Class of WINHAC Monte Carlo generator. It is used to model charged Drell-Yan process type.
 */
class WINHAC: public VINHAC::Manager {

private:
	void PrintLogo();

public:
	/**
	 * A WINHAC constructor.
	 * A more elaborate description of the constructor.
	 */
	WINHAC();
	~WINHAC();
};
}

#endif
