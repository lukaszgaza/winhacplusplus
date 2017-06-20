/*
 * SeedTable2.h
 *
 *  Created on: 2008-11-23
 *      Author: kamil
 */

#ifndef SEEDTABLE2_H_
#define SEEDTABLE2_H_


//#include "CLHEP/Random/defs.h"

namespace CLHEP {

class SeedTable2{
public:

static void getTheTableSeeds(long* seeds, int index);

static const long seedTable[215][2];
};




}  // namespace CLHEP

#ifdef ENABLE_BACKWARDS_COMPATIBILITY
//  backwards compatibility will be enabled ONLY in CLHEP 1.9
using namespace CLHEP;
#endif




#endif /* SEEDTABLE2_H_ */
