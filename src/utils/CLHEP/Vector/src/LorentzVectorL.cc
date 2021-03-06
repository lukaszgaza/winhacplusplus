// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is part of the implementation of the HepLorentzVector class:
// Those methods which might, if coded in other modules, force loading
// of the LorentzRotation.cc code module.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

//#include "CLHEP/Vector/defs.h"
#include "../Vector/LorentzVector.h"
#include "../Vector/LorentzRotation.h"

namespace CLHEP  {

HepLorentzVector &
HepLorentzVector::operator *= (const HepLorentzRotation & m) {
  return *this = m.vectorMultiplication(*this);
}

HepLorentzVector &
HepLorentzVector::transform(const HepLorentzRotation & m){
  return *this = m.vectorMultiplication(*this);
}

}  // namespace CLHEP
