
#ifndef _BEAMPARTICLE_H_
#define _BEAMPARTICLE_H_
#include "Particle.h"
#include "../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include "PDF.h"
#include <vector>

using namespace std;

namespace VINHAC {

class BeamParticle;
class Particle;
}

namespace VINHAC {

//!  Class representing beam particle in BeamHandler
/*!
 * @see BeamHandler
*/
class BeamParticle: public Particle {
private:
	PDF *pdf; //!< Parton density function of the beam particle;
public:
	/**
	* \brief Constructor of BeamParticle
	*/
	BeamParticle();

	/**
	* \brief Destructor of BeamParticle
	*/
	~BeamParticle();

	/**
	* \brief Constructor of BeamParticle
	*
	* @see PDF
	* @param PDGid PDG id of particle
	* @param pdf pointer to PDF object for this particle
	*/
	BeamParticle(int PDGid, PDF* pdf);

	/**
	* \brief Constructor of BeamParticle
	*
	* @see PDF
	* @see CLHEP::HepLorentzVector
	* @param PDGid PDG id of particle
	* @param momentum momentum of particle
	* @param pdf pointer to PDF object for this particle
	*/
	BeamParticle(int PDGid, CLHEP::HepLorentzVector momentum, PDF* pdf);

	/**
	* \brief sets PDF object for this particle
	*
	* @see PDF
	* @param pdf pointer to PDF object for this particle
	*/
	inline void setPDF(PDF* pdf) {
		this->pdf = pdf;
	}

	/**
	* \brief returns PDF object for this particle
	*
	* @see PDF
	* @return pdf pointer to PDF object for this particle
	*/
	inline PDF* getPDF() {
		return pdf;
	}

};
}
#endif /* BEAMPARTICLE_H_ */
