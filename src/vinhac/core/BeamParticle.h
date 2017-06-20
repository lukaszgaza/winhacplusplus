
#ifndef _BEAMPARTICLE_H_
#define _BEAMPARTICLE_H_
#include <vector>
#include "Particle.h"
#include "../pdf/PDF.h"
using namespace std;

namespace VINHAC {

class BeamParticle;
class Particle;
class PDF;
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
	BeamParticle(const int& PDGid, PDF* pdf);

	/**
	* \brief Constructor of BeamParticle
	*
	* @see PDF
	* @see CLHEP::HepLorentzVector
	* @param PDGid PDG id of particle
	* @param momentum momentum of particle
	* @param pdf pointer to PDF object for this particle
	*/
	BeamParticle(const int& PDGid,const  CLHEP::HepLorentzVector& momentum, PDF* pdf);

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
