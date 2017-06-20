/*
 * \file ppWintegrand.h
 * \brief
 * \date Created on: 2009-06-11
 * \author     Author: siudmy
 */

#ifndef PPWINTEGRAND_H_
#define PPWINTEGRAND_H_
#include <iostream>
#include "BeamHandler.h"
#include "../utils/FOAM/TFOAM_INTEGRAND.h"
#include "DataSource.h"
#include "BeamParticle.h"
#include "VinhacException.h"

namespace VINHAC {
class DataSource;
class BeamHandler;
}

namespace VINHAC {

/**
 * \brief Class used by FOAM to generate Bjorken's xs.
 *
 * @see Density
 */
class ppWintegrand: public TFOAM_INTEGRAND {
private:

	DataSource *ds; //!< pointer to object with input parameters
	BeamHandler *beamHandler;

	double s;
	double x_max; //!< Max value of x for PDF
	double x_min; //!< Min value of x for PDF
	double q2minCut;
	double xminimal;
	double a;
	double b;
	double tau_min;
	double tau_max;
	double deltau;

	BeamParticle* p_beamParticleA;
	BeamParticle* p_beamParticleB;

	vector<int> dQuarks;
	vector<int> dAntiQuarks;
	vector<int> uQuarks;
	vector<int> uAntiQuarks;

	bool isWplus;
	bool isWminus;

public:
	/**
	 * \brief constructor
	 *
	 * @param ds pointer to object with input data
	 * @param beamHandler pointer to BeamHandler object
	 */
	ppWintegrand(DataSource *ds, BeamHandler *beamHandler);
	~ppWintegrand() {
	}

	/*!
	 * \brief Two-dimensional function used for generation of quarks Bjorken
	 * x's: x1 of quark and x2 of antiquark.
	 * This is a "density" function used by the general purpose MC
	 * sampler Foam to generate variables \f$(\eta,\zeta)\f$, where \f$0<\eta,\zeta<1\f$.
	 * \param Xarg - array of arguments from the range (0,1);
	 *         Xarg[0] = eta, Xarg[1] = zet
	 * \param ndim dimesion of Xarg array
	 *
	 * The following mapping of integration variables has been applied:
	 *       \li \f$x_1 = x_{min}  (\frac{x_{max}}{x_{min}})^\zeta\f$,
	 *       \li \f$x_2 = \frac{\tau}{x1}\f$,\n
	 *  where:
	 *  	 - \f$tau  = a + b\tan(\alpha)\f$,\n
	 *  and:
	 *       - \f$\alpha = alpha_{min} + (\alpha_{max} - \alpha_{min})\eta\f$,
	 *       - \f$\alpha_{min} = \arctan( (\frac{x_{min} x_1 - a}{b} )\f$,
	 *       - \f$\alpha_{max} = \arctan( (\frac{x_{max} x_1 - a}{b} )\f$,\n
	 *  with:
	 *  	 - \f$a = \frac{m_W^2}{s}\f$,
	 *  	 - \f$b = \frac{m_W \Gamma_W}{s}\f$,
	 *  	 - \f$s = E_{cm}^2\f$.
	 *       - \f$y = \frac{1 + \cos(\theta)}{2}\f$.
	 *
	 * @return function value
	 */
	double Density(int ndim, double* Xarg);

	//setters

	/**
	 * \brief sets Mandelstam variable s
	 *
	 * @param s to set
	 */
	void setS(double s);

	/**
	 * \brief sets pointer to particle from beam A
	 *
	 * @param p_beamParticle pointer to particle from beam A
	 */
	inline void setBeamParticleA(BeamParticle* p_beamParticle) {
		p_beamParticleA = p_beamParticle;
	}

	/**
	 * \brief sets pointer to particle from beam B
	 *
	 * @param p_beamParticle pointer to particle from beam B
	 */
	inline void setBeamParticleB(BeamParticle* p_beamParticle) {
		p_beamParticleB = p_beamParticle;
	}

	/**
	 * @return minimal x
	 */
	inline double getXminimal() {
		return xminimal;
	}

	/**
	 * @return \f$x_{min}\f$
	 */
	inline double getX_min() {
		return x_min;
	}

	/**
	 * @return \f$x_{max}\f$
	 */
	inline double getX_max() {
		return x_max;
	}

	/**
	 * @return a
	 */
	inline double getA() {
		return a;
	}

	/**
	 * @return b
	 */
	inline double getB() {
		return b;
	}

};
}
#endif /* PPWINTEGRAND_H_ */
