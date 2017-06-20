/*
 * \file ppWintegrand.cpp
 * \brief
 * \date Created on: 2009-06-11
 * \author     Author: siudmy
 */

#include "ppWintegrand.h"
#include "../utils/FortranFunction.h"

VINHAC::ppWintegrand::ppWintegrand(DataSource *ds, BeamHandler *beamHandler) :
	ds(ds), beamHandler(beamHandler), isWplus(false), isWminus(false) {
	x_min = atof(ds->beamA["xMin"].c_str());
	x_max = atof(ds->beamA["xMax"].c_str());
	q2minCut = atof(ds->beamA["Q2min"].c_str());

	xminimal = log(x_max / x_min);

	for (vector<int>::iterator it = ds->initialParticles.begin(); it
			!= ds->initialParticles.end(); ++it) {
		if ((*it) > 0) {
			if ((*it) % 2 == 1) {
				dQuarks.push_back((*it));
			} else {
				uQuarks.push_back((*it));
			}

		} else {
			if ((-(*it)) % 2 == 1) {
				dAntiQuarks.push_back((*it));
			} else {
				uAntiQuarks.push_back((*it));
			}

		}
	}

	for (vector<int>::iterator it = ds->intermediateParticles.begin(); it
			!= ds->intermediateParticles.end(); ++it) {
		if ((*it) == 24)
			isWplus = true;
		if ((*it) == -24)
			isWminus = true;

	}

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: ppWintegrand::ppWintegrand() # xmin xmax q2min "<<x_min<<" "<<x_max<<" "<<q2minCut<<std::endl;
	std::cout<<"DEBUG: ppWintegrand::ppWintegrand() # isW+ isW- "<<isWplus<<" "<<isWminus<<" "<<std::endl;
#endif
}

void VINHAC::ppWintegrand::setS(double energy) {
	s = energy;
	a = ds->MW * ds->MW / s;
	b = ds->MW * ds->GW / s;
	tau_min = max(q2minCut / s, x_min * x_min);
	tau_max = min(1.0, x_max * x_max);
	deltau = tau_max - tau_min;
}

double VINHAC::ppWintegrand::Density(int nDim, double *Xarg) {
#ifdef DEBUG_FLAG

	//std::cout<<"DEBUG: ppWintegrand::Density()"<<std::endl;

#endif
	//return fotranDensity(Xarg);
	double FunDis = 0.0;

	// Auxiliary variables
	double eta = Xarg[0];
	double zet = Xarg[1];
	double x = x_min * exp(xminimal * zet);
	double alf_min = atan((x_min * x - a) / b);
	double alf_max = atan((x_max * x - a) / b);
	double delalf = alf_max - alf_min;
	double alf = alf_min + delalf * eta;
	double tau = a + b * tan(alf);

	if (tau < tau_min)
		return FunDis;

	double x1 = min(1.0, x);
	if (x1 < x_min || x1 > x_max)
		return FunDis;
	double x2 = min(1.0, tau / x1);
	if (x2 < x_min || x2 > x_max)
		return FunDis;

	double Q2;
	if (ds->keyFac == 0) {
		Q2 = tau * s;
	} else {
		Q2 = ds->MW * ds->MW;
	}

	double Q = sqrt(Q2);

	std::vector<double> x1pdf = p_beamParticleA->getPDF()->xfx(x1, Q);
	std::vector<double> x2pdf = p_beamParticleB->getPDF()->xfx(x2, Q);

	beamHandler->setXPDFa(x1pdf);
	beamHandler->setXPDFb(x2pdf);

	if (ds->keyCol == 0) {

		if (isWminus) {
			for (vector<int>::iterator it_d = dQuarks.begin(); it_d
					!= dQuarks.end(); ++it_d) {
				for (vector<int>::iterator it_u = uAntiQuarks.begin(); it_u
						!= uAntiQuarks.end(); ++it_u) {

					FunDis += x1pdf[(*it_d) + 6] * x2pdf[(*it_u) + 6]
							* ds->ckm2[(-(*it_u)) / 2 - 1][(*it_d) / 2];
				}

			}
		}

		if (isWplus) {
			for (vector<int>::iterator it_d = dAntiQuarks.begin(); it_d
					!= dAntiQuarks.end(); ++it_d) {
				for (vector<int>::iterator it_u = uQuarks.begin(); it_u
						!= uQuarks.end(); ++it_u) {

					FunDis += x1pdf[(*it_u) + 6] * x2pdf[(*it_d) + 6]
							* ds->ckm2[(*it_u) / 2 - 1][(-(*it_d)) / 2];
				}

			}
		}
		FunDis *= 2.0;

	} else {
		throw VinhacException(
				"ppWintegrand::Density # not implemented beam particles");
	}

	// Parton-level crude x-section divided by the CKM matrix element
	FunDis *= beamHandler->totalXsCrude(ds->ckm[0][0], Q2) / ds->ckm2[0][0];

	// Jacobian factor
	FunDis *= xminimal * delalf * ((Q2 - ds->MW * ds->MW) * (Q2 - ds->MW
			* ds->MW) + ds->MW * ds->GW * ds->MW * ds->GW) / (ds->MW * ds->GW)
			/ Q2;

	if (FunDis < 0.0) {
		FunDis = 0.0;
		//todo dorobic printout
	}
#ifdef DEBUG_FORTRAN_INTEGRAND
	double xarg2[2];
	//xarg2[0]=Xarg[1];
	//xarg2[1]=Xarg[0];
	double test = wh_funppwtest_(Xarg);
	//if(FunDis>1E-4)
	//std::cout<<"DEBUG: ppWintegrand::Density()"<<test/FunDis<<std::endl;

#endif

	return FunDis;

}


