/*
 * \file ppWintegrand.cpp
 * \brief
 * \date Created on: 2009-06-11
 * \author     Author: siudmy
 */

#include "ppWintegrand.h"
#include "../../utils/FortranFunction.h"
#include "../input/DataSource.h"
#include "BeamHandler.h"
#include "../core/BeamParticle.h"
#include "../core/VinhacException.h"
#include <cstdlib>
#include <cmath>

using namespace std;

VINHAC::ppWintegrand::ppWintegrand(BeamHandler *beamHandler) :
	beamHandler(beamHandler), isWplus(false), isWminus(false) {
	x_min = DataSource::get().pdfXMin;
	x_max = DataSource::get().pdfXMax;
	q2minCut = DataSource::get().pdfQ2min;

	xminimal = log(x_max / x_min);

	for (vector<int>::iterator it = DataSource::get().initialParticles.begin(); it
			!= DataSource::get().initialParticles.end(); ++it) {
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

	for (vector<int>::iterator it = DataSource::get().intermediateParticles.begin(); it
			!= DataSource::get().intermediateParticles.end(); ++it) {
		if ((*it) == 24)
			isWplus = true;
		if ((*it) == -24)
			isWminus = true;

	}

#ifdef DEBUG_FLAG

	cout<<"DEBUG: ppWintegrand::ppWintegrand() # xmin xmax q2min "<<x_min<<" "<<x_max<<" "<<q2minCut<<endl;
	cout<<"DEBUG: ppWintegrand::ppWintegrand() # isW+ isW- "<<isWplus<<" "<<isWminus<<" "<<endl;
#endif
}

void VINHAC::ppWintegrand::setS(const double& energy) {
	s = energy;
	a = DataSource::get().MW * DataSource::get().MW / s;
	b = DataSource::get().MW * DataSource::get().GW / s;
	tau_min = max(q2minCut / s, x_min * x_min);
	tau_max = min(1.0, x_max * x_max);
	deltau = tau_max - tau_min;
}

double VINHAC::ppWintegrand::Density(int nDim, double *Xarg) {
#ifdef DEBUG_FLAG

	//cout<<"DEBUG: ppWintegrand::Density()"<<endl;

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
	if (DataSource::get().keyFac == 0) {
		Q2 = tau * s;
	} else {
		Q2 = DataSource::get().MW * DataSource::get().MW;
	}

	double Q = sqrt(Q2);

	vector<double> x1pdf = p_beamParticleA->getPDF()->xfx(x1, Q);
	vector<double> x2pdf = p_beamParticleB->getPDF()->xfx(x2, Q);

	beamHandler->setScalePDF(Q);
	beamHandler->setXPDFa(x1pdf);
	beamHandler->setXPDFb(x2pdf);

	if (DataSource::get().keyCol == 0) {

		if (isWminus) {
			for (vector<int>::iterator it_d = dQuarks.begin(); it_d
					!= dQuarks.end(); ++it_d) {
				for (vector<int>::iterator it_u = uAntiQuarks.begin(); it_u
						!= uAntiQuarks.end(); ++it_u) {

					FunDis += x1pdf[(*it_d) + 6] * x2pdf[(*it_u) + 6]
							* DataSource::get().ckm2[(-(*it_u)) / 2 - 1][(*it_d) / 2];
				}

			}
		}

		if (isWplus) {
			for (vector<int>::iterator it_d = dAntiQuarks.begin(); it_d
					!= dAntiQuarks.end(); ++it_d) {
				for (vector<int>::iterator it_u = uQuarks.begin(); it_u
						!= uQuarks.end(); ++it_u) {

					FunDis += x1pdf[(*it_u) + 6] * x2pdf[(*it_d) + 6]
							* DataSource::get().ckm2[(*it_u) / 2 - 1][(-(*it_d)) / 2];
				}

			}
		}
		FunDis *= 2.0;

	} else {
		throw VinhacException(
				"ppWintegrand::Density # not implemented beam particles");
	}

	// Parton-level crude x-section divided by the CKM matrix element
	FunDis *= beamHandler->totalXsCrude(DataSource::get().ckm[0][0], Q2) / DataSource::get().ckm2[0][0];

	// Jacobian factor
	FunDis *= xminimal * delalf * ((Q2 - DataSource::get().MW * DataSource::get().MW) * (Q2 - DataSource::get().MW
			* DataSource::get().MW) + DataSource::get().MW * DataSource::get().GW * DataSource::get().MW * DataSource::get().GW) / (DataSource::get().MW * DataSource::get().GW)
			/ Q2;

	if (FunDis < 0.0) {
		FunDis = 0.0;

		if(DataSource::get().printOutLevel == PrintOutLevel::Detailed){
			if(FunDis < -1e-7){
		            cout<<">>>>>>>>>>>>>>>>> ppWintegrand: <<<<<<<<<<<<<<<<<<<<<<"<<endl
					<<"x1,x2,Q2="<<x1<<" "<<x2<<" "<<Q2<<endl
					<<"FunDis,SigPar="<<FunDis<<" "<<endl;
					for(unsigned i = 0 ; i <x1pdf.size(); ++i ){
						cout<<"x1pdf["<<i<<"]="<<x1pdf[i]<<endl;
					}
					for(unsigned i = 0 ; i <x2pdf.size(); ++i ){
						cout<<"x2pdf["<<i<<"]="<<x2pdf[i]<<endl;
					}
		            cout<<"==================================================="<<endl;
			}
		}
	}
#ifdef DEBUG_FORTRAN_INTEGRAND
	double xarg2[2];
	//xarg2[0]=Xarg[1];
	//xarg2[1]=Xarg[0];
	double test = wh_funppwtest_(Xarg);
	//if(FunDis>1E-4)
	//cout<<"DEBUG: ppWintegrand::Density()"<<test/FunDis<<endl;

#endif

	return FunDis;

}


