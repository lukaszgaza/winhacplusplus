/*
 * YFSUtil.cpp
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#include "YFSUtil.h"
#include "../../input/DataSource.h"
#include <iostream>

namespace VINHAC {

YFSUtil::YFSUtil() {

}

YFSUtil::~YFSUtil() {
}


double VINHAC::YFSUtil::YWldec(const CLHEP::HepLorentzVector& Q,
		const CLHEP::HepLorentzVector& ql, const double& aKmax) {

	// Dummy photon mass
	double amg = 1.0E-20;
	//!
	double aM = Q.mag();
	double aml = ql.mag();
	double t = (Q - ql).mag2();
	// YFS IR factor Y, i.e. sum of virtual and real photon IR functions
	double Bvirt = ReBtug(t, aM, aml, amg);
	double Breal = Btilde(Q, ql, aM, aml, aKmax, amg);

	return Bvirt + Breal;
}

double VINHAC::YFSUtil::ReBtug(const double& t, const double& xm1,
		const double& xm2, const double& amg) {

	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / alfinv / pi;
	double epst = 1.E-6;
	double am1, am2;

	if (xm1 > xm2) {
		am1 = xm1;
		am2 = xm2;
	} else {
		am1 = xm2;
		am2 = xm1;
	}
	double am1s = pow(am1, 2);
	double am2s = pow(am2, 2);
	double am12 = am1 * am2;
	double amds = (am1 - am2) * (am1 + am2);
	// IR regulator LOG
	double RegLog = -2 * log(sqrt(am12) / amg);
	double ReB2pi;
	//--------------------------------------------------------------------
	// Non-zero t case:
	if (fabs(t) > epst) {
		double xnu = (am1s + am2s - t) / 2;
		double xla = sqrt((xnu + am12) * (xnu - am12));
		double xlanu = xla + xnu;
		// Function A(p1,p2)
		double A = log(xlanu / am12) / xla;
		// Function A1(p1,p2)
		double A1 = amds / t * log(am1 / am2) - 2 * pow(xla, 2) / t * A - 2;
		// Function A3(p1,p2)
		double tix = t / (2 * xla * (xlanu - am2s));
		double eta = am2s * tix;
		double zet = xlanu * tix;
		double xlg1 = log(xlanu / am1s);
		double xlg2 = log(xlanu / am2s);
		double xar1 = fabs((xla - xnu + am1s) / t);
		double xar2 = fabs((xlanu - am2s) / am2s);
		double xlet = log(fabs(eta)) * log(1 + eta);
		double xlze = log(fabs(zet)) * log(1 + zet);
		double A3 = log(2 * xla / am12) * A + (0.25 * (xlg1 + 2 * log(xar1))
				* xlg1 + 0.25 * (xlg2 - 2 * log(xar2)) * xlg2 + 0.5 * (xlet
				- xlze) + DILOGY(-eta) - DILOGY(-zet)) / xla;
		//	Function ReB(s)
		ReB2pi = (xnu * A - 1) * RegLog + 0.5 * A1 - xnu * A3;
		//--------------------------------------------------------------------
		// Zero t case:
	} else {
		ReB2pi = ((am1s + am2s) / amds * log(am1 / am2) - 1) * (RegLog + 0.5);
	}
	return alfpi * ReB2pi;
}

double VINHAC::YFSUtil::DILOGY(const double& X) {
	//C-------------------------------------------- REMARKS ---------------
	//C DILOGARITHM FUNCTION: DILOG(X)=INT( -LN(1-Z)/Z ) , 0 < Z < X .
	//C THIS IS THE CERNLIB VERSION.
	//C--------------------------------------------------------------------
	//IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	double Z = -1.644934066848226;
	double T, S;

	if (X < -1.0) {
		T = 1.0 / X;
		S = -0.5;
		Z = Z - 0.5 * pow(log(fabs(X)), 2);
	} else if (X <= 0.5) {
		T = X;
		S = 0.5;
		Z = 0.0;
	} else if (X == 1.0) {
		return 1.644934066848226;
	} else if (X <= 2.0) {
		T = 1.0 - X;
		S = -0.5;
		Z = 1.644934066848226 - log(X) * log(fabs(T));
	} else {
		Z = 3.289868133696453;
		T = 1.0 / X;
		S = -0.5;
		Z = Z - 0.5 * pow(log(fabs(X)), 2);
	}
	double Y = 2.666666666666667 * T + 0.666666666666667;
	double B = 0.000000000000001;
	double A = Y * B + 0.000000000000004;
	B = Y * A - B + 0.000000000000011;
	A = Y * B - A + 0.000000000000037;
	B = Y * A - B + 0.000000000000121;
	A = Y * B - A + 0.000000000000398;
	B = Y * A - B + 0.000000000001312;
	A = Y * B - A + 0.000000000004342;
	B = Y * A - B + 0.000000000014437;
	A = Y * B - A + 0.000000000048274;
	B = Y * A - B + 0.000000000162421;
	A = Y * B - A + 0.000000000550291;
	B = Y * A - B + 0.000000001879117;
	A = Y * B - A + 0.000000006474338;
	B = Y * A - B + 0.000000022536705;
	A = Y * B - A + 0.000000079387055;
	B = Y * A - B + 0.000000283575385;
	A = Y * B - A + 0.000001029904264;
	B = Y * A - B + 0.000003816329463;
	A = Y * B - A + 0.000014496300557;
	B = Y * A - B + 0.000056817822718;
	A = Y * B - A + 0.000232002196094;
	B = Y * A - B + 0.001001627496164;
	A = Y * B - A + 0.004686361959447;
	B = Y * A - B + 0.024879322924228;
	A = Y * B - A + 0.166073032927855;
	A = Y * A - B + 1.935064300869969;
	return S * T * (A - B) + Z;
}

double VINHAC::YFSUtil::Btilde(const CLHEP::HepLorentzVector& p1,
		const CLHEP::HepLorentzVector& p2, const double& am1,
		const double& am2, const double& aKmax, const double& amgam) {

	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / alfinv / pi;
	double epsb = 1.0E-10;

	double E1 = p1.e();
	double E2 = p2.e();
	double am12 = am1 * am2;
	double p1p2 = p1 * p2;
	if (p1p2 - am12 < 1.0E-10)
		return 0;
	double xlam = sqrt((p1p2 - am12) * (p1p2 + am12));
	// Function A(p1,p2)
	double A = 1.0 / xlam * log((p1p2 + xlam) / am12);
	// betas and big logs - limit beta->0 has to be treated carefully!
	double bt1s = fabs((1.0 - am1 / E1) * (1.0 + am1 / E1));
	double bet1 = 0, b1ln = 0;
	if (bt1s > epsb) {
		bet1 = sqrt(bt1s);
		b1ln = log((1.0 + bet1) * E1 / am1) / bet1;
	} else {
		b1ln = 1.0;
	}
	double bt2s = fabs((1.0 - am2 / E2) * (1.0 + am2 / E2));
	double bet2 = 0, b2ln = 0;
	if (bt2s > epsb) {
		bet2 = sqrt(bt2s);
		b2ln = log((1 + bet2) * E2 / am2) / bet2;
	} else {
		b2ln = 1.0;
	}
	// Function A4(p1,p2)
	double A4 = YFSUtil::A4anal(p1p2, E1, E2, am1, am2);

	// B-tilde(p1,p2;aKmax,amgam)
	double Btian = (p1p2 * A - 1) * 2 * log(2 * aKmax / amgam) + b1ln + b2ln
			+ p1p2 * A4;
	return alfpi * Btian;
}

double VINHAC::YFSUtil::A4anal(const double& p1p2, const double& E1,
		const double& E2, const double& am1, const double& am2) {

	bool RF1, RF2;
	double A4 = 0.0, A4p, A4z;
	//
	double Q2 = 2 * p1p2 - pow(am1, 2) - pow(am2, 2);
	// Check if momenta in rest frame of p1 or p2
	RF1 = (1 - am1 / E1) > 1.0E-10;
	RF2 = (1 - am2 / E2) > 1.0E-10;
	// Q2 <> 0
	if (fabs(Q2) > 1.0E-10) {
		//     ==========================
		//... NOT in p1/p2 rest frame - usual case
		if (RF1 && RF2) {
			// Q2 > 0 (usual case)
			if (Q2 > 0.0) {
				A4 = YFSUtil::A4Q2pos(p1p2, E1, E2, am1, am2);
				// Q2 < 0 (e.g. radiative leptonic W-deays):
			} else if (Q2 < 0.0) {
				A4 = YFSUtil::A4Q2neg(p1p2, E1, E2, am1, am2);
				A4p = YFSUtil::A4Q2pos(p1p2, E1, E2, am1, am2);
				A4z = YFSUtil::A4Q2zer(E1, E2, am1, am2);
			}
			//... In p1 rest frame
		} else if (RF2) {
			A4 = YFSUtil::A4Q2fRF(E2, am2, am1);
			// In p2 rest frame
		} else if (RF1) {
			A4 = YFSUtil::A4Q2fRF(E1, am1, am2);
		} else {
			A4 = 0.0;
		}
		// Q2=0 case (e.g. leptonic W-decays):
	} else {
		//     ====
		if (RF1 && RF2) {
			A4 = YFSUtil::A4Q2zer(E1, E2, am1, am2);
			//... In p1 rest frame
		} else if (RF2) {
			A4 = YFSUtil::A4Q2zRF(am2, am1);
			//... In p2 rest frame
		} else if (RF1) {
			A4 = YFSUtil::A4Q2zRF(am1, am2);
		} else {
			A4 = 0.0;
		}
	}
	//     =====
	return A4;
}

double VINHAC::YFSUtil::A4Q2pos(const double& p1p2, double E1, double E2,
		double am1, double am2) {

	double p1s = pow(E1, 2) - pow(am1, 2);
	double p2s = pow(E2, 2) - pow(am2, 2);
	if (p1s < p2s) {
		double tmp = am1;
		am1 = am2;
		am2 = tmp;
		tmp = E1;
		E1 = E2;
		E2 = tmp;
	}
	double Ep = E1 + E2;
	double Em = E1 - E2;
	double sm = am1 + am2;
	double dm = am1 - am2;
	double Q2 = 2 * p1p2 - pow(am1, 2) - pow(am2, 2);
	double xl = sqrt((Q2 + pow(sm, 2)) * (Q2 + pow(dm, 2)));
	double xq = sqrt(Q2 + pow(Em, 2));
	double qp = xq + Em;
	double qm = xq - Em;
	double et0 = sqrt(pow(E2, 2) - pow(am2, 2));
	if (p1p2 > E1 * E2)
		et0 = -et0;
	double et1 = sqrt(pow(E1, 2) - pow(am1, 2)) + xq;
	double y1 = 0.5 * ((xq - Ep) + (sm * dm + xl) / qp);
	double y2 = y1 - xl / qp;
	double y3 = 0.5 * ((xq + Ep) + (sm * dm + xl) / qm);
	double y4 = y3 - xl / qm;
	// Some auxiliary functions
	double Eln;
	if (fabs(Em) > 1.0E-10) {
		Eln = log(fabs(qm / qp)) * (log(fabs((et1 - y1) * (et1 - y4) / (et1 - y2)
				/ (et1 - y3))) - log(fabs((et0 - y1) * (et0 - y4) / (et0 - y2)
				/ (et0 - y3))));
	} else {
		Eln = 0.0;
	}
	double Vet0 = YFSUtil::Yijeta(y1, y4, et0) + YFSUtil::Yijeta(y2,
			y1, et0) + YFSUtil::Yijeta(y3, y2, et0)
			- YFSUtil::Yijeta(y3, y4, et0) + 0.5 * etaln(y1, y2, y3, y4,
			et0) * etaln(y2, y3, y1, y4, et0);
	double Vet1 = YFSUtil::Yijeta(y1, y4, et1) + YFSUtil::Yijeta(y2,
			y1, et1) + YFSUtil::Yijeta(y3, y2, et1)
			- YFSUtil::Yijeta(y3, y4, et1) + 0.5 * etaln(y1, y2, y3, y4,
			et1) * etaln(y2, y3, y1, y4, et1);
	// Function A4(p1,p2)
	return (Eln + Vet1 - Vet0) / xl;
}

double VINHAC::YFSUtil::Yijeta(const double& yi, const double& yj,
		const double& eta) {

	return 2 * DILOGY((yj - yi) / (eta - yi)) + 0.5 * pow(log(fabs((eta - yi)
			/ (eta - yj))), 2);
}

double VINHAC::YFSUtil::A4Q2neg(const double& p1p2, double E1, double E2,
		double am1, double am2) {

	// Some auxiliary variables
	if (am1 < am2) {
		double tmp = am1;
		am1 = am2;
		am2 = tmp;
		tmp = E1;
		E1 = E2;
		E2 = tmp;
	}
	double p1m = sqrt((E1 - am1) * (E1 + am1));
	double p2m = sqrt((E2 - am2) * (E2 + am2));

	double Ep = E1 + E2;
	double Em = E1 - E2;
	double sm = am1 + am2;
	double dm = am1 - am2;
	double Q2 = 2 * p1p2 - pow(am1, 2) - pow(am2, 2);
	double xl = sqrt((Q2 + pow(sm, 2)) * (Q2 + pow(dm, 2)));
	double xq = sqrt(Q2 + pow(Em, 2));
	double qp = xq + Em;
	double et0 = p2m / xq;
	double et1 = p1m / xq - 1;
	double y1 = -(xq + Ep + (sm * dm + xl) / Q2 * qp) / 2 / xq;
	double y2 = y1 + xl / Q2 * qp / xq;
	double y3 = -(xq - Ep + (sm * dm + xl) / qp) / 2 / xq;
	double y4 = y3 + xl / qp / xq;
	// Some auxiliary functions
	double Eln;
	if (fabs(Em) > 1.0E-10) {
		Eln = log(fabs(pow(qp, 2) / Q2)) * (etaln(y1, y4, y2, y3, et1) - etaln(
				y1, y4, y2, y3, et0));
	} else {
		Eln = 0.0;
	}
	double Vet0 = Yijeta(y1, y4, et0) + Yijeta(y2, y1, et0) + Yijeta(y3, y2,
			et0) - Yijeta(y3, y4, et0) + 0.5 * etaln(y1, y2, y3, y4, et0)
			* etaln(y2, y3, y1, y4, et0);
	double Vet1 = Yijeta(y1, y4, et1) + Yijeta(y2, y1, et1) + Yijeta(y3, y2,
			et1) - Yijeta(y3, y4, et1) + 0.5 * etaln(y1, y2, y3, y4, et1)
			* etaln(y2, y3, y1, y4, et1);
	// Function A4(p1,p2)
	double A4 = (Eln + Vet1 - Vet0) / xl;
	if (A4 < 1.0E6) {

	} else {
		if(DataSource::get().printOutLevel == PrintOutLevel::Detailed){
		std::cout<<"A4Q2neg"<<std::endl
		<<"E1,am1="<<E1<<""<<am1<<std::endl
		<<"E2,am2="<<E2<<""<<am2<<std::endl
		<<"p1p2,E1E2="<<p1p2<<""<<E1*E2<<std::endl
		<<"et0,et1="<<et0<<""<<et1<<std::endl
		<<"y1,y2  ="<<y1<<""<<y2<<std::endl
		<<"y3,y4  ="<<y3<<""<<y4<<std::endl
		<<"Eln,Vet0,Vet1="<<Eln<<""<<Vet0<<""<<Vet1<<std::endl
		<<"A4="<<A4<<std::endl;
		}
	}
	return A4;
}

double VINHAC::YFSUtil::A4Q2zer(double E1, double E2, double am1,
		double am2) {

	// Some auxiliary variables

	if (am1 < am2) {
		double tmp = am1;
		am1 = am2;
		am2 = tmp;
		tmp = E1;
		E1 = E2;
		E2 = tmp;
	}

	double am2s = pow(am2, 2);
	double Em = E1 - E2;
	double dm = (am1 + am2) * (am1 - am2);
	double xem = dm / 2 / pow(Em, 2);
	double p1m = sqrt((E1 - am1) * (E1 + am1));
	double p2m = sqrt((E2 - am2) * (E2 + am2));
	double et0 = p2m / Em;
	double et1 = p1m / Em - 1;
	double y1 = E2 / Em;
	double y2 = y1 - xem;
	double y3 = -y1 + 2 * am2s / dm;
	// eta_i - y_j, y_j - y_i
	double e0y1 = -am2s / Em / (E2 + p2m);
	double e0y2 = e0y1 + xem;
	double e0y3 = et0 - y3;
	double e1y1 = et1 - y1;
	double e1y2 = et1 - y2;
	double e1y3 = et1 - y3;
	double y2y1 = -xem;
	double y3y2 = y3 - y2;
	// Some auxiliary functions
	double Eln = -log(xem) * log(fabs(e1y1 / e0y1 * e0y2 / e1y2 * e0y3 / e1y3));
	double Vet0 = 0.5 * pow(log(fabs(e0y1 * e0y2 / e0y3)), 2) + log(fabs(e0y1))
			* log(fabs(e0y1) / pow(e0y2, 2)) + 2 * DILOGY(y2y1 / e0y1) + 2
			* DILOGY(y3y2 / e0y2);
	double Vet1 = 0.5 * pow(log(fabs(e1y1 * e1y2 / e1y3)), 2) + log(fabs(e1y1))
			* log(fabs(e1y1) / pow(e1y2, 2)) + 2 * DILOGY(y2y1 / e1y1) + 2
			* DILOGY(y3y2 / e1y2);
	// Function A4(p1,p2)
	double A4 = (Eln + Vet1 - Vet0) / dm;
	if (A4 < 1.0E5) {

	} else {
		if(DataSource::get().printOutLevel == PrintOutLevel::Detailed){
			std::cout<<""<<std::endl
			 <<"am1,am2="<<am1<<" "<<am2<<std::endl
			 <<"E1,E2="<<E1<<" "<<E2<<std::endl
			 <<"xem,et0,et1="<<xem<<" "<<et0<<" "<<et1<<std::endl
			 <<"y1,y2,y3="<<y1<<" "<<y2<<" "<<y3<<std::endl
			 <<"e0y1,e0y2,e0y3="<<e0y1<<" "<<e0y2<<" "<<e0y3<<std::endl
			 <<"e1y1,e1y2,e1y3="<<e1y1<<" "<<e1y2<<" "<<e1y3<<std::endl
			 <<"y2y1,y3y2="<<y2y1<<" "<<y3y2<<std::endl
			 <<"Eln,Vet0,Vet1="<<Eln<<" "<<Vet0<<" "<<Vet1<<std::endl
			 <<"A4="<<A4<<std::endl;
		}
	}
	return A4;
}

double VINHAC::YFSUtil::A4Q2fRF(const double& E1, const double& am1,
		const double& am2) {

	// Some auxiliary variables
	double am1s = pow(am1, 2);
	double p1 = sqrt((E1 - am1) * (E1 + am1));
	double Epp = E1 + p1;
	double Emp = am1s / Epp;
	double xmm = am2 - Emp;
	double xmp = (am2 * (am2 - 2 * E1) + am1s) / xmm;
	double twp = 2 * p1;
	double tmp = twp * am2;
	double ymp = am2 * Epp - am1s;
	double ymm = am2 * Emp - am1s;
	// Function A4(p1,p2)

	double A4 = (log(fabs(xmm / xmp)) * log(Epp / am2) - 2 * log(twp * xmm / am2
			/ am1) * log(Epp / am1) + 2 * DILOGY(Emp / am2) - 2 * DILOGY(Epp
			/ am2) + DILOGY(-xmp / twp) - DILOGY(xmm / twp) + DILOGY(ymp / tmp)
			- DILOGY(-ymm / tmp)) / tmp;
	return A4;
}

double VINHAC::YFSUtil::A4Q2zRF(const double& am1, const double& am2) {

	//
	// Some auxiliary variables
	double xm12 = (am1 - am2) * (am1 + am2);
	// Function A4(p1,p2)
	double A4 = -2 * (pow(log(am1 / am2), 2) + DILOGY(xm12 / pow(am1, 2)))
			/ xm12;
	return A4;
}

double VINHAC::YFSUtil::deltaQED(const double& amW, const double& aml) {

	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / pi / alfinv;
	//
	//	 Based on: Bardin & Passarino, "The standard model ...", Oxford, 1999.
	//	*WP      del = alfpi * ( LOG(amW/aml) - 1.5 )
	//	 Based on: Marciano & Sirlin, Phys. Rev. D8 (1973) 3612.
	double del = alfpi * (log(amW / aml) + 0.5);
	return del;
}

double VINHAC::YFSUtil::StPair(const double& am1s, const double& am2s,
		const CLHEP::HepLorentzVector& p1, const CLHEP::HepLorentzVector& p2,
		const CLHEP::HepLorentzVector& pk) {

	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / pi / alfinv;

	double p1p2 = p1 * p2;
	double p1k = p1 * pk;
	double p2k = p2 * pk;
	double SoFa = (p1p2 / p2k - am1s / p1k) / p1k + (p1p2 / p1k - am2s / p2k)
			/ p2k;

	return alfpi / pi / 4.0 * SoFa;
}

double VINHAC::YFSUtil::VirSofQED(const double& amW, const double& aml,
		const double& aKmax) {

	double pi = 3.1415926535897932, alfinv = 137.03599911;
	double alfpi = 1.0 / pi / alfinv;

	double Blog = log(amW / aml);
	double viso = 2.0 * (Blog - 1.0) * log(2.0 * aKmax / amW) + 1.5 * Blog
			- pow(pi, 2.0) / 6.0 + 1.0;
	return alfpi * viso;
}


}
