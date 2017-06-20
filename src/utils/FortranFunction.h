/*
 * FortranFunction.h
 *
 *  Created on: 2009-11-26
 *      Author: kamil
 */

#ifndef FORTRANFUNCTION_H_
#define FORTRANFUNCTION_H_
extern "C" {
typedef struct {
	double r, i;
} fcomplex;
double wh_mat0_(int *keyWid, double *amW, double *GaW, double *sW2,
		double *Uij, double p1[4], double p2[4], double Q[4], double q1[4],
		double q2[4]);
double wh_matfsr1_(int *keyWid, double *amW, double *GaW, double *sW2,
		double *Uij, double *QW, double p1[4], double p2[4], double Q[4],
		double q1[4], double q2[4], double pk[4]);
void wprobor_(double *sW2, double *Uij, double p1[4], double p2[4],
		double q[4], fcomplex amp[2][2][3]);
void wdecbor_(double *sW2, double *Cof, double *Uij, double q[4], double p1[4],
		double p2[4], fcomplex awdb[3][2][2]);
void wdecrad_(double *sW2, double *Cof, double *Uij,double *QW, double *Qf1, double q[4],
		double p1[4], double p2[4], double pk[4], fcomplex awdb[3][2][2][2]);
double wh_funppwtest_(double *Xarg);
void wh_initpdf_();

double rebtug_(double *t,double *xm1,double *xm2,double *amg);
double btilde_(double p1[4],double p2[4],double *am1,double *am2,double *aKmax,double *amgam);
double a4anal_(double *p1p2, double *E1, double *E2, double *am1,double *am2);
double a4q2frf_(double *E2, double *am2,double *am1);
double dilogy_(double *x);

}

#endif /* FORTRANFUNCTION_H_ */
