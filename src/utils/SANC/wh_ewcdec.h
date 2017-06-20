/*
 * wh_ewcdec.h
 *
 *  Created on: 2010-03-21
 *      Author: kamil
 */

#ifndef WH_EWCDEC_H_
#define WH_EWCDEC_H_

extern "C" {

	double wh_delewdec_(double *aM,int *Idle,int *KeyGmu,
			double *aMW,double *aMZ, double *aMH,double Xmfe[20],
			double Xckm[3][3]);

}


#endif /* WH_EWCDEC_H_ */
