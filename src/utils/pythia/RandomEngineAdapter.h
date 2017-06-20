/*
 * RandomEngineAdapter.h
 *
 *  Created on: 03-08-2011
 *      Author: sobol
 */

#ifndef RANDOMENGINEADAPTER_H_
#define RANDOMENGINEADAPTER_H_

#include "Pythia.h"

class RandomEngineAdapter: public Pythia8::RndmEngine{
private:
	TRND& generator;
public:
	RandomEngineAdapter(TRND& target):generator(target){}
	double flat(){
		return generator.Flat();
	}
};

#endif /* RANDOMENGINEADAPTER_H_ */
