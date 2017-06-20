#include "CLHEP/Random/Random/RanluxEngine.h"

#define RANDGENERATOR_H_

namespace VINHAC
{
	class Generator
	{
	public:
			Generator(){};
			virtual double Next() = 0 ;
			virtual ~Generator(){};
	};

	class RandGenerator:public Generator
	{
	private:
		double a,b;
		bool correct;
		CLHEP::RanluxEngine *gen;

	public:
		RandGenerator();
		bool SetSeed(double A,double B);
		double Next();
		~RandGenerator();
	};
}
