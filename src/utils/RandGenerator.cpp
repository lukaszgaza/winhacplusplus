#include<cmath>
using namespace std;

#include "RandGenerator.h"

VINHAC::RandGenerator::RandGenerator(){
	a=.0,b=.0;correct=true;
	gen=new CLHEP::RanluxEngine(300,3);
}

VINHAC::RandGenerator::~RandGenerator(){
	delete gen;
}

bool VINHAC::RandGenerator::SetSeed(double A,double B){

	if((2.0*A-B)*(2.0*A+B)<0)
	{
		cout<<"Niepoprawna inicjalizacja!!!"<<endl;
		correct=false;
		return false;
	}
	a=A,b=B;
}

double VINHAC::RandGenerator::Next(){
	if(!correct) {cout<<"Niepoprawna inicjalizacja!!!"<<endl;return -1;}

	int i=0;
	if((2.0*a+b)>0) i=0;else i=1;

	switch(i)
	{
	case 0:

		while(1)
		{
		double x=.0;
		if(gen->flat()<=(4.0*a)/(6.0*a-b)) 	x=pow(3.0*gen->flat()/(6.0*a-b),1.0/3.0);
		else  x=gen->flat()/(6.0*a-b);
		double u=gen->flat();
		if(u*(4.0*a*x*x+2.0*a-b)<=(4.0*a*x*x-2.0*(2.0*a-b)*x+2.0*a-b)) return 2.0*x-1;
		}

	case 1:

		while(1)
		{
		double x=.0;
		if(gen->flat()<=(2.0*a)/(b)) 	x=pow(3.0*gen->flat()/(2.0*b),1.0/3.0);
		else  x=pow(2.0*gen->flat()/(2.0*b),1.0/2.0);
		double u=gen->flat();
		if(u*(4.0*a*x*x-2.0*(2.0*a-b)*x)<=(4.0*a*x*x-2.0*(2.0*a-b)*x+2.0*a-b)) return 2.0*x-1;
		}


	default:
		return -1;
	}

	return -1;
}
