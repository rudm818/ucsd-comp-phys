#ifndef _TEMPLATE_FUNCTIONS_
#define _TEMPLATE_FUNCTIONS_

#include "../complexObject/complex.h"

//functions used in hw1_template.cpp

template <class floatType>
void fillRange(floatType* Axis, floatType LOW, floatType HI, int DIM){
	
	floatType delta = (HI-LOW)/float(DIM-1); //note the '-1' is required to fill whole range
	for(int i=0; i<DIM; i++){
		Axis[i] = LOW + floatType(i)*delta;
	}
}


//here we are using a funtion as a parameter
//described here: http://stackoverflow.com/questions/9410/how-do-you-pass-a-function-as-a-parameter-in-c
template <class floatType>
void fillPropagator1d(complex<floatType>* K,floatType* Xpos,floatType dX, floatType (*V)(floatType),int DIM){
	//some propagator code...
	//floatType a = (*V)(Xpos[i])
}

//fills and normalizes wave function
template <class floatType>
fillWaveFxn(complex<floatType>* psi,int DIM){
	//do it yourself
	//also normalize
}

#endif