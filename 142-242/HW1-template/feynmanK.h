// AUTHOR: MICHAEL FOLKERTS - http://bit.ly/folkerts
// April 2012

#ifndef _FEYNMAN_H_
#define _FEYNMAN_H_

#include "../complexObject/complex.h"

//functions used in hw1_template.cpp

//fills axis with range of values
template <class floatType>
floatType fillRange(floatType* Axis, floatType LOW, floatType HI, int DIM){

	floatType delta = (HI-LOW)/float(DIM-1); //note the '-1' is required to fill whole range
	for(int i=0; i<DIM; i++){
		Axis[i] = LOW + floatType(i)*delta;
	}
	return delta;
}


//here we are using a function (the potential) as a parameter
//described here: http://stackoverflow.com/questions/9410/how-do-you-pass-a-function-as-a-parameter-in-c
template <class floatType>
void fillPropagator1D(complex<floatType>* K,floatType dX, floatType dT, floatType h, floatType m, floatType (*V_Fxn)(floatType),floatType* x, int DIM){
  
  //Normalization
  complex<floatType> A;
	//NEEDS CORRECT VALUE
	A.a = 0.0f; //real part
	A.b = 1.0f; //complex part
	printf("MESSAGE from Mike on line %i of %s: You may want to define the correct A value!\n",__LINE__,__FILE__);

  for (int i=0; i<DIM; i++) {
		for (int j=0; j<DIM; j++) {
			floatType argument = dT/h*(  m/2.0 * pow( (x[j]-x[i])/dT , 2.0 ) - (*V_Fxn)( (x[j]+x[i])/2.0 )  ) ;
			
      //compute complex number exp[i * arg] = cos(arg) + i*sin(arg)
			complex<floatType> Kelement;
			complex<floatType> preFactor = cDiv(dX,A);//complex division of real dX by complex A (dX/A)
			
      Kelement.a = preFactor.a*cos(argument);// real part
			Kelement.b = preFactor.b*sin(argument);// complex part
			
			K[i*DIM + j]= Kelement; //note dX has already been multiplied by each element (see lecture notes...)
		}
	}
}


//computes normalized initial wave function
template <class floatType>
void fillGaussianFxn(complex<floatType>* psi,floatType norm, floatType x0, floatType sigma, floatType* xAxis, int DIM){
	
	floatType sigmaSquared = sigma*sigma; //whatever you want really
	for(int i=0; i<DIM; i++){
		psi[i].a = norm * exp( pow(xAxis[i]-x0,2.0)/sigmaSquared );
		psi[i].b = 0.0;
	}
}

//computes the expectatino value of x
template <class floatType>
floatType expectationX(complex<floatType>* psi,floatType* xAxis,floatType dX, int DIM){
	
  //notes: <X> = integrate[conj(psi) * X * psi]
	floatType sum = 0.0f;
	
	for (int i=0; i<DIM; i++) {
		sum += dX * xAxis[i] * cMagSq(psi[i]);//x[i]*Conj(psi)*psi = x[i]*cMag(psi)^2
	}
	return sum;
}


#endif