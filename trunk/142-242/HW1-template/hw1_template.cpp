/*
 *  hw1_template.cpp
 *  
 *
 *  Created by Michael M Folkerts on 4/23/12.
 *
 */

//cuda
//#include <cuda.h>

//blas
//#include <cblas.h>

#include <iostream>
using std::cout;
using std::endl;

#include<fstream>
using std::ofstream;

#include <algorithm>
using std::fill;

//allow complex numbers
#include "../complexObject/complex.h"

//include function definitions
#include "templateFunctions.cpp"

//change below to switch between double and float
#define FP float

//__PARAMETERS__

//defines length of wave-funtion potential
int X_DIM = 10;

//wave function parameters
FP SIGMA = 1.0f;

//defines x-axis from - X_RANGE_LOW to + X_RANGE_HI
FP X_RANGE_LOW = -1.0f;
FP X_RANGE_HI = 1.0f;


//integration factors and limits
FP DELTA_X = 0.01f;
int K_POWER = 4;

MAX_K_OPS = 100;


//define the potentail in a function
FP potential(FP x){
	return x*x;	
}

FP waveFxn(FP position){
	FP sigmaSquared = SIGMA*SIGMA; //whatever you want really
	return exp((x-X0)*(x-X0)/sigmaSquared);
}

int main(int argc, char* argc[]){

	//create and fill descretized x-axis:
	FP* xAxis = new FP[X_DIM];
	fillRange(xAxis, X_RANGE_LOW, X_RANGE_HI, X_DIM);


	//create and fill K-epsilon matrix
	complex<FP>* Ke = new complex<FP>[X_DIM][X_DIM];
	fillPropagator1D(Ke,xAxis,DELTA_X,potential,X_DIM); //potential should be a FUNCTION that returns the potential as a function of position
	
	//create an aggregate K matrix
	complex<FP>* KeN = new complex<FP>[X_DIM][X_DIM];
	cSingleSqMatMatMult(KeN, Ke, Ke, X_DIM);
	for (int n=1; n<K_POWER; n++) {
		cSingleSqMatMatMult(KeN, KeN, Ke, X_DIM);
	}
	
	//create and fill waveFxn
	complex<FP>* psi = new complex<FP>[X_DIM];
	fillWaveFxn(psi,waveFxn,X_DIM); //fills and normalizes wave function
			
	for (int timeStep = 0; timeStep<MAX_K_OPS; timeStep++) {
		//multiply wave function by KeN (propegate)
		cSingleSqSymMatVectMult(psi,KeN,psi,X_DIM);
		
		//then compute expectation values of x,v,potE,kinE,totE ...
		//save data to file...
	}
	
	//clean up memory
	delete xAxis;
	delete Ke;
	delete KeN;
	delete psi;
}
