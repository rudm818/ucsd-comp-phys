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
#include <math.h>
#include <iostream>
using std::cout;
using std::endl;

#include<fstream>
using std::ofstream;

#include <algorithm>
using std::fill;

//allow complex numbers (from "our" library)
#include "../complexObject/complex.h"

//include function definitions
#include "feynmanK.h"

//change below to switch between double and float (double is SLOWER)
#define FLOAT float

//__CONSTANTS__
FLOAT PI = acos(-1.0);
FLOAT PLANCKS = 1.0;//Planck's constant
FLOAT MASS   = 1.0;//mass of particle


//__PARAMETERS__

//defines length of wave-funtion potential
const int X_DIM = 600;

//wave function parameters
const FLOAT SIGMA = 1.0f;
const FLOAT X0 = 0.75f;
const FLOAT NORM = pow(2.0/PI,0.25);

//define potential function (used for Legrangian = KE - PotentialFxn)
FLOAT PotentialFxn( FLOAT x){
	return pow(x,2.0f)/2.0f; //1D harmonic oscilator
}

//defines x-axis from - X_RANGE_LOW to + X_RANGE_HI
const FLOAT X_RANGE_LOW = -4.0f;
const FLOAT X_RANGE_HI  =  4.0f;

//integration factors and limits
const FLOAT DELTA_T = 0.01f;
const int K_POWER =  4; //the power of the aggregate KeN = (Kepsilon)^K_POWER
const int MAX_K_OPS = 100; // totalTime = K_POWER*DELTA_T * MAX_K_OPS


//__main function__
int main(int argc, char* argv[]){

	//sanity check
	if( 2*sizeof(FLOAT) != sizeof(complex<FLOAT>) ){
		cout << "ERROR, complex structure is bad!\n";
		exit(1);
	}

	//create and fill descretized x-axis:
	FLOAT* xAxis = new FLOAT[X_DIM];
	FLOAT deltaX = fillRange(xAxis, X_RANGE_LOW, X_RANGE_HI, X_DIM);

	cout << "Filling K-epsilon matrix..." << endl;
	//create and fill K-epsilon matrix
	complex<FLOAT>* Ke = (complex<FLOAT>*)malloc(sizeof(complex<FLOAT>)*X_DIM*X_DIM); //malloc();
	fillPropagator1D(Ke,deltaX,DELTA_T,PLANCKS,MASS,PotentialFxn,xAxis,X_DIM); //potential should be a FUNCTION that returns the potential as a function of position
	
	
	cout << "Computing aggregate K matrix, Ke^" << K_POWER << "..."<<endl;
	//create an aggregate K matri
	complex<FLOAT>* KeN = (complex<FLOAT>*)malloc(sizeof(complex<FLOAT>)*X_DIM*X_DIM);
	cSymSqMatMatMult(KeN, Ke, Ke, X_DIM);//Ke^2
	for (int n=2; n<K_POWER; n++) { //Ke^K_POWER (dX^N term is already added)
		cSymSqMatMatMult(KeN, KeN, Ke, X_DIM);
	}
	
	cout << "Generating wave function..." << endl;
	//create and fill waveFxn
	complex<FLOAT>* psi = (complex<FLOAT>*)malloc(sizeof(complex<FLOAT>)*X_DIM); //empty
	fillGaussianFxn(psi,NORM,X0,SIGMA,xAxis,X_DIM); //fills initial wave function
		
	for (int timeStep = 0; timeStep<MAX_K_OPS; timeStep++) {
		//multiply wave function by KeN (propegate)
		cSqMatVectMult(psi,KeN,psi,X_DIM);
		
		//then compute expectation values of x,v,potE,kinE,totE ...
		//computeExpectation of x
		cout << "<x> : " << expectationX(psi,xAxis,deltaX,X_DIM) << endl;	
		//save data to file...
	}
	
	//clean up memory
	delete xAxis;
	free(Ke);
	free(KeN);
	free(psi);
}
