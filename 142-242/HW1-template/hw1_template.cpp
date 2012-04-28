/*
 *  hw1_template.cpp
 *  
 *
 *  Created by Michael M Folkerts on 4/23/12.
 *  Last Modified 4/25/12 -> fixed complex number issue in cblas, now using c++'s stl complex object
 *
 */

#include <iostream>
using std::cout;
using std::endl;

#include<fstream>
using std::ofstream;

#include <algorithm>
using std::fill;

//allow complex numbers (from "our" library)
//#include "../complexObject/complex.h"

//c++ complex number class, documented here: http://www.cplusplus.com/reference/std/complex/
#include <complex>
using std::complex;

//include function definitions
#include "feynmanK.h"

//change below to switch between double and float (double is SLOWER)
#define FLOAT double


//__CONSTANTS__
FLOAT PI      = acos(-1.0);
FLOAT PLANCKS = 1.0f;//Planck's constant
FLOAT MASS    = 1.0f;//mass of particle


//__PARAMETERS__

//defines length of wave-funtion potential
const int X_DIM = 601;

//Gaussian wave function parameters, see class notes
const FLOAT ALPHA = 2.0;
const FLOAT    X0 = 0.75;

//define potential function (used for Legrangian = KE - PotentialFxn)
//can be anything really!
FLOAT PotentialFxn( FLOAT x){
	return x*x/2.0; //1D harmonic oscilator with K=1
}

//defines x-axis from - X_RANGE_LOW to + X_RANGE_HI
const FLOAT X_RANGE_LOW = -4.0;
const FLOAT  X_RANGE_HI =  4.0;

//integration factors and limits
const int     N_STEPS = 128;    // totalTime = K_POWER*EPSILON_T * N_STEPS
const int     K_POWER = 4;     // the power of the aggregate KeN = (Ke)^K_POWER
const FLOAT EPSILON_T = 2.0*PI/FLOAT(N_STEPS);

//__main function__
int main(int argc, char* argv[]){

	cout << endl << argv[0] << " Running... " << endl << endl;

	//create and fill descretized x-axis:
	FLOAT* xAxis = new FLOAT[X_DIM];
	FLOAT deltaX = fillRange(xAxis, X_RANGE_LOW, X_RANGE_HI, X_DIM);

	cout << "Filling K-epsilon matrix... with dX=" << deltaX << endl;
	//create and fill K-epsilon matrix
	complex<FLOAT>* Ke = new complex<FLOAT>[X_DIM*X_DIM];
	fillPropagator1D(Ke,deltaX,EPSILON_T,PLANCKS,MASS,PotentialFxn,xAxis,X_DIM);
	//NOTE: in above 'PotentialFxn' should be a FUNCTION that returns the value of the potential with position as input
	
	
	cout << "Computing aggregate K matrix, Ke^" << K_POWER << "..."<<endl;
	//create an aggregate K matri
	complex<FLOAT>* KeN = new complex<FLOAT>[X_DIM*X_DIM];
	cSymSqMatMatMult(KeN, Ke, Ke, X_DIM);//Ke^2
	for (int n=2; n<K_POWER; n++) { //Ke^K_POWER (dX^N term is already added)
		cSymSqMatMatMult(KeN, KeN, Ke, X_DIM);
	}
	
	cout << "Generating wave function..." << endl;
	//create and fill waveFxn
	complex<FLOAT>* psi = new complex<FLOAT>[X_DIM]; //empty
	fillGaussianFxn(psi,X0,ALPHA,xAxis,X_DIM); //fills and normalizes initial wave function

	
	cout << "INITIAL <psi|psi>:" << sqWaveFxn(psi,deltaX,X_DIM) << endl;
	cout << "INITIAL <x> : " << expectationX(psi,xAxis,deltaX,X_DIM) << endl;
	
	//propagate
	for (int timeStep = 0; timeStep<N_STEPS; timeStep++) {

		//multiply wave function by KeN (propagate deltaT = K_POWER * EPSILON_T)
		cSqMatVectMult(psi,KeN,psi,X_DIM);
		
		//then compute expectation values of x,v,potE,kinE,totE ...
		
		//computeExpectation of x
		cout << "<x> : " << expectationX(psi,xAxis,deltaX,X_DIM) << endl;		
		
		//check normalization:
		cout << "<psi|psi>:" << sqWaveFxn(psi,deltaX,X_DIM) << endl;
			
		//save data to file...
	}
	
	
	//clean up memory
	delete[] xAxis;
	delete[] Ke;
	delete[] KeN;
	delete[] psi;
	
	cout << endl << argv[0] << " Done! " << endl << endl;
}
