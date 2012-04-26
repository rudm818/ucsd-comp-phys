// AUTHOR: MICHAEL FOLKERTS - http://bit.ly/folkerts
// April 2012

#ifndef _FEYNMAN_H_
#define _FEYNMAN_H_

//c++ complex number class, documented here: http://www.cplusplus.com/reference/std/complex/
#include <complex>
using std::complex;

extern "C" {

 #include <cblas.h>

};


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
  
	floatType pi = acos(-1.0);
	
	//Normalization
	floatType temp = 1.0; //NEEDS CORRECT VALUE
	complex<floatType> A(temp,temp); // constructor (real part, imag. part)
	printf("MESSAGE from Mike on line %i of %s: You may want to define the correct A value!\n",__LINE__,__FILE__);
	

	complex<floatType> preFactor = dX/A;
	cout << "preFactor dX/A: " << preFactor <<endl;
	
	for (int i=0; i<DIM; i++) {
		for (int j=0; j<DIM; j++) {
			floatType argument = dT/h*(  m/2.0 * pow( (x[j]-x[i])/dT , 2.0 ) - (*V_Fxn)( (x[j]+x[i])/2.0 )  ) ;
			
			//compute complex number exp[i * arg] = cos(arg) + i*sin(arg):
			complex<floatType> Kelement(preFactor.real()*cos(argument),   //real part
										preFactor.imag()*sin(argument));  //imag. part
			
			K[i*DIM + j] = Kelement; //note dX has already been multiplied by each element (see preFactor)
		}
	}
}



//computes and normalizes initial wave function
template <class floatType>
void fillGaussianFxn(complex<floatType>* psi, floatType x0, floatType alpha, floatType* xAxis, int DIM){

	floatType pi = acos(-1.0);
	
	for(int i=0; i<DIM; i++){
		psi[i] = complex<floatType>( pow(alpha/pi,0.25) * exp( - (alpha/2)*pow(xAxis[i]-x0,2.0)),
									 0.0);
	}
	
	normalizeWaveFxn(psi,DIM);
	
}

//normalizes given wave fxn
template <class floatType>
void normalizeWaveFxn(complex<floatType>* psi, int DIM){
	//now normalize!!
	floatType normFactor = 1.0/sqrt(sqWaveFxn(psi, DIM));

	for(int i=0; i<DIM; i++){
		psi[i] *= normFactor;
	}

}


//computes the expectatino value of x
template <class floatType>
floatType expectationX(complex<floatType>* psi,floatType* xAxis,floatType dX, int DIM){
	
  //notes: <X> = integrate[conj(psi) * X * psi]
	floatType sum = 0.0f;
	
	for (int i=0; i<DIM; i++) {
		sum += xAxis[i] * norm(psi[i]);// norm is abs(Z)^2
	}
	return sum;
}


template <class floatType>
floatType sqWaveFxn(complex<floatType>* wave, int dim){
	floatType sum = 0;
	
	for(int i=0; i<dim; i++){
		sum += norm(wave[i]);
	}
	return sum;
}

//this was just to check that cblas was working properly... it's slow
template <class floatType>
void myMatVectTest(complex<floatType>* result, complex<floatType>* K, complex<floatType>* psi, int dim){
	
	complex<floatType>* temp = new complex<floatType>[dim];
	
	for(int i=0; i<dim; i++){ //rows of K
		temp[i] = complex<floatType>(0.0,0.0);
		for(int j=0; j<dim; j++){ //cols of K
			temp[i] = temp[i] + ( K[i*dim+j] * psi[j]);
		}
	}
	
	//copy result back
	for(int i=0; i<dim; i++){
		result[i]=temp[i];
	}	
	
	delete[] temp;
}



//___cblas___
//cublas ref: http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi?cmd=getdoc&coll=0650&db=man&fname=3%20INTRO_CBLAS
//cublas doc: http://www.prism.gatech.edu/~ndantam3/cblas-doc/doc/html/main.html

//complex single and double symmetric square matrix-matrix multiplication
template <class floatType>
void cSymSqMatMatMult(complex<floatType>* result, complex<floatType>* matA, complex<floatType>* matB, const int dim){
	
	complex<floatType>* tempResult = new complex<floatType>[dim*dim];
	
	int M = dim; //rows of A and C
	int N = dim; //cols of B and C
	//int K = dim; //cols of A and rows of B
	int iDimA = dim; // first dim of A
	int iDimB = dim; // first dim of B
	int iDimC = dim; // fisrt dim of C
	
	complex<floatType> alpha = 1.0f; //THESE MUST BE COMPLEX!
	complex<floatType> beta = 0.0f;  //THESE MUST BE COMPLEX!
	
	//check precision dynamically
	if( sizeof(floatType) == sizeof(float)){ //assume single precision
		
		// update C = alpha*A*B + beta*C where A = Transpose(A) (symmetric)
		cblas_csymm(CblasRowMajor,CblasLeft,CblasUpper,
					M,N,
					&alpha, matA, iDimA, matB, iDimB,
					&beta,  tempResult,iDimC);
		
	}else if(sizeof(floatType) == sizeof(double)){ //assume double precision
		
		// update C = alpha*A*B + beta*C where A = Transpose(A) (symmetric)
		cblas_zsymm(CblasRowMajor,CblasLeft,CblasUpper,
					M,N,
					&alpha, matA, iDimA, matB, iDimB,
					&beta,  tempResult,iDimC);
	}
	
	//copy result back
	//maybe use memcpy here
	for(int i=0; i<dim*dim; i++){
		result[i]=tempResult[i];
	}
	
	//clean up
	delete[] tempResult;
	
}

//complex single and double precision square matrix-vector multiplication
template <class floatType>
void cSqMatVectMult(complex<floatType>* resultY, complex<floatType>* matA, complex<floatType>* vectX, int dim){
	
	complex<floatType>* tempY = new complex<floatType>[dim];
	int M = dim; //rows of A
	int N = dim; //cols of A
	int iDimA = dim; // first dim of A
	int incX = 1;    // uses every incX'th value in X
	int incY = 1;    // uses every incY'th value in Y
	
	
	complex<floatType> alpha = 1.0f; //THESE MUST BE COMPLEX!
	complex<floatType> beta = 0.0f;  //THESE MUST BE COMPLEX!
	
	if( sizeof(floatType) == sizeof(float)){ //assume single precision complex number
		
		// update Y = alpha*A*X + beta*Y
		cblas_cgemv(CblasRowMajor,CblasNoTrans,
					M,N,
					&alpha, matA,  iDimA, vectX, incX,
					&beta,  tempY, incY);
		
	}else if(sizeof(floatType) == sizeof(double)){ //assume double precision complex number
		
		// update Y = alpha*A*X + beta*Y
		cblas_zgemv(CblasRowMajor,CblasNoTrans,
					M,N,
					&alpha, matA,  iDimA, vectX, incX,
					&beta,  tempY, incY);	
		
	}
	
	//copy result back to user's vector
	for(int i=0; i<dim; i++){
		resultY[i]=tempY[i];
	}
	
	//clean up
	delete[] tempY;
	
}

#endif