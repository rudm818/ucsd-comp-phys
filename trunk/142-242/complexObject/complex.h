#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include <cbals.h>

template <class FP>
struct complex{
	FP a;
	FP b;
};

//C=A+B
template <class FP>
inline void cAdd(complex<FP>& C, const complex<FP>& A, const complex<FP>& B){
	C.a = A.a + B.a;
	C.b = A.b + B.b;
}

//C=A-B
template <class FP>
inline void cSub(complex<FP>& C, const complex<FP>& A, const complex<FP>& B){
	C.a = A.a - B.a;
	C.b = A.b - B.b;
}

//A=B
template <class FP>
inline void cCpy(complex<FP>& A, const complex<FP>& B){
	A.a = B.a;
	A.b = B.b;
}

//C=A*B
template <class FP>
inline void cMult(complex<FP>& C, const complex<FP>& A, const complex<FP>& B){
	printf("you can code multiplication yourself\n");
}


//__cblas__
//cublas ref: http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi?cmd=getdoc&coll=0650&db=man&fname=3%20INTRO_CBLAS
//cublas doc: http://www.prism.gatech.edu/~ndantam3/cblas-doc/doc/html/main.html

//single square matrix matrix multiplication
void cSingleSqMatMatMult(complex<float>* result, complex<float>* matA, complex<float>* matB, int dim){
	//look up blas rutine
	complex<float>* tempResult = new complex<float>[dim][dim];
	
	int M = dim; //rows of A and C
	int N = dim; //cols of B and C
	int K = dim; //cols of A and rows of B
	int iDimA = dim; // first dim of A
	int iDimB = dim; // first dim of B
	int iDimC = dim; // fisrt dim of C
	
	float alpha = 1.0f;
	float beta = 0.0f;
	
	// update C = alpha*A*B + beta*C
	cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
				M,N,K,
				alpha, (float*)matA,iDimA, (float*)matB, iDimB,
				beta,  (float*)tempResult,iDimC);
	
	//copy result back
	for(int i=0; i<dim*dim; i++){
		cCpy(result[i],tempResult[i]);
	}
	
	//clean up
	delete tempResult;
	
}

//single square matrix vector multiplication
void cSingleSqSymMatVectMult(complex<float>* result, complex<float>* matA, complex<float>* vectX, int dim){
	complex<float>* tempVectY = new complex<float>[dim];
	
	int M = dim; //rows of A
	int N = dim; //cols of A
	int iDimA = dim; // first dim of A
	int incX = 1;    // uses every incX'th value in X
	int incY = 1;    // uses every incY'th value in Y

	
	float alpha = 1.0f;
	float beta = 0.0f;
	
	// update Y = alpha*A*X + beta*X
	cblas_sgemv(CblasRowMajor,CblasNoTrans,
				M,N,
				alpha, (float*)matA,iDimA, (float*)vectX, incX,
				beta,  (float*)tempVectY,incY);
	
	//copy result back
	for(int i=0; i<dim*dim; i++){
		cCpy(result[i],tempVectY[i]); //complex number copy
	}
	
	//clean up
	delete tempVectY;
}

#endif