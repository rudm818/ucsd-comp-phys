#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include <cblas.h>

template <class floatType>
struct complex{
	floatType a;
	floatType b;
};

//C=A+B
template <class floatType>
inline void cAdd(complex<floatType>& C, const complex<floatType>& A, const complex<floatType>& B){
	C.a = A.a + B.a;
	C.b = A.b + B.b;
}

//C=A-B
template <class floatType>
inline void cSub(complex<floatType>& C, const complex<floatType>& A, const complex<floatType>& B){
	C.a = A.a - B.a;
	C.b = A.b - B.b;
}

//A=B
template <class floatType>
inline void cCpy(complex<floatType>& A, const complex<floatType>& B){
	A.a = B.a;
	A.b = B.b;
}

//C=A*B
template <class floatType>
inline void cMult(complex<floatType>& C, const complex<floatType>& A, const complex<floatType>& B){
	printf("you can code multiplication yourself\n");
}

// Z = c/B where c is real
template <class floatType>
inline complex<floatType> cDiv(floatType c, const complex<floatType>& B){
	complex<floatType> temp;
	
	temp.a = c*B.a/cMagSq(B); //real
	temp.b = -c*B.b/cMagSq(B);//complex
	
	return temp;
}

template <class floatType>
inline floatType cMagSq(complex<floatType> Z){
	return Z.a*Z.a+Z.b*Z.b;
}


//__cblas__
//cublas ref: http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi?cmd=getdoc&coll=0650&db=man&fname=3%20INTRO_CBLAS
//cublas doc: http://www.prism.gatech.edu/~ndantam3/cblas-doc/doc/html/main.html

//single square matrix matrix multiplication
void cSingleSqMatMatMult(complex<float>* result, complex<float>* matA, complex<float>* matB, const int dim){
	//look up blas rutine
	complex<float>* tempResult = (complex<float>*)malloc(sizeof(complex<float>)*dim*dim);
	
	int M = dim; //rows of A and C
	int N = dim; //cols of B and C
	int K = dim; //cols of A and rows of B
	int iDimA = dim; // first dim of A
	int iDimB = dim; // first dim of B
	int iDimC = dim; // fisrt dim of C
	
	float alpha = 1.0f;
	float beta = 0.0f;
	
	// update C = alpha*A*B + beta*C
	cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
				M,N,K,
				&alpha, (float*)matA,iDimA, (float*)matB, iDimB,
				&beta,  (float*)tempResult,iDimC);
	
	//copy result back
	for(int i=0; i<dim*dim; i++){
		cCpy(result[i],tempResult[i]);
	}
	
	//clean up
	free(tempResult);
	
}

//single square matrix vector multiplication
void cSingleSqSymMatVectMult(complex<float>* result, complex<float>* matA, complex<float>* vectX, int dim){
	complex<float>* tempVectY = (complex<float>*)malloc(sizeof(complex<float>)*dim);
	
	int M = dim; //rows of A
	int N = dim; //cols of A
	int iDimA = dim; // first dim of A
	int incX = 1;    // uses every incX'th value in X
	int incY = 1;    // uses every incY'th value in Y

	
	float alpha = 1.0f;
	float beta = 0.0f;
	
	// update Y = alpha*A*X + beta*Y
	cblas_cgemv(CblasRowMajor,CblasNoTrans,
				M,N,
				&alpha, (float*)matA,iDimA, (float*)vectX, incX,
				&beta,  (float*)tempVectY,incY);
	
	//copy result back
	for(int i=0; i<dim*dim; i++){
		cCpy(result[i],tempVectY[i]); //complex number copy
	}
	
	//clean up
	free(tempVectY);
}

#endif