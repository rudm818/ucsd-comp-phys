#ifndef _COMPLEX_H_
#define _COMPLEX_H_


extern "C" {
	#include <cblas.h>
};

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

// |Z|^2
template <class floatType>
inline floatType cMagSq(complex<floatType> Z){
	return Z.a*Z.a+Z.b*Z.b;
}


//__cblas__
//cublas ref: http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi?cmd=getdoc&coll=0650&db=man&fname=3%20INTRO_CBLAS
//cublas doc: http://www.prism.gatech.edu/~ndantam3/cblas-doc/doc/html/main.html

//single symmetric square matrix-matrix multiplication
template <class floatType>
void cSymSqMatMatMult(complex<floatType>* result, complex<floatType>* matA, complex<floatType>* matB, const int dim){
	//look up blas rutine
	complex<floatType>* tempResult = (complex<floatType>*)malloc(sizeof(complex<floatType>)*dim*dim);
	
	int M = dim; //rows of A and C
	int N = dim; //cols of B and C
	//int K = dim; //cols of A and rows of B
	int iDimA = dim; // first dim of A
	int iDimB = dim; // first dim of B
	int iDimC = dim; // fisrt dim of C
	
	floatType alpha = 1.0f;
	floatType beta = 0.0f;
	
	//check precision dynamically
	if( sizeof(floatType) == sizeof(float)){ //assume single precision
		
		// update C = alpha*A*B + beta*C where A = Transpose(A) (symmetric)
		cblas_csymm(CblasRowMajor,CblasLeft,CblasUpper,
					M,N,
					&alpha, (float*)matA,iDimA, (float*)matB, iDimB,
					&beta,  (float*)tempResult,iDimC);
		
	}else if(sizeof(floatType) == sizeof(double)){ //assume double precision
		
		// update C = alpha*A*B + beta*C where A = Transpose(A) (symmetric)
		cblas_zsymm(CblasRowMajor,CblasLeft,CblasUpper,
					M,N,
					&alpha, (double*)matA,iDimA, (double*)matB, iDimB,
					&beta,  (double*)tempResult,iDimC);
	}
	
	//copy result back
	for(int i=0; i<dim*dim; i++){
		cCpy(result[i],tempResult[i]);
	}
	
	//clean up
	free(tempResult);
	
}

//single square matrix-vector multiplication
template <class floatType>
void cSqMatVectMult(complex<floatType>* resultY, complex<floatType>* matA, complex<floatType>* vectX, int dim){
		
	int M = dim; //rows of A
	int N = dim; //cols of A
	int iDimA = dim; // first dim of A
	int incX = 1;    // uses every incX'th value in X
	int incY = 1;    // uses every incY'th value in Y

	
	floatType alpha = 1.0f;
	floatType beta = 0.0f;
	
	if( sizeof(floatType) == sizeof(float)){ //assume single precision complex number
		
		// update Y = alpha*A*X + beta*Y
		cblas_cgemv(CblasRowMajor,CblasNoTrans,
					M,N,
					&alpha, (float*)matA,iDimA, (float*)vectX, incX,
					&beta,  (float*)resultY,incY);
		
	}else if(sizeof(floatType) == sizeof(double)){ //assume double precision complex number
		
		// update Y = alpha*A*X + beta*Y
		cblas_zgemv(CblasRowMajor,CblasNoTrans,
					M,N,
					&alpha, (double*)matA,iDimA, (double*)vectX, incX,
					&beta,  (double*)resultY,incY);	
	
	}

}

#endif