/*
 *  complex.cpp
 *  
 *
 *  Created by Michael M Folkerts on 4/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>

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



int main(int argc, char* argv[]){
//	complex<float> x,y,z;
	complex<double> x,y,z;

	x.a = 1.0;
	x.b = 2.0;
	
	cCpy(y,x);   // y = x;
	
	cAdd(z,x,y); //z = x + y
	printf("Z = %f + i*%f\n", z.a, z.b);
	
	cSub(z,x,y); //z = x - y
	printf("Z = %f + i*%f\n", z.a, z.b);
	
	cMult(z,x,y);
	printf("Z = %f + i*%f\n", z.a, z.b);
	
	//test Indexing
	printf("Testing typcasting and indexing\n");
	complex<float> A[10];
	
	for(int i; i<10; i++){
		A[i].a = float(i);
		A[i].b = float(i)*float(i);
	}
	
	float* cArray = (float*)A;
	
	for(int i; i<2*10; i+=2){
		printf("Z = %f + i*%f\n", cArray[i], cArray[i+1]);
	}
	
}