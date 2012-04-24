/*
 *  complex.cpp
 *  
 *
 *  Created by Michael M Folkerts on 4/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include "complex.h"

int main(int argc, char* argv[]){
//	complex<float> x,y,z;
	complex<double> x,y,z;

	x.a = 1.0;
	x.b = 2.0;
	
	//test copy
	cCpy(y,x);   // y = x;
	
	//test add
	cAdd(z,x,y); //z = x + y
	printf("Z = %f + i*%f\n", z.a, z.b);
	
	//test subtraction
	cSub(z,x,y); //z = x - y
	printf("Z = %f + i*%f\n", z.a, z.b);
	
	//test mulitplication
	cMult(z,x,y);
	printf("Z = %f + i*%f\n", z.a, z.b);
	
	//test Indexing into an array
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