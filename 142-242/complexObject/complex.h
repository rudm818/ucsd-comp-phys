#ifndef _COMPLEX_H_
#define _COMPLEX_H_

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

#endif