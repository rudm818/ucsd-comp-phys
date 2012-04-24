#include <cuda.h>

#include <iostream>
using std::cout;
using std::endl;

#include<fstream>
using std::ofstream;

#include <algorithm>
using std::fill;


/***** ERROR CHECKING MACRO *****/
cudaError_t _TempErrorCode;
#define CHECK_CUDA_ERROR() _TempErrorCode = cudaGetLastError(); if(_TempErrorCode) fprintf(stderr,"!!CUDA ERROR in %s at line %d : %s\n",__FILE__,__LINE__,cudaGetErrorString(_TempErrorCode));



/***** CUSTOM COMMAND LINE ARGUMENT PARSING *****/
//list of global variables (with default values)
int NumberOfArgs = 1; //how many constants are listed below

// you can add your own global variables to be parsed here
// (I start with underline to distiguish that it is a global variable):
int   _ArraySize = 1024;
float _IncrementValue = 1.0f;
char _OutputFile[] = "output.txt";


//this will display the global variable values before program starts running
void displayGlobals(void){
	cout<<"Setting ArraysSize to " << _ArraySize << endl;
	cout<<"Setting IncrementValue to " << _IncrementValue << endl;
}


//this parses the command line arguments
void parseArguments(int arg_count, char* args[]){
	
	//the first argument is always the program
	cout << "Running (" << args[0] << ")" << endl;
	
	if(arg_count > NumberOfArgs){
		// add your string to whatever parsing here
		_ArraySize = atoi(args[1]);
		_IncrementValue = atof(args[2]);
		//for strings just copy the pointer? (address):
		//OutputFile = args[3];
		displayGlobals();
		
	}else{
		//output usage
		cout << "Usage: "<< args[0] << " <ArraySize> <IncrementValue> " << endl;// <OutputFile>" << endl;
		//show default values
		displayGlobals();
	}

}

/***** A DEVICE FUNCTION *****/ 
__device__ float AddNum(float a, float b){
	return a + b;
}

/***** CUDA KERNEL ******/

/** this function increments the inArray by increment for all indicies less than MaxIndex **/
__global__ void incrementKernel(float* outArray,float* inArray, int MaxIndex, float increment){
	
	//the objects (gridDim,blockIdx,blockDim,threadIdx) are already defined:
	int threadIndex = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(threadIndex < MaxIndex){ //keep it safe
		outArray[threadIndex] = AddNum(inArray[threadIndex], increment);
	}
	
}


/***** MAIN *****/
int main(int argc, char* argv[]){

	parseArguments(argc, argv); //this will set the global variables

	//Device array pointers
	float* inArray_dev; //set to zero to avoid compile warnings
	float* outArray_dev;
	
	//Host array pointers
	float* inArray_host;
	float* outArray_host;
	
	//initialize arrays on host (using c++)
	inArray_host = new float[_ArraySize]; //equiv. to (float*)malloc(sizeof(float)*ARRAY_SIZE);
	outArray_host = new float[_ArraySize];
	//fill
	fill(inArray_host, inArray_host+_ArraySize, 1.0f); //fill with ones
	fill(outArray_host, outArray_host+_ArraySize, 0.0f); //fill with zeros
	
	//initialize arrays on device (GPU)
	cudaMalloc((void**)&inArray_dev, sizeof(float)*_ArraySize);
	CHECK_CUDA_ERROR();
	cudaMalloc((void**)&outArray_dev, sizeof(float)*_ArraySize);
	CHECK_CUDA_ERROR();
	
	//fill
	cudaMemset(inArray_dev, sizeof(float)*_ArraySize, 0); //set input array to zero (must be a byte value)
	CHECK_CUDA_ERROR();
	cudaMemset(outArray_dev, sizeof(float)*_ArraySize, 0); //set output array to zero (must be a byte value)
	CHECK_CUDA_ERROR();
	
	//copy input array to device
	//cudaMemcpy(DestinationPointer, SourcePointer, NumberOfBytes, cudaMemcpy[Host|Device]To[Host|Device]);
	cudaMemcpy(inArray_dev, inArray_host, sizeof(float)*_ArraySize, cudaMemcpyHostToDevice);
	
	//__LAUNCH KERNEL__
	//in general this geometry can be 3D, but for now we are just indexing a linear array
	int threadsPerBlock = 512; //this is typically the max for most GPUs except Fermi
	int blockCount;
	
	//special case for small array size:
	if(_ArraySize <= threadsPerBlock){
		blockCount = 1;
	}else{
		blockCount = _ArraySize/threadsPerBlock + 1; //max block size
	}
	
	incrementKernel <<< blockCount,threadsPerBlock >>> (outArray_dev, inArray_dev, _ArraySize, _IncrementValue);
	cudaThreadSynchronize();
	CHECK_CUDA_ERROR();

	// copy back results
	cudaMemcpy(outArray_host, outArray_dev, sizeof(float)*_ArraySize, cudaMemcpyDeviceToHost);
	CHECK_CUDA_ERROR();

	ofstream outFileStream;
	outFileStream.open(_OutputFile);
	
	//print output
	for(int i = 0; i<_ArraySize; i++){
		outFileStream << outArray_host[i] << endl;
	}

	cout << "Data saved in file " << _OutputFile << endl;

	//cleanup
	outFileStream.close();
	delete inArray_host; //like free()
	delete outArray_host;
	cudaFree(inArray_dev);
	cudaFree(outArray_dev);
	
}
