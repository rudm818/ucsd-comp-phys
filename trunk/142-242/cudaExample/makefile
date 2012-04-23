
#set up cuda environment
CUDA_ROOT = /software/common/cuda/

all:
	nvcc cudaExample.cu -o cudaExample.x
	chmod +x cudaExample.x
run:all
	./cudaExample.x 128 4.0
	
clean:
	rm *.x