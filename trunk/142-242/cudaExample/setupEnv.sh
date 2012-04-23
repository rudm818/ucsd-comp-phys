#!/bin/bash

# this stuff should be put in your ~/.bashrc file eventually (quick fix)
export PATH=$PATH:/software/common/cuda/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/common/cuda/lib/
echo "Path set: $PATH"
echo "LD Path set: $LD_LIBRARY_PATH"