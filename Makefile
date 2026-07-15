# dsa110-bbproc — offline M8 voltage-dump processing
# Builds on h23 (2x RTX 2080 Ti, sm_75, CUDA 11.1).

NVCC     ?= /usr/local/cuda/bin/nvcc
ARCH     ?= sm_75
NVFLAGS  = -arch=$(ARCH) -O3 -std=c++14 -Xcompiler="-Wall -pthread" \
           -I src -L/usr/local/cuda/lib64 -lcudart -lm

all: toolkit fake_voltages

toolkit: src/toolkit.cu src/bbproc.h
	$(NVCC) $(NVFLAGS) -o $@ src/toolkit.cu

fake_voltages: src/fake_voltages.cu src/bbproc.h
	$(NVCC) $(NVFLAGS) -o $@ src/fake_voltages.cu

.PHONY: all clean test

test: toolkit fake_voltages
	bash tests/roundtrip.sh

clean:
	rm -f toolkit fake_voltages
