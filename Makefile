beamformer: beamformer.c
	gcc -o $@ $^ -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran

dsacorr: dsacorr.c
	gcc -o $@ $^ -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran

dsacorr_little: dsacorr_little.c
	gcc -o $@ $^ -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran

rotate_file: rotate_file.c
	gcc -o $@ $^ -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran

extract_antennas: extract_antennas.c
	gcc -o $@ $^ -g -O3 -Wall -pthread -march=native -I/usr/local/include -I/usr/local/include/src -I/usr/local/cfitsio-3.47/include/ -L/usr/local/lib -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran -lm  -lsigproc

single_baseline: single_baseline.c
	gcc -o $@ $^ -g -O3 -Wall -pthread -march=native -I/usr/local/include -I/usr/local/include/src -I/usr/local/cfitsio-3.47/include/ -L/usr/local/lib -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran -lm  -lsigproc

toolkit: toolkit.cu
	/usr/local/cuda/bin/nvcc -D CUDA -ccbin=g++ -o $@ $^ -I/usr/local/include -I/usr/local/include/src -I/usr/local/cfitsio-3.47/include -arch=sm_75 -O3 -Xcompiler="-pthread" -DMATRIX_ORDER_TRIANGULAR -std=c++14 -L/usr/local/lib -lpsrdada -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran -L/usr/local/cuda/lib64 -lcudart -lm -L/usr/local/cfitsio-3.47/lib -lcfitsio

toolkit_dev: toolkit_dev.cu
	/usr/local/cuda/bin/nvcc -D CUDA -ccbin=g++ -o $@ $^ -I/usr/local/include -I/usr/local/include/src -I/usr/local/cfitsio-3.47/include -arch=sm_75 -O3 -Xcompiler="-pthread" -DMATRIX_ORDER_TRIANGULAR -std=c++14 -L/usr/local/lib -lpsrdada -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran -L/usr/local/cuda/lib64 -lcudart -lm -L/usr/local/cfitsio-3.47/lib -lcfitsio

.PHONY: clean all

clean:
	rm -f toolkit_dev toolkit single_baseline extract_antennas rotate_file dsacorr dsacorr_little beamformer

all: toolkit_dev toolkit single_baseline extract_antennas rotate_file dsacorr dsacorr_little beamformer

