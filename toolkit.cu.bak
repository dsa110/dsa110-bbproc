// -*- c++ -*-
/*

Strategy is to operate on a single voltage file, and produce a heap of stuff. 
 - read in and send to GPU
 - simple promote
 - optionally calibrate voltages
 - correlate and write out
 - optionally remove delays from visibilities
 - optionally average visibilities in frequency
 - rotate visibilities to particular beam (later can be RA/DEC)
 - write out beamformed filterbank

*/

#include <iostream>
#include <algorithm>
using std::cout;
using std::cerr;
using std::endl;
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <unistd.h>
#include <netdb.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <time.h>
#include <syslog.h>
#include <pthread.h>

#include <mma.h>
#include <cuda.h>
#include "cuda_fp16.h"
//#include "dada_cuda.h"
#include "dada_client.h"
#include "dada_def.h"
#include "dada_hdu.h"
#include "multilog.h"
#include "ipcio.h"
#include "ipcbuf.h"
#include "dada_affinity.h"
#include "ascii_header.h"
#include <thrust/device_ptr.h>
#include <thrust/fill.h>

#include <cuda_runtime_api.h>
using namespace nvcuda;

#define NANT 63
#define NCHAN 384
#define NT 30720
#define NBASE 2016
#define NPTR 8 // pols, small times, r/i
#define sep 1.0 // arcmin
#define AV 8
#define PI 3.141592653589793238
#define CVAC 299792458.0


// promoter to fp32
// run with NANT*NCHAN*NPTR/2/32 blocks of 32 threads
__global__ void promoter(char *input, float *output) {

  int bidx = blockIdx.x; // assume 16*48*NANT
  int tidx = threadIdx.x; // assume 32
  int iidx = bidx*32+tidx;
  
  output[2*iidx] = (float)((char)(((unsigned char)(input[iidx]) & (unsigned char)(15)) << 4) >> 4); //r
  output[2*iidx+1] = (float)((char)(((unsigned char)(input[iidx]) & (unsigned char)(240))) >> 4); //i

}

// correlator
// input is two packed time ints for all antennas
// also input antenna 1 and 2 indices for each baseline
// output is [2x time, baseline, freq, pols, r/i]
// run with NBASE*NCHAN/32 blocks of 32 threads
__global__ void correlator(float *input, float *output, int *a1, int *a2, float scfac, float *weights) {

  int bidx = blockIdx.x; // assume 16*48*NANT                                                              
  int tidx = threadIdx.x; // assume 32                                                                     
  int iidx = bidx*32+tidx;
  int basel = (int)(iidx/NCHAN); // baseline number
  int chgidx = (int)(bidx % (NCHAN/32)); // index of 32-channel group for this block
  int ch = (int)(iidx % NCHAN); // channel number
  
  // each block operates on 32 channels (one per thread)
  __shared__ float d1[32*NPTR];
  __shared__ float d2[32*NPTR];
  // start indices for each antenna from input
  int idx0_1 = a1[basel]*NCHAN*NPTR + chgidx*32*NPTR;
  int idx0_2 = a2[basel]*NCHAN*NPTR + chgidx*32*NPTR;

  // pull data into shared mem, for each antenna
  int ii = tidx*NPTR;
  for (int i=idx0_1+tidx*NPTR; i<idx0_1+(tidx+1)*NPTR; i++) {
    d1[ii] = input[i];
    ii++;
  }
  ii=tidx*NPTR;
  for (int i=idx0_2+tidx*NPTR; i<idx0_2+(tidx+1)*NPTR; i++) {
    d2[ii] = input[i];
    ii++;
  }

  // get weights for a1 and a2;
  float w_a1[4], w_a2[4];
  for (int i=0;i<4;i++) {
    w_a1[i] = weights[a1[basel]*192 + (int)(ch/8)*4 + i];
    w_a2[i] = weights[a2[basel]*192 + (int)(ch/8)*4 + i];
  }
  
  // now each thread can happily operate on a single channel
  // order is [time, pol, R/I]
  // make two separate arrays, each with [X*X / X*Y / Y*X / Y*Y, complexity]
  float output_tims[2][8], a1r, a1i, a2r, a2i;
  float w1r, w1i, w2r, w2i;
  // loop over times
  for (int ti=0;ti<2;ti++) {
    // loop over pols
    ii=0;
    for (int p1=0;p1<2;p1++) {
      for (int p2=0;p2<2;p2++) {
	
	a1r = d1[tidx*NPTR + ti*4 + p1*2];
	a1i = d1[tidx*NPTR + ti*4 + p1*2 + 1];
	a2r = d2[tidx*NPTR + ti*4 + p2*2];
	a2i = d2[tidx*NPTR + ti*4 + p2*2 + 1];

	w1r = a1r*w_a1[2*p1] - a1i*w_a1[2*p1+1];
	w2r = a2r*w_a2[2*p2] - a2i*w_a2[2*p2+1]; 
	w1i = a1r*w_a1[2*p1+1] + a1i*w_a1[2*p1];
	w2i = a2r*w_a2[2*p2+1] + a2i*w_a2[2*p2]; 
	
	output_tims[ti][2*ii] = w1r*w2r + w1i*w2i;
	output_tims[ti][2*ii+1] = w1r*w2i - w1i*w2r;

	ii++;
	
      }
    }
  }
  
  // write to output
  ii = basel*NCHAN*8 + ch*8;
  for (int i=0;i<8;i++) output[ii+i] += output_tims[0][i]*scfac;
  ii += NBASE*NCHAN*8;
  for (int i=0;i<8;i++) output[ii+i] += output_tims[1][i]*scfac;
  
  
}

// input has shape NBASE*NCHAN*8
// reduce to stokes I along NBASE axis using shared memory
// run with NCHAN blocks of 512 threads - will add 2016 baselines
__global__ void reduce_corrs(float *input, unsigned char *output, float scfac, int *a1, int *a2) {

  int bidx = blockIdx.x; // assume NCHAN
  int tidx = threadIdx.x; // assume 512                                                                  
  int iidx = bidx*512+tidx;

  volatile __shared__ float summer[512];

  // add into shared memory
  summer[tidx] = 0.;
  if (tidx<504) {
    if (a1[tidx]!=a2[tidx])
      summer[tidx] += input[tidx*NCHAN*8 + bidx*8] + input[tidx*NCHAN*8 + bidx*8 + 6];
    if (a1[tidx+504]!=a2[tidx+504])
      summer[tidx] += input[(tidx+1*504)*NCHAN*8 + bidx*8] + input[(tidx+1*504)*NCHAN*8 + bidx*8 + 6];
    if (a1[tidx+2*504]!=a2[tidx+2*504])
      summer[tidx] += input[(tidx+2*504)*NCHAN*8 + bidx*8] + input[(tidx+2*504)*NCHAN*8 + bidx*8 + 6];
    if (a1[tidx+3*504]!=a2[tidx+3*504])
      summer[tidx] += input[(tidx+3*504)*NCHAN*8 + bidx*8] + input[(tidx+3*504)*NCHAN*8 + bidx*8 + 6];
  }

  __syncthreads();

  // now reduce in shared memory
  if (tidx<256) {
    summer[tidx] += summer[tidx+256];
    __syncthreads();
    summer[tidx] += summer[tidx+128];
    __syncthreads();
    summer[tidx] += summer[tidx+64];
    __syncthreads();
    summer[tidx] += summer[tidx+32];
    __syncthreads();
    summer[tidx] += summer[tidx+16];
    __syncthreads();
    summer[tidx] += summer[tidx+8];
    __syncthreads();
    summer[tidx] += summer[tidx+4];
    __syncthreads();
    summer[tidx] += summer[tidx+2];
    __syncthreads();
    summer[tidx] += summer[tidx+1];
  }

  __syncthreads();

  if (tidx==0) output[bidx] = (unsigned char)(summer[0]*scfac);

}

// this kernel removes baseline delays by multiplying by exp(-2*pi*i*nu*tau)
// run with NBASE*NCHAN*4/32 blocks of 32 threads
// delays in ns
__global__ void delayer(float *input, float *freqs, float *delays) {

  int bidx = blockIdx.x; // assume 16*48*NANT
  int tidx = threadIdx.x; // assume 32
  int iidx = bidx*32+tidx;
  int bci = (int)(iidx/4);
  int basel = (int)(bci / NCHAN);
  int ch = (int)(bci % NCHAN);

  float vr, vi, arg=-2.*PI*freqs[ch]*delays[basel]*1e-9;
  vr = input[2*iidx]*cosf(arg) - input[2*iidx+1]*sinf(arg);
  vi = input[2*iidx]*sinf(arg) + input[2*iidx+1]*cosf(arg);

  input[2*iidx] = vr;
  input[2*iidx+1] = vi;
  
}

// kernel to enable frequency averaging of visibility output
// will only output XX and YY pols
// run with NBASE*NCHAN*4/AV/32 blocks of 32 threads
__global__ void fscrunch(float *input, float *output) {

  int bidx = blockIdx.x; 
  int tidx = threadIdx.x; // assume 32
  int iidx = bidx*32+tidx;
  int bcli = (int)(iidx/4);
  int poli = (int)(iidx % 4);
  int basel = (int)(bcli / (NCHAN/AV));
  int lch = (int)(bcli % (NCHAN/AV));

  int sumss[4];
  sumss[0] = 0;
  sumss[1] = 1;
  sumss[2] = 6;
  sumss[3] = 7;

  output[iidx] = 0.;
  for (int i=0;i<AV;i++) 
    output[iidx] += input[basel*NCHAN*8 + (AV*lch+i)*8 + sumss[poli]];
    
}

// really simple - adds the two times in correlator output
// run with NBASE*NCHAN*8/32 blocks of 32 threads
__global__ void adder(float *input, float *output) {

  int bidx = blockIdx.x; // assume 16*48*NANT
  int tidx = threadIdx.x; // assume 32
  int iidx = bidx*32+tidx;
  
  output[iidx] = input[iidx] + input[NBASE*NCHAN*8 + iidx];

}

// really simple - zeros correlator output
// run with 2*NBASE*NCHAN*8/32 blocks of 32 threads
__global__ void zeroer(float *input) {

  int bidx = blockIdx.x; // assume 16*48*NANT                                                              
  int tidx = threadIdx.x; // assume 32                                                                     
  int iidx = bidx*32+tidx;

  input[iidx] = 0.;

}


// CPU functions
int init_weights(char *wnam, float *antpos, float *weights, char *flagnam, int weight, int doflag);
// loads in weights
int init_weights(char * wnam, float *antpos, float *weights, char *flagnam, int weight, int doflag) {

  // assumes 64 antennas
  // antpos: takes only easting
  // weights: takes [ant, NW==48] 

  FILE *fin;
  FILE *fants;
  float wnorm;

  if (weight) {
    if (!(fin=fopen(wnam,"rb"))) {
      printf("Couldn't open weights file %s\n",wnam);
      return 1;
    }

    fread(antpos,64*sizeof(float),1,fin);
    fread(weights,64*48*2*2*sizeof(float),1,fin);

    for (int i=0;i<64*48*2;i++) {
      wnorm = sqrt(weights[2*i]*weights[2*i] + weights[2*i+1]*weights[2*i+1]);
      if (wnorm!=0.0) {
	weights[2*i] /= wnorm*wnorm;
	weights[2*i+1] /= wnorm*wnorm;
      }
    }           

    fclose(fin);
  }
  else {

    for (int i=0;i<64*48*2;i++) {
      weights[2*i] = 1.;
      weights[2*i+1] = 0.;
    }

  }
 

  int ant;
  if (doflag) {
    if (!(fants=fopen(flagnam,"r"))) {
      printf("Couldn't open flag ants file %s\n",flagnam);
      return 1;
    }
    
    while (!feof(fants)) {
      fscanf(fants,"%d\n",&ant);
      for (int j=0;j<48*2*2;j++) {
	weights[ant*48*2*2+j] = 0.0;
      }
    }

    fclose(fants);
    
  }
  
  printf("Loaded antenna positions and weights\n");
  return 0;

}

void calc_voltage_weights(float *antpos, float *weights, float *freqs, float *bfweights, int nBeamNum);
void calc_voltage_weights(float *antpos, float *weights, float *freqs, float *bfweights, int nBeamNum) {

  float theta, afac, twr, twi;
  theta = sep*(127.-(float)nBeamNum)*PI/10800.; // radians
  for(int nAnt=0;nAnt<64;nAnt++){
    for(int nChan=0;nChan<48;nChan++){
      for(int nPol=0;nPol<2;nPol++){
	afac = -2.*PI*freqs[nChan*8+4]*theta/CVAC; // factor for rotate
	twr = cos(afac*antpos[nAnt]);
	twi = sin(afac*antpos[nAnt]);

	bfweights[nAnt*(48*2*2)+nChan*2*2+nPol*2] = (twr*weights[(nAnt*(48*2)+nChan*2+nPol)*2] - twi*weights[(nAnt*(48*2)+nChan*2+nPol)*2+1]);
	bfweights[nAnt*(48*2*2)+nChan*2*2+nPol*2+1] = (twi*weights[(nAnt*(48*2)+nChan*2+nPol)*2] + twr*weights[(nAnt*(48*2)+nChan*2+nPol)*2+1]);
      }
    }
  }

}

// only does Stokes I for now
void sum_to_filterbank(float *corrout, unsigned char *filout) {

  float val;
  for (int iCh=0; iCh<NCHAN; iCh++) {
    val = 0.;
    for (int b=0; b<NBASE; b++) 
      val += (0.5*(corrout[b*NCHAN*8 + iCh*8] + corrout[b*NCHAN*8 + iCh*8 + 6]));
    filout[iCh]	= (unsigned char)(val);
  }

}

void usage()
{
  fprintf (stdout,
	   "toolkit [options]\n"
	   " -i input filename [no default]\n"
	   " -o output filename [no default - will not write if none given]\n"
	   " -t number of time integrations x 32.768us [default 8]\n"
	   " -w optional weights file\n"
	   " -f optional antenna flags file\n"
	   " -c set frequency of first channel in MHz [default 1530.0]\n"
	   " -b optionally set in 0-255 to rotate voltages to beam\n"
	   " -p coherent philterbank writing file [no default - will not write if none given]\n"
	   " -d file with NBASE delays to remove from baselines [optional - no default]\n"
	   " -a average visibilities by 8x in frequency\n"
	   " -h print usage\n");
}


// MAIN

int main (int argc, char *argv[]) {

  // use cuda device 1
  printf("Using GPU 1\n");
  cudaSetDevice(1);

  // command line arguments
  int arg = 0;
  char * finnam;
  finnam=(char *)malloc(sizeof(char)*100);
  char * foutnam;
  foutnam=(char *)malloc(sizeof(char)*100);
  int writing=0;
  int tint = 8;
  char * wnam;
  wnam=(char *)malloc(sizeof(char)*100);
  char * flagnam;
  flagnam=(char *)malloc(sizeof(char)*100);
  int weight=0, doflag=0;
  int beamn = -1;
  float fch1 = 1530.;
  char * filnam;
  filnam = (char *)malloc(sizeof(char)*100);
  int philwriting = 0;
  char * delnam;
  delnam = (char *)malloc(sizeof(char)*100);
  int delaying = 0;
  int averaging = 0;

  while ((arg=getopt(argc,argv,"i:o:t:w:f:c:b:p:d:ah")) != -1)
    {
      switch (arg)
	{
	case 'i':
	  if (optarg)
	    {
	      strcpy(finnam,optarg);
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-i flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'o':
	  if (optarg)
	    {
	      strcpy(foutnam,optarg);
	      writing=1;
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-o flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'p':
	  if (optarg)
	    {
	      strcpy(filnam,optarg);
	      philwriting=1;
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-p flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'd':
	  if (optarg)
	    {
	      strcpy(delnam,optarg);
	      delaying=1;
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-d flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'w':
	  if (optarg)
	    {
	      strcpy(wnam,optarg);
	      weight=1;
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-w flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'f':
	  if (optarg)
	    {
	      strcpy(flagnam,optarg);
	      doflag=1;
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-f flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'c':
	  if (optarg)
	    {
	      fch1 = atof(optarg);
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-c flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 't':
	  if (optarg)
	    {
	      tint = atoi(optarg);
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-t flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'b':
	  if (optarg)
	    {
	      beamn = atoi(optarg);
	      break;
	    }
	  else
	    {
	      syslog(LOG_ERR,"-b flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'a':
	  averaging=1;
	  break;
	case 'h':
	  usage();
	  return EXIT_SUCCESS;
	}
    }

  if (writing) printf("Reading from %s, writing to %s\n",finnam,foutnam);
  else printf("Reading from %s, no visibilities written\n",finnam);
  if (philwriting) printf("Will write coherent filterbank to %s\n",filnam);
  printf("Integrating by %d ints - check that this is power of 2\n",tint);
  printf("Assuming fch1 %f MHz\n",fch1);
  if (weight) printf("Will weight voltages using %s\n",wnam);
  if (doflag) printf("Will flag antennas using %s\n",flagnam);
  if (beamn>=0 && beamn <=255) printf("Will rotate voltages to beam %d\n",beamn);
  else
    printf("Not rotating voltages with beamn %d\n",beamn);
  if (averaging) printf("Will average visibilities by 8x in frequency\n");
  if (delaying) printf("Will apply baseline delays from %s\n",delnam);

  
  // allocate all memory

  // CPU
  char *indata = (char *)malloc(sizeof(char)*NANT*NCHAN*NPTR/2);
  float *outdata = (float *)malloc(sizeof(float)*NBASE*NCHAN*8);
  unsigned char *filout = (unsigned char *)malloc(sizeof(unsigned char)*NCHAN);
  int *h_a1 = (int *)malloc(sizeof(int)*NBASE);
  int *h_a2 = (int *)malloc(sizeof(int)*NBASE);
  // GPU
  char *d_indata;
  float *d_promoted, *d_corrout, *d_finalout, *d_avout;
  int *d_a1, *d_a2;
  unsigned char *d_filout;
  cudaMalloc((void **)&d_indata, NANT*NCHAN*(NPTR/2)*sizeof(char));
  cudaMalloc((void **)&d_promoted, NANT*NCHAN*NPTR*sizeof(float));
  cudaMalloc((void **)&d_corrout, 2*NBASE*NCHAN*8*sizeof(float));
  cudaMalloc((void **)&d_finalout, NBASE*NCHAN*8*sizeof(float));
  cudaMalloc((void **)&d_avout, NBASE*(NCHAN/AV)*4*sizeof(float));
  cudaMalloc((void **)&d_a1, NBASE*sizeof(int));
  cudaMalloc((void **)&d_a2, NBASE*sizeof(int));
  cudaMalloc((void **)&d_filout, NCHAN*sizeof(unsigned char));

  // load in delays
  float * h_delays = (float *)malloc(sizeof(float)*NBASE);
  float * d_delays;
  cudaMalloc((void **)&d_delays, NBASE*sizeof(float));
  FILE *fdel;
  if (delaying) {
    if (!(fdel=fopen(delnam,"r"))) {
      printf("could not open delay file %s\n",delnam);
      return(1);
    }
    for (int i=0;i<NBASE;i++)
      fscanf(fdel,"%f\n",&h_delays[i]);
    fclose(fdel);
    cudaMemcpy(d_delays,h_delays,NBASE*sizeof(float),cudaMemcpyHostToDevice);
  }
  
  // load in weights and antpos
  float * antpos = (float *)malloc(sizeof(float)*64); // easting
  float * weights = (float *)malloc(sizeof(float)*64*48*2*2); // complex weights [ant, NW, pol, r/i]
  float * bfweights = (float *)malloc(sizeof(float)*64*48*2*2); // complex weights [ant, NW, pol, r/i]
  float * freqs = (float *)malloc(sizeof(float)*NCHAN); // freq
  float * d_freqs;
  cudaMalloc((void **)&d_freqs, NCHAN*sizeof(float));
  for (int i=0;i<NCHAN;i++) freqs[i] = (fch1 - i*250./8192.)*1e6;
  cudaMemcpy(d_freqs,freqs,NCHAN*sizeof(float),cudaMemcpyHostToDevice);
  init_weights(wnam,antpos,weights,flagnam,weight,doflag);
  if (beamn>=0 && beamn<=255)
    calc_voltage_weights(antpos,weights,freqs,bfweights,beamn);
  float *d_weights;
  cudaMalloc((void **)&d_weights, 64*48*2*2*sizeof(float));
  if (beamn>=0 && beamn<=255)
    cudaMemcpy(d_weights,bfweights,64*48*2*2*sizeof(float),cudaMemcpyHostToDevice);
  else
    cudaMemcpy(d_weights,weights,64*48*2*2*sizeof(float),cudaMemcpyHostToDevice);
  
  // set up a1 and a2
  int ctr=0;
  for (int i=0;i<63;i++) {
    for (int j=i;j<63;j++) {
      h_a1[ctr] = i;
      h_a2[ctr] = j;
      ctr++;
    }
  }
  cudaMemcpy(d_a1,h_a1,NBASE*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(d_a2,h_a2,NBASE*sizeof(int),cudaMemcpyHostToDevice);

  // open input and output files
  FILE *fin, *fout, *flout;
  if (!(fin=fopen(finnam,"rb"))) 
    printf("could not open input file\n");
  if (writing) {
    if (!(fout=fopen(foutnam,"wb"))) 
      printf("could not open output file\n");
  }
  if (philwriting) {
    if (!(flout=fopen(filnam,"wb"))) 
      printf("could not open filterbank output file\n");
  }

  // loop over input

  printf("starting loop stay patient... :-)\n");
  
  int timi=0;
  ctr = 0;
  for(int bigI=0;bigI<NT;bigI++) {

    // read data, send to GPU, promote
    fread(indata, sizeof(char), NANT*NCHAN*NPTR/2, fin);
    cudaMemcpy(d_indata, indata, (NANT*NCHAN*NPTR/2)*sizeof(char), cudaMemcpyHostToDevice);
    promoter<<<NANT*NCHAN*NPTR/2/32,32>>>(d_indata, d_promoted);
    
    // deal with time integration
    if (timi==0) 
      zeroer<<<2*NBASE*NCHAN*8/32,32>>>(d_corrout);

    // correlate
    correlator<<<NBASE*NCHAN/32,32>>>(d_promoted,d_corrout,d_a1,d_a2,(1./(1.*tint)),d_weights);
    timi+=2;

    // deal with time integration
    if (timi>=tint) {

      // don't add up
      if (tint==1) {
		
	if (writing) {
	  if (delaying)
	    delayer<<<NBASE*NCHAN*4/32, 32>>>(d_corrout, d_freqs, d_delays);
	  if (averaging) {
	    fscrunch<<<NBASE*NCHAN*4/AV/32, 32>>>(d_corrout, d_avout);
	    cudaMemcpy(outdata, d_avout, NBASE*(NCHAN/AV)*4*sizeof(float), cudaMemcpyDeviceToHost);
	    fwrite(outdata,sizeof(float),NBASE*(NCHAN/AV)*4,fout);
	  }
	  else {
	    cudaMemcpy(outdata, d_corrout, NBASE*NCHAN*8*sizeof(float), cudaMemcpyDeviceToHost);
	    fwrite(outdata,sizeof(float),NBASE*NCHAN*8,fout);
	  }
	}
	if (philwriting) {
	  reduce_corrs<<<NCHAN,512>>>(d_corrout, d_filout, 0.25, d_a1, d_a2);
	  cudaMemcpy(filout, d_filout, NCHAN*sizeof(unsigned char), cudaMemcpyDeviceToHost);
	  fwrite(filout,sizeof(unsigned char),NCHAN,flout);
	}
		
	if (writing) {
	  if (delaying)
	    delayer<<<NBASE*NCHAN*4/32, 32>>>(d_corrout + NBASE*NCHAN*8, d_freqs, d_delays);
	  if (averaging) {
	    fscrunch<<<NBASE*NCHAN*4/AV/32, 32>>>(d_corrout + NBASE*NCHAN*8, d_avout);
	    cudaMemcpy(outdata, d_avout, NBASE*(NCHAN/AV)*4*sizeof(float), cudaMemcpyDeviceToHost);
	    fwrite(outdata,sizeof(float),NBASE*(NCHAN/AV)*4,fout);
	  }
	  else {
	    cudaMemcpy(outdata, d_corrout + NBASE*NCHAN*8, NBASE*NCHAN*8*sizeof(float), cudaMemcpyDeviceToHost);
	    fwrite(outdata,sizeof(float),NBASE*NCHAN*8,fout);
	  }
	}
	if (philwriting) {
	  reduce_corrs<<<NCHAN,512>>>(d_corrout + NBASE*NCHAN*8, d_filout, 0.25, d_a1, d_a2);
	  cudaMemcpy(filout, d_filout, NCHAN*sizeof(unsigned char), cudaMemcpyDeviceToHost);
	  fwrite(filout,sizeof(unsigned char),NCHAN,flout);
	}
	
      }

      // add up
      else {

	adder<<<NBASE*NCHAN*8/32,32>>>(d_corrout,d_finalout);	
	if (writing) {
	  if (delaying)
	    delayer<<<NBASE*NCHAN*4/32, 32>>>(d_finalout, d_freqs, d_delays);
	  if (averaging) {
	    fscrunch<<<NBASE*NCHAN*4/AV/32, 32>>>(d_finalout, d_avout);
	    cudaMemcpy(outdata, d_avout, NBASE*(NCHAN/AV)*4*sizeof(float), cudaMemcpyDeviceToHost);
	    fwrite(outdata,sizeof(float),NBASE*(NCHAN/AV)*4,fout);
	  }
	  else {
	    cudaMemcpy(outdata, d_finalout, NBASE*NCHAN*8*sizeof(float), cudaMemcpyDeviceToHost);
	    fwrite(outdata,sizeof(float),NBASE*NCHAN*8,fout);
	  }
	}
	if (philwriting) {
	  reduce_corrs<<<NCHAN,512>>>(d_finalout, d_filout, 0.125, d_a1, d_a2);
	  cudaMemcpy(filout, d_filout, NCHAN*sizeof(unsigned char), cudaMemcpyDeviceToHost);	  
	  fwrite(filout,sizeof(unsigned char),NCHAN,flout);
	}

      }

      ctr++;
      //printf("done with integration %d of %d\n",ctr,NT*2/tint);
      timi = 0;
    }
      

  }

  fclose(fin);
  if (writing) fclose(fout);
  if (philwriting) fclose(flout);
  
  cudaFree(d_indata);
  cudaFree(d_corrout);
  cudaFree(d_finalout);
  cudaFree(d_promoted);
  cudaFree(d_a1);
  cudaFree(d_a2);
  cudaFree(d_weights);
  cudaFree(d_filout);
  cudaFree(d_avout);
  cudaFree(d_delays);
  cudaFree(d_freqs);
  free(filout);
  free(antpos);
  free(weights);
  free(freqs);
  free(wnam);
  free(flagnam);
  free(indata);
  free(outdata);
  free(h_a1);
  free(h_a2);
  free(finnam);
  free(foutnam);
  free(filnam);
  free(delnam);
  free(h_delays);

}
