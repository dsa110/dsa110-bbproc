/*gcc -o beamformer beamformer.c -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran
python beamformer was too slow, decided to use python to write up header etc but do actual beamforming in C.
This code should take 8 parameters:
* data file name
* calibration file name
* number of antennas in voltage file
* number of antennas to use in beamforming
* start frequency
* separation
* beam number
* output file name
assumes 48 channels for beamformer (weights for 8 data channels), ONLY 1 beam
greg hellbourg
ghellbourg@astro.caltech.edu
*/


#include "stdio.h"
#include "stdlib.h"
#include "sys/types.h"
#include "sys/socket.h"
#include "string.h"
#include "netinet/in.h"
#include "netdb.h"
#include <unistd.h>
#include <pthread.h>
#include <arpa/inet.h>
#include <math.h>

int NW = 48;    // number of channels for the beamformer
float PI = 3.141592653589793238;
float CVAC = 299792458.0;



int init_weights(char * fnam, char *flagants, float *antpos, float *weights, int nPols) {

        // assumes 64 antennas
        // antpos: takes only easting
        // weights: takes [ant, NW==48]

        FILE *fin;
	FILE *fflag;
        FILE *fants;
        int rd;

	int flags[64], nflag=0;
        fflag = fopen(flagants,"r");
	while (!feof(fflag)) {
	  fscanf(fflag,"%d\n",&flags[nflag]);
	  nflag++;
	}
	fclose(fflag);	 
	
        fin=fopen(fnam,"rb");

        rd = fread(antpos,64*sizeof(float),1,fin);
        rd = fread(weights,64*NW*nPols*2*sizeof(float),1,fin);
        float wnorm;
	int i;
	for (int ii=0;ii<64;ii++) {
	  for (int jj=0;jj<NW*nPols;jj++) {
	    i = ii*NW*nPols+jj;
	    wnorm = sqrt(weights[2*i]*weights[2*i] + weights[2*i+1]*weights[2*i+1]);
	    for (int kk=0;kk<nflag;kk++) {
	      if (flags[kk]==ii) {
		weights[2*i] = 0.;
		weights[2*i+1] = 0.;
	      }
	    }
	    if (wnorm!=0.0) {
	      weights[2*i] /= wnorm;
	      weights[2*i+1] /= wnorm;
	    }
	  }
	}

        fclose(fin);
        return 0;

}

void calc_weights(float *antpos, float *weights, float *freqs, float *wr, float *wi, float sep, int nBeamNum, int nPols) {


        float theta, afac, twr, twi;

        theta = sep*(127.-(float)nBeamNum)*PI/10800.; // radians
        for(int nAnt=0;nAnt<64;nAnt++){
                for(int nChan=0;nChan<48;nChan++){
                        for(int nPol=0;nPol<nPols;nPol++){
                                afac = -2.*PI*freqs[nChan*8+4]*theta/CVAC; // factor for rotate
                                twr = cos(afac*antpos[nAnt]);
                                twi = sin(afac*antpos[nAnt]);

                                wr[nAnt*(48*nPols)+nChan*nPols+nPol] = (twr*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2] - twi*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2+1]);
                                wi[nAnt*(48*nPols)+nChan*nPols+nPol] = (twi*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2] + twr*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2+1]);
                        }
                }
        }

}

void beamformer(char *input, float *wr, float *wi, float *output, int nChans, int nAnts, int nTimes, int nPols, int ri, int pol) {

  float inr_x, ini_x, inr_y, ini_y;
  float wrx, wix, wry, wiy;
  float rx, ix, ry, iy;
  float tmprealX, tmpimagX, tmprealY, tmpimagY;
  char v;
  
  for(int nTime=0;nTime<nTimes;nTime++){
    for(int nChan=0;nChan<48;nChan++){
      for(int i=0;i<8;i++){
	rx = 0;
	ix = 0;
	ry = 0;
	iy = 0;
	for(int nAnt=0;nAnt<nAnts;nAnt++){
	  v = input[nAnt*(nChans*nPols*nTimes)+(nChan*8+i)*(nPols*nTimes)+nTime*2];				  
	  inr_x = (float)((char)(((unsigned char)(v) & (unsigned char)(15)) << 4) >> 4);
	  //inr_x = (float)(((char)((v & 15) << 4)) >> 4);
	  v = input[nAnt*(nChans*nPols*nTimes)+(nChan*8+i)*(nPols*nTimes)+nTime*2];
	  ini_x = (float)((char)(((unsigned char)(v) & (unsigned char)(240))) >> 4);
	  //ini_x = (float)(((char)((v & 240))) >> 4);
	  v = input[nAnt*(nChans*nPols*nTimes)+(nChan*8+i)*(nPols*nTimes)+nTime*2+1];
	  inr_y = (float)((char)(((unsigned char)(v) & (unsigned char)(15)) << 4) >> 4);
	  //inr_y = (float)(((char)((v & 15) << 4)) >> 4);
	  v = input[nAnt*(nChans*nPols*nTimes)+(nChan*8+i)*(nPols*nTimes)+nTime*2+1];
	  ini_y = (float)((char)(((unsigned char)(v) & (unsigned char)(240))) >> 4);
	  //ini_y = (float)(((char)((v & 240))) >> 4);
	  
	  wrx = wr[nAnt*(48*nPols)+nChan*nPols];
	  wix = wi[nAnt*(48*nPols)+nChan*nPols];
	  wry = wr[nAnt*(48*nPols)+nChan*nPols+1];
	  wiy = wi[nAnt*(48*nPols)+nChan*nPols+1];
	  
	  rx += inr_x*wrx - ini_x*wix;
	  ix += inr_x*wix + ini_x*wrx;
	  ry += inr_y*wry - ini_y*wiy;
	  iy += inr_y*wiy + ini_y*wry;
	  
	}
	if (ri==0 && pol==0) 
	  output[nTime*nChans+nChan*8+i] = rx;
	if (ri==1 && pol==0)
	  output[nTime*nChans+nChan*8+i] = ix;
	if (ri==0 && pol==1) 
	  output[nTime*nChans+nChan*8+i] = ry;
	if (ri==1 && pol==1)
	  output[nTime*nChans+nChan*8+i] = iy;

      }
    }
  }
}

void usage()
{
  fprintf (stdout,
           "beamformer_volts [options]\n"
           " -d voltage data file name [no default]\n"
           " -f calibration file name [no default]\n"
           " -o output file name [no default]\n"
           " -z fch1 in MHz [default 1530]\n"
           " -s interbeam separation in arcmin [default 1.4]\n"
           " -n beam number [0 -- 255, default 127]\n"
	   " -q flagants file [no default]\n"
	   " -p pol [default B]\n"
	   " -c complexity [default real]\n"
           " -h print usage\n");
}


int main (int argc, char *argv[]) {


  int nChans = 384;
  int nAnts = 63;
  int nPols = 2;
  int nTimes = 2;
  int opol = 0;
  int ori = 0;

  // read params : fch1, fnam, fdataname, sep
  int arg = 0;
  float fch1 = 1530.0;
  float sep = 1.0;
  int nBeamNum = 127;
  char * fnam;
  fnam=(char *)malloc(sizeof(char)*200);
  char * fflag;
  fflag=(char *)malloc(sizeof(char)*200);
  sprintf(fnam,"nofile");
  char * fdata;
  fdata=(char *)malloc(sizeof(char)*200);
  sprintf(fdata,"nofile");
  char * fout;
  fout=(char *)malloc(sizeof(char)*200);
  sprintf(fout,"nofile");
  
  while ((arg=getopt(argc,argv,"d:f:o:u:z:s:n:q:p:c:ih")) != -1)
    {
      switch (arg)
	{
	case 'd':
	  if (optarg)
	    {
	      strcpy(fdata,optarg);
	      break;
	    }
	  else
	    {
	      printf("-d flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'f':
	  if (optarg)
	    {
	      strcpy(fnam,optarg);
	      break;
	    }
	  else
	    {
	      printf("-f flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'q':
	  if (optarg)
	    {
	      strcpy(fflag,optarg);
	      break;
	    }
	  else
	    {
	      printf("-q flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'o':
	  if (optarg)
	    {
	      strcpy(fout,optarg);
	      break;
	    }
	  else
	    {
	      printf("-o flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'p':
	  if (optarg)
	    {
	      opol = atoi(optarg);
	      break;
	    }
	  else
	    {
	      printf("-p flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'c':
	  if (optarg)
	    {
	      ori = atoi(optarg);
	      break;
	    }
	  else
	    {
	      printf("-c flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'z':
	  if (optarg)
	    {
	      fch1 = atof(optarg);
	      break;
	    }
	  else
	    {
	      printf("-z flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 's':
	  if (optarg)
	    {
	      sep = atof(optarg);
	      break;
	    }
	  else
	    {
	      printf("-s flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'n':
	  if (optarg)
	    {
	      nBeamNum = atoi(optarg);
	      break;
	    }
	  else
	    {
	      printf("-n flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'h':
	  usage();
	  return EXIT_SUCCESS;
	}
    }


  // compute beamformer weights
  float * output = (float *)malloc(sizeof(float)*nChans*nTimes);
  unsigned char * input = (char *)malloc(sizeof(char)*nAnts*nChans*nTimes*nPols);
  float * antpos = (float *)malloc(sizeof(float)*64); // easting
  float * weights = (float *)malloc(sizeof(float)*64*NW*nPols*2); // complex weights [ant, NW, pol, r/i]
  float * wr = (float *)malloc(sizeof(float)*64*NW*nPols); // complex weights [ant, NW, pol]
  float * wi = (float *)malloc(sizeof(float)*64*NW*nPols); // complex weights [ant, NW, pol]
  float * freqs = (float *)malloc(sizeof(float)*nChans); // freq
  for (int i=0;i<nChans;i++) freqs[i] = (fch1 - i*250./8192.)*1e6;
  init_weights(fnam,fflag,antpos,weights,nPols);
  calc_weights(antpos,weights,freqs,wr,wi,sep,nBeamNum,nPols);
  
  FILE *ptr;
  FILE *write_ptr;
  ptr = fopen(fdata,"rb");  // r for read, b for binary
  write_ptr = fopen(fout,"wb");  // w for write, b for binary
  
  long int sz;
  fseek(ptr, 0L, SEEK_END);
  sz = ftell(ptr);
  rewind(ptr);
  int nTotSam = 32768;

  fseek(ptr, 677376000, SEEK_SET);
  int rd;
  for(int nSam = 7000; nSam < 7000+4096; nSam++) {
    
    rd = fread(input,nAnts*nChans*nTimes*nPols,1,ptr);
    
    beamformer(input,wr,wi,output,nChans,nAnts,nTimes,nPols,ori,opol);
    
    fwrite(output,sizeof(float),nChans*nTimes,write_ptr);
    
  }
  fclose(ptr);
  fclose(write_ptr);

}
