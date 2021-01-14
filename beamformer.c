/*
gcc -o beamformer beamformer.c -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran
python beamformer was too slow, decided to use python to write up header etc but do actual beamforming in C.
This code should take 5 parameters:
* data file name
* calibration file name
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

int NW = 48;	// number of channels for the beamformer
float PI = 3.141592653589793238;
float CVAC = 299792458;



int init_weights(char * fnam, float *antpos, float *weights) {

	// assumes 64 antennas
	// antpos: takes only easting
	// weights: takes [ant, NW==48] 

	FILE *fin;
	FILE *fants;
	int rd;

	fin=fopen(fnam,"rb");

	rd = fread(antpos,64*sizeof(float),1,fin);
	rd = fread(weights,64*NW*2*2*sizeof(float),1,fin);
	float wnorm;
	for (int i=0;i<64*NW*2;i++) {
		wnorm = sqrt(weights[2*i]*weights[2*i] + weights[2*i+1]*weights[2*i+1]);
		if (wnorm!=0.0) {
			weights[2*i] /= wnorm;
			weights[2*i+1] /= wnorm;
		}
	}

	fclose(fin);
	return 0;

}

void calc_weights(float *antpos, float *weights, float *freqs, float *wr, float *wi, float sep, int nBeamNum) {


	float theta, afac, twr, twi;
	
	theta = sep*(127.-(float)nBeamNum)*PI/10800.; // radians
	for(int nAnt=0;nAnt<64;nAnt++){
		for(int nChan=0;nChan<48;nChan++){
			for(int nPol=0;nPol<2;nPol++){
				afac = -2.*PI*freqs[nChan*8+4]*theta/CVAC; // factor for rotate
				twr = cos(afac*antpos[nAnt]);
				twi = sin(afac*antpos[nAnt]);

				wr[nAnt*(48*2)+nChan*2+nPol] = (twr*weights[(nAnt*(48*2)+nChan*2+nPol)*2] - twi*weights[(nAnt*(48*2)+nChan*2+nPol)*2+1]);
				wi[nAnt*(48*2)+nChan*2+nPol] = (twi*weights[(nAnt*(48*2)+nChan*2+nPol)*2] + twr*weights[(nAnt*(48*2)+nChan*2+nPol)*2+1]);
			}
		}
	}


}

void beamformer(char *input, float *wr, float *wi, unsigned char *output) {

	float inr_x, ini_x, inr_y, ini_y;
	float wrx, wix, wry, wiy;
	float rx, ix, ry, iy;
	float tmp;
	
	for(int nTime=0;nTime<2;nTime++){
		for(int nChan=0;nChan<48;nChan++){
			for(int i=0;i<8;i++){
				tmp = 0;
				for(int nAnt=0;nAnt<24;nAnt++){
					inr_x = (float)(((char)((input[nAnt*(384*2*2)+(nChan*8+i)*(2*2)+nTime*2] & 15) << 4)) >> 4);
					ini_x = (float)(((char)((input[nAnt*(384*2*2)+(nChan*8+i)*(2*2)+nTime*2] & 240))) >> 4);
					inr_y = (float)(((char)((input[nAnt*(384*2*2)+(nChan*8+i)*(2*2)+nTime*2+1] & 15) << 4)) >> 4);
					ini_y = (float)(((char)((input[nAnt*(384*2*2)+(nChan*8+i)*(2*2)+nTime*2+1] & 240))) >> 4);
					
					wrx = wr[nAnt*(48*2)+nChan*2];
					wix = wi[nAnt*(48*2)+nChan*2];
					wry = wr[nAnt*(48*2)+nChan*2+1];
					wiy = wi[nAnt*(48*2)+nChan*2+1];
					
					rx = inr_x*wrx + ini_x*wix;
					ix = -inr_x*wix + ini_x*wrx;
					ry = inr_y*wry + ini_y*wiy;
					iy = -inr_y*wiy + ini_y*wry;
					
					tmp += (rx*rx+ix*ix+ry*ry+iy*iy) / 24;
				}
				output[nTime*384+nChan*8+i] = (unsigned char)tmp;
			}
		}
	}
	
	
}

void usage()
{
  fprintf (stdout,
	   "t3_beamformer [options]\n"
	   " -d voltage data file name [no default]\n"
	   " -f calibration file name[no default]\n"
	   " -o output file name[no default]\n"
	   " -z fch1 in MHz [default 1530]\n"
	   " -s interbeam separation in arcmin [default 1.4]\n"
	   " -n beam number [0 -- 255, default 127]\n"
	   " -h print usage\n");
}


int main (int argc, char *argv[]) {
	
	
	// read params : fch1, fnam, fdataname, sep
	int arg = 0;
	float fch1 = 1530.0;
	float sep = 1.4;
	int nBeamNum = 127;
	char * fnam;
	fnam=(char *)malloc(sizeof(char)*200);
	sprintf(fnam,"nofile");
	char * fdata;
	fdata=(char *)malloc(sizeof(char)*200);
	sprintf(fdata,"nofile");
	char * fout;
	fout=(char *)malloc(sizeof(char)*200);
	sprintf(fout,"nofile");
	
	while ((arg=getopt(argc,argv,"d:f:o:z:s:n:h")) != -1)
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
	unsigned char * output = (char *)malloc(sizeof(char)*384*2);
	unsigned char * input = (char *)malloc(sizeof(char)*24*384*2*2);
	float * antpos = (float *)malloc(sizeof(float)*64); // easting
	float * weights = (float *)malloc(sizeof(float)*64*NW*2*2); // complex weights [ant, NW, pol, r/i]
	float * wr = (float *)malloc(sizeof(float)*64*NW*2); // complex weights [ant, NW, pol]
	float * wi = (float *)malloc(sizeof(float)*64*NW*2); // complex weights [ant, NW, pol]
	float * freqs = (float *)malloc(sizeof(float)*384); // freq
	for (int i=0;i<384;i++) freqs[i] = (fch1 - i*250./8192.)*1e6;
	init_weights(fnam,antpos,weights);
	calc_weights(antpos,weights,freqs,wr,wi,sep,nBeamNum);
 
	FILE *ptr;
	FILE *write_ptr;
	ptr = fopen(fdata,"rb");  // r for read, b for binary
	write_ptr = fopen(fout,"wb");  // w for write, b for binary
	int rd;
	for(int nSam = 0; nSam < 61440; nSam++) {
		
		rd = fread(input,24*384*2*2,1,ptr);

		beamformer(input,wr,wi,output);

		fwrite(output,384*2,1,write_ptr);
		
	}
	fclose(ptr);
	fclose(write_ptr);


}
