/*
gcc -o incoh_beamformer incoh_beamformer.c -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran
inoherent beamformer
This code should take 5 parameters:
* data file name
* output file name
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


void beamformer(char *input, int nAnts, int nChans, int nTimes, int nPols, unsigned char *output) {

	float inr_x, ini_x, inr_y, ini_y;
	float tmp;
	
	for(int nt=0;nt<nTimes;nt++){
		for(int nc=0;nc<nChans;nc++){
			tmp = 0;
			for(int na=0;na<nAnts;na++){
				inr_x = (float)(((char)((input[na*(nChans*nTimes*nPols)+nc*(nTimes*nPols)+nt*nPols] & 15) << 4)) >> 4);
				ini_x = (float)(((char)((input[na*(nChans*nTimes*nPols)+nc*(nTimes*nPols)+nt*nPols] & 240))) >> 4);
				inr_y = (float)(((char)((input[na*(nChans*nTimes*nPols)+nc*(nTimes*nPols)+nt*nPols+1] & 15) << 4)) >> 4);
				ini_y = (float)(((char)((input[na*(nChans*nTimes*nPols)+nc*(nTimes*nPols)+nt*nPols+1] & 240))) >> 4);
				
				tmp += (inr_x*inr_x+ini_x*ini_x+inr_y*inr_y+ini_y*ini_y);
			}
			tmp /= nAnts;
			output[nt*nChans+nc] = (unsigned char)tmp;
		}
	}
	
	
}

void usage()
{
  fprintf (stdout,
	   "t3_INCOHERENT_beamformer [options]\n"
	   " -d voltage data file name [no default]\n"
	   " -o output file name[no default]\n"
	   " -a number of antennas in voltage file[default : 24]\n"
	   " -u number of antennas used in beamformer [default : 24]\n"
	   " -h print usage\n");
}


int main (int argc, char *argv[]) {
	
	// data file constants
	int nAnts = 24;
	int nUsedAnts = 24;
	int nChans = 384;
	int nTimes = 2;
	int nPols = 2;
	
	// read params : fch1, fnam, fdataname, sep
	int arg = 0;
	char * fdata;
	fdata=(char *)malloc(sizeof(char)*200);
	sprintf(fdata,"nofile");
	char * fout;
	fout=(char *)malloc(sizeof(char)*200);
	sprintf(fout,"nofile");
	
	while ((arg=getopt(argc,argv,"d:o:a:u:h")) != -1)
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
			case 'a':
			if (optarg)
			{
				nAnts = atoi(optarg);
				break;
			}
			else
			{
				printf("-a flag requires argument");
				usage();
				return EXIT_FAILURE;
			}
			case 'u':
			if (optarg)
			{
				nUsedAnts = atoi(optarg);
				break;
			}
			else
			{
				printf("-u flag requires argument");
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
	unsigned char * input = (char *)malloc(sizeof(char)*nAnts*384*2*2);

	FILE *ptr;
	FILE *write_ptr;
	ptr = fopen(fdata,"rb");  // r for read, b for binary
	write_ptr = fopen(fout,"wb");  // w for write, b for binary
	int rd;
	
	long int sz;
	fseek(ptr, 0L, SEEK_END);
	sz = ftell(ptr);
	rewind(ptr);
	int nTotSam = (int)(floor(sz / (nAnts*nChans*nTimes*nPols)));
	
	
	for(int nSam = 0; nSam < nTotSam; nSam++) {
		
		rd = fread(input,nAnts*nChans*nTimes*nPols,1,ptr);

		beamformer(input,nUsedAnts,nChans,nTimes,nPols,output);

		fwrite(output,nChans*nTimes,1,write_ptr);
		
	}
	fclose(ptr);
	fclose(write_ptr);


}
