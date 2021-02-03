/*
gcc -o dsacorr corr.c -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran
correlates voltage data

* data file name
* output file name
* time integration in ms
* frequency integration in MHz

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

void usage()
{
  fprintf (stdout,
	   "correlator [options]\n"
	   " -d voltage data file name [no default]\n"
	   " -o output file name[no default]\n"
	   " -t time integration in ms [default : 0.03]\n"
	   " -f frequency resolution in MHz [default : 0.03]\n"
	   " -h print usage\n");
}


int main (int argc, char *argv[]) {
	
	// data file constants
	int nAnts = 24;
	int nChans = 384;
	int nTimes = 2;
	int nPols = 2;
	float fDF = 0.03051757812; // frequency bandwidth per channel in MHz
	
	// read params : fInt, fRes
	int arg = 0;
	float fInt = 1./fDF/1000.;	// time integration in ms
	float fRes = fDF;	// frequency resolution in MHz
	
	char * fdata;
	fdata=(char *)malloc(sizeof(char)*100);
	sprintf(fdata,"nofile");
	char * fout;
	fout=(char *)malloc(sizeof(char)*100);
	sprintf(fout,"nofile");
	
	while ((arg=getopt(argc,argv,"d:o:t:f:h")) != -1)
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
			case 't':
			if (optarg)
			{
				fInt = atof(optarg);
				break;
			}
			case 'f':
			if (optarg)
			{
				fRes = atof(optarg);
				break;
			}
			case 'h':
			usage();
			return EXIT_SUCCESS;
		}
	}
	
	FILE *ptr;
	FILE *write_ptr;
	ptr = fopen(fdata,"rb");
	write_ptr = fopen(fout,"wb");
	int rd;
	
	long int sz;
	fseek(ptr, 0L, SEEK_END);
	sz = ftell(ptr);
	rewind(ptr);
	int nTotSam = (int)(floor(sz / (nAnts*nChans*nTimes*nPols)));
	
	float * corrs = (float *)malloc(sizeof(float)*(int)(nChans*nAnts*(nAnts-1)/2*nPols*2));
	float * autos = (float *)malloc(sizeof(float)*(int)(nChans*nAnts*nPols));
	float * crosspols = (float *)malloc(sizeof(float)*(int)(nChans*nAnts*nAnts*2));
	
	
	unsigned char * input = (unsigned char *)malloc(sizeof(char)*nAnts*nChans*nTimes*nPols);
	
	/* autos : [time, chan, antenna, pol] */
	/* corr_real : [time, chan, antenna1, antenna2, pol] */
	
	int nTimeSam = (int)ceil((double)(fInt / (2./fDF/1000.)));	// number of time integration
	if (nTimeSam > nTotSam) nTimeSam = nTotSam;
	if (nTimeSam <= 0) nTimeSam = 1;
	int nFreqsSam = (int)floor((double)(fRes / fDF));	// number of frequency integration
	if (nFreqsSam > nChans) nFreqsSam = nChans;
	if (nFreqsSam <= 0) nFreqsSam = 1;
	
	// find averages over freq channels so that all channels are used
	int idx = 0;
	while(nChans%(nFreqsSam+idx) != 0) idx++;
	nFreqsSam = nFreqsSam+idx;
	
	// arrays to be written
	float * corrswrite = (float *)malloc(sizeof(float)*(int)(nChans/nFreqsSam*nAnts*(nAnts-1)/2*nPols*2));
	float * autoswrite = (float *)malloc(sizeof(float)*(int)(nChans/nFreqsSam*nAnts*nPols));
	float * crosspolswrite = (float *)malloc(sizeof(float)*(int)(nChans/nFreqsSam*nAnts*nAnts*2));
	
	float in1rx, in1ix, in1ry, in1iy;
	float in2rx, in2ix, in2ry, in2iy;
	int i, j, nSam, nt, nfr, npl, nant1, nant2, timeslo;
	
	printf("\n");
	printf("processing file %s.\n", fdata);
	printf("found %d time samples.\n", nTotSam);
	printf("writing to file %s.\n", fout);
	printf("time resolution : %f ms.\n", nTimeSam * 2./fDF/1000.);
	printf("averaging %d time samples together.\n", (int)(nTimeSam*2));
	printf("frequency resolution : %f MHz.\n", nFreqsSam * fDF);
	printf("averaging %d frequency channels together.\n", nFreqsSam);
	printf("autocorrelations float encoded over %d bytes.\n", (int)sizeof(float));
	printf("crosscorrelations float encoded over %d bytes.\n", (int)sizeof(float));
	printf("\n");
	
	
	for(timeslo = 0; timeslo < (int)(nTotSam/nTimeSam); timeslo++) {
		
//		printf("writing time sample %d / %d\n", timeslo+1, (int)(nTotSam/nTimeSam));
		
		memset(autos,0,(nChans*nAnts*nPols)*sizeof(float));
		memset(corrs,0,(nChans*nAnts*(nAnts-1)/2*nPols*2)*sizeof(float));
		memset(crosspols,0,(nChans*nAnts*nAnts*2)*sizeof(float));
		for(nSam = 0; nSam < nTimeSam; nSam++) {	// avg sets of 2 time samples
			rd = fread(input,nAnts*nChans*nTimes*nPols,1,ptr);
			for(nt = 0; nt < nTimes; nt++) {
				for(nfr = 0; nfr < nChans; nfr++) {
					for(nant1 = 0; nant1 < nAnts; nant1++) {
						
						in1rx = (float)(((char)((input[nant1*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols] & 15) << 4)) >> 4);
						in1ix = (float)(((char)((input[nant1*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols] & 240))) >> 4);
						in1ry = (float)(((char)((input[nant1*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols+1] & 15) << 4)) >> 4);
						in1iy = (float)(((char)((input[nant1*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols+1] & 240))) >> 4);
						autos[nfr*(nAnts*nPols)+nant1*nPols] += (in1rx*in1rx + in1ix*in1ix);
						autos[nfr*(nAnts*nPols)+nant1*nPols+1] += (in1ry*in1ry + in1iy*in1iy);
						
						for(nant2 = nant1+1; nant2 < nAnts; nant2++) {
							in2rx = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols] & 15) << 4)) >> 4);
							in2ix = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols] & 240))) >> 4);
							in2ry = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols+1] & 15) << 4)) >> 4);
							in2iy = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols+1] & 240))) >> 4);
							
							corrs[ nfr*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2    ] += (in1rx*in2rx + in1ix*in2ix);
							corrs[ nfr*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 1] += (-in1rx*in2ix + in1ix*in2rx);
							corrs[ nfr*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 2] += (in1ry*in2ry + in1iy*in2iy);
							corrs[ nfr*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 3] += (-in1ry*in2iy + in1iy*in2ry);
							
							
						}
						for(nant2 = 0; nant2 < nAnts; nant2++) {
							in2rx = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols] & 15) << 4)) >> 4);
							in2ix = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+nfr*(nTimes*nPols)+nt*nPols] & 240))) >> 4);
							crosspols[nfr*nAnts*nAnts*2 + nant1*nAnts*2 + nant2*2] += (in1ry*in2rx + in1iy*in2ix);
							crosspols[nfr*nAnts*nAnts*2 + nant1*nAnts*2 + nant2*2+1] += (-in1ry*in2ix + in1iy*in2rx);
						}
					}
				}
			}
		}
		// average over frequencies
		memset(autoswrite,0,(nChans/nFreqsSam*nAnts*nPols)*sizeof(float));
		memset(corrswrite,0,(nChans/nFreqsSam*nAnts*(nAnts-1)/2*nPols*2)*sizeof(float));
		memset(crosspolswrite,0,(nChans/nFreqsSam*nAnts*nAnts*2)*sizeof(float));
		for (i = 0; i < (int)(nChans/nFreqsSam); i++){
			for (j = 0; j < nFreqsSam; j++){
				for (nant1 = 0; nant1 < nAnts; nant1++) {
					autoswrite[i*nAnts*nPols+nant1*nPols] += autos[(i*nFreqsSam+j)*nAnts*nPols+nant1*nPols];
					autoswrite[i*nAnts*nPols+nant1*nPols+1] += autos[(i*nFreqsSam+j)*nAnts*nPols+nant1*nPols+1];
					for (nant2 = nant1+1; nant2 < nAnts; nant2++) {
						corrswrite[i*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2    ] += corrs[ (i*nFreqsSam+j)*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2    ];
						corrswrite[i*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 1] += corrs[ (i*nFreqsSam+j)*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 1];
						corrswrite[i*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 2] += corrs[ (i*nFreqsSam+j)*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 2];
						corrswrite[i*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 3] += corrs[ (i*nFreqsSam+j)*(nAnts*(nAnts-1)/2*nPols*2) + ((nAnts-1)*nAnts/2 - (nAnts-1-nant1)*(nAnts-nant1)/2 + nant2 - (nant1+1))*nPols*2 + 3];
					}
					for (nant2 = 0; nant2 < nAnts; nant2++) {
						crosspolswrite[i*nAnts*nAnts*2 + nant1*nAnts*2 + nant2*2    ] += crosspols[(i*nFreqsSam+j)*nAnts*nAnts*2 + nant1*nAnts*2 + nant2*2    ];
						crosspolswrite[i*nAnts*nAnts*2 + nant1*nAnts*2 + nant2*2 + 1] += crosspols[(i*nFreqsSam+j)*nAnts*nAnts*2 + nant1*nAnts*2 + nant2*2 + 1];
					}
				}
			}
		}
		
		for (i = 0; i < nChans/nFreqsSam*nAnts*nPols; i++) autoswrite[i] = (autoswrite[i] / nTimes / nTimeSam / nFreqsSam);
		for (i = 0; i < nChans/nFreqsSam*nAnts*(nAnts-1)/2*nPols*2; i++) corrswrite[i] = (corrswrite[i] / nTimes / nTimeSam / nFreqsSam);
		for (i = 0; i < nChans/nFreqsSam*nAnts*nAnts*2; i++) crosspolswrite[i] = (crosspolswrite[i] / nTimes / nTimeSam / nFreqsSam);
		for (i = 0; i < nChans/nFreqsSam; i++) {
			fwrite(&autoswrite[i*nAnts*nPols],sizeof(float),nAnts*nPols,write_ptr);
			fwrite(&corrswrite[i*nAnts*(nAnts-1)/2*nPols*2],sizeof(float),nAnts*(nAnts-1)/2*nPols*2,write_ptr);
			fwrite(&crosspolswrite[i*nAnts*nAnts*2],sizeof(float),nAnts*nAnts*2,write_ptr);
		}
	}


	fclose(ptr);
	fclose(write_ptr);


}
