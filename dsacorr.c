/*
gcc -o dsacorr corr.c -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran
correlates voltage data
* data file name
* output file name
* time integration in ms
* frequency integration in MHz
* number of antennas
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
	   " -a number of antennas [default : 24]\n"
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
	
	while ((arg=getopt(argc,argv,"d:o:t:f:a:h")) != -1)
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
			case 'a':
			if (optarg)
			{
				nAnts = atoi(optarg);
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
	
	unsigned char * input = (unsigned char *)malloc(sizeof(char)*nAnts*nChans*nTimes*nPols);
	
	int nTimeSam = (int)ceil((double)(fInt / (2./fDF/1000.)));	// number of time integration
	if (nTimeSam > nTotSam) nTimeSam = nTotSam;
	if (nTimeSam <= 1) nTimeSam = 1;
	int nFreqsSam = (int)floor((double)(fRes / fDF));	// number of frequency integration
	if (nFreqsSam > nChans) nFreqsSam = nChans;
	if (nFreqsSam <= 1) nFreqsSam = 1;
	
	// find averages over freq channels so that all channels are used
	int idx = 0;
	while(nChans%(nFreqsSam+idx) != 0) idx++;
	nFreqsSam = nFreqsSam+idx;
	
	int nBLines = (int)(nAnts*(nAnts+1)/2);
	int nSubChans = (int)floor((double)(nChans/nFreqsSam));
	int nAllSams = (int)floor((double)(nTotSam/nTimeSam));
	// arrays to be written
	float * corrs = (float *)malloc(sizeof(float)*(int)(nSubChans*nBLines*nPols*2*2)); // nChans/nFreqsSam channels, nAnts*(nAnts+1)/2 baselines, 4 pols (XX, YY, XY, YX), 2 R/I
	
	float in1rx, in1ix, in1ry, in1iy;
	float in2rx, in2ix, in2ry, in2iy;
	int i, j, nSam, nt, nc, npl, nant1, nant2, timeslo;
	
	printf("\n");
	printf("processing file %s.\n", fdata);
	printf("found %d time samples.\n", nTotSam*2);
	printf("writing to file %s.\n", fout);
	printf("time resolution : %f ms.\n", nTimeSam * 2./fDF/1000.);
	printf("averaging %d time samples together.\n", (int)(nTimeSam*2));
	printf("frequency resolution : %f MHz.\n", nFreqsSam * fDF);
	printf("averaging %d frequency channels together.\n", nFreqsSam);
	printf("correlations float encoded over %d bytes.\n", (int)sizeof(float));
	printf("producing %d time samples.\n", nAllSams);
	printf("producing %d frequency channels.\n", nSubChans);
	printf("computing %d baselines.\n", nBLines);
	printf("\n");
	
	size_t numwri;
	
	for(timeslo = 0; timeslo < nAllSams; timeslo++) {
		
//		printf("writing time sample %d / %d\n", timeslo+1, (int)(nTotSam/nTimeSam));
		
		memset(corrs,0,sizeof(float)*(int)(nSubChans*nBLines*nPols*2*2));
		
		for(nSam = 0; nSam < nTimeSam; nSam++) {	// avg sets of 2 time samples
			rd = fread(input,nAnts*nChans*nTimes*nPols,1,ptr);
			for(nt = 0; nt < nTimes; nt++) {
				for(nc = 0; nc < nSubChans; nc++) {
					for (i = 0; i < nFreqsSam; i++) {

						for(nant1 = 0; nant1 < nAnts; nant1++) {
							
							in1rx = (float)(((char)((input[nant1*(nChans*nTimes*nPols)+(nc*nFreqsSam+i)*(nTimes*nPols)+nt*nPols] & 15) << 4)) >> 4);
							in1ix = (float)(((char)((input[nant1*(nChans*nTimes*nPols)+(nc*nFreqsSam+i)*(nTimes*nPols)+nt*nPols] & 240))) >> 4);
							in1ry = (float)(((char)((input[nant1*(nChans*nTimes*nPols)+(nc*nFreqsSam+i)*(nTimes*nPols)+nt*nPols+1] & 15) << 4)) >> 4);
							in1iy = (float)(((char)((input[nant1*(nChans*nTimes*nPols)+(nc*nFreqsSam+i)*(nTimes*nPols)+nt*nPols+1] & 240))) >> 4);
							
							for(nant2 = nant1; nant2 < nAnts; nant2++) {
								in2rx = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+(nc*nFreqsSam+i)*(nTimes*nPols)+nt*nPols] & 15) << 4)) >> 4);
								in2ix = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+(nc*nFreqsSam+i)*(nTimes*nPols)+nt*nPols] & 240))) >> 4);
								in2ry = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+(nc*nFreqsSam+i)*(nTimes*nPols)+nt*nPols+1] & 15) << 4)) >> 4);
								in2iy = (float)(((char)((input[nant2*(nChans*nTimes*nPols)+(nc*nFreqsSam+i)*(nTimes*nPols)+nt*nPols+1] & 240))) >> 4);
								
								corrs[ (nBLines- (nAnts-nant1)*(nAnts-nant1+1)/2 + (nant2-nant1))*(nSubChans*nPols*2*2) + nc*(nPols*2*2)   ] += (float)(in1rx*in2rx + in1ix*in2ix) / 2. / nFreqsSam / nTimeSam; // real XX
								corrs[ (nBLines- (nAnts-nant1)*(nAnts-nant1+1)/2 + (nant2-nant1))*(nSubChans*nPols*2*2) + nc*(nPols*2*2) +1] += (float)(in1rx*in2ix - in1ix*in2rx) / 2. / nFreqsSam / nTimeSam; // imag XX
								corrs[ (nBLines- (nAnts-nant1)*(nAnts-nant1+1)/2 + (nant2-nant1))*(nSubChans*nPols*2*2) + nc*(nPols*2*2) +2] += (float)(in1rx*in2ry + in1ix*in2iy) / 2. / nFreqsSam / nTimeSam; // real XY
								corrs[ (nBLines- (nAnts-nant1)*(nAnts-nant1+1)/2 + (nant2-nant1))*(nSubChans*nPols*2*2) + nc*(nPols*2*2) +3] += (float)(in1rx*in2iy - in1ix*in2ry) / 2. / nFreqsSam / nTimeSam; // imag XY
								corrs[ (nBLines- (nAnts-nant1)*(nAnts-nant1+1)/2 + (nant2-nant1))*(nSubChans*nPols*2*2) + nc*(nPols*2*2) +4] += (float)(in1ry*in2rx + in1iy*in2ix) / 2. / nFreqsSam / nTimeSam; // real YX
								corrs[ (nBLines- (nAnts-nant1)*(nAnts-nant1+1)/2 + (nant2-nant1))*(nSubChans*nPols*2*2) + nc*(nPols*2*2) +5] += (float)(in1ry*in2ix - in1iy*in2rx) / 2. / nFreqsSam / nTimeSam; // imag YX
								corrs[ (nBLines- (nAnts-nant1)*(nAnts-nant1+1)/2 + (nant2-nant1))*(nSubChans*nPols*2*2) + nc*(nPols*2*2) +6] += (float)(in1ry*in2ry + in1iy*in2iy) / 2. / nFreqsSam / nTimeSam; // real YY
								corrs[ (nBLines- (nAnts-nant1)*(nAnts-nant1+1)/2 + (nant2-nant1))*(nSubChans*nPols*2*2) + nc*(nPols*2*2) +7] += (float)(in1ry*in2iy - in1iy*in2ry) / 2. / nFreqsSam / nTimeSam; // imag YY
							}
						}
					}
				}
			}
		}
		numwri = fwrite(corrs,sizeof(float),(int)(nSubChans*nBLines*nPols*2*2),write_ptr);
//		printf("sampl # %d -- written : %d\n", timeslo, (int)(numwri));
	}

	fclose(ptr);
	fclose(write_ptr);
}
