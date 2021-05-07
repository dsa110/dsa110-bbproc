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
	
/*
0 0
0 1
0 2
0 3
0 4
0 5
0 6
0 7
0 8
0 9
0 10
0 11
0 12
0 13
0 14
0 15
0 16
0 17
0 18
0 19
0 20
0 21
0 22
0 23
1 1
1 2
1 3
1 4
1 5
1 6
1 7
1 8
1 9
1 10
1 11
1 12
1 13
1 14
1 15
1 16
1 17
1 18
1 19
1 20
1 21
1 22
1 23
2 2
2 3
2 4
2 5
2 6
2 7
2 8
2 9
2 10
2 11
2 12
2 13
2 14
2 15
2 16
2 17
2 18
2 19
2 20
2 21
2 22
2 23
3 3
3 4
3 5
3 6
3 7
3 8
3 9
3 10
3 11
3 12
3 13
3 14
3 15
3 16
3 17
3 18
3 19
3 20
3 21
3 22
3 23
4 4
4 5
4 6
4 7
4 8
4 9
4 10
4 11
4 12
4 13
4 14
4 15
4 16
4 17
4 18
4 19
4 20
4 21
4 22
4 23
5 5
5 6
5 7
5 8
5 9
5 10
5 11
5 12
5 13
5 14
5 15
5 16
5 17
5 18
5 19
5 20
5 21
5 22
5 23
6 6
6 7
6 8
6 9
6 10
6 11
6 12
6 13
6 14
6 15
6 16
6 17
6 18
6 19
6 20
6 21
6 22
6 23
7 7
7 8
7 9
7 10
7 11
7 12
7 13
7 14
7 15
7 16
7 17
7 18
7 19
7 20
7 21
7 22
7 23
8 8
8 9
8 10
8 11
8 12
8 13
8 14
8 15
8 16
8 17
8 18
8 19
8 20
8 21
8 22
8 23
9 9
9 10
9 11
9 12
9 13
9 14
9 15
9 16
9 17
9 18
9 19
9 20
9 21
9 22
9 23
10 10
10 11
10 12
10 13
10 14
10 15
10 16
10 17
10 18
10 19
10 20
10 21
10 22
10 23
11 11
11 12
11 13
11 14
11 15
11 16
11 17
11 18
11 19
11 20
11 21
11 22
11 23
12 12
12 13
12 14
12 15
12 16
12 17
12 18
12 19
12 20
12 21
12 22
12 23
13 13
13 14
13 15
13 16
13 17
13 18
13 19
13 20
13 21
13 22
13 23
14 14
14 15
14 16
14 17
14 18
14 19
14 20
14 21
14 22
14 23
15 15
15 16
15 17
15 18
15 19
15 20
15 21
15 22
15 23
16 16
16 17
16 18
16 19
16 20
16 21
16 22
16 23
17 17
17 18
17 19
17 20
17 21
17 22
17 23
18 18
18 19
18 20
18 21
18 22
18 23
19 19
19 20
19 21
19 22
19 23
20 20
20 21
20 22
20 23
21 21
21 22
21 23
22 22
22 23
23 23
*/
