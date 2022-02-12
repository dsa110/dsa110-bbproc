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

#define NINT 30
#define NCHAN 384
#define NPACKET 2048
#define NANT 63
#define FILE_BYTES 5945425920

int antennas[64] = {24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 20, 19, 18, 17, 16, 15, 14, 13, 50, 51, 102, 116, 103, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 44, 45, 46, 115, 36, 37, 38, 39, 40, 41, 42, 43, 47, 48, 49, 100};

// output has shape [NCHAN, NINT, POL, R/I]
// input is in NINT-sized gulps. provide nT
void correlator(char *input, float *output, int a1, int a2, int nT) {

  float inr1, ini1;
  float inr2, ini2;
  float vr, vi;
  int idx1, idx2;

  for (int nF=0;nF<NCHAN;nF++) {
    for (int pol=0;pol<2;pol++) {
      
      vr=0.;
      vi=0.;
      
      for (int nP=0;nP<NPACKET;nP++) {
	for (int ntim=0;ntim<2;ntim++) {
	  
	  idx1 = nP*(NANT*NCHAN*2*2) + a1*NCHAN*2*2 + nF*2*2 + ntim*2 + pol;
	  idx2 = nP*(NANT*NCHAN*2*2) + a2*NCHAN*2*2 + nF*2*2 + ntim*2 + pol;
	  inr1 = (float)(((char)((input[idx1] & 15) << 4)) >> 4);
	  inr2 = (float)(((char)((input[idx2] & 15) << 4)) >> 4);
	  ini1 = (float)(((char)((input[idx1] & 240))) >> 4);
	  ini2 = (float)(((char)((input[idx2] & 240))) >> 4);
	  
	  vr += inr1*inr2 + ini1*ini2;
	  vi += inr1*ini2 - ini1*inr2;
	  
	}
      }
      
      output[nF*(NINT*2*2) + nT*2*2 + pol*2] = vr;
      output[nF*(NINT*2*2) + nT*2*2 + pol*2 + 1] = vi;
            
    }
  }

}

void usage()
{
  fprintf (stdout,
           "single_baseline [options]\n"
           " -d voltage data file name [no default]\n"
	   " -a antenna 1 (field ordering) [default 24]\n"
	   " -b antenna 2 (field ordering) [default 25]\n"
           " -h print usage\n");
}


int main (int argc, char *argv[]) {

  // read params : fch1, fnam, fdataname, sep
  int arg = 0;
  char * fnam=(char *)malloc(sizeof(char)*200);
  sprintf(fnam,"nofile");
  int a1 = 24;
  int a2 = 25;

  while ((arg=getopt(argc,argv,"d:a:b:h")) != -1)
    {
      switch (arg)
	{
	case 'd':
	  if (optarg)
	    {
	      strcpy(fnam,optarg);
	      break;
	    }
	  else
	    {
	      printf("-d flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'a':
	  if (optarg)
	    {
	      a1 = atoi(optarg);
	      break;
	    }
	  else
	    {
	      printf("-a flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }
	case 'b':
	  if (optarg)
	    {
	      a2 = atoi(optarg);
	      break;
	    }
	  else
	    {
	      printf("-b flag requires argument");
	      usage();
	      return EXIT_FAILURE;
	    }	  
	case 'h':
	  usage();
	  return EXIT_SUCCESS;
        }
    }

  // get antennas
  int i1, i2, aa;
  for (aa=0;aa<64;aa++) {
    if (antennas[aa]==a1) i1=aa;
    if (antennas[aa]==a2) i2=aa;
  }

  //printf("Antennas %d %d\n",i1,i2);
  
  // data
  char * input = (char *)malloc(sizeof(char)*FILE_BYTES/NINT);
  float *output = (float *)malloc(sizeof(float)*NINT*NCHAN*2*2);
  FILE *ptr;
  ptr = fopen(fnam,"rb");  // r for read, b for binary

  // calculate correlaton
  //printf("Correlating... ");
  for (int i=0;i<NINT;i++) {    
  
    fread(input,sizeof(char),FILE_BYTES/NINT,ptr);
    correlator(input, output, i1, i2, i);
    //printf("%d ",i+1);

  }
  //printf("\n");

  for (int i=0;i<NCHAN;i++) {
    for (int j=0;j<NINT;j++) {
      printf("%f %f %f %f\n",output[i*NINT*2*2 + j*2*2],output[i*NINT*2*2 + j*2*2 + 1],output[i*NINT*2*2 + j*2*2 + 2],output[i*NINT*2*2 + j*2*2 + 3]);
    }
  }

  fclose(ptr);
  free(input);
  free(output);
  free(fnam);

}
