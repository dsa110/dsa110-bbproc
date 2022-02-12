/*-o beamformer beamformer.c -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran
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

#include <src/sigproc.h>
#include <src/header.h>

#define nsamps 16384
#define offset 16384
#define nchan 768

FILE *output;

void send_string(char *string) /* includefile */
{
  int len;
  len=strlen(string);
  fwrite(&len, sizeof(int), 1, output);
  fwrite(string, sizeof(char), len, output);
}

void send_float(char *name,float floating_point) /* includefile */
{
  send_string(name);
  fwrite(&floating_point,sizeof(float),1,output);
}

void send_double (char *name, double double_precision) /* includefile */
{
  send_string(name);
  fwrite(&double_precision,sizeof(double),1,output);
}

void send_int(char *name, int integer) /* includefile */
{
  send_string(name);
  fwrite(&integer,sizeof(int),1,output);
}

void send_char(char *name, char integer) /* includefile */
{
  send_string(name);
  fwrite(&integer,sizeof(char),1,output);
}


void send_long(char *name, long integer) /* includefile */
{
  send_string(name);
  fwrite(&integer,sizeof(long),1,output);
}

void send_coords(double raj, double dej, double az, double za) /*includefile*/
{
  if ((raj != 0.0) || (raj != -1.0)) send_double("src_raj",raj);
  if ((dej != 0.0) || (dej != -1.0)) send_double("src_dej",dej);
  if ((az != 0.0)  || (az != -1.0))  send_double("az_start",az);
  if ((za != 0.0)  || (za != -1.0))  send_double("za_start",za);
}


// assume this runs on a single [nAnts, 384, 2times, 2pols] sample
// returns concatenated Tee-antenna spectra at 48 chan res. 

void extractor(char *input, unsigned char *opt, int *antennas, int nAnts) {

  float inr_x, ini_x;
  float rx;
  int ant;

  for (int nAnt=0;nAnt<nAnts;nAnt++) {
    ant = antennas[nAnt];
    for(int nChan=0;nChan<48;nChan++){

      rx = 0;
      
      for(int i=0;i<32;i++){
	inr_x = (float)(((char)((input[ant*48*32 + nChan*32 + i] & 15) << 4)) >> 4);
	ini_x = (float)(((char)((input[ant*48*32 + nChan*32 + i] & 240))) >> 4);
	rx += inr_x*inr_x + ini_x*ini_x; 
      }

      opt[nAnt*48+nChan] = (unsigned char)(rx/32.);
      
    }
  }

  
}

void usage()
{
  fprintf (stdout,
           "extract_antennas [options]\n"
           " -d voltage data file names x 16 [no default]\n"
	   " -a antennas file [no default]\n"
           " -h print usage\n");
}


int main (int argc, char *argv[]) {


  // read file names from command line
  FILE *fant;
  FILE **fin;
  fin = malloc( 16 * sizeof(FILE*) );

  for (int i=1;i<argc;i++) {

    if (strcmp(argv[i],"-a")==0) {
      fant=fopen(argv[i+1],"r");
    }

    if (strcmp(argv[i],"-d")==0) {
      for (int j=0;j<16;j++) fin[j] = fopen(argv[i+1+j],"rb");
    }

    if (strcmp(argv[i],"-h")==0) {
      usage();
      exit(0);
    }
    
  }

  // get antennas from file
  int *antennas = (int *)malloc(sizeof(int)*64);
  int nAnts=0;
  while (!feof(fant)) {
    fscanf(fant,"%d\n",&antennas[nAnts]);
    nAnts++;
  }
  fclose(fant);

  printf("read command-line args\n");

  // set up arrays
  unsigned char **opt = (unsigned char **)malloc(sizeof(unsigned char *)*nAnts);
  for (int i=0;i<nAnts;i++) opt[i] = (unsigned char *)malloc(sizeof(unsigned char)*nsamps*nchan);
  char *input = (char *)malloc(sizeof(char)*63*384*2*2);
  unsigned char *tout = (unsigned char *)malloc(sizeof(unsigned char)*nAnts*48);

  printf("allocated memory\n");

  // advance file pointers
  for (int i=0;i<16;i++) fseek(fin[i], offset*63*384*2*2, SEEK_SET);
  
  // fill output array

  for (int samp=0;samp<nsamps;samp++) {

    printf("sample %d of %d\n",samp+1,nsamps+1);
    
    for (int chgroup=0;chgroup<16;chgroup++) {

      // read data
      fread(input,sizeof(char),63*384*2*2,fin[chgroup]);

      // extract antennas into tout - [nAnts, spec]
      extractor(input,tout,antennas,nAnts);

      // place in output
      for (int ant=0;ant<nAnts;ant++) 
	memcpy(opt[ant] + samp*nchan + chgroup*48, tout + ant*48, 48);

    }
  }
  for (int j=0;j<16;j++) fclose(fin[j]);
  printf("extracted\n");

  // write to filterbanks
  char foutnam[200];
  for (int i=0;i<nAnts;i++) {
    sprintf(foutnam,"/home/ubuntu/tmp/output_%d.fil",i);
    output=fopen(foutnam,"wb");
    
    send_string("HEADER_START");
    send_string("source_name");
    send_string("TEST");
    send_int("machine_id",1);
    send_int("telescope_id",82);
    send_int("data_type",1); // filterbank data
    send_double("fch1",1498.75); // THIS IS CHANNEL 0 :)
    send_double("foff",-0.244140625);
    send_int("nchans",768);
    send_int("nbits",8);
    send_double("tstart",55000.0);
    send_double("tsamp",8.192e-6*8.);
    send_int("nifs",1);
    send_string("HEADER_END");

    fwrite(opt[i],sizeof(unsigned char),nsamps*nchan,output);
    fclose(output);
  }

  free(fin);
  free(antennas);
  free(opt);
  free(input);
  free(tout);
    
}
