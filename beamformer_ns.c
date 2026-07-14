/*gcc -o beamformer_ns beamformer_ns.c -I/usr/local/include -L/usr/local/lib -lm -g -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran
python beamformer was too slow, decided to use python to write up header etc but do actual beamforming in C.
This code should take parameters:
* data file name
* calibration file name
* number of antennas in voltage file
* number of antennas to use in beamforming
* start frequency
* separation (E/W)
* separation (N/S)
* beam number (0-255 for E/W, 256-511 for N/S)
* DEC (declination in degrees, required for N/S beams)
* output file name
assumes 48 channels for beamformer (weights for 8 data channels), ONLY 1 beam

Extended from beamformer.c to support N/S beams (256-511) using the second 48 antennas.
N/S beam weights calculated as in dsaX_bfCorr.cu.

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
float DSA_LAT = 37.23;  // DSA-110 latitude in degrees
int NANT = 96;


int init_weights(char * fnam, char *flagants, float *antpos, float *weights, int nPols) {

        // assumes NANT antennas
        // antpos: takes eastings (first NANT) and northings (second NANT)
        // weights: takes [ant, NW==48]

        FILE *fin;
	FILE *fflag;
        FILE *fants;
        int rd;

	int flags[NANT], nflag=0;
        fflag = fopen(flagants,"r");
	while (!feof(fflag)) {
	  fscanf(fflag,"%d\n",&flags[nflag]);
	  nflag++;
	}
	fclose(fflag);	 
	
        fin=fopen(fnam,"rb");

        rd = fread(antpos,2*NANT*sizeof(float),1,fin);
        rd = fread(weights,NANT*NW*nPols*2*sizeof(float),1,fin);
        float wnorm;
	int i;
	for (int ii=0;ii<NANT;ii++) {
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

void calc_weights(float *antpos, float *weights, float *freqs, float *wr, float *wi, float sep, float sep_ns, float nBeamNum, float dec, int nPols) {

        float theta, afac, twr, twi;
        int is_ns = (nBeamNum >= 256.0);  // N/S beam if >= 256
        float internal_beam = is_ns ? (nBeamNum - 256.0) : nBeamNum;
        int ant_offset = is_ns ? (NANT/2) : 0;  // Use antennas 48-95 for N/S, 0-47 for E/W

        if (is_ns) {
                // N/S beam: theta includes DEC offset, use sin(theta) in afac
                theta = sep_ns*(127.-internal_beam)*PI/10800. - (PI/180.)*(DSA_LAT - dec); // radians
                for(int nAnt=0;nAnt<NANT;nAnt++){
                        for(int nChan=0;nChan<48;nChan++){
                                for(int nPol=0;nPol<nPols;nPol++){
                                        afac = -2.*PI*freqs[nChan*8+4]*sin(theta)/CVAC; // factor for rotate (sin(theta) for N/S)
                                        // Use northing position for antennas 48-95: antpos[nAnt + NANT]
                                        // Note: antpos layout is [96 eastings, 96 northings]
                                        // So antpos[nAnt + NANT] = northing of antenna nAnt
                                        twr = cos(afac*antpos[nAnt + NANT]);
                                        twi = sin(afac*antpos[nAnt + NANT]);

                                        wr[nAnt*(48*nPols)+nChan*nPols+nPol] = (twr*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2] - twi*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2+1]);
                                        wi[nAnt*(48*nPols)+nChan*nPols+nPol] = (twi*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2] + twr*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2+1]);
                                }
                        }
                }
        } else {
                // E/W beam: original calculation
                theta = sep*(127.-internal_beam)*PI/10800.; // radians
                for(int nAnt=0;nAnt<NANT;nAnt++){
                        for(int nChan=0;nChan<48;nChan++){
                                for(int nPol=0;nPol<nPols;nPol++){
                                        afac = -2.*PI*freqs[nChan*8+4]*theta/CVAC; // factor for rotate
                                        // Use easting position: antpos[nAnt]
                                        twr = cos(afac*antpos[nAnt]);
                                        twi = sin(afac*antpos[nAnt]);

                                        wr[nAnt*(48*nPols)+nChan*nPols+nPol] = (twr*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2] - twi*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2+1]);
                                        wi[nAnt*(48*nPols)+nChan*nPols+nPol] = (twi*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2] + twr*weights[(nAnt*(48*nPols)+nChan*nPols+nPol)*2+1]);
                                }
                        }
                }
        }
}

void beamformer(char *input, float *wr, float *wi, unsigned char *output, int nChans, int nTimes, int nPols, int incoh, float nBeamNum) {

        float inr_x, ini_x, inr_y, ini_y;
        float wrx, wix, wry, wiy;
        float rx, ix, ry, iy;
        float tmprealX, tmpimagX, tmprealY, tmpimagY;
	char v;

        int is_ns = (nBeamNum >= 256.0);  // N/S beam if >= 256
        int ant_start = is_ns ? (NANT/2) : 0;  // Start at antenna 48 for N/S, 0 for E/W
        int ant_end = is_ns ? NANT : (NANT/2);  // End at antenna 96 for N/S, 48 for E/W

        for(int nTime=0;nTime<nTimes;nTime++){
                for(int nChan=0;nChan<48;nChan++){
                        for(int i=0;i<8;i++){
                                rx = 0;
                                ix = 0;
                                ry = 0;
                                iy = 0;
                                for(int nAnt=ant_start;nAnt<ant_end;nAnt++){
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

                                        /***********/
                                        /*towrite[nAnt*(nChans*nTimes*nPols*2)+(nChan*8+i)*(nPols*nTimes*2)+nTime*nPols*2  ] = (int)inr_x;
                                        towrite[nAnt*(nChans*nTimes*nPols*2)+(nChan*8+i)*(nPols*nTimes*2)+nTime*nPols*2+1] = (int)ini_x;
                                        towrite[nAnt*(nChans*nTimes*nPols*2)+(nChan*8+i)*(nPols*nTimes*2)+nTime*nPols*2+2] = (int)inr_y;
                                        towrite[nAnt*(nChans*nTimes*nPols*2)+(nChan*8+i)*(nPols*nTimes*2)+nTime*nPols*2+3] = (int)ini_y;*/
                                        /***********/

                                        wrx = wr[nAnt*(48*nPols)+nChan*nPols];
                                        wix = wi[nAnt*(48*nPols)+nChan*nPols];
                                        wry = wr[nAnt*(48*nPols)+nChan*nPols+1];
                                        wiy = wi[nAnt*(48*nPols)+nChan*nPols+1];

					if (!incoh) {
					  rx += inr_x*wrx - ini_x*wix;
					  ix += inr_x*wix + ini_x*wrx;
					  ry += inr_y*wry - ini_y*wiy;
					  iy += inr_y*wiy + ini_y*wry;
					}
					else {
					  rx += inr_x*inr_x*wrx*wrx + inr_y*inr_y*wry*wry + ini_x*ini_x*wix*wix + ini_y*ini_y*wiy*wiy;
					}

                                }
				if (!incoh)
				  output[nTime*nChans+nChan*8+i] = (unsigned char)(128.*16.*(rx*rx + ix*ix + ry*ry + iy*iy) / NANT / NANT);
				else
				  output[nTime*nChans+nChan*8+i] = (unsigned char)(64.*64.*rx/NANT/NANT);
                        }
                }
        }
}

void usage()
{
  fprintf (stdout,
           "beamformer_ns [options]\n"
           " -d voltage data file name [no default]\n"
           " -f calibration file name [no default]\n"
           " -o output file name [no default]\n"
           " -a number of antennas in file [default 30]\n"
           " -u number of antennas to be used [default 24]\n"
           " -z fch1 in MHz [default 1530]\n"
           " -s E/W interbeam separation in arcmin [default 1.4]\n"
           " -S N/S interbeam separation in arcmin [default 1.0]\n"
           " -n beam number [0-255 for E/W, 256-511 for N/S, default 127]\n"
           " -g DEC in degrees [required for N/S beams 256-511]\n"
           " -i incoherent beamforming\n"
           " -q flagants file [no default]\n"
           " -h print usage\n");
}


int main (int argc, char *argv[]) {


        int nChans = 384;
        int nPols = 2;
        int nTimes = 2;
	int incoh = 0;

        // read params : fch1, fnam, fdataname, sep
        int arg = 0;
        float fch1 = 1530.0;
        float sep = 1;
        float sep_ns = 1.0;  // N/S beam separation
        float nBeamNum = 127.;
        float dec = 0.0;  // DEC in degrees
        int dec_provided = 0;  // flag to check if DEC was provided
	int nUsedAnts;
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

        while ((arg=getopt(argc,argv,"d:f:o:a:u:z:s:S:n:g:q:ih")) != -1)
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
                        case 'S':
                        if (optarg)
                        {
                                sep_ns = atof(optarg);
                                break;
                        }
                        else
                        {
                                printf("-S flag requires argument");
                                usage();
                                return EXIT_FAILURE;
                        }
                        case 'n':
                        if (optarg)
                        {
                                nBeamNum = atof(optarg);
                                break;
                        }
                        else
                        {
                                printf("-n flag requires argument");
                                usage();
                                return EXIT_FAILURE;
                        }
                        case 'g':
                        if (optarg)
                        {
                                dec = atof(optarg);
                                dec_provided = 1;
                                break;
                        }
                        else
                        {
                                printf("-g flag requires argument");
                                usage();
                                return EXIT_FAILURE;
                        }
                        case 'h':
                        usage();
                        return EXIT_SUCCESS;
                        case 'i':
			  incoh=1;
			  break;
                }
        }

        // Validate: DEC is required for N/S beams (256-511)
        if (nBeamNum >= 256.0 && !dec_provided) {
                printf("Error: DEC (-g) is required for N/S beams (256-511)\n");
                usage();
                return EXIT_FAILURE;
        }

        // Validate beam number range
        if (nBeamNum < 0 || nBeamNum > 511) {
                printf("Error: Beam number must be between 0 and 511\n");
                usage();
                return EXIT_FAILURE;
        }

        // Print info about which array arm is being used
        if (nBeamNum >= 256.0) {
                printf("Using N/S array (antennas 48-95), beam %d (internal beam %d), DEC=%.2f deg\n", 
                       (int)nBeamNum, (int)(nBeamNum - 256), dec);
        } else {
                printf("Using E/W array (antennas 0-47), beam %d\n", (int)nBeamNum);
        }


        // compute beamformer weights
        //unsigned char * output = (char *)malloc(sizeof(char)*nChans*nTimes);
        unsigned char * output = (unsigned char *)malloc(sizeof(unsigned char)*nChans*nTimes);
        unsigned char * input = (char *)malloc(sizeof(char)*NANT*nChans*nTimes*nPols);
        float * antpos = (float *)malloc(sizeof(float)*NANT*2); // easting and northing
        float * weights = (float *)malloc(sizeof(float)*NANT*NW*nPols*2); // complex weights [ant, NW, pol, r/i]
        float * wr = (float *)malloc(sizeof(float)*NANT*NW*nPols); // complex weights [ant, NW, pol]
        float * wi = (float *)malloc(sizeof(float)*NANT*NW*nPols); // complex weights [ant, NW, pol]
        float * freqs = (float *)malloc(sizeof(float)*nChans); // freq
        for (int i=0;i<nChans;i++) freqs[i] = (fch1 - i*250./8192.)*1e6;
        init_weights(fnam,fflag,antpos,weights,nPols);
        calc_weights(antpos,weights,freqs,wr,wi,sep,sep_ns,nBeamNum,dec,nPols);

        FILE *ptr;
        FILE *write_ptr;
        ptr = fopen(fdata,"rb");  // r for read, b for binary
        write_ptr = fopen(fout,"wb");  // w for write, b for binary

        //long int sz;
        //fseek(ptr, 0L, SEEK_END);
        //sz = ftell(ptr);
        //rewind(ptr);
	//fseek(ptr, 369169920L, SEEK_SET);
        //int nTotSam = (int)(floor(sz / (nAnts*nChans*nTimes*nPols)));
	int nTotSam = 32768;
	
        int rd;
        for(int nSam = 0; nSam < nTotSam; nSam++) {

                rd = fread(input,NANT*nChans*nTimes*nPols,1,ptr);

                beamformer(input,wr,wi,output,nChans,nTimes,nPols,incoh,nBeamNum);

                fwrite(output,sizeof(unsigned char),nChans*nTimes,write_ptr);

        }
        fclose(ptr);
        fclose(write_ptr);


}
