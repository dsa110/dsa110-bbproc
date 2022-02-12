#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define blocksize 94371840
#define nblocks 30

// arguments are: voltage file, blocks to rotate by

int main(int argc, char *argv[]) {

  FILE *fin, *fout;
  fin=fopen(argv[1],"rb");
  fout=fopen("output.dat","wb");

  long long int nbl = (long long int)(atoi(argv[2]));
  char * buf2 = (char *)malloc(sizeof(char)*nbl*blocksize);
  char * buf = (char *)malloc(sizeof(char)*(nblocks-nbl)*blocksize);

  printf("%lld %lld\n",nbl*blocksize,(nblocks-nbl)*blocksize);
  
  fread(buf2,sizeof(char),nbl*blocksize,fin);
  fread(buf,sizeof(char),(nblocks-nbl)*blocksize,fin);

  fwrite(buf,sizeof(char),(nblocks-nbl)*blocksize,fout);
  fwrite(buf2,sizeof(char),nbl*blocksize,fout);
  
  free(buf);
  free(buf2);
  fclose(fin);
  fclose(fout);

}
