#include <math.h>
#define NTAB 32
double ran1(void){
  int IA=16807,IM=2147483647,IQ=127773,IR=2836,NDIV=1+(IM-1)/NTAB,kcas,jcas;
  double casuale,AM,EPS,ENMX,RNMX;
  EPS=1.7*pow(10.,-7.);
  RNMX=1.-EPS;
  AM=1./IM;
  static int iv[NTAB];
  static int iy=0;
  static int seme = 2351653;
  if(seme<=0 || iy==0){
    seme = -seme>1 ? -seme : 1;
    for(jcas=NTAB+8;jcas>=1;jcas--){
      kcas=seme/IQ;
      seme=IA*(seme-kcas*IQ)-IR*kcas;
      if(seme<=0) seme += IM;
      if(jcas<=NTAB) iv[jcas]=seme;
    }
    iy = iv[1];
  }
  kcas=seme/IQ;
  seme=IA*(seme-kcas*IQ)-IR*kcas;
  if(seme<=0) seme+=IM;
  jcas=1+iy/NDIV;
  iy=iv[jcas];
  iv[jcas]=seme;
  casuale=AM*iy < RNMX ? AM*iy : RNMX;  
  return casuale;
}
