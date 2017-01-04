#include "../include/Matematica.h"
#include <time.h>
void SigErr(int IsWrong,const char * s, ...){
#ifdef DEBUG
  if(IsWrong){
    va_list args;
    va_start(args, s);
    fprintf(stderr, "Error: ");
    vfprintf(stderr, s, args);
    fprintf(stderr, "exiting the program\n");
    va_end(args);
    abort();
  }
#else  
  return;
#endif
}
void abort_(const char * s, ...)
{
#ifdef DEBUG
  va_list args;
  va_start(args, s);
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
  abort();
#else
  return;
#endif
}
void Message(const char * s, ...)
{
#ifdef DEBUG
  va_list args;
  va_start(args, s);
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
#else
  return;
#endif
}
Matematica::Matematica(){
  time_t TempoSeme;
  (void)time(&TempoSeme);
  int seme = (int)(TempoSeme);
  //seme = 2351653;
  NPassi = (int)1000;
  PrecMinimo=1e-6;
#ifdef __GSL__
  rng = (gsl_rng *)malloc(sizeof(gsl_rng));
  rng = gsl_rng_alloc(gsl_rng_mt19937);
 //rng = gsl_rng_alloc(gsl_rng_tt800);
  gsl_rng_set(rng, seme);
#else
  time_t Now;
  time( &Now );
  unsigned long Seed = (unsigned long)Now;
  init_genrand( Seed ) ; 
  seme = (int)TempoSeme;
  //  printf("seme %d\n",seme);
#endif
  //SeGaussiano = 0;
  //fp = new double;
  Func = &Matematica::ContactAngle;
  //  ElabFunc = &Matematica::Derivata;
  Elab = &Matematica::Derivata;
  SpMem.a2 = 0.;
  SpMem.a1 = 0.;
  // NTemp=0;
  // Temp=0.;
  // casuale=0.;
  // segno=1;
  // DeltaGauss=0.;
  // SeGaussiano=0;
  //  matrix <double> F(4,4);
}
Matematica::~Matematica(){
#ifdef __GSL__
  gsl_rng_free(rng);
#endif
}
int Matematica::Factorial(int times){
  if(times <0)
    printf("the number is negative %d\n",times);
  int Resp=1;
  for(int t=1;t<times;t++){
    Resp *= t;
  }
  return Resp;
}
double Matematica::Binomial(int times,int n){
  double Resp=1.;
  Resp = (double)Factorial(times)/(double)(Factorial(n)*Factorial(times-n));
  return Resp;
}
int NumeroStringa(const char *NumeroC){
  int Cifra=0;
  int Numero=0;
  int FattoreDieci=1;
  int Potenza = strlen(NumeroC);
  for(int j=0;j<Potenza;j++){
    switch (NumeroC[j]){
    case '0':
      Cifra=0;
      break;
    case '1':
      Cifra=1;
      break;
    case '2':
      Cifra=2;
      break;
    case '3':
      Cifra=3;
      break;
    case '4':
      Cifra=4;
      break;
    case '5':
      Cifra=5;
      break;
    case '6':
      Cifra=6;
      break;
    case '7':
      Cifra=7;
      break;
    case '8':
      Cifra=8;
      break;
    case '9':
      Cifra=9;
      break;
    }
    //    cout << "Cifra " << Cifra << " Numero " << NumeroC[j] << " Potenza " << Potenza << endl;
    FattoreDieci=1;
    for(int i=Potenza-j;i>1;i--)
      FattoreDieci *=10;
    Numero+=Cifra*FattoreDieci;
  }
  return Numero;
}
