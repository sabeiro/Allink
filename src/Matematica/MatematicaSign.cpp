//#include "../include/Matematica.h"
#include <Matematica.h>

#ifdef __GSL__
#include <gsl/gsl_randist.h>
double Matematica::Casuale(){
  return gsl_rng_uniform(rng);
}
bool Matematica::InizializzaGaussiano(double Scarto,int N){
  return 0;
}
double Matematica::Gaussiano(double Media,double Scarto){//Gia inizializzato
  return Media + gsl_ran_gaussian(rng,Scarto);
}
#else
bool Matematica::InizializzaGaussiano(double Scarto,int N){
  zigset(86947731);
  return 0;
}
/**
   Mersenne Twister + Box-Muller transform
 */
double Matematica::Gaussiano(double Media,double Scarto){
  double Rand = genrand_real1() ;
  double Length = -log( Rand ) ;
  double Theta = 2. * M_PI * genrand_real2() ;
  double kdX = sqrt( 2. * Length ) * cos( Theta ) ;
  return   Scarto * kdX + Media ;
  return  RandomGaussian(Media,Scarto);
}
double Matematica::Casuale(){
  return genrand_real2();
  return RandomUniform();
  int IA=16807,IM=2147483647,IQ=127773,IR=2836,NDIV=1+(IM-1)/NTAB,kcas,jcas;
  double AM,EPS,RNMX;
  EPS=1.7*pow(10.,-7.);
  RNMX=1.-EPS;
  AM=1./IM;
  static int iv[NTAB];
  static int iy=10;
  static int seme = 2351653;
  //  printf("%d,%d\n",seme,iy);
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
  double casuale=AM*iy < RNMX ? AM*iy : RNMX;
  //  fprintf(CONTROLLA,"%f\n",casuale);
  return casuale;
}
#endif
double Matematica::RandDiscrProb(double *Prob,int NBin){
  double Ran = Casuale();
  double InvNBin = 1./(double)NBin;
  int j = 0;
  for(int i=0;i<NBin-1;i++){
    if(Ran > Prob[i] && Ran <= Prob[i+1]){
      j = i;
      break;
    }
  }
  double xi = Ran - floor(Ran*NBin)*InvNBin;
  double m = (Prob[j+1] - Prob[j])*NBin;
  double y = j*InvNBin + m*xi;
  return y;
}
#ifdef USE_FFTW
void Matematica::Spettro(double *st,double *sw,int NMax){//N=int^2
  //printf("using fftw\n");
  fftw_complex *uscita = (fftw_complex *) fftw_malloc( sizeof(fftw_complex)*(NMax) );
  fftw_complex *entrata = (fftw_complex *) fftw_malloc( sizeof(fftw_complex)*(NMax) );
  fftw_plan p;
  for(int i=0;i<NMax;i++){
    entrata[0][i] = st[i];
  }
  //for(int i=0;i<NMax;i++)printf("%d %lf \n",i,st[i]);  
  //p = fftw_plan_dft_r2c_1d(NMax,st,uscita,FFTW_FORWARD);
  p = fftw_plan_dft_1d(NMax,entrata,uscita,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p);
  for(int i=0;i<NMax;i++){
    sw[i] = QUAD(uscita[i][0])+QUAD(uscita[i][1]);
  }
  //for(int i=0;i<NMax;i++) printf("%d %lf %lf\n",i,st[i],sw[i]);  
  fftw_destroy_plan(p);
  fftw_free(entrata);
  fftw_free(uscita);
}
void Matematica::Spettro2d(double *st,double *sw,int NMax){//N=int^2
  //printf("using fftw\n");
  fftw_complex *uscita = (fftw_complex *) fftw_malloc( sizeof(fftw_complex)*(NMax*NMax) );
  fftw_plan p = fftw_plan_dft_r2c_2d(NMax,NMax,st,uscita,FFTW_MEASURE);
  fftw_execute(p);
  for(int i=0;i<NMax;i++){
    sw[i] = 0.;
    for(int j=0;j<NMax;j++){
      sw[i*NMax+j] = SQR(uscita[i*NMax+j][0])+SQR(uscita[i*NMax+j][0]);
    }
  }
  fftw_destroy_plan(p);
  fftw_free(uscita);
}
void Matematica::Spettro2d(double *st,double **sw,int NMax){//N=int^2
  //printf("using fftw\n");
  double NMaxInv = 1./(double)NMax;
  int NHalf = (int)(NMax/2.);
  fftw_plan p;
  //fftw_complex *uscita = (fftw_complex *) fftw_malloc( sizeof(fftw_complex)*(NMax*NMax) );
  fftw_complex *uscita = (fftw_complex *)calloc(NMax*NMax,sizeof(fftw_complex));
  //for(int i=0;i<NMax;i++)for(int j=0;j<NMax;j++)printf("%d %d %lf \n",i,j,st[i*NMax+j]);
  p = fftw_plan_dft_r2c_2d(NMax,NMax,st,uscita,FFTW_MEASURE);
  fftw_execute(p);
  for(int i=0;i<NMax;i++){
    for(int j=0;j<NMax;j++){
      sw[i][j] = QUAD(uscita[0][i*NMax+j])+QUAD(uscita[1][i*NMax+j]);
      sw[i][j] *= QUAD(NMaxInv);
      //printf("%d %d %lf %lf %lf\n",i,j,st[i*NMax+j],uscita[0][i*NMax+j],uscita[1][i*NMax+j]);
    }
  }
  fftw_destroy_plan(p);
  //fftw_free(uscita);
  free(uscita);
}
#else
void Matematica::Spettro(double *st,double *sw,int NMax){//N=int^2
  double Re1=0.,Im1=0.,Re2=0.,Im2=0.;
  double FMass=(double)NMax/DUE_PI;
  double dNMax = 1./(double)NMax;
  for(int j=0;j<NMax;j++){
    Re1=0.;Re2=0.;Im1=0.;Im2=0.;
    for(int i=0;i<NMax;i++){
      Re1+=st[i]*cos(DUE_PI*dNMax*(i)*j);
      Im1-=st[i]*sin(DUE_PI*dNMax*(i)*j);
    }
    sw[j]= sqrt(QUAD((Re1+Re2))*dNMax+QUAD((Im1+Im2))*dNMax);
  }
  // double Parseval1 = 0.;
  // double Parseval2 = 0.;
  // for(int j=0;j<NMax;j++){
  //   Parseval1 += QUAD(st[j]);
  //   Parseval2 += sw[j];
  // }
  // printf("Parseval %lf=%lf\n",Parseval1,Parseval2);
}
void Matematica::Spettro2d(double *st,double *sw,int NMax){
  printf("dft\n");
  double dNMax = 1./(double)NMax;
  int NHalf = (int)(NMax/2.);
  for(int kx=-NMax/2;kx<NMax/2;kx++){
    double qx = kx*dNMax;
    for(int ky=-NMax/2;ky<NMax/2;ky++){
      double qy = ky*dNMax;
      double Re2=0.,Im2=0.;
      double Re1=0.,Im1=0.;
      for(int lx=0;lx<NMax;lx++){
	for(int ly=0;ly<NMax;ly++){
	  double Arg = 2.*M_PI*(lx*kx + ly*ky);
	  double cy = cos(Arg);
	  double sy = sin(Arg);
	  Re1 += st[lx*NMax+ly]*cy;
	  Im1 += st[lx*NMax+ly]*sy;
	}
      }
      int kkx = kx + NMax/2;
      int kky = ky + NMax/2;
      sw[kkx*NMax+kky] = SQR(Re1*dNMax) + SQR(Im1*dNMax);
    }
  }
}
// void Matematica::Spettro2d(double *st,double *sw,int NMax){
//   double dNMax = 1./(double)NMax;
//   int NHalf = (int)(NMax/2.);
//   for(int kx=0;kx<NMax;kx++){
//     for(int ky=0;ky<NMax;ky++){
//       double Re2=0.,Im2=0.;
//       for(int lx=0;lx<NMax;lx++){
// 	double cx = cos(kx*lx*dNMax*DUE_PI);
// 	double sx = sin(kx*lx*dNMax*DUE_PI);
// 	double Re1=0.,Im1=0.;
// 	for(int ly=0;ly<NMax;ly++){
// 	  double cy = cos(ky*ly*dNMax*DUE_PI);
// 	  double sy = sin(ky*ly*dNMax*DUE_PI);
// 	  Re1 += st[lx*NMax + ly]*cy;
// 	  Im1 += st[lx*NMax + ly]*sy;
// 	}
// 	Re2 += cx*Re1 - sx*Im1;
// 	Im2 -= sx*Re1 + cx*Im1;
//       }
//       sw[kx*NMax+ky] = SQR(Re2*dNMax) + SQR(Im2*dNMax);
//     }
//   }
// }
void Matematica::Spettro2d(double *st,double **sw,int NMax){
  printf("No fftw\n");
  double dNMax = 1./(double)NMax;
  int NHalf = (int)(NMax/2.);
  for(int l=0;l<NMax;l++){
    for(int k=0;k<NMax;k++){
      double Re1=0.,Im1=0.,Re2=0.,Im2=0.;
      for(int n=0;n<NMax;n++){
	double Cosn = cos(DUE_PI*l*n*dNMax);
	double Sinn = sin(DUE_PI*l*n*dNMax);
	for(int m=0;m<NMax;m++){
	  double Arg = DUE_PI*((l-NHalf)*n*dNMax + (k-NHalf)*m*dNMax);
	  // Re1 += cos(Arg)*st[n*NMax + m];
	  // Im2 += sin(Arg)*st[n*NMax + m];
	  Re1 += Cosn*cos(DUE_PI*k*m*dNMax)*st[n*NMax + m];
	  Re2 -= Sinn*sin(DUE_PI*k*m*dNMax)*st[n*NMax + m];
	  Im1 += Cosn*sin(DUE_PI*k*m*dNMax)*st[n*NMax + m];
	  Im2 += Sinn*cos(DUE_PI*k*m*dNMax)*st[n*NMax + m];
	}
      }
      sw[l][k] = QUAD((Re1+Re2)*dNMax*dNMax) + QUAD((Im1+Im2)*dNMax*dNMax);
    }
  }
  // double Parseval1 = 0.;
  // double Parseval2 = 0.;
  //  for(int v=0;v<NMax;v++){
  //   for(int vv=0;vv<NMax;vv++){
  //     Parseval1 += QUAD(sw[v][vv]);
  //     Parseval2 += QUAD(st[v*NMax+vv]);
  //   }
  // }
  //  printf("Parseval %lf=%lf\n",Parseval1,Parseval2);
}
#endif
void Matematica::Radice(double *st,double *sw,int N){
  for(int j=0;j<N;j++){
    if(st[j] > 0.)
      sw[j] = sqrt(st[j]);
    else 
      sw[j] = sqrt(-st[j]);
  }
}
void Matematica::Autocor(double *st,double *sAutocor,int NMax){
  for(int i=0;i<NMax;i++){
    sAutocor[i] = 0.;
    for(int j=0,k=0;j<NMax;j++){
      k = i-j;
      if(k<0)
	k = NMax - i + j;
      sAutocor[i] += (double)(st[k]*st[j]);
    }
    sAutocor[i] /= (double)(NMax);
  }
}
void Matematica::MediaMobile(double *st,int NMax,double *sw,int Parti){
  if(Parti <= 1){
    return;
  }
  double PartiInv = 1./(double)Parti;
  int Nw = (int) (NMax/(double)Parti);
  for(int i=0;i<Nw;i++)
    sw[i] = 0.;
  for(int i=0,j=0;i<NMax;i++){
    sw[j] += st[i];
    if( (i%Parti)==0 ){
      sw[j] *= PartiInv;
      j++;
    }
  }
}
int Matematica::MediaMobile(double *st,int NMax,double *sw,double *sErr,int NParti){
  if(NParti <= 0 ){
    return NParti;
  }
  int Nw = (int) (NMax/(double)NParti);
  double InvNParti = 1./(double)NParti;
  for(int i=0;i<Nw;i++){
    sw[i] = 0.;
    sErr[i] = 0.;
  }
  for(int i=0,j=0;i<NMax;i++){
    sw[j] += st[i];
    sErr[j] += SQR((st[i]));
    if( ((i+1)%NParti)==0 ){
      sw[j] *= InvNParti;
      sErr[j] = sqrt( (sErr[j] - SQR(sw[j])*NParti)*InvNParti );
      j++;
      if(j==Nw) return Nw;
    }
  }
  return Nw;
}
int Matematica::CorrelaDuePunti(double *st,int NMax,double *sw,int Punti){
  int PuntiMass = NMax - Punti;
  for(int i=0;i<PuntiMass;i++){
    sw[i] = (st[i] + st[i+Punti]) / 2.;
  }
  return PuntiMass;
}
void Matematica::Autosimilarita(double *st,int NMax,double *sw,int Potenze){
  for(int j=0;j<Potenze;j++){
    sw[j] = 0.;
  }
  for(int i=0;i<NMax;i++){
    double Temp = 1;
    for(int j=0;j<Potenze;j++){
      Temp *= st[i];
      sw[j] += Temp;
    }
  }
  for(int j = 0;j<Potenze;j++){
    sw[j] = log10(POS(sw[j]))/log10((double)NMax);
    //    printf("sw[%d] %f\n",j,sw[j]);
  }
}
int Matematica::Normalizza(double *st,double *sw,int NMax){
  double yMass=st[0];
  double yMin = st[0];
  double Media=0.;
  int Zeri=0;
  int Campioni=0;
  for(int i = 1;i<NMax;i++){
    if(yMin > st[i])
      yMin = st[i];
    if(yMass< st[i])
      yMass = st[i];
    Media += st[i];
  }
  Media /= (double)NMax;
  for(int i=0;i<NMax;i++){
    sw[i] = (st[i]-yMin)/(yMass-yMin);
    if(st[i] > Media && st[i-1] < Media)
      Zeri++;
  }
  Campioni = (int)(NMax/(double)Zeri);
  return Campioni;
}
int Matematica::Normalizza(double *st,int NMax){
  double yMass=st[0];
  double yMin = st[0];
  double Media=0.;
  int Zeri=0;
  int Campioni=0;
  for(int i = 1;i<NMax;i++){
    if(yMin > st[i])
      yMin = st[i];
    if(yMass< st[i])
      yMass = st[i];
    Media += st[i];
  }
  Media /= (double)NMax;
  for(int i=0;i<NMax;i++){
    st[i] = (st[i]-yMin)/(yMass-yMin);
    if(st[i] > Media && st[i-1] < Media)
      Zeri++;
  }
  Campioni = (int)(NMax/(double)Zeri);
  return Campioni;
}
int Matematica::NormalizeArea(double *st,int NMax){
  double Area = 0.;
  for(int i = 0;i<NMax;i++){
    Area += st[i];
  }
  for(int i=0;i<NMax;i++){
    st[i] /= Area;
  }
  return 0;
}
MOMENTI Matematica::Distribuzione(const double *st,int NMax){
  MOMENTI m1; m1.Uno=0.; m1.Due=0.; m1.Tre=0.;m1.Delta=0.;m1.Num=0;
  m1.Min = st[0];
  m1.Max = st[0];
  for(int i=0;i<NMax;i++){
    m1.Uno += st[i];
    if(st[i] < m1.Min) m1.Min = st[i];
    if(st[i] > m1.Max) m1.Max = st[i];
  }
  m1.Uno /= (double)(NMax);
  m1.Delta = (m1.Max-m1.Min);
  for(int i=0;i<NMax;i++){
    m1.Due+=QUAD((st[i]-m1.Uno));
    m1.Tre+=QUAD((st[i]-m1.Uno))*(st[i]-m1.Uno);
    m1.Num++;
  }
  m1.Due=sqrt(m1.Due/(double)(NMax-1));
  m1.Tre = pow(m1.Tre / (double)(NMax - 2),.33333);
  return m1;
}
MOMENTI Matematica::DistrErr(const double *st,int NMax,double *Distr,double *Err,int NBin,double *Confine,int IfNorm){
  MOMENTI m1 = Distribuzione(st,NMax,Distr,NBin,Confine,IfNorm);
  for(int i=0;i<NBin;i++)
    Err[i] = 0.;
  double Add = 1./(double)NMax;
  for(int i=0;i<NMax;i++){
    int j = (int)((st[i]-Confine[0])/(Confine[1]-Confine[0])*NBin);
    if(j >= NBin || j < 0) continue;
    Err[j] += Add;
  }
  for(int i=0;i<NBin;i++){
    Err[i] *= Distr[i];
  }
  return m1;
}
MOMENTI Matematica::Distribuzione(const double *st,int NMax,double *Distr,int NBin,double *Confine,int IfNorm){
  MOMENTI m1; m1.Uno=0.; m1.Due=0.; m1.Tre=0.;m1.Delta=0.;m1.Num=0;
  m1.Min = Confine[0];
  m1.Max = Confine[1];
  for(int i=0;i<NBin;i++) Distr[i]=0.;
  for(int i=0;i<NMax;i++) m1.Uno += st[i];
  m1.Uno/=(double)(NMax);
  m1.Delta = (m1.Max - m1.Min)/(double)NBin;
  double Add = !IfNorm ? 1. : 1./(double)NMax;
  for(int i=0;i<NMax;i++){
    m1.Due += QUAD(st[i]-m1.Uno);
    m1.Tre += m1.Due*(st[i]-m1.Uno);
    int j = (int)((st[i]-Confine[0])/(Confine[1]-Confine[0])*NBin);
    if(j >= NBin || j < 0) continue;
    Distr[j] += Add;
    m1.Num++;
  }
  m1.yMin = Distr[0];
  m1.yMax = Distr[0];
  for(int v=0;v<NBin;v++){
    if(m1.yMin > Distr[v]) m1.yMin = Distr[v];
    if(m1.yMax < Distr[v]) m1.yMax = Distr[v];
  }
  m1.Due=sqrt(m1.Due/((double)(NMax-1)));
  m1.Tre = pow(m1.Tre,.33333) / (double)(NMax - 2);
  //   printf("[-]    Massimo %g m1.Minimo %g Media %g          [-] \n[-]    Scarto %g Terzo %g Deltax/3 %g      [-]\n",m1.Massimo,m1.Minimo,m1.Uno,m1.Due,m1.Tre,(m1.Massimo-m1.Minimo)/6.);
  return m1;
}
MOMENTI Matematica::Distribuzione(const double *st,int NMax,double *Distr,int NBin,int IfNorm){
  MOMENTI m1; m1.Uno=0.; m1.Due=0.; m1.Tre=0.;m1.Delta=0.;m1.Num=0;
  m1.Min = st[0];
  m1.Max = st[0];
  double Sum2 = 0.;
  double Sum3 = 0.;
  for(int i=0;i<NBin;i++) Distr[i] = 0.;
  for(int i=0;i<NMax;i++){
    if(st[i] < m1.Min) m1.Min = st[i];
    if(st[i] > m1.Max) m1.Max = st[i];
    m1.Uno += st[i];
    Sum2 += QUAD(st[i]);
    Sum3 += st[i]*st[i]*st[i];
  }
  double Add = !IfNorm ? 1. : 1./(double)NMax;
  m1.Delta = (m1.Max-m1.Min)/(double)NBin;
  for(int i=0;i<NMax;i++){
    int j = (int)((st[i]-m1.Min)/(m1.Max-m1.Min)*NBin);
    if(j >= NBin) continue;
    Distr[j] += Add;
    m1.Num++;
  }
  m1.yMin = Distr[0];
  m1.yMax = Distr[0];
  for(int v=0;v<NBin;v++){
    if(m1.yMin > Distr[v]) m1.yMin = Distr[v];
    if(m1.yMax < Distr[v]) m1.yMax = Distr[v];
  }
  m1.Uno/=(double)(NMax);
  m1.Due = sqrt((Sum2 - QUAD(m1.Uno)*NMax)/((double)NMax-1.));
  m1.Tre = Sum3-3.*Sum2*m1.Uno+3.*NMax*QUAD(m1.Uno)-NMax*m1.Uno;
  m1.Tre = pow( m1.Tre , 1./3.)/(double)(NMax-1);
  //   printf("[-]    Massimo %g m1.Minimo %g Media %g          [-] \n[-]    Scarto %g Terzo %g Deltax/3 %g      [-]\n",m1.Massimo,m1.Minimo,m1.Uno,m1.Due,m1.Tre,(m1.Massimo-m1.Minimo)/6.);
  return m1;
}
void Matematica::DistrSample(double *Px,double *Py,int NMax,double **Distr,int NBin,const int NSample,int IfNorm,double *yBorder){
  double xBorder[2];
  double yDelta;
  double Norm[NSample];
  xBorder[1] = Px[0];
  xBorder[0] = Px[0];
  for(int n=0;n<NMax;n++){
    if(xBorder[1] < Px[n]) xBorder[1] = Px[n];
    if(xBorder[0] > Px[n]) xBorder[0] = Px[n];
  }
  double Dx = (xBorder[1] - xBorder[0]) > 0. ? 1./(xBorder[1] - xBorder[0]) : 1.;
  if(1==0){//search for borders
    for(int s=0;s<NSample;s++){
      yBorder[0] = Py[s];
      yBorder[1] = Py[s];
    }
    for(int n=0;n<NMax;n++){
      int sx = (int)((Px[n]-xBorder[0])*Dx*NSample);
      if(sx < 0 || sx >= NSample) continue;
      if(yBorder[0] > Py[n]) yBorder[0] = Py[n];
      if(yBorder[1] < Py[n]) yBorder[1] = Py[n];    
    }
  }
  yDelta = yBorder[1] - yBorder[0] > 0. ? 1./(yBorder[1] - yBorder[0]) : 1.;
  for(int n=0;n<NMax;n++){
    int sx = (int)((Px[n]-xBorder[0])*Dx*NSample);
    if(sx < 0 || sx >= NSample) continue;
    int by = (int)((Py[n]-yBorder[0])*yDelta*NBin);
    if(by < 0 || by >= NBin) continue;
    Distr[sx][by] += 1.;
  }
  if(IfNorm){
    for(int sx=0;sx<NSample;sx++){
      Norm[sx] = 0.;
      for(int by=0;by<NBin;by++){
	Norm[sx] += Distr[sx][by];
      }
    }
    for(int sx=0;sx<NSample;sx++){
      Norm[sx] = Norm[sx] > 1. ? 1./Norm[sx] : 1.;
      for(int by=0;by<NBin;by++){
	Distr[sx][by] *= Norm[sx];
      }
    }
  }
}
MOMENTI Matematica::DistribuzioneGauss(const double *st,int NMax,double *Distr,double *dInt,int NBin,int IfNorm){
  MOMENTI m1;
  m1 = Distribuzione(st,NMax,Distr,NBin,IfNorm);
  m1.Chi = 0.;
  for(int v=0;v<NBin;v++){
    dInt[v] = NMax*Gauss( m1.Uno , m1.Due , (v)*m1.Delta + m1.Min)*m1.Delta;
    m1.Chi += QUAD( (Distr[v] - dInt[v])/m1.Due);
  }
  m1.Chi /= NBin-2;
  return m1;
}
MOMENTI Matematica::DistribuzioneMaxwell(const double *st,int NMax,double *Distr,double *dInt,int NBin,int IfNorm){
  MOMENTI m1;
  m1 = Distribuzione(st,NMax,Distr,NBin,IfNorm);
  m1.Chi = 0.;
  for(int i=0;i<NBin;i++){
    dInt[i] =  2*DUE_PI*QUAD((i*m1.Delta + m1.Min))*NMax*Gauss( m1.Uno , m1.Due , i*m1.Delta + m1.Min)*m1.Delta;
    m1.Chi += QUAD( Distr[i] - dInt[i]);
  }
  m1.Chi /= NBin-2;
  return m1;
}
MOMENTI Matematica::WeightAverage(const double *sx,const double *sy,int NMax){
  MOMENTI m1; m1.Uno=0.; m1.Due=0.; m1.Tre=0.;m1.Delta=0.;m1.Num=0;
  m1.Min = sx[0];
  m1.Max = sx[0];
  double TotWei = 0.;
  for(int i=0;i<NMax;i++){
    if(sx[i] < m1.Min) m1.Min = sx[i];
    if(sx[i] > m1.Max) m1.Max = sx[i];
    m1.Uno += sx[i]*sy[i];
    TotWei += sy[i];
    m1.Due += sy[i] > 0. ? 1./sy[i] : 0.;
  }
  m1.Uno /= TotWei;
  m1.Due = sqrt(1./m1.Due);
  m1.Delta = (m1.Max-m1.Min);
  return m1;
}
void Matematica::WeightHisto(double **hist,double *Border,int NBin,int NHisto,double tolerance,double *OrPos,double *kSpring){
  if(1==0) {
    cout<<"Usage:\n\n";
    cout<<"do_wham_c <TEMPERATURE> <MIN> <MAX> <BINS> <TOLERANCE> <OUTPUT> <INPUT1> <INPUT2> ....\n\n";
    cout<<"TEMPERATURE\t\tTemperature that simulations were performed at\n";
    cout<<"MIN\t\t\tMinimum value used in histogram. All histograms must have been computed with matching value.\n";
    cout<<"MAX\t\t\tMaximum value used in histogram. All histograms must have been computed with matching value.\n";
    cout<<"NBIN\t\t\tNumber of NBin in histogram. All histograms must have been computed with matching value.\n";
    cout<<"TOLERANCE\t\tTolerance value used to determine convergence. Value of 0.00001 seems to work fine.\n";
    cout<<"OUTPUT\t\t\tName of output file.\n";
    cout<<"INPUT\t\t\tNames of input files. These are histogram files produced from the .pdo files using make_histo.pl.\n";
    cout<<endl;
  }
  double temperature = 1.;
  double min         = Border[0];
  double max         = Border[1];
  double RT = temperature * 8.314472e-3;
  ofstream output("WeightHistoAnal.dat");
  double* F      = new double[NHisto];
  double* old_F  = new double[NHisto];
  double* result = new double[NBin];
  double*   N    = new double[NHisto];
  double* energy = new double[NBin];
  double * numerator = new double[NBin];
  double ** expU  = new double *[NHisto];
  double ** NexpU = new double *[NHisto];
  //filling the exponential and normalizing factor
  for(int h = 0;h < NHisto;++h) {
    F[h] = 1.0;
    old_F[h] = 0.0;
    N[h] = 0.0;
    expU[h]  = new double[NBin];
    NexpU[h] = new double[NBin];
    for(int i = 0; i < NBin; ++i) {
      N[h] += hist[h][i];
      double pos = (i + 0.5)/NBin * (max - min) + min;
      expU[h][i] = exp(-(0.5 * kSpring[h] * (pos - OrPos[h]) * (pos - OrPos[h]))/RT);
    }
    for(int i=0; i < NBin; ++i) {
      NexpU[h][i] = N[h] * expU[h][i];
    }
  }
  for(int bin = 0; bin < NBin; ++bin) {
    result[bin] = 0.0;
  }
  double d_max = 0.;
  for(int bin = 0; bin < NBin; ++bin) {
    numerator[bin] = 0.0;
    for(int file = 0; file < NHisto; ++file) {
      numerator[bin] += hist[file][bin];
    }
  }
  cerr<<"Done loading histograms. Doing WHAM."<<endl;
  int count = 0;
  int IfContinue = 1;
  do {
    d_max = 0.;
    for(int bin=0; bin < NBin; ++bin) {
      double denom = 0.0;
      for(int file = 0; file < NHisto; ++file) {
	denom += NexpU[file][bin] / F[file];
      }
      //printf("%lf %lf %lf\n",denom,numerator[bin],result[bin]);
      if(denom > 1.e9) IfContinue = 0;
      denom = denom <= 0.0000001 ? 1. : denom;
      result[bin] = numerator[bin] / denom;
    }
    double norm = 0.0;
    for(int bin=0; bin < NBin; ++bin) norm += result[bin];
    norm = norm > 0. ? norm : 1.;
    for(int bin=0; bin < NBin; ++bin) result[bin] /= norm;
    for(int file = 0; file < NHisto; ++file) {
      double Z = 0.;
      for (int bin = 0; bin < NBin; ++bin){
	Z += result[bin] * expU[file][bin];
      }
      if(isnan(Z)){IfContinue = 0;break;}
      if(isinf(Z)){IfContinue = 0;break;}
      if(Z <= 1.e-9){IfContinue = 0;break;}
      double dZ = (max - min)/(double)NBin;
      Z = Z * dZ;
      old_F[file] = F[file];
      //printf("%d %d %lf %lf\n",count,file,F[file],Z);
      F[file] = Z;
    }
    if(!IfContinue) break;
    double temp = F[0];
    for(int file = 0; file < NHisto; ++file) {
      temp = RT * fabs(log(F[file]) - log(old_F[file]));
      if(temp > d_max)
	d_max = temp;
    }
    if(!(count % 1000)) 
      cerr<<"Step: "<<count<<"\tTolerance: "<<d_max<<endl;
    ++count;
    if(count == 2) break;
  } while(d_max > tolerance);
  for(int bin = 0; bin < NBin; ++bin) {
    energy[bin] = 0.0;
  }
  for(int bin = 0; bin < NBin; ++bin) {
    double res = result[bin] > 0. ? log(result[bin]) : 0.;
    energy[bin] = -RT*res;
    double pos = (double)bin / NBin * (max - min) + min;
    output << pos << "  " << result[bin] << "  " << energy[bin] << endl;
  }
  output.close();
  delete [] F;
  delete [] old_F;
  delete [] result;
  delete [] N;
  delete [] energy;
  delete [] numerator;
  for(int current = 0;current < NHisto;++current) {
    delete [] expU[current];
    delete [] NexpU[current];
  }
  delete [] expU;
  delete [] NexpU;
}
void Matematica::Sort(double *Sign,int NMax){
  for(int i=0;i<NMax;i++){
    for(int j=i;j>0;j--){
      if(Sign[j] < Sign[j-1])
	Swap(j,j-1,Sign);
      else 
	break;
    }
  }
}
void Matematica::Swap(int i, int j,double *Sign){
  double Temp = Sign[j];
  Sign[j] = Sign[i];
  Sign[i] = Temp;
}
void Matematica::Swap(double *s,int si,double *t,int ti,const int NDim){
  double Temp[NDim];
  for(int d=0;d<NDim;d++){
    Temp[d] = s[si*NDim+d];
    s[si*NDim+d] = t[ti*NDim+d];
    t[ti*NDim+d] = Temp[d];
  }
}
void Matematica::Sort(int *Sign,int NMax){
  for(int i=0;i<NMax;i++){
    for(int j=i;j>0;j--){
      if(Sign[j] < Sign[j-1])
	Swap(j,j-1,Sign);
      else 
	break;
    }
  }
}
void Matematica::Swap(int i, int j,int *Sign){
  int Temp = Sign[j];
  Sign[j] = Sign[i];
  Sign[i] = Temp;
}
void Matematica::ConvWeight(double *st,int NMax,double *sw,int *WIndex,int NWeight){
  double *st1 = (double *)calloc(NMax,sizeof(double));
  for(int n=0;n<NMax;n++){
    st1[n] = st[n];
  }
  for(int n=0;n<NMax;n++){
    st[n] = 0;
    if(n+WIndex[0] < 0 || n+WIndex[NWeight-1] >= NMax){
      st[n] = st1[n];
      continue;
    }
    for(int w=0;w<NWeight;w++){
      st[n] += sw[w]*st1[n+WIndex[w]];
    }
  }
  free(st1);
}
void Matematica::FillWeightGauss(double *st,int *WIndex,int NWeight,double CutOff,double Sigma){
  int Half = NWeight/2;
  double Norm = 0.;
  for(int w=0;w<NWeight;w++){
    WIndex[w] = w-Half;
    double x = CutOff*(w-Half)/(double)NWeight;
    double r2 = SQR(x);
    double Gauss = exp(-r2*.5/SQR(Sigma));
    st[w] = Gauss;
    Norm += Gauss;
  }
  for(int w=0;w<NWeight;w++){
    st[w] /= Norm;
  }
}
void Matematica::ExecCommand(double *st,double *sw,int NMax,char *cmd){
  for(int i=0;i<NMax;i++){
    sw[i] = ExecFormula(st[i],1.,cmd);
  }
}
double Matematica::ExecFormula(double x,double y,char *cmd){
  if(strlen(cmd) == 0) return x;
  const int NOp = 20;
  int PosOp[NOp];
  int PosPar[NOp];
  int PosNum[2*NOp];
  PosNum[0] = 0;
  int nOp = 1;
  int Comm[NOp];
  double Numb;
  char *cPart = cmd;
  for(int i=0;i<NOp;i++){
    cPart = strpbrk(cPart+1,"*/^+-");
    if(cPart == NULL) break;
    Comm[i] = (int)*cPart;
    PosOp[i] = cPart-cmd+1;
    PosNum[nOp] = PosOp[i];
    nOp++;
  }
  Sort(PosNum,nOp);
  PosNum[nOp++] = strlen(cmd);
  double Res = 0.;
  for(int i=0;i<nOp-1;i++){
    if( *(cmd+PosNum[i]) == 'x'){
      Numb = x;
    }
    if( *(cmd+PosNum[i]) == 'y'){
      Numb = y;
    }
    else
      sscanf(cmd+PosNum[i],"%lf",&Numb);
    if(i == 0){ 
      Res = Numb;
      continue;
    }
    switch(Comm[i-1]){
    case '+':
      Res += Numb;
      break;
    case '-':
      Res -= Numb;
      break;
    case '/':
      Res /= Numb;
      break;
    case '*':
      Res *= Numb;
      break;
    }
    //printf("%d) %s %d %lf %lf %lf %c %lf\n",i,cmd+PosNum[i],PosNum[i],Numb,x,y,Comm[i],Res);
  }
  return Res;
  // cPart = cmd;
  // for(int i=0;i<NOp;i++){
  //   cPart = strpbrk(cPart+1,"()");
  //   if(cPart == NULL) break;
  //   PosPar[i] = cPart-cmd+1;
  //   PosNum[nOp] = PosOp[i];
  //   nOp++;
  // }
  // 
}
double Matematica::ExecFormula(double **st,int n,char *cmd){
  if(strlen(cmd) == 0) return st[0][n];
  const int NOp = 20;
  int PosOp[NOp];
  int PosPar[NOp];
  int PosNum[2*NOp];
  PosNum[0] = 0;
  int nOp = 1;
  int Comm[NOp];
  double Numb;
  char *cPart = cmd;
  for(int i=0;i<NOp;i++){
    cPart = strpbrk(cPart+1,"*/^+-");
    if(cPart == NULL) break;
    Comm[i] = (int)*cPart;
    PosOp[i] = cPart-cmd+1;
    PosNum[nOp] = PosOp[i];
    nOp++;
  }
  Sort(PosNum,nOp);
  PosNum[nOp++] = strlen(cmd);
  double Res = 0.;
  int NVar = 0;
  for(int i=0;i<nOp-1;i++){
    if( *(cmd+PosNum[i]) == 'x'){
      sscanf(cmd+PosNum[i]+1,"%d",&NVar);
      Numb = st[NVar][n];
    }
    else
      sscanf(cmd+PosNum[i],"%lf",&Numb);
    //printf("%d) %s %d %lf %lf %c\n",i,cmd+PosNum[i],PosNum[i],Numb,x,Comm[i]);
    if(i == 0){ 
      Res = Numb;
      continue;
    }
    switch(Comm[i-1]){
    case '+':
      Res += Numb;
      break;
    case '-':
      Res -= Numb;
      break;
    case '/':
      Res /= Numb;
      break;
    case '*':
      Res *= Numb;
      break;
    }
  }
  return Res;
}
