#include "../include/Matematica.h"

double Matematica::ContactAngle(double x){
  double Ris;
  //Ris = PreFact*pow(3./(2.+3.*sin(x)-CUB(sin(x)) ), 4./3. )*QUAD(QUAD(cos(x))) - Ypsilon;
  Ris = PreFact*pow(1-cos(x), -4./3. )*pow(sin(x),2.) - Ypsilon;
  return Ris;
}
double Matematica::fProva(double x){
  double Ris;
  //Ris = pow(x,5.)/((exp(x)-1)*(1-exp(-x)));
  // Ris = pow(x,5.)*exp(x)/QUAD(exp(x)-1);
  // Ris = 1.*pow(x-150,4.)-1.*pow(x-150,5.)+x;
  // Ris = pow(x-150,2.)+x-150;
  Ris = 1./(13.9332529115775*x*x - 4.*13.9332529115775*0.13710248135043*x + 7.86306085419497) + .5/(-0.0805855688278211*x + 10.6925178203198);
  //Ris = 1./(12.5342*x*x) + .5/(-0.207807999999993*x + 17.4283);
  return Ris;
}
double Matematica::F(double T,double TD){
  return pow((T/TD),5.)*Integrazione(.00001,TD/T);
}
double Matematica::Integrazione(double a,double b){
  double Delta = (b-a)/(double)NPassi;
  double Risp = 0.;
  for(int i=0;a+((double)i)*Delta<b;i+=3){//+3 se si sovrappongono
    Risp += 3.*3.*Delta*( Evalx(a+Delta*i)+3.*Evalx(a+Delta*(i+1))+3.*Evalx(a+Delta*(i+2))+Evalx(a+Delta*(i+3)) )/8.;
  }
  return Risp;
}
double Matematica::Integrazione(double *st,double *sw,int NMass){
  double NMassInv = 1/((double)NMass);
  sw[0] = st[0];
  for(int i=1;i<NMass;i++){
    sw[i] = st[i] + sw[i-1];
  }
  return sw[NMass-1];
}
double Matematica::Df(double x,double Delta){
  return ( Evalx(x+Delta)-Evalx(x) )/Delta;
}
void Matematica::Derivata(double *st,double *sw,int NMass){
  for(int i=0;i<NMass-1;i++){
    sw[i] = (st[i+1]-st[i])/2.;
  }
}
void Matematica::DerO4(double *st,double *sw,int NMass){
  sw[0] = 0.;
  sw[NMass-1] = 0.;
  for(int i=1;i<NMass-1;i++){
    if(i<2)
      sw[i] = (st[i+1] - st[i-1])/2.;
    else if (i < NMass - 2)
      sw[i] = st[i-2] - 8.*st[i-1] + 8.*st[i+1] - st[i+2];
    else if (i < NMass - 1)
      sw[i] = (st[i+1] - st[i-1])/2.;
  }
}
int Matematica::Zeri(double a,double b,double *Radici,int NRadici){
  double Uno;double Due;double Delta;
  int rad=0;
  RADICE Rad;
  if( a > b ){
    Uno = b ; Due = a;
  }
  else{
    Uno = a ; Due = b;
  }
  Delta = (Due - Uno)/2.;
  //  for(int p=0;p<4;p++){
  for(int i=0;i<NRadici*2;i++){
      Delta = (Due - Uno)/(2*NRadici);
      Rad = RegulaFalsi(Uno+i*Delta,Due-(2*NRadici-i-1)*Delta);
      //printf("%d %d) %g|%g %g|%g %g|%g\n",NRadici,i,Uno+i*Delta,Rad.iLim,Due-(2*NRadici-i-1)*Delta,Rad.sLim,Delta,Rad.sLim-Rad.iLim);
      if(Rad.IfRis == 1){
	printf("Found a root in %lf\n",Rad.Zero);
	Radici[rad] = Rad.Zero;
	rad++;
      }
    }
    //  }
  return rad;
}
double Matematica::Estremo(double a,double b){
  double Uno;double Due;double Tre;
  //  a<b : Delta=(b-a)/NPassi ? Delta=(a-b)/NPassi;
  double Delta=0.;
  int NLim = 1000;
  Delta=(b-a)/(double)NLim;
  Uno=a;Due=b;Tre=0.;
  for(int i=0;i<NLim;i++){
    if( ASS((Evalx(Tre)-0.)) < PrecMinimo){
      break;
    }
    Tre = Due - (Due - Uno)/(Df(Due,Delta)-Df(Uno,Delta))*Df(Due,Delta);
    Due = Tre;
    Uno = Due;
  }
  return Tre;
}
RADICE Matematica::RegulaFalsi(double a,double b){
  RADICE Rad; 
  double Uno=a;double Due=b;double Tre=0.;
  double dIncr = 100.;
  double Delta = (b-a)/dIncr;
  Uno=a;Due=b;Tre=0.;
  FILE *CONTROLLA;
  CONTROLLA = fopen("RegulaFalsi.dat","w");
  for(int i=0;i<NPassi;i++){
    if( ASS(Evalx(Tre)) < PrecMinimo){
      break;
    }
    else if( ASS(Evalx(Due)) < PrecMinimo){
      Tre = Due;
      break;
    }
    else if( ASS(Evalx(Uno)) < PrecMinimo){
      Tre = Uno;
      break;
    }
    if( Evalx(Due) < 0. && Evalx(Uno) > 0.){
      if( Evalx(Due) < 0. &&  Evalx(Uno + Delta) > 0.){
	Uno += Delta;
      }
      if(Evalx(Due - Delta) < 0. && Evalx(Uno) > 0.){
	Due -= Delta;
      }
      else {
	dIncr *= 10.;
      }
    }
    else if( Evalx(Due) > 0. && Evalx(Uno) < 0.){
      if( Evalx(Due) > 0. && Evalx(Uno+Delta) < 0.){
	Uno += Delta;
      }
      if(Evalx(Due - Delta) > 0.  && Evalx(Uno) < 0.){
	Due -= Delta;
      }
      else {
	dIncr *= 10.;
      }
    }
    else{ 
      Uno += (b-a) / dIncr;
    }
    if(Uno > Due){
      Uno = a;
      dIncr *= 10.;
    }
    Delta = (b-a)/dIncr;
    Tre = (Due + Uno)/2.;
    //printf("%d) Uno %g|%g  Due %g|%g  Evalx(Tre) %g|%g Delta %g Incr %g\n",i,Uno,Evalx(Uno),Due,Evalx(Due),Tre,Evalx(Tre),Delta,dIncr);
    //fprintf(CONTROLLA,"%lf %lf\n",Tre,Evalx(Tre));
  }
  if( !(ASS((Evalx(Tre)-0.)) < PrecMinimo)){
    Rad.IfRis = 0;
  }
  else 
    Rad.IfRis = 1;
  Rad.Zero = Tre;
  Rad.iLim = Uno;
  Rad.sLim = Due;
  fclose(CONTROLLA);
  return Rad;
}
RADICE Matematica::Newton(double a){
  RADICE Rad;
  double Uno=a,Due =0.,Tre=0.;
  double m=0.,q=0.;
  FILE *CONTROLLA;
  CONTROLLA = fopen("Newton.dat","w");
  for(int i=0;i<NPassi;i++){
    if( ASS((Evalx(Tre))) < PrecMinimo){
      break;
    }
    Due = Uno + 1e-5;
    m = (Evalx(Due) - Evalx(Uno))/(Due - Uno);
    q = Evalx(Uno) - m*Uno;
    Tre = -q/m;
    printf("%d) %lf  %lf  %lf\n",i,Uno,Due,Evalx(Tre));
    fprintf(CONTROLLA,"%lf %lf\n",Tre,Evalx(Tre));
    Uno = Tre;
  }
  if( !(ASS((Evalx(Tre)-0.)) < PrecMinimo)){
    printf("Calculation failed\n");
    Rad.IfRis = 0;
  }
  else 
    Rad.IfRis = 1;
  Rad.Zero = Tre;
  Rad.iLim = Uno;
  Rad.sLim = Due;
  fclose(CONTROLLA);
  return Rad;
}
double Matematica::Gauss(double Media,double Scarto,double x){
  return 1./(Scarto*sqrt(DUE_PI))*exp(- .5*SQR((x-Media)/Scarto) );
}
double Matematica::IntegrazioneGauss(double a,double b,double Scarto){
  double Delta = (b-a)/(double)NPassi;
  double Risp = 0.;
  for(int i=0;a+((double)i)*Delta<b;i+=3){//+3 se si sovrappongono
    Risp += 3*Delta*( Gauss(0.,Scarto,a+Delta*i)+3*Gauss(0.,Scarto,a+Delta*(i+1))+3*Gauss(0.,Scarto,a+Delta*(i+2))+Gauss(0.,Scarto,a+Delta*(i+3)) )/8;
  }
  return Risp;
}
void Matematica::SquareGradient(double *st,double *sw,int NMass){
  DerO4(st,sw,NMass);
}
void Matematica::NormalizeVect(double *st,int NMass){
  double Norm = 0.;
  for(int n=0;n<NMass;n++){
    Norm += SQR(st[n]);
  }
  Norm = 1./sqrt(Norm);
  for(int n=0;n<NMass;n++){
    st[n] = st[n]*Norm;
  }  
}
double Matematica::Norm(double *st,int NMass){
  double Norm = 0.;
  for(int n=0;n<NMass;n++){
    Norm += SQR(st[n]);
  }
  return sqrt(Norm);
}
void Matematica::Modulo(double *st,double *sw,int NMass){
  for(int n=0;n<NMass;n++){
    sw[n] = POS(st[n]);
  }
}
double Matematica::Fattoriale(int n){
  if(n < 0) { printf("Il fattoriale di numeri negativi non ha senso\n"); return 0;}
  if( n == 0 ) return 1;
  int Ris=1;
  for(int i=n;i>0;i--){
    //    printf("Fatt %d\n",i);
    Ris *= i;
  }
  return (double) Ris;
}
double Matematica::Gamma(int n){
  return Fattoriale(n-1);
}
double Matematica::Elevato(double x,int Volte){
  double Risp=1.;  
  double Moltx = Volte >= 0 ? x : 1./x;
  for(int v=0;v<Volte;v++)
    Risp *= Moltx;
  return Risp;
}
double Matematica::Bessel(double Val,int Ord){
  int NOrd = POS(Ord);
  int NMax = 10;
  double Risp=0.;
  for(int n=0;n<NMax;n++){
    //    printf("  %d\n",n);
    Risp += Elevato(-1.,n)/(Fattoriale(n)*Gamma(NOrd+n+1))*Elevato(.5*Val,2*n+NOrd);
  }
  return Ord>=0 ? Risp : Elevato(-1.,Ord)*Risp;
}
double Matematica::Neumann(double Val,int Ord){
  double Angolo = Ord*.5*DUE_PI+.0001;
  return Bessel(Val,Ord)*(cos(Angolo) - Elevato(-1.,Ord))/sin(Angolo);
}
double Matematica::QuasiBessel(double Val,int Ord){
  if( Val > (double)Ord +1.)
    return sqrt((2./(PI*Val))) *cos( Val- (2.*Ord+1.)*PI*.25);
  else
    return Elevato(.5*Val,Ord)/(Fattoriale(Ord));
}
double Matematica::QuasiNeumann(double Val,int Ord){
  if(Ord == 0)
    return 2./PI*log(Val);
  if(Val > (double) Ord+1.)
    return sqrt(2./(PI*Val))*sin(Val-(2.*Ord+1.)*PI*.25);
  else
    return Elevato(2,Ord)*Fattoriale(Ord-1)/PI*Elevato(Val,-Ord);
}
double Matematica::Segno(int n){
  return (n%2)==1 ? -1. : 1.;
}
double Matematica::WeightFunction(double x,double a){
  double Risp=2.*x*x*x - 3.*(a+1.)*x*x-3.*a*a+1.;
  Risp /= CUB((1.-a));
  return Risp;
}
double Matematica::WeightFunction2(double x,double a){
  double Risp= -a*x*x + 1. + 4./3.*PI*a;
  Risp /= CUB((1.-a));
  return Risp;
}
double Matematica::LJHamaker(double Rad,double RadNp,double Theta){
  //return 1./CUB(Rad-RadNp);
  double Num = 2.*DUE_PI*QUAD(RadNp)*sin(Theta);
  double Den = QUAD(RadNp) + QUAD(Rad+RadNp) - 2.*(Rad+RadNp)*RadNp*cos(Theta);
  return Num/CUB(Den);// - Num/CUB(CUB(Den));
}
double Matematica::LJ39(double r,double r_np){
  return pow(1./(r-r_np),9.) - pow(1./(r-r_np),3.);
}
double Matematica::LJHamaker(double Rad,double RadNp){
  double ThetaMax = PI;
  double ThetaMin = 0.;
  double ThetaDelta = (ThetaMax-ThetaMin)/100.;
  double Risp = 0.;
  for(double Theta = ThetaMin;Theta<ThetaMax;Theta+=ThetaDelta){
    Risp += ThetaDelta*.5*(LJHamaker(Rad,RadNp,Theta) + LJHamaker(Rad,RadNp,Theta+ThetaDelta) );
  }
  return Risp;
}
double Matematica::LJHamakerCum(double Rad,double RadNpMin,double RadNpMax){
  double RadNpDelta = (RadNpMax-RadNpMin)/100.;
  double Risp = 0.;
  for(double RadNp = RadNpMin;RadNp<RadNpMax;RadNp+=RadNpDelta){
    //Risp += RadNpDelta*LJHamaker(Rad,RadNp+RadNpDelta*0.5);
    Risp += RadNpDelta*.5*(LJHamaker(Rad,RadNp) + LJHamaker(Rad,RadNp+RadNpDelta) );
    //Risp += 3.*RadNpDelta*(LJHamaker(Rad,RadNp) + 3.*LJHamaker(Rad,RadNp+RadNpDelta) + 3.*LJHamaker(Rad,RadNp+2.*RadNpDelta) + 3.*LJHamaker(Rad,RadNp+3.*RadNpDelta) )/8.;
  }
  return Risp;
}
double Potenziale(double Dist, double RadNp){
  double Rad  = Dist + RadNp;
  double Pre1 = 2./(12.*Rad);
  double Post1 = RadNp / CUB(Rad+RadNp) + RadNp/CUB(Rad-RadNp);
  double Pre2 = 1./(12.*Rad);
  double Post2 = 1./QUAD(Rad+RadNp) - 1./QUAD(Rad+RadNp);
  return PI*(Pre1*Post1 + Pre2*Post2);
}
double Potenziale2(double Dist, double RadNp){
  double Rad = Dist + RadNp;
  double Pre1 = -( QUAD(RadNp) - QUAD(Rad) )/(4.*Rad);
  double Post1 = 1./QUAD(QUAD(Rad+RadNp)) - 1./QUAD(QUAD(Rad-RadNp));
  double Pre2  = -2./(3.);
  double Post2 = 1./CUB(Rad+RadNp) - 1./CUB(Rad+RadNp);
  double Pre3  = 1./(2.*Rad);
  double Post3 = 1./QUAD(Rad+RadNp) - 1./QUAD(Rad+RadNp);
  return PI*(Pre1*Post1 + Pre2*Post2 + Pre3*Post3);
}
void Matematica::IntegraA3(){
  double RadNpMin = 0.;
  double RadNpMax = .001;
  double RadNpDelta = (RadNpMax - RadNpMin)/20.;
  double RadMin = 0.1;
  double RadMax = 3.;
  double RadDelta = (RadMax-RadMin)/100.;
  char *FileName = (char *)calloc(60,sizeof(char));
  double RadNp = RadNpMax;
  //for(double RadNp = RadNpMin;RadNp <= RadNpMax;RadNp += RadNpDelta)
    {
    sprintf(FileName,"Potential%0.2f.dat",RadNp);
    FILE *POT = fopen(FileName,"w");
    for(double Rad = RadMin;Rad < RadMax;Rad += RadDelta){
      fprintf(POT,"%lf %g %g %g \n",Rad,LJHamaker(Rad,RadNp),LJHamakerCum(Rad,RadNpMin,RadNp),Potenziale(Rad,RadNp) );
    }
    fclose(POT);
  }
}
