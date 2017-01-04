#include "VarSegnali.h"
Variabili::Variabili(int nSegn,int nPeriodi,double xMediaSu,double xScartoSu,double xMediaGiu,double xScartoGiu,double xPasso,double xIncr){
  NSegn = nSegn;
  Periodi = nPeriodi;
  Passo = xPasso;
  MediaSu = xMediaSu;
  ScartoSu = xScartoSu;
  MediaGiu = xMediaGiu;
  ScartoGiu = xScartoGiu;
  Incr = xIncr;
  NMass = (int)((double)Periodi/(double)Passo);
  sMedia = new double[NMass];
  st = new bool[NMass];
  mSu = new Matematica();
  mGiu = new Matematica();
}
bool Variabili::Simula(){
  double Fase=0.;
  bool Livello=0;
  double Delta=0;
  double IncrSu=0.;
  double IncrGiu=0.;
  //------------------Per le volte successive------------
  if(NMass != (int)((double)Periodi/(double)Passo)){
    NMass = (int)((double)Periodi/(double)Passo);
    delete sMedia;
    delete st;
    sMedia = new double[NMass];
    st = new bool[NMass];
  }
  printf("Inizializzo %d Segnali da %d Punti (%d Periodi con Passo %.3f)\n               con livello 1 %.2f +- %.2f con livello 0 %.2f +- %.2f\n               incrementando di %.2f\n",
	 NSegn,NMass,Periodi,Passo,MediaSu,ScartoSu,MediaGiu,ScartoGiu,Incr);
  if(!mSu->InizializzaGaussiano(ScartoSu,(int)(MediaSu/Passo)) )
    return 0;
  if(!mGiu->InizializzaGaussiano(ScartoGiu,(int)(MediaGiu/Passo)) )
    return 0;
  double Temp=0.;
  for(int i=0;i<NMass;i++){
    sMedia[i]=0.;
  }
  //---------------Inizia-------------------
  for(int i=0;i<NSegn;i++){
    int m=0;//indice di passo
    if( (i%1000)==0 )
      printf("Processo il segnale numero %d di %d\n",i,NSegn);
    if(Casuale()<=(MediaSu)/(MediaGiu)){
      Livello = 1;//Livello di partenza
      Temp = POS(mSu->Gaussiano((MediaSu+IncrSu),(ScartoSu+IncrSu)));//Fase iniziale (m=0)
      Fase=Casuale()*(MediaSu+IncrSu);
    }
    else{
      Livello=0;
      Temp = POS(mGiu->Gaussiano(MediaGiu+IncrGiu,ScartoGiu+IncrGiu));
      Fase=Casuale()*(MediaGiu+IncrGiu);
    }// Primo periodo
    for(int l=0;Fase+Passo*(double)l<Temp;l++){
      if(i == 0) st[m] = Livello;
      sMedia[m]+=(double)Livello;
      m++;
    }//inizia l'altro periodo
    for(int j=1;m<NMass;j++){
      Livello = !Livello;
      if(Livello){
	Temp= POS(mSu->Gaussiano(MediaSu+IncrSu,ScartoSu+IncrSu));
      }
      else if(!Livello){
	Temp = POS(mGiu->Gaussiano(MediaGiu+IncrGiu,ScartoGiu+IncrGiu));
      }
      for(int l=0;Passo*(double)l<Temp;l++){
	//	printf("Livello %d\n",Livello);
	if(i == 0) st[m] = Livello;
	sMedia[m]+=(double)Livello;
	m++;
	if(m>=NMass)break;
      }
    }
    IncrSu += Incr*(MediaSu)/(MediaGiu);
    IncrGiu += Incr;
  }
  char *Comando;
  char *Nome;
  char *Cartella;
  Comando = new char[50];
  Cartella = new char[50];
  Nome = new char[50];
  //   sprintf(Cartella,"%.0fk%dp%di%.0fm%.0fs",(double)(NSegn/1000),Periodi,Inter,Media,Scarto);
  //   sprintf(Comando,"mkdir %s",Cartella);
  //   system(Comando);
  sprintf(Cartella,".");
  sprintf(Nome,"%s/Segnale.dat",Cartella);
  if((SEGNALE = fopen(Nome,"w"))==0)
    printf("Non s'apre %s, Permessi?\n",Nome);
  sprintf(Nome,"%s/MediaSegnali.dat",Cartella);
  if((MEDIA = fopen(Nome,"w"))==0)
    printf("Non s'apre %s, Permessi?\n",Nome);
  for(int i=0;i<NMass;i++){
    fprintf(SEGNALE,"%d\n",st[i]);
    if(i>(int)(MediaGiu/Passo)) fprintf(MEDIA,"%g\n",sMedia[i]);
  }
  fclose(SEGNALE);
  fclose(MEDIA);
  return 1;
}
Variabili::~Variabili(){
  delete st;
  delete sMedia;
}
void Variabili::Cambia(char *cosa,int n){
  if (n > 0){
    if(!strcmp(cosa,"num"))
      NSegn = n;
    else if(!strcmp(cosa,"per")){
      Periodi = n;
    }
  }
}
void Variabili::Cambia(char *cosa,double x){
  if(x>0){
    if(!strcmp(cosa,"mSu") )
      MediaSu = x;
    else if(!strcmp(cosa,"sSu") )
      ScartoSu = x;
    else if(!strcmp(cosa,"mGiu") )
      MediaGiu = x;
    else if(!strcmp(cosa,"sGiu") )
      ScartoGiu = x;
    else if(!strcmp(cosa,"pas") )
      Passo = x;
  }
  if(!strcmp(cosa,"incr") )
    Incr = x;
}
