#include "../../include/Matematica.h"
#include "../../include/VarDatFile.h"

#include <iostream>
using namespace std;

int main(int argc,char **argv){
  FILE *SEGNALE;
  char Uscita[256];
  char Comando[50];
  char Executex[60];
  char Executey[60];
  double *Dis = new double;
  int CoordX = 0;
  int NBin = 200;
  int NMax=0;
  int NVar =0;
  int CoordY = 1;
  int CoordY1 = 1;
  int CoordX1 = 1;
  int NMedia=10;
  int LogLog = 0;
  int IfNorm = 0;
  int IfBorder = 0;
  int NSample = 10;
  double *Intervalli;
  double *dInter;
  //vector <double> st;
  double *sp;
  double *sw;
  double Border[4] = {0.,0.,0.,0.};
  double Param = 1.;
  VarDatFile *v1;
  Matematica *m1 = new Matematica();
  MOMENTI Mom;
  sprintf(Uscita,"Risposta.dat");
  //  for(int i=0;i<argc;i++){
  int NFile = 0;
  int *ArgFile = (int *)calloc(argc,sizeof(int));
  int IfUser = 1;
  int NElMin = 0;
  int NElMax = 0;
  for(int i=1;i<argc;i++){
    if(!strcmp(*(argv+i),"--col")){
      if(i<argc-1){
	i+=1;
      }	
    }
    else if(!strcmp(*(argv+i),"--exec")){
      if(i<argc-1){
	sprintf(Executex,"%s",*(argv+i+1));
	sprintf(Executey,"%s",*(argv+i+2));
	i+=2;
	sprintf(Comando,"exec");
	IfUser = 0;
      }	
    }
    else if(!strncmp(*(argv+i),"--Border",8)){
      if(argc < i){printf("How many values?\n");return 0;}
      sscanf(argv[i+1],"%lf",Border);
      sscanf(argv[i+2],"%lf",Border+1);
      IfBorder=1;
      i+=2;
    }
    else if(!strncmp(*(argv+i),"-l",2)){
      if(argc < i){printf("Whic log scale?\n");return 0;}
      char LogScale[12];
      sscanf(argv[i+1],"%s",LogScale);
      LogLog = 0;
      if(!strcmp(LogScale,"logx"))
	DIS_ADD_TYPE(LogLog,DIS_LOGX);
      else if(!strcmp(LogScale,"logy"))
	DIS_ADD_TYPE(LogLog,DIS_LOGY);
      else if(!strcmp(LogScale,"loglog"))
	DIS_ADD_TYPE(LogLog,DIS_LOGLOG);
      else{
	printf("Log scale not recognized\n");
      }
      i+=1;
    }
    else if(!strncmp(*(argv+i),"--NElMin",8)){
      if(argc < i){printf("How many values?\n");return 0;}
      sscanf(argv[i+1],"%d",&NElMin);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"--NElMax",8)){
      if(argc < i){printf("How many values?\n");return 0;}
      sscanf(argv[i+1],"%d",&NElMax);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-m",2)){
      if(argc < i){printf("which interval for the running average?\n");return 0;}
      sscanf(argv[i+1],"%d",&NMedia);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"--norm",6)){
      IfNorm = 1;
    }
    else if(!strncmp(*(argv+i),"--NSample",8)){
      if(argc < i){printf("How many samples?\n");return 0;}
      sscanf(argv[i+1],"%d",&NSample);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-o",2)){
      if(argc < i){printf("which output file?\n");return 0;}
      sscanf(argv[i+1],"%s",Uscita);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-p",2)){
      if(argc < i){printf("what parameter?\n");return 0;}
      sscanf(argv[i+1],"%lf",&Param);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-v",2)){
      if(argc < i){printf("How many values?\n");return 0;}
      sscanf(argv[i+1],"%d",&NBin);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-x",2)){
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-y",2)){
      i+=1;
    }
    else if(!strncmp(*(argv+i),"--",2)){
      sprintf(Comando,"%s",argv[i]+2);
      IfUser = 0;
    }
    else if(strncmp(*(argv+i),"-",1)){
      ArgFile[NFile] = i;
      NFile++;
    }
  }
  if(NFile == 0){
    printf("File arguments are missing\n");
    //return 0;
  }
  v1 = new VarDatFile(argv,ArgFile,NFile,NBin);
  for(int i=1;i<argc;i++){
    if(!strncmp(*(argv+i),"-a",2)){
      v1->ImpSequence();
      }
    if(!strncmp(*(argv+i),"--col",5)){
      if(argc < i){printf("Which column?\n");return 0;}
      sscanf(argv[i+1],"%d",&CoordY);
      v1->Punta(CoordY);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-l",2)){
      if(argc < i){printf("Which coord for the  log scale?\n");return 0;}
      char LogScale[12];
      sscanf(argv[i+1],"%s",LogScale);
      LogLog = 0;
      if(!strcmp(LogScale,"logx"))
  	DIS_ADD_TYPE(LogLog,DIS_LOGX);
      else if(!strcmp(LogScale,"logy"))
  	DIS_ADD_TYPE(LogLog,DIS_LOGY);
      else if(!strcmp(LogScale,"loglog"))
  	DIS_ADD_TYPE(LogLog,DIS_LOGLOG);
      else{
  	printf("Log scale not recognized\n");
      }
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-y1",3)){
      if(argc < i){printf("Which column?\n");return 0;}
      sscanf(argv[i+1],"%d",&CoordY1);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-x1",3)){
      if(argc < i){printf("Which column?\n");return 0;}
      sscanf(argv[i+1],"%d",&CoordX1);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-x",2)){
      if(argc < i){printf("which x coord?\n");return 0;}
      sscanf(argv[i+1],"%d",&CoordX);
      v1->ImpCoordX(CoordX);
      i+=1;
    }
    else if(!strncmp(*(argv+i),"-y",2)){
      if(argc < i){printf("Which column?\n");return 0;}
      sscanf(argv[i+1],"%d",&CoordY);
      v1->ImpCoordY(CoordY);
      i+=1;
    }
  }
  NMax = v1->pNMax();
  if(NElMax == 0 || NElMax > NMax) NElMax = NMax;
  if(NElMin > NMax || NElMin < 0) NElMin = 0;
  NVar = v1->pNVar();
  sp = new double;
  sw = new double[NMax];
  Intervalli = new double[NBin];
  dInter = new double[NBin];
  printf("Loaded a file of %d points %d cols from %d files: %s ... to write on %s\n",NMax,NVar,NFile,argv[ArgFile[0]],Uscita);
  //goto vet;
  v1->pMinMaxGlob(NElMin,NElMax);
  while(strcmp(Comando,"q") ){
    if(IfUser){
      printf("Matematica> ");
      scanf("%s",Comando);
    }
    else if(!strcmp(Comando,"agg")){
      char NomeFile[120];
      printf("Nome file :");
      scanf("%s",NomeFile);
      v1->Aggiungi(NomeFile);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"create")){
      FILE *Ciccia = fopen(Uscita,"w");
      if(Ciccia == NULL){printf("%s could not be created\n",Uscita);exit(0);}
      int NPassi = 100;
      double InvNPassi = 1./(double)NPassi;
      double Freq1 = 32.;
      double Freq2 = 4.;
      for(int px=0;px<NPassi;px++){
	double x = px*InvNPassi*DUE_PI;
	fprintf(Ciccia,"%lf %lf\n",x,sin(x*Freq1));
      }
    }
    else if(!strcmp(Comando,"create2d")){
      FILE *Ciccia = fopen(Uscita,"w");
      int NPassi = 100;
      double InvNPassi = 1./(double)NPassi;
      double Freq1 = 32.;
      double Freq2 = 4.;
      int IfOctave = 1;
      if(!IfOctave){
	for(int px=0;px<NPassi;px++){
	  for(int py=0;py<NPassi;py++){
	    double x = px*InvNPassi*DUE_PI;
	    double y = py*InvNPassi*DUE_PI;
	    fprintf(Ciccia,"%lf %lf %lf\n",x,y,sin(x*y*Freq1));
	  }
	}
      }
      else{
	fprintf(Ciccia,"# Created by ElPoly\n");
	fprintf(Ciccia,"# name: tz\n");
	fprintf(Ciccia,"# type: matrix\n");
	fprintf(Ciccia,"# rows: %d\n",NPassi);	
	fprintf(Ciccia,"# columns: %d\n",NPassi);	
	for(int px=0;px<NPassi;px++){
	  for(int py=0;py<NPassi;py++){
	    double x = px*InvNPassi*DUE_PI;
	    double y = py*InvNPassi*DUE_PI;
	    fprintf(Ciccia,"%lf ",sin(x*y*Freq2));
	  }
	  fprintf(Ciccia,"\n");
	}
      }
      fclose(Ciccia);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"dis")){
      Mom = v1->DistrGaussSegnale(NElMin,NElMax,IfNorm);
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      }
      fprintf(USCITA,"# %lf %lf %lf\n",Mom.Uno,Mom.Due,Mom.Tre);
      for(int i=0;i<NBin;i++){
	double x = Mom.Delta*i+Mom.Min;
	fprintf(USCITA,"%lf %g\n",x,v1->pInter(i));
      }
      fclose(USCITA);
      printf("Media %lf Scarto %lf N %d\n",Mom.Uno,Mom.Due,Mom.Num);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"disErr")){
      Mom = v1->DistrSignErr(NElMin,NElMax,Border,IfNorm);
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      }
      fprintf(USCITA,"# %lf %lf %lf\n",Mom.Uno,Mom.Due,Mom.Tre);
      for(int i=0;i<NBin;i++){
	double x = Mom.Delta*i+Mom.Min;
	fprintf(USCITA,"%lf %g %g\n",x,v1->pInter(i),v1->pError(i));
      }
      fclose(USCITA);
      printf("Media %lf Scarto %lf N %d\n",Mom.Uno,Mom.Due,Mom.Num);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"disB")){
      Mom = v1->DistrSegnale(NElMin,NElMax,Border,IfNorm);
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n","Distribuzione.dat");
      }
      fprintf(USCITA,"# %lf %lf %lf\n",Mom.Uno,Mom.Due,Mom.Tre);
      for(int i=0;i<NBin;i++){
	double x = Mom.Delta*i+Mom.Min;
	fprintf(USCITA,"%lf %g\n",x,v1->pInter(i));
      }
      fclose(USCITA);
      printf("Media %lf Scarto %lf N %d\n",Mom.Uno,Mom.Due,Mom.Num);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"disLog")){
      Mom = v1->DistrLogSegnale(NElMin,NElMax,IfNorm);
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n","Distribuzione.dat");
      }
      fprintf(USCITA,"# %lf %lf %lf\n",Mom.Uno,Mom.Due,Mom.Tre);
      for(int i=0;i<NBin;i++){
	double x = Mom.Delta*i+Mom.Min;
	fprintf(USCITA,"%lf %g\n",x,v1->pInter(i));
      }
      fclose(USCITA);
      printf("Media %lf Scarto %lf N %d\n",Mom.Uno,Mom.Due,Mom.Num);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"disExp")){
      Mom = v1->DistrExpSegnale(NElMin,NElMax,IfNorm);
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n","Distribuzione.dat");
      }
      fprintf(USCITA,"# %lf %lf %lf\n",Mom.Uno,Mom.Due,Mom.Tre);
      for(int i=0;i<NBin;i++){
	double x = Mom.Delta*i+Mom.Min;
	fprintf(USCITA,"%lf %g\n",x,v1->pInter(i));
      }
      fclose(USCITA);
      printf("Media %lf Scarto %lf N %d\n",Mom.Uno,Mom.Due,Mom.Num);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"disSample")){
      double xBound[2] = {Border[0],Border[1]};
      double **Distr;
      Distr = new double*[NSample];
      for(int s=0;s<NSample;s++){
	Distr[s] = new double[NBin];
      }
      v1->DistrSignSample(NElMin,NElMax,Distr,NSample,IfNorm,xBound);
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n","Distribuzione.dat");
      }
      fprintf(USCITA,"# distribution per sample\n");
      for(int b=0;b<NBin;b++){
	double x = b/(double)NBin*(xBound[1]-xBound[0]) + xBound[0];
	fprintf(USCITA,"%lf ",x);
	for(int s=0;s<NSample;s++){
	  fprintf(USCITA,"%lf ",Distr[s][b]);
	}
	fprintf(USCITA,"\n");
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"DrawEq")){
      FILE *FWrite = fopen(Uscita,"w");
      if(!IfBorder){
	Border[0] = v1->pxMinGlob(NElMin,NElMax);
	Border[1] = v1->pxMaxGlob(NElMin,NElMax);
      }
      double Dx = (Border[1]-Border[0])/(double)NBin;
      for(int i=0;i<NBin;i++){
	double x = i*Dx+Border[0];
	fprintf(FWrite,"%lf %lf\n",x,v1->fProva(x));
      }
      fclose(FWrite);
      if(!IfUser) return 0;
    }    
    else if(!strcmp(Comando,"coord")){
      printf("Number of coordinate: ");
      int NCoord = 0;
      scanf("%d",&NCoord);
      if(NCoord >= 0 && NCoord<NVar)
	v1->Punta(NCoord);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"MatchCoord")){
      double *st = (double *)calloc(NMax,sizeof(double));
      double tInit = v1->Val(CoordX1,0);
      int mOld = 0;
      for(int n=NElMin,m=0;n<NElMax;n++){
	mOld = m-1;
	do {
	  double t = tInit + v1->Val(CoordX,m);
	  //printf("%d %d %lf <= %lf < %lf -> %lf\n",n,m,v1->Val(CoordX1,n),t,v1->Val(CoordX1,n+1),st[n]);
	  if(v1->Val(CoordX1,n) <= t && t < v1->Val(CoordX1,n+1))
	    break;
	  m++;
	  if(m >= NMax){
	    m = mOld;
	    break;
	  }
	} while(1==1);
	st[n] = v1->Val(CoordY,m);
	//printf("%d %d %lf <= %lf < %lf -> %lf\n",n,m,v1->Val(CoordX1,n),v1->Val(CoordX,m)+tInit,v1->Val(CoordX1,n+1),st[n]);
	m++;
      }
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0)
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      for(int i=NElMin;i<NElMax;i++){
	fprintf(USCITA,"%lf %.5g %.5g\n",st[i],v1->Val(CoordY1,i),v1->Val(CoordY1+1,i));
      }
      fclose(USCITA);
      free(st);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"norm")){
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0)
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      v1->NormalizzaSegnale(NElMin,NElMax);
      for(int i=NElMin;i<NElMax;i++){
	fprintf(USCITA,"%lf %.5g\n",v1->Abscissa(CoordY,i),v1->Val(CoordY,i));
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"sin")){
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      }
      for(int i=0;i<100;i++){
	double x = i*DUE_PI*.01;
	fprintf(USCITA,"%g\n",sin(x));
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"sum")){
      double Sum = 0.;
      for(int i=NElMin;i<NElMax;i++){
	Sum += v1->Val(CoordY,i);
      }
      printf("Sum: %lf\n",Sum);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"SpD")){
      printf("Calcolo lo spettro \n");
      m1->Spettro(sp,sw,NMax);
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n",Uscita);
	return 0;
      }
      for(int i=0;i<NMax/2;i++){
	fprintf(USCITA,"%g\n",sw[i]);
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"sort")){
      v1->Punta(0);
      v1->Sort();
      v1->ScriviTutto(Uscita,0,NElMin,NElMax);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"Cas")){
      printf("Proviamo il generatore di numeri casuali \n");
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      }
      for(int i=0;i<10000;i++){
	fprintf(USCITA,"%lf\n",m1->Casuale());
      }
      fclose(USCITA);
       if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"Func")){
      FILE *USCITA;
      printf("Proviamo una funzione \n");
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      }
      for(double x=0.,Dx=.001,Max=5.;x<Max;x+=Dx){
	printf("%lf\n",x/Max);
	fprintf(USCITA,"%lf %g %g %g\n",x,m1->QuasiBessel(x,0),m1->QuasiBessel(x,1),m1->Bessel(x,0));
	//fprintf(USCITA,"%lf %g %g %g\n",x,m1->Neumann(x,0),m1->Neumann(x,1),m1->Neumann(x,-1));
	x += Dx;
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"SpeLine")){
      int NBin = NElMax-NElMin;
      double InvNVar = .5/(double)(NVar);
      double InvNBin = 1./(double)(NMax);
      double *Spe = (double *)calloc(NBin,sizeof(double));
      v1->SpeLine(NElMin,NElMax,NBin,Spe);
      FILE *FWrite = fopen("StalkLineSpectrum.dat","w");
      fprintf(FWrite,"# q=2\\pi n/N <h(q)^2>\n");
      char FName[120];
      for(int n=0;n<NBin;n++){
	fprintf(FWrite,"%lf %g \n",DUE_PI*InvNBin*n,Spe[n]*InvNVar*InvNBin);
      }
      fclose(FWrite);
      free(Spe);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"spe") ){
      FILE *USCITA;
      m1->Spettro(v1->sp,sw,NMax);
      if((USCITA = fopen(Uscita,"w"))==0){
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      }
      for(int i=0;i<NMax;i++){
	fprintf(USCITA,"%g\n",sw[i]);
      }
      fclose(USCITA);
      break;
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"RadTens") ){
      if(NFile != 4){
	printf("I need the three components of the radial pressure\n");
	if(!IfUser) return 0;	
      }
      FILE *RadTens = fopen(Uscita,"w");
      double Vol1 = 0.;
      double Vol2 = 0.;
      for(int n=0;n<NMax-1;n++){
	Vol1 = DUE_PI*SQR(v1->Abscissa(CoordY,n+1));
	double Prr = sqrt( SQR(v1->Val(1,n)) + SQR(v1->Val(6,n)));
	double Pzz = v1->Val(11,n);
	double PTrace = v1->Val(1,n)+v1->Val(5,n)+v1->Val(11,n);
	double Ptt = .5*(PTrace - Prr);
	double SurfTens = (Prr - .5*(Ptt+Pzz))/(Vol1-Vol2);
	fprintf(RadTens,"%lf %lf\n",v1->Abscissa(CoordY,n),SurfTens);
	Vol2 = Vol1;
      }
      fclose(RadTens);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"perm") ){
      int NMax = 16;
      PERMUTE *Perm = (PERMUTE *)calloc(NMax,sizeof(PERMUTE));
      m1->PermuteRandomAll(Perm,NMax);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"solve") ){
      Matrice *A = new Matrice(4,4);
//       Matrice Pasquale(4);
//       Matrice Romualdo(4);
//       Matrice Adalgiso(4);
//       Matrice Goffredo(4);
//       Pasquale.RandomFill(1.);
//       Romualdo.RandomFill(1.);
//       //Pasquale = Romualdo;
//       Adalgiso.Mult(Pasquale,Romualdo);
//       Adalgiso.Print();
//       exit(0);
      A->Set(0,0,1.);A->Set(0,1,6.);A->Set(0,2,5.);A->Set(0,3,8.);
      A->Set(1,0,5.);A->Set(1,1,2.);A->Set(1,2,5.);A->Set(1,3,5.);
      A->Set(2,0,6.);A->Set(2,1,3.);A->Set(2,2,7.);A->Set(2,3,1.);
      A->Set(3,0,8.);A->Set(3,1,7.);A->Set(3,2,2.);A->Set(3,3,6.);

//       A->Set(0,0,1.);A->Set(0,1,0.);A->Set(0,2,0.);A->Set(0,3,1.);
//       A->Set(1,0,0.);A->Set(1,1,1.);A->Set(1,2,0.);A->Set(1,3,0.);
//       A->Set(2,0,0.);A->Set(2,1,0.);A->Set(2,2,1.);A->Set(2,3,0.);
//       A->Set(3,0,2.);A->Set(3,1,0.);A->Set(3,2,0.);A->Set(3,3,1.);
      //      A->RandomFill(20.);
      A->Print();
      double *Known = (double *)calloc(4,sizeof(double));
      Known[0] = 1.; Known[1] = 4.; Known[2] = 3.; Known[3] = 2.;
      double *UnKnown = (double *)calloc(4,sizeof(double));
      A->Solve(Known,UnKnown);
      for(int i=0;i<4;i++){
	printf("%lf %lf\n",Known[i],UnKnown[i]);
      }
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"zeri") ){
      double a=0.,b=0.;
      int rad;
//       printf("Minimo valore: ");
//       scanf("%lf",&a);
//       printf("Massimo valore: ");
//       scanf("%lf",&b);
//%% 2.5 
//%% 2.0
//%% 1.5 23.376860  3817.700000
//%% 1.0 21.643605  3817.700000
      int NRadici = 4;
      double *Radici;
      m1->Ypsilon = 23.37;
      m1->PreFact = 3817.7;;
      m1->Func = &Matematica::ContactAngle;
      Radici = (double *)malloc(NRadici*sizeof(double));
      rad = m1->Zeri(0.,DUE_PI/2.,Radici,NRadici);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"fal") ){
      double a=0.,b=0.;
      RADICE Zero;
//       printf("Minimo valore: ");
//       scanf("%lf",&a);
//       printf("Massimo valore: ");
//       scanf("%lf",&b);
      Zero = m1->RegulaFalsi(.1,3.);
      printf("Lo Zero e in x %lf\n",Zero.Zero);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"new") ){
      double a=0.;
      RADICE Zero;
      printf("Punto iniziale: ");
      scanf("%lf",&a);
      Zero = m1->Newton(a);
      printf("Lo Zero e in x %lf\n",Zero.Zero);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"fun") ){
      FILE *USCITA;
      USCITA = fopen("Funzione.dat","w");
      double A = 0.;
      double B = DUE_PI;
      int NPunti = 10000;
      double Delta = (B-A)/(double)NPunti; 
      m1->Ypsilon = 3.814;
      m1->PreFact = 0.169370;
      m1->Func = &Matematica::ContactAngle;
      for(int i=0;i<NPunti;i++){
	fprintf(USCITA,"%g  %g\n",Delta*(double)i,m1->Evalx(Delta*(double)i));
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"FitLin") ){
      FILE *USCITA = fopen(Uscita,"w");
      RETTA r1 = v1->InterRettSegnale(CoordY,NElMin,NElMax,LogLog);
      printf("y = (%lf +- %lf)x + (%lf +- %lf)\n",r1.m,r1.ErrM,r1.q,r1.ErrQ);
      fprintf(USCITA,"#y = (%lf +- %lf)x + (%lf +- %lf)\n",r1.m,r1.ErrM,r1.q,r1.ErrQ);
      if(!LogLog){
	for(int i=0;i<NMax;i++){
	  fprintf(USCITA,"%lf %lf %lf\n",v1->Abscissa(CoordY,i),v1->Val(CoordY,i),v1->Abscissa(CoordY,i)*r1.m+r1.q);
	}
      }
      else{
	for(int i=0;i<NMax;i++){
	  if(v1->Abscissa(CoordY,i) <= 0.) continue;
	  if(v1->Val(CoordY,i) <= 0.) continue;
	  fprintf(USCITA,"%lf %lf %lf\n",log10(v1->Abscissa(CoordY,i)),log10(v1->Val(CoordY,i)),log10(v1->Abscissa(CoordY,i))*r1.m+r1.q);
	}
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"FitPara") ){
      FILE *USCITA = fopen(Uscita,"w");
      PARABOLA p1 = v1->ParabolaSegnale(CoordY,NElMin,NElMax,LogLog);
      printf("%.2g + %.2g x + %.2g x^2\n",p1.a0,p1.a1,p1.a2);
      printf("Minimum %.2g \n",p1.Minimo);
      fprintf(USCITA,"#%.2g + %.2g x + %.2g x^2\n",p1.a0,p1.a1,p1.a2);
      fprintf(USCITA,"#Minimum %.2g \n",p1.Minimo);
      for(int i=0;i<NMax;i++){
	fprintf(USCITA,"%lf %lf %lf\n",v1->Abscissa(CoordY,i),v1->Val(CoordY,i),v1->Abscissa(CoordY,i)*p1.a1+v1->Abscissa(CoordY,i)*v1->Abscissa(CoordY,i)*p1.a2+p1.a0);
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"FitExp") ){
      FILE *USCITA = fopen(Uscita,"w");
      RETTA r1 = v1->InterExpSegnale(CoordY,NElMin,NElMax,LogLog);
      printf("y=(%.3g+-%.3g)exp(x/(%.3g+-%.3g))\n",log10(r1.q),log10(r1.ErrQ),1./r1.m,r1.ErrM/SQR(r1.m));
      printf("decay t_d = %.4g +- %.4g\n",1./r1.m,r1.ErrM/SQR(r1.m));
      fprintf(USCITA,"#y=(%.3g+-%.3g)exp(x/(%.3g+-%.3g))\n",log10(r1.q),log10(r1.ErrQ),1./r1.m,r1.ErrM/SQR(r1.m));
      fprintf(USCITA,"#decay t_d = %.4g +- %.4g\n",1./r1.m,r1.ErrM/SQR(r1.m));
      for(int i=0;i<NMax;i++){
	fprintf(USCITA,"%lf %lf %lf\n",v1->Abscissa(CoordY,i),v1->Val(CoordY,i),exp(r1.q+r1.m*v1->Abscissa(CoordY,i)));
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"Prob2Nrg") ){
      FILE *File2Write = fopen(Uscita,"w");
      if(File2Write == NULL){
	printf("Could not open %s\n",Uscita);
	return 1;
      }
      MOMENTI Mom = v1->DistrSegnale(NElMin,NElMax,IfNorm);
      v1->NormalizzaInter();
      for(int v=0;v<NBin;v++){
	double x = Mom.Min + v*Mom.Delta;
	if(v1->pInter(v) > 0.)
	  fprintf(File2Write,"%lf %lf\n",x,-log(v1->pInter(v)));
      }
      fclose(File2Write);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"min") ){
      double Min = v1->Val(CoordY,0);
      int nMin = 0;
      for(int n=0;n<NMax;n++){
	if(Min > v1->Val(CoordY,n)){
	  Min = v1->Val(CoordY,n);
	  nMin = n;
	}
      }
      printf("%lf %lf\n",v1->Val(CoordX,nMin),v1->Val(CoordY,nMin));
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"mm") ){
      printf("Number of points per sample: %d\n",NMedia);
      int NSample = v1->MediaMobSegnale(NMedia);
      double Deltax = 1.;//(v1->pxMax()-v1->pxMin())/(double)NSample;
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0)
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      for(int i=0;i<NSample;i++){
	fprintf(USCITA,"%lf %lf %lf\n",v1->Abscissa(CoordY,i*NMedia)+Deltax,v1->pPunti(i),v1->pPuntiErr(i));
      }
      fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"WeightHisto") ){
      v1->WeightHistoSign(NFile);
      // FILE *USCITA;
      // if((USCITA = fopen(Uscita,"w"))==0)
      // 	printf("Non s'apre %s\n Permessi?\n",Uscita);
      // for(int i=0;i<NSample;i++){
      // 	fprintf(USCITA,"%lf %lf %lf\n",v1->Abscissa(CoordY,i*NMedia)+Deltax,v1->pPunti(i),v1->pPuntiErr(i));
      // }
      // fclose(USCITA);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"int")){
      double Int = v1->IntSegnale();
      FILE *USCITA;
      if((USCITA = fopen(Uscita,"w"))==0)
	printf("Non s'apre %s\n Permessi?\n",Uscita);
      printf("%d %d %d\n",NElMin,NElMax,v1->pNMax());
      for(int i=NElMin;i<NElMax;i++){
	fprintf(USCITA,"%lf %lf\n",v1->Abscissa(CoordY,i),v1->pPunti(i));
      }
      printf("%lf\n",Int);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"mom")){
      for(int v=0;v<NVar;v++){
	v1->Punta(v);
	Mom = v1->DistrGaussSegnale(NElMin,NElMax,IfNorm);
	printf("%d) %lf pm %lf\n",v,Mom.Uno,Mom.Due);
      }
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"exec")){
      if(IfUser){
	printf("Inserire la formula per x: ");
	scanf("%s",Executex);
	printf("Inserire la formula per y: ");
	scanf("%s",Executey);
      }
      FILE *Ciccia = fopen(Uscita,"w");
      for(int i=NElMin;i<NElMax;i++){
	for(int v=0;v<NVar;v++){
	  double Val = v1->Val(v,i);
	  if(v==CoordX) 
	    Val = v1->ExecFormula(v1->Val(CoordX,i),v1->Val(CoordY,i),Executex);
	  else if(v==CoordY) 
	    Val = v1->ExecFormula(v1->Val(CoordX,i),v1->Val(CoordY,i),Executey);
	  fprintf(Ciccia,"%.5g ",Val);
	}
	fprintf(Ciccia,"\n");
      }
      fclose(Ciccia);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"widom")){
      FILE *Widom = fopen(Uscita,"w");
      for(int i=NElMin;i<NElMax;i++){
	if(v1->Val(3,i)*v1->Val(1,i) > 0.){
	  double Arg = v1->Val(3,i)/(v1->Val(1,i));
	  double ChemPot = v1->Val(0,i) + log(Arg);
	  fprintf(Widom,"%lf %lf\n",v1->Val(0,i),ChemPot);
	}
      }
      fclose(Widom);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"DecayTime")){
      double *st = new double[NVar];
      double *sw = new double[NVar];
      FILE *DecayTime = fopen(Uscita,"w");
      for(int f=0;f<NVar;f++){
	for(int i=NElMin;i<NElMax;i++){
	  st[f] = v1->Val(0,i)*v1->Val(f,i);
	}
	st[f] /= (double)v1->pNMax();
      }
      m1->Spettro(st,sw,NVar);
      for(int f=0;f<NVar;f++){
	fprintf(DecayTime,"%d %g %g\n",f,sw[f],st[f]);
      }
      fclose(DecayTime);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"reverse") ){
      v1->Reverse();
      v1->ScriviTutto(Uscita,0,NElMin,NElMax);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"DoubleDist") ){
      v1->DoubleDistFluct();
      v1->ScriviPunti(Uscita);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"smooth") ){
      v1->Smooth(Param,CoordY,NElMin,NElMax);
      v1->Smooth(Param,CoordY,NElMin,NElMax);
      v1->Smooth(Param,CoordY,NElMin,NElMax);
      v1->ScriviPunti(Uscita);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"smoothGauss") ){
      v1->SmoothGauss(Param,CoordY,NElMin,NElMax);
      v1->ScriviTutto(Uscita,0,NElMin,NElMax);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"salva") ){
      v1->ScriviTutto(Uscita,0,NElMin,NElMax);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"tecplot") ){
      v1->TecPlot(Uscita);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"txvl") ){
      v1->ExportTxvl(Uscita,NElMin,NElMax);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"ResBulk") ){
      v1->RescaleToBulk(Uscita);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"conv") ){
      FILE *FWrite = fopen(Uscita,"w");
      double x = 0.;
      for(int n=0;n<NMax;n++){
	x += .1;
	fprintf(FWrite,"%lf %lf\n",x,v1->Val(0,n)*1e12);
      }
      fclose(FWrite);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"salvacol") ){
      v1->ScriviFile(Uscita,CoordY,0,NElMin,NElMax);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"col") ){
      printf("Number of column: ");
      int NCol=0;
      scanf("%d",&NCol);
      v1->Punta(NCol);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"print") ){
      v1->Print();
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"vet")){
      Vettore v(1.,0.,1.);
      Vettore u(0.,1.,0.);
      Vettore w(0.,0.,1.);
      Vettore n(3);
      v.VetV(&u);
      //printf("%lf %lf %lf\n",n.Angle(v,u)*RAD_DEG,n.SinAngle(v,u),n.CosAngle(v,u));
      n.Print();
      break;
       if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"OrdAv")){
      double *Distr = new double[NBin];
      if(!IfBorder){
	v1->pMinMaxGlob(NElMin,NElMax);
	v1->pGlobBorder(Border,Border+1,Border+2,Border+3);
	printf("%lf %lf\n",Border[0],Border[1]);
	Border[1] = Border[1] - Border[0];
	Border[1] = 1./(1.001*Border[1]);
      }
      else{
	Border[1] = Border[1] - Border[0];
	Border[1] = 1./(1.001*Border[1]);
      }
      printf("%lf %lf\n",Border[0],Border[1]);
      v1->AverageOrdinate(NElMin,NElMax,Distr,Border);
      FILE *WFile = fopen(Uscita,"w");
      for(int b=0;b<NBin;b++){
	double x = b/(Border[1]*NBin) + Border[0];
	fprintf(WFile,"%lf %lf\n",x,Distr[b]);
      }
      free(Distr);
      fclose(WFile);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"WAv")){
      MOMENTI m1 = v1->WeightAverageSet(CoordY,NElMin,NElMax);
      printf("Mean %lf +- %lf\n",m1.Uno,m1.Due);
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"q") ){}
    else{
      printf("Comando sbagliato, riprovare\n");
      IfUser = 1;
    }
  }
  printf("Te se qe ve be te ne?\n");
  delete m1;
  return 0;
}
