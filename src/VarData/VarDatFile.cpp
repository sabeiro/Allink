#include "../include/VarDatFile.h"
#include <stdarg.h>

VarDatFile::VarDatFile(int ExtNMax,int ExtNVar,int ExtNBin){
  Init();
  NMax = ExtNMax;
  NMaxPunti = NMax;
  NVar = ExtNVar;
  NBin = ExtNBin;
  Allocate();
  sp = st[0];
}
VarDatFile::VarDatFile(double **ExtSt,int ExtNMax,int ExtNVar,int ExtNBin){
  Init();
  NMax = ExtNMax;
  NMaxPunti = NMax;
  NVar = ExtNVar;
  NBin = ExtNBin;
  Allocate();
  st = ExtSt;
  sp = st[0];
}
// VarDatFile::VarDatFile(char *file,int NNBin){
//   NMax = 0;
//   NVisMin = 0;
//   NVar=0;
//   NBin = NNBin;
//   FILE *InFile;
//   if((InFile = fopen(file,"r"))==0){
//     printf("Non s'apre %s\n",file);
//     exit(0);
//   }
//   Init();
//   int NCol = 0;
//   NVisMax = NMax = FileInfo(InFile,&NCol);
//   NVar = NCol;
//   Allocate();
//   char cLine[STRLEN];
//   double *dLine = (double *)calloc(NVar,sizeof(double));
//   for(int v=0;v<NVar;v++) SetNMax[v] = NMax;
//   for(int l=0; l<NMax;l++){
//     fgets(cLine,512,InFile);
//     if(ReadLine(cLine,dLine)){l--;continue;}
//     for(int v=0;v<NVar;v++){
//       st[v][l] = dLine[v];
//       if(l==0){
// 	sMin[v] = st[v][l];
// 	sMax[v] = st[v][l];
//       }
//       if(sMin[v] > st[v][l])
// 	sMin[v] = st[v][l];
//       if(sMax[v] < st[v][l])
// 	sMax[v] = st[v][l];
//     }
//   }
//   free(dLine);
//   if(NVar > 1){
//     CoordY = 1;
//     CoordX = 0;
//     RefAbsc[CoordX] = REF_ASC;
//     for(int v=1;v<NVar;v++)
//       RefAbsc[v] = CoordX;
//   }
//   else{
//     CoordY = 0;
//     RefAbsc[CoordY] = REF_SEQ;
//   }
//   fclose(InFile);
//   sp = st[0];
// }
VarDatFile::VarDatFile(char **FileList,int *Pos,int NFile,int NNBin){
  NMax=0;
  NVar=0;
  Init();
  NBin = NNBin;
  int *NMaxFile = (int *)calloc(NFile,sizeof(int));
  int *NColFile = (int *)calloc(NFile+1,sizeof(int));
  NColFile[0] = 0;
  for(int f=0;f<NFile;f++){
    FILE *FileIn = fopen(FileList[Pos[f]],"r");
    if(FileIn == 0){
      printf("File %s not found\n",FileList[Pos[f]]);
      return;
    }
    int NCol = 0;
    NMaxFile[f] = FileInfo(FileIn,&NCol);
    NColFile[f+1] = NCol + NColFile[f];
    if(NMaxFile[f] > NMax) NMax = NMaxFile[f];
    NMaxPunti = NMax;
    NVar += NCol;
    fclose(FileIn);
  }
  Allocate();
  char cLine[STRLEN];
  double *dLine = (double *)calloc(NVar,sizeof(double));
  for(int f=0;f<NFile;f++){
    for(int v=NColFile[f];v<NColFile[f+1];v++)
      SetNMax[v] = NMaxFile[f];
    FILE *FileIn = fopen(FileList[Pos[f]],"r");
    for(int l=0;l<NMaxFile[f];l++){
      if(!fgets(cLine,512,FileIn))break;
      if(ReadLine(cLine,dLine)){l--;continue;}
      for(int v=NColFile[f],vv=0;v<NColFile[f+1];v++,vv++){
	st[v][l] = dLine[vv];
      }
    }
    fclose(FileIn);
    if(NColFile[f+1]-NColFile[f] == 1){
      RefAbsc[NColFile[f]] = NVar;
    }
    else {
      for(int v=NColFile[f];v<NColFile[f+1];v++){
	RefAbsc[v] = NColFile[f];
      }
    }
  }
  if(NVar == 1) Punta(0);
  else Punta(1);
  free(dLine);
  free(NMaxFile);
  free(NColFile);
}
void VarDatFile::Allocate(){
  dInter = (double *)calloc(NBin,sizeof(*dInter));
  dError = (double *)calloc(NBin,sizeof(*dError));
  dInter1 = (double *)calloc(NBin,sizeof(*dInter));
  sp = (double *)malloc(sizeof(double));
  sType = (int *)calloc(NVar,sizeof(int));
  SetNMax = (int *)calloc(NVar,sizeof(int));
  sColor = (int *)calloc(NVar,sizeof(int));
  sMin = (double *)calloc(NVar,sizeof(double));
  sMax = (double *)calloc(NVar,sizeof(double));
  sColor = (int *)calloc(NVar,sizeof(int));
  RefAbsc = (int *)calloc(NVar,sizeof(int));
  st = (double **)calloc(NVar+1,sizeof(double));
  if(st == NULL){
    printf("Could not allocate st \n");
    exit(0);
  }
  //if( !st){printf("st non s'alloca\n");return ;}
  for(int v=0;v<NVar+1;v++){
    st[v] = (double *)calloc(NMax,sizeof(double));
    if(st[v] == NULL){
      printf("Could not allocate st[%d] \n",v);
      exit(0);
    }
    //if( !st[v]){printf("st non s'alloca\n");return ;}
  }
  for(int n=0;n<NMax;n++){
    st[NVar][n] = (double)n;
  }
  //sc = NULL;
  //sc = (char **)calloc(NVar,sizeof(char));
  //if(!sc){printf("sc non s'alloca\n");return ;}
  for(int v=0;v<NVar;v++){
    //sc[v] = (char *)calloc(120,sizeof(char));
    //if(!sc[v]){printf("st non s'alloca\n");return ;}
  }
}
VarDatFile::~VarDatFile(){
  for(int v=0;v<NVar;v++){
    free(st[v]);
    //if(st[v]){printf("st %d is not freed\n",v);}
  }
  free(st);
  //if(st){printf("st is not freed\n");}
  for(int v=0;v<NVar;v++){
    //free(sc[v]);
    //if(sc[v] != NULL){printf("sc %d is not freed\n",v);}
  }
  //free(sc);//if(sc != NULL){printf("sc is not freed\n");}
  free(sMin);
  free(sMax);
  free(sType);
  free(SetNMax);
  free(sColor);
  free(RefAbsc);
  sp = NULL;
  free(sp);
  free(dInter);
  free(dInter1);
  PuntiFree();
}
void VarDatFile::Init(){
  IfPuntiAlloc = 1;
  IfPuntiFree = 0;
}
//FIXME: there is still color, type... to be allocated
int VarDatFile::Aggiungi(char *file){
  FILE *InFile;
  if((InFile = fopen(file,"r"))==NULL){
    printf("Non s'apre %s\n",file);
    return 1;
  }
  char cLine[512];
  double **stTemp = (double **)calloc(NVar,sizeof(double));
  if(stTemp == NULL){printf("stTemp non s'alloca\n");return 0;}
  for(int v=0;v<NVar;v++){
    stTemp[v] = (double *)calloc(NMax,sizeof(double));
    if( stTemp[v] == NULL){printf("stTemp non s'alloca\n");return 0;}
  }
  for(int v=0;v<NVar;v++){
    for(int n=0;n<NMax;n++){
      stTemp[v][n] = st[v][n];//sign[v*NMax+n];
    }
  }
  int NVarOld = NVar;
  int NCol = 0;
  int nPoint = FileInfo(InFile,&NCol);
  NVar += NCol;
  if(nPoint > NMax) NMax = nPoint;
  NMaxPunti = NMax;
  {
    for(int v=0;v<NVarOld;v++){
      free(st[v]);
      //free(sc[v]);
    }
    free(st);
    //free(sc);
    st = (double **)calloc(NVar,sizeof(double));
    if( st  == NULL){printf("st non s'alloca\n");return 0;}
    for(int v=0;v<NVar;v++){
      st[v] = (double *)calloc(NMax,sizeof(double));
      if( st[v] == NULL){printf("st non s'alloca\n");return 0;}
    }
  }
  // {
  //   st = (double **)realloc(st,NVar*sizeof(double));
  //   for(int v=0;v<NVar;v++){
  //     st[v] = (double *)realloc(st[v],NMax*sizeof(double));
  //   }
  // }
  double *dLine = (double *)calloc(NVar,sizeof(double));
  for(int l=0; l<NMax;l++){
    fgets(cLine,512,InFile);
    if(ReadLine(cLine,dLine))continue;
    for(int v=NVarOld;v<NVar;v++){
      st[v][l] = dLine[v-NVarOld];
      if(sMin[v] > st[v][l])
	sMin[v] = st[v][l];
      if(sMax[v] < st[v][l])
	sMax[v] = st[v][l];
    }
  }
  GlobMax = sMax[0];
  GlobMin = sMin[0];
  for(int v=0;v<NVar;v++){
    if(GlobMax < sMax[v])
      GlobMax = sMax[v];
    if(GlobMin > sMin[v])
      GlobMin = sMin[v];
  }
  for(int n=0;n<NMax;n++)
    for(int v=0;v<NVarOld;v++)
      st[v][n] = stTemp[v][n];
  for(int v=0;v<NVarOld;v++)
    free(stTemp[v]);
  free(stTemp);
  free(dLine);
  fclose(InFile);
  PuntiFree();
  //printf("Uscita\n");
  return 0;
}
int VarDatFile::FileInfo(FILE *InFile,int *NCol){
  int nPoint = 0;
  char cLine[STRSIZE];
  for(int k=0;!(fgets(cLine,STRSIZE,InFile)==NULL);k++){
    if(k == 0 && cLine[0] == '#'){
      sprintf(Header,"%s",cLine);
      continue;
    }
    int iLen = (int) (strlen(cLine));
    if(cLine[0] != '#' && cLine[0] != '\n')
      nPoint++;
    for(int i=0,var=1;i<iLen;i++){
      if(cLine[i] == '#'){
	break;//Comment
      }
      if(cLine[i] == '\n') break;
      if(i>0){
	if( (cLine[i-1] == ' ' || cLine[i-1] == '\t' ) &&
	    (cLine[i] != ' ' && cLine[i] != '\t' ))
	  var++;
      }
      if(*NCol < var) *NCol = var;
    }
  }
  rewind(InFile);
  return nPoint;
}
int VarDatFile::ReadLine(char *cLine,double *sLine){
  int iLen = (int) (strlen(cLine));
  if(cLine[0] == ' ' || cLine[0] == '\t'){
    for(int i=0;i<iLen;i++){
      if(cLine[i] != ' ' && cLine[i] != '\t'){
  	cLine += i;
  	break;
      }
    }
  }
  sscanf(cLine,"%lf",sLine);
  for(int i=0,v=1;i<iLen;i++){
    if(cLine[i]=='#') return 1;
    if(v>=NVar) return 0;
    if(cLine[i]==' ' || cLine[i] =='\n' || cLine[i]=='\t'){
      sscanf(cLine+i,"%lf",sLine+v);
      //printf("%d) %d %s %lf\n",v,i,cLine+i,st[v]);
      i++;
      v++;
    }
  }
  return 0;
}
int VarDatFile::ReadLine(FILE *InFile, const char *Format, ...){
  va_list ap;
  int NVariable;
  char line[256];
  va_start(ap, Format);
  fgets(line, 256, InFile);
  NVariable = vsscanf(line, Format, ap);
  va_end(ap);
  return NVariable;
}
void VarDatFile::AverageOrdinate(int ElMin, int ElMax,double *Distr,double *xBound){
  double *Count = (double *)calloc(NBin,sizeof(double));
  for(int n=0;n<NMax;n++){
    for(int v=0;v<NVar;v++){
      if(IsAbscissa(v)) continue;
      if(isnan(st[v][n]))continue;
      int bx = (int)((Abscissa(v,n) - xBound[0])*xBound[1]*NBin);
      if(bx < 0 || bx >= NBin) continue;
      Distr[bx] += st[v][n];
      Count[bx] += 1.;
    }
  }
  for(int bx=0;bx<NBin;bx++){
    Distr[bx] /= Count[bx] > 0. ? Count[bx] : 1.;
  }
  free(Count);
}
MOMENTI VarDatFile::DistrSegnale(int ElMin, int ElMax,int IfNorm){
  return Distribuzione(sp+ElMin,ElMax-ElMin,dInter,NBin,IfNorm);
}
MOMENTI VarDatFile::DistrLogSegnale(int ElMin, int ElMax,int IfNorm){
  PuntiAlloc();
  for(int n=0;n<NMax;n++){
    Punti[n] = log10(sp[n]);
  }
  return Distribuzione(Punti+ElMin,ElMax-ElMin,dInter,NBin,IfNorm);
}
MOMENTI VarDatFile::DistrExpSegnale(int ElMin, int ElMax,int IfNorm){
  PuntiAlloc();
  for(int n=0;n<NMax;n++){
    Punti[n] = exp(sp[n]);
  }
  return Distribuzione(Punti+ElMin,ElMax-ElMin,dInter,NBin,IfNorm);
}
MOMENTI VarDatFile::DistrSegnale(int ElMin, int ElMax,double *Border,int IfNorm){
  return Distribuzione(sp+ElMin,ElMax-ElMin,dInter,NBin,Border,IfNorm);
}
MOMENTI VarDatFile::DistrSignErr(int ElMin, int ElMax,double *Border,int IfNorm){
  return DistrErr(sp+ElMin,ElMax-ElMin,dInter,dError,NBin,Border,IfNorm);
}
void VarDatFile::DistrSignSample(int ElMin, int ElMax,double **Distr,int NSample,int IfNorm,double *xBound){
  DistrSample(st[0]+ElMin,sp+ElMin,ElMax-ElMin,Distr,NBin,NSample,IfNorm,xBound);
}
MOMENTI VarDatFile::DistrGaussSegnale(int ElMin,int ElMax,int IfNorm){
  return DistribuzioneGauss(sp+ElMin,ElMax-ElMin,dInter,dInter1,NBin,IfNorm);
}
MOMENTI VarDatFile::WeightAverageSet(int CoordY,int ElMin, int ElMax){
  int CoordX = RefAbsc[CoordY];
  MOMENTI m1 = WeightAverage(st[CoordX]+ElMin,sp+ElMin,ElMax-ElMin);
  return m1;
}
int VarDatFile::ElabSegnale(int ElMin,int ElMax){
  PuntiAlloc();
  //  Elabora(sp+ElMin,Punti,ElMax-ElMin);
  ElabSt(sp+ElMin,Punti,ElMax-ElMin);
  return ElMax-ElMin;
}
int VarDatFile::AutocorSegnale(int ElMin,int ElMax){
  PuntiAlloc();
  Autocor(sp+ElMin,Punti,ElMax-ElMin);
  NMaxPunti = NMax/2;
  return NMaxPunti;
}
int VarDatFile::SpettroSegnale(int ElMin,int ElMax){
  PuntiAlloc();
  Spettro(sp+ElMin,Punti,ElMax-ElMin);
  NMaxPunti = NMax;
  return (int)( (NMax));
}
void VarDatFile::SpeLine(int ElMin,int ElMax,int NBin,double *Spe){
#ifdef USE_FFTW
  double InvNBin = 1./(double)(NBin);
  fftw_complex *out = (fftw_complex *)fftw_malloc(NBin*sizeof(fftw_complex));
  fftw_complex *in = (fftw_complex *)fftw_malloc(NBin*sizeof(fftw_complex));
  fftw_plan plan = fftw_plan_dft_1d(NBin,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan plan2 = fftw_plan_dft_1d(NBin,out,in,FFTW_BACKWARD,FFTW_ESTIMATE);
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v))continue;
    for(int n=0;n<NBin;n++) in[n][0] = st[v][n];
    fftw_execute(plan);
    for(int n=0;n<NBin;n++){
      Spe[n] += SQR(out[n][0]) + SQR(out[n][1]);
    }
  }
  fftw_free(out);
  fftw_free(in);
#else
  printf("fftw not present\n");
#endif
}
int VarDatFile::NormalizzaInter(){
  Normalizza(dInter,NBin);
}
int VarDatFile::NormalizzaSegnale(int ElMin,int ElMax){
  int Campioni = NormalizeArea(sp+ElMin,ElMax-ElMin);
  return Campioni;
}
bool VarDatFile::RadiceSegnale(){
  PuntiAlloc();
  Radice(sp,Punti,NMax);
  return 1;
}
double VarDatFile::IntSegnale(){
  PuntiAlloc();
  return Integrazione(sp,Punti,NMax);
}
double VarDatFile::SumSegnale(int CoordY,int ElMin,int ElMax){
  double Resp = 0.;
  int CoordX = RefAbsc[CoordY];
  double Dx = (st[CoordX][ElMax]-st[CoordX][ElMin]);
  double NumInv = 1./(double)(ElMax-ElMin);
  for(int i=ElMin;i<ElMax;i++){
    Resp += st[CoordY][ElMax];
  }
  return Resp*Dx*NumInv;
}
void VarDatFile::SommaSegnali(){
  PuntiAlloc();
  for(int nm=0;nm<NMax;nm++){ 
    Punti[nm] = 0.;
    for(int nv=0;nv<NVar;nv++){
      if(IsAbscissa(nv)) continue;
      Punti[nm] += st[nv][nm];
    }
  }
}
double VarDatFile::VarieSegnale(int ElMin,int ElMax){
  if(NVar < 2) return 0.;
  PuntiAlloc();
  double *temp1 = (double *)malloc(NMax*sizeof(double));
  double *temp2 = (double *)malloc(NMax*sizeof(double));
  double VolA=0.,VolB=0.;
  double VolFracA=0.,VolFracB=0.;
  double SqrGradCoeff=0.;
  for(int n=0;n<NMax;n++){
    temp1[n] = 0.;
    temp2[n] = 0.;
  }
  double Coeff[4] = {1.,0.,0.,1.};
  double MaxA=0.,MaxB=0.;
  for(int n=0;n<NMax;n++){
    if(MaxA < st[1][n]) MaxA = st[1][n];
    if(MaxB < st[2][n]) MaxB = st[2][n];
  }
  for(int n=0;n<NMax;n++){
    if(st[1][n] > MaxA*.5)
      VolA += 1.;
    if(st[2][n] > MaxB*.5)
      VolB += 1.;
  }
  VolFracA = VolA / (VolA + VolB);
  VolFracB = VolB / (VolA + VolB);
  SqrGradCoeff = 19.13/(24.*VolFracA*VolA)+3.04/(24.*VolFracB*VolB);
  sp = st[1];
  DerO4(sp+ElMin,temp1,ElMax-ElMin);
  sp = st[2];
  DerO4(sp+ElMin,temp2,ElMax-ElMin);
  for(int i=0;i<4;i++){
    Coeff[i] *= SqrGradCoeff;
  }
  for(int n=0;n<NMax;n++){
    Punti[n] = Coeff[0]*temp1[n]*temp1[n];
    Punti[n] += Coeff[1]*temp1[n]*temp2[n];
    Punti[n] += Coeff[2]*temp2[n]*temp1[n];
    Punti[n] += Coeff[3]*temp2[n]*temp2[n];
    Punti[n] /= 2.;
  }
  double Risp=0.;
  for(int n=75;n<125;n++){
    Risp += Punti[n];
  }
  free(temp1);
  free(temp2);
  return Risp;
}
RETTA VarDatFile::InterRettSegnale(int CoordY,int ElMin,int ElMax,int LogLog){
  RETTA r1;
  double *Puntix;
  double *Puntiy;
  int CoordX = RefAbsc[CoordY];
  if(DIS_IF_TYPE(LogLog,DIS_LOGX)){
    Puntix = (double *)calloc(ElMax-ElMin,sizeof(*Puntix));
    for(int i=ElMin;i<ElMax;i++){
      if(st[CoordX][i] > 0.)
	Puntix[i-ElMin] = log10(st[CoordX][i]);
    }
  }
  else{
    Puntix = st[CoordX]+ElMin;
  }
  if(DIS_IF_TYPE(LogLog,DIS_LOGY)){
    Puntiy = (double *)calloc(ElMax-ElMin,sizeof(*Puntiy));
    for(int i=ElMin;i<ElMax;i++){
      if(st[CoordY][i] > 0.)
	Puntiy[i-ElMin] = log10(st[CoordY][i]);
    }
  }
  else{
    Puntiy = st[CoordY]+ElMin;
  }
  r1 = InterRett(Puntix,Puntiy,ElMax-ElMin);
  if(DIS_IF_TYPE(LogLog,DIS_LOGX)) free(Puntix);
  if(DIS_IF_TYPE(LogLog,DIS_LOGY)) free(Puntiy);
  return r1;
}
RETTA VarDatFile::InterExpSegnale(int CoordY,int ElMin,int ElMax,int LogLog){
  RETTA r1;
  int CoordX = RefAbsc[CoordY];
  r1 = InterExp(st[CoordX]+ElMin,st[CoordY]+ElMin,ElMax-ElMin);
  return r1;
}
MOMENTI VarDatFile::InterGaussSegnale(int CoordY,int ElMin,int ElMax,int LogLog){
  int CoordX = RefAbsc[CoordY];
  MOMENTI m1 = InterGauss(st[CoordX]+ElMin,st[CoordY]+ElMin,ElMax-ElMin);
  return m1;
}
PARABOLA VarDatFile::ParabolaSegnale(int CoordY,int ElMin,int ElMax,int LogLog){
  PARABOLA p1;
  double *Puntix;
  double *Puntiy;
  int CoordX = RefAbsc[CoordY];
  if(DIS_IF_TYPE(LogLog,DIS_LOGX)){
    Puntix = (double *)calloc(ElMax-ElMin,sizeof(*Puntix));
    for(int i=ElMin;i<ElMax;i++){
      if(st[CoordX][i] > 0.)
	Puntix[i] = log(st[CoordX][i]);
    }
  }
  else{
    Puntix = st[CoordX];
  }
  if(DIS_IF_TYPE(LogLog,DIS_LOGY)){
    Puntiy = (double *)calloc(ElMax-ElMin,sizeof(*Puntix));
    for(int i=ElMin;i<ElMax;i++){
      if(st[CoordY][i] > 0.)
	Puntiy[i] = log10(st[CoordY][i]);
    }
  }
  else{
    Puntiy = st[CoordY];
  }
  p1 = MinimoParabola(Puntix+ElMin,Puntiy+ElMin,ElMax-ElMin);
  if(DIS_IF_TYPE(LogLog,DIS_LOGX))
    free(Puntix);
  if(DIS_IF_TYPE(LogLog,DIS_LOGY))
    free(Puntiy);
  return p1;
}
int VarDatFile::MediaMobSegnale(int n){
  PuntiAlloc();
  NMaxPunti = MediaMobile(sp,NMax,Punti,PuntiErr,n);
  return NMaxPunti;
}
int VarDatFile::WeightHistoSign(int NHisto){
  int NBin = NMax;
  double **Histo = (double **)calloc(NHisto,sizeof(double));
  for(int h=0;h<NHisto;h++){
    Histo[h] = (double *)calloc(NBin,sizeof(double));
  }
  for(int v=0,h=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    for(int b=0;b<NBin;b++){
      Histo[h][b] = st[v][b];
    }
    h++;
  }
  /// Find the borders globally
  pMinMaxGlob(0,NBin);  
  double Border[2] = {xGlobMin,xGlobMax};
  double tolerance = 0.00001;
  double *NanoDist = (double *)calloc(NHisto,sizeof(double));
  double *kSpring = (double *)calloc(NHisto,sizeof(double));
  for(int h=0;h<NHisto;h++){
    kSpring[h] = Histo[h][NBin-2];
    NanoDist[h] = Histo[h][NBin-1];
  }
  // NanoDist[1]=3.5;
  // NanoDist[2]=4.0;
  // NanoDist[3]=5.07;                      
  // NanoDist[4]=5.6;                 
  // NanoDist[5]=6.4;
  // NanoDist[6]=7.2;                     
  // NanoDist[7]=8;
  // NanoDist[8]=8.8;                          
  // NanoDist[9]=9.6;
  // NanoDist[10]=10.4;                        
  // NanoDist[11]=10.76;
  // NanoDist[12]=11.56;               
  // NanoDist[13]=12.8;
  WeightHisto(Histo,Border,NBin-2,NHisto,tolerance,NanoDist,kSpring);
  for(int h=0;h<NHisto;h++)
    free(Histo[h]);
  free(Histo);
  free(NanoDist);
  free(kSpring);
}
int VarDatFile::CorrelaADuePunti(int n){
  PuntiAlloc();
  NMaxPunti = CorrelaDuePunti(sp,NMax,Punti,n);
  return NMaxPunti;
}
void VarDatFile::AutosimilaritaSegnale(int n){
  PuntiAlloc();
  Autosimilarita(sp,NMax,Punti,n);
  NMaxPunti = NMax;
}
void VarDatFile::CambiaNBin(int n){
  if(n < 0)
    return;
  NBin = n;
  dInter = (double *)realloc(dInter,NBin*sizeof(double));
  dInter1 = (double *)realloc(dInter1,NBin*sizeof(double));
}
void VarDatFile::Punta(int n){
  if(n < 0 || n > NVar) 
    return;
  sp = st[n];
}
void VarDatFile::Punta(double *spExt,int n){
  if(n < 0 || n >= NVar) 
    return;
  spExt = st[n];
}
void VarDatFile::PuntaInt(double *spExt){
  sp = spExt;
}
void VarDatFile::Punta(double **ExtSt,int n){
  if(n <0 || n>NVar) 
    return;
  st = ExtSt;
  sp = ExtSt[n];
}
double VarDatFile::Val(int CoordY,int n){
#ifdef FILE_DEBUG
  if(n<0 || n>NMax-1 ){printf("Wrong row\n");return 0.;}
  if(CoordY<0 || v>NBin-1){printf("Wrong col\n"); return 0.;}
#endif
  return st[CoordY][n];
}
int VarDatFile::pNRow(int CoordY){
  return SetNMax[CoordY];
}
double VarDatFile::Abscissa(int CoordY,int n){
#ifdef FILE_DEBUG
  if(CoordY<0 || n>NMax-1){ printf("Wrong row\n");return 0.;}
#endif
  return st[RefAbsc[CoordY]][n];
}
double VarDatFile::pPunti(int n){
#ifdef FILE_DEBUG
  if(n<0 || n>NMax-1){ printf("Wrong row\n");return 0.;}
  if(IfPuntiAlloc) return 0.;
#endif
  return Punti[n];
}
double VarDatFile::pPuntiErr(int n){
#ifdef FILE_DEBUG
  if(n<0 || n>NMax-1){ printf("Wrong row\n");return 0.;}
  if(IfPuntiAlloc) return 0.;
#endif
  return PuntiErr[n];
}
char *VarDatFile::PrintHeader(){
  return Header;
}
double VarDatFile::PuntiMin(){
#ifdef FILE_DEBUG
  if(IfPuntiAlloc) return 0.;
#endif
  double Min = Punti[0];
  for(int p=0;p<NMaxPunti;p++)
    if(Min > Punti[p])
      Min = Punti[p];
  return Min;
}
double VarDatFile::PuntiMax(){
#ifdef FILE_DEBUG
  if(IfPuntiAlloc) return 0.;
#endif
  double Max = Punti[0];
  for(int p=0;p<NMaxPunti;p++)
    if(Max < Punti[p])
      Max = Punti[p];
  return Max;
}
void VarDatFile::CambiaPunti(){
  if(!IfPuntiFree)
    return;
  Punti1 = (double *)realloc(Punti1,NMaxPunti*sizeof(double));
  for(int n=0;n<NMax;n++){
    Punti1[n] = Punti[n];
  }
  sp = Punti1;
}
void VarDatFile::Print(){
  for(int c=0;c<NVar;c++){
    printf("%d ",RefAbsc[c]);
  }
  printf("\n--------------\n");
  for(int r=0;r<NMax;r++){
    printf("%d) ",r);
    for(int c=0;c<NVar;c++){
      if(IsAbscissa(r)) printf("%.5g ",Abscissa(r,c));
      printf("%.5g ",st[r][c]);
    }
    printf("\n");
  }
}
void VarDatFile::Sort(){
  for(int i=0;i<NMax;i++){
    for(int j=i;j>0;j--){
      if(sp[j] < sp[j-1])
	for(int n=0;n<NVar;n++){
	  Swap(j,j-1,st[n]);
	}
      else
	break;
    }
  }
}
void VarDatFile::ScriviPunti(char *file){
  FILE *SALVA;
  if((SALVA = fopen(file,"w"))==0)
     printf("Non s'apre %s, Permessi?\n",file);
  for(int n=0;n<NMaxPunti;n++)
    fprintf(SALVA,"%lf %lf\n",Punti1[n],Punti[n]);
  fclose(SALVA);
}
void VarDatFile::ScriviFile(char *file,int CoordY,int LogLog,int NVisMin,int NVisMax){
  FILE *SALVA;
  int CoordX = RefAbsc[CoordY];
  if((SALVA = fopen(file,"w"))==0)
    printf("Non s'apre %s, Permessi?\n",file);
  for(int n=NVisMin;n<NVisMax;n++){
    if(DIS_IF_TYPE(LogLog,DIS_LOGLOG))
      fprintf(SALVA,"%lf %lf \n",log10(st[CoordX][n]),log10(st[CoordY][n]));
    else 
      fprintf(SALVA,"%lf %lf \n",st[CoordX][n],st[CoordY][n]);
  }
  fclose(SALVA);
}
void VarDatFile::ScriviTutto(char *file,int LogLog,int NVisMin,int NVisMax){
  FILE *SALVA;
  if((SALVA = fopen(file,"w"))==0)
     printf("Non s'apre %s, Permessi?\n",file);
  if(DIS_IF_TYPE(LogLog,DIS_LOGLOG)){
    for(int n=NVisMin;n<NVisMax;n++){
      for(int v=0;v<NVar;v++){
	fprintf(SALVA,"%.4g ",log10(st[v][n]));
      }
      fprintf(SALVA,"\n");
    }
  }
  else {
    for(int n=NVisMin;n<NVisMax;n++){
      for(int v=0;v<NVar;v++){
	fprintf(SALVA,"%.4g ",st[v][n]);
      }
      fprintf(SALVA,"\n");
    }
  }
  fclose(SALVA);
}
void VarDatFile::RescaleToBulk(char *FName){
  double *BulkVal = (double *)calloc(NVar,sizeof(double));
  double Count = 0;
  int NCol = 0;
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    NCol++;
  }
  int NThick = (int)(NMax/(double)NCol);
  int Nx = NMax;
  for(int v=0;v<NVar;v++){
    for(int nx=NMax/2;nx<3*(NMax/4);nx++){
      if(isnan(st[v][nx]))continue;
      BulkVal[v] += st[v][nx];
      Count += 1.;
    }
  }
  for(int v=0;v<NVar;v++){
    BulkVal[v] /= Count > 0. ? Count/(double)NVar : 1.;
    BulkVal[v] = fabs(BulkVal[v]) > 0. ? 1./BulkVal[v] : 1.;
  }
  FILE *FWrite = fopen(FName,"w");
  for(int nx=0;nx<NMax;nx++){
    fprintf(FWrite,"%lf ",Abscissa(1,nx));
    for(int v=1;v<NVar;v++){
      fprintf(FWrite,"%lf ",st[v][nx]*BulkVal[v]);
    }
    fprintf(FWrite,"\n");
  }
  fclose(FWrite);
  free(BulkVal);
}
void VarDatFile::ExportTxvl(char *FName,int NElMin,int NElMax){
  FILE *FOut = fopen(FName,"w");
  int p=0;
  pMinMaxGlob(0,NMax);
  double Dx = 1./(xGlobMax-xGlobMin);
  Dx = 0.05/1.;
  double Dz = 1./(GlobMax-GlobMin);
  fprintf(FOut,"# l(1 1 1) v[%d] d[part] \n",NMax/2);
  double *Cm = (double *)calloc(NVar,sizeof(double));
  double CmTot = 0.;
  double *Bulk = (double *)calloc(NVar,sizeof(double));
  int NBulk = NElMax-NElMin;
  for(int v=0;v<NVar;v++){
    for(int n=0;n<NMax;n++){
      if(n >= NBulk) Cm[v] += Val(v,n);
    }
  }
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    Cm[v] /= (double)(NMax-NBulk);
    CmTot += Cm[v];
  }
  CmTot = CmTot/(double)(NVar-1);
  for(int v=0;v<NVar;v++){
    Bulk[v] = .1/fabs(Cm[v]-CmTot);
  }
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    double OffSet = .6 - .2*(v-1);
    for(int n=0;n<NMax;n++){
      double Zed = (Val(v,n)-Cm[v])*Bulk[v]+OffSet;
      if(Abscissa(v,n)*Dx > 1.) continue;
      fprintf(FOut,"{t[%d %d %d] x(%lf %lf %lf)}\n",p++,2,2,Abscissa(v,n)*Dx,.5,Zed);
    }
  }
  free(Cm);
  free(Bulk);
  fclose(FOut);
}
void VarDatFile::TecPlot(char *FName){
  FILE *TecPlot = fopen(FName,"w");
  double *BulkVal = (double *)calloc(NVar,sizeof(double));
  double Count = 0;
  int NCol = 0;
  pMinMaxGlob(0,NMax);
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    NCol++;
  }
  int NThick = (int)(NMax/(double)NCol);
  int Nx = 2*(int)(NMax/3.);
  int Ny = Nx;
  for(int v=0;v<NVar;v++){
    for(int nx=NMax/2;nx<3*(NMax/4);nx++){
      if(isnan(st[v][nx]))continue;
      BulkVal[v] += st[v][nx];
      Count += 1.;
    }
  }
  for(int v=0;v<NVar;v++){
    BulkVal[v] /= Count > 0. ? Count/(double)NVar : 1.;
    printf("%d %lf\n",v,BulkVal[v]);
    BulkVal[v] = fabs(BulkVal[v]) > 0. ? 1./BulkVal[v] : 1.;
  }
  fprintf(TecPlot,"VARIABLES = \"R\", \"Z\", \"percentage\"\n");
  fprintf(TecPlot,"ZONE J=%d, K=%d, F=POINT\n",Nx,Ny);
  double Int = .5*(xGlobMax-xGlobMin)/(double)(Ny);
  for(int nx=0,v=0;nx<Nx;nx++){
    for(int ny=0;ny<Ny;ny++){
      if((ny%(NThick/2))==0) v++;
      if(v >= NVar) v=0;
      if(IsAbscissa(v)) v++;
      double z = Val(v,nx)*BulkVal[v];
      if(isnan(z)) z = 1.;
      fprintf(TecPlot,"%lf %lf %lf\n",Abscissa(v,nx),Int*ny,z);
    }
  }
  fclose(TecPlot);
  free(BulkVal);
}
void VarDatFile::PuntiAlloc(){
  if(IfPuntiAlloc){
    Punti = (double *)calloc(NMaxPunti,sizeof(*Punti));
    if( Punti  == NULL){printf("Punti non s'alloca\n");return ;}
    Punti1 = (double *)calloc(NMaxPunti,sizeof(*Punti));
    if( Punti1  == NULL){printf("Punti non s'alloca\n");return ;}
    PuntiErr = (double *)calloc(NMaxPunti,sizeof(*Punti));
    if( PuntiErr  == NULL){printf("Punti non s'alloca\n");return ;}
    IfPuntiAlloc = 0;
    IfPuntiFree = 0;
  }  
}
void VarDatFile::ImpSequence(int v){
  RefAbsc[v] = NVar;
}
void VarDatFile::ImpSequence(){
  for(int v=0;v<NVar;v++){
    RefAbsc[v] = NVar;
  }
}
void VarDatFile::ImpCoordX(int vSet,int vAbs){
  if(vAbs > NVar) return;
  if(IsAbscissa(vSet)) return;
  // one set refers to the sequence
  if(NVar == 1) RefAbsc[vSet] = NVar;
  // set to the sequence set
  else if(vAbs == REF_SEQ) RefAbsc[vSet] = NVar;
  else RefAbsc[vSet] = vAbs;
}
void VarDatFile::ImpCoordX(int vAbs){
  if(vAbs > NVar || vAbs < REF_ASC) return;
  // one set refers to the sequence
  if(NVar == 1){ RefAbsc[0] = NVar;return;}
  for(int v=0;v<NVar;v++){
    if(vAbs == REF_SEQ) RefAbsc[v] = NVar;
    else RefAbsc[v] = vAbs;
  }
}
void VarDatFile::ImpCoordY(int Ext){
  if(Ext < NVar && Ext >= 0){
    sp = st[Ext];
  }
}
void VarDatFile::PuntiFree(){
  if(IfPuntiFree){
    free(Punti);
    IfPuntiAlloc = 1;
  }
}
double VarDatFile::pMaxGlob(int NVisMin,int NVisMax){
  double Resp = st[0][0];
  // set the initial value
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v))continue;
    Resp = st[v][NVisMin];
  }
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v))continue;
    for(int n=NVisMin;n<NVisMax;n++){
      if(n > pNRow(v)) continue;
      if(isnan(st[v][n]))continue;
      if(Resp < st[v][n])
	Resp = st[v][n];
    }
  }
  GlobMax = Resp;
  return Resp;
}
double VarDatFile::pMinGlob(int NVisMin,int NVisMax){
  double Resp = st[0][0];
  // set the initial value
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v))continue;
    Resp = st[v][NVisMin];
  }
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v))continue;
    for(int n=NVisMin;n<NVisMax;n++){
      if(n > pNRow(v)) continue;
      if(isnan(st[v][n])) continue;
      if(Resp > st[v][n])
	Resp = st[v][n];
    }
  }
  return Resp;  
}
void VarDatFile::pMinMaxGlob(int NVisMin,int NVisMax){
  double Border[4];
  // set the initial value
  for(int v=0;v<NVar;v++){
    if(!IsAbscissa(v)){
      Border[0] = Abscissa(v,NVisMin);
      Border[1] = Abscissa(v,NVisMin);
      Border[2] = st[v][NVisMin];
      Border[3] = st[v][NVisMin];
      break;
    }
  }
  // find the borders
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    for(int n=NVisMin;n<NVisMax;n++){
      if(n >= pNRow(v)) continue;
      if(isnan(st[v][n])) continue;
      if(Border[0] > Abscissa(v,n)) Border[0] = Abscissa(v,n);
      if(Border[1] < Abscissa(v,n)) Border[1] = Abscissa(v,n);
      if(Border[2] > st[v][n]) Border[2] = st[v][n];
      if(Border[3] < st[v][n]) Border[3] = st[v][n];
    }
  }
  xGlobMin = Border[0];
  xGlobMax = Border[1];
  GlobMin = Border[2];
  GlobMax = Border[3];
}
double VarDatFile::pxMaxGlob(int NVisMin,int NVisMax){
  double Resp = st[0][NVisMin];
  // set the initial value
  for(int v=0;v<NVar;v++){
    if(RefAbsc[v] == NVar) return NVisMin;
    if(IsAbscissa(v)){
      Resp = st[v][NVisMin];
      break;
    }
  }
  for(int v=0;v<NVar;v++){
    if(!IsAbscissa(v))continue;
    for(int n=NVisMin;n<NVisMax;n++){
      if(n > pNRow(v)) continue;
      if(isnan(st[v][n]))continue;
      if(isinf(st[v][n]))continue;
      if(Resp < st[v][n])
	Resp = st[v][n];
    }
  }
  xGlobMax = Resp;
  printf("%lf\n",Resp);
  return Resp;
}
double VarDatFile::pxMinGlob(int NVisMin,int NVisMax){
  double Resp = st[0][NVisMin];
  // set the initial value
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)){
      Resp = st[v][NVisMin];
      break;
    }
  }
  for(int v=0;v<NVar;v++){
    if(!IsAbscissa(v)) continue;
    for(int n=NVisMin;n<NVisMax;n++){
      if(n >= pNRow(v)) continue;
      if(Resp > st[v][n]){
	Resp = st[v][n];
      }
    }
  }
  xGlobMin = Resp;
  return Resp;  
}
double VarDatFile::pMaxGlobLog(int NVisMin,int NVisMax){
  double Resp = st[0][NVisMin];
  int IfContinue = 1;
  // set the initial value
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    for(int n=0;n<NMax;n++){
      if(st[v][n] <= 0.) continue;
      Resp = st[v][n];
      IfContinue = 0;
      break;
    }
    if(IfContinue) break;
  }
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    for(int n=NVisMin;n<NVisMax;n++){
      if(st[v][n] <= 0.) continue;
      if(n >= pNRow(v)) continue;
      if(Resp < st[v][n]) Resp = st[v][n];
    }
  }
  return log10(Resp > 0. ? Resp : 1.);
}
double VarDatFile::pMinGlobLog(int NVisMin,int NVisMax){
  double Resp = st[0][NVisMin];
  int IfContinue = 1;
  // set the initial value
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    for(int n=0;n<NMax;n++){
      if(st[v][n] <= 0.) continue;
      Resp = st[v][n];
      IfContinue = 0;
      break;
    }
    if(!IfContinue) break;
  }
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    for(int n=NVisMin;n<NVisMax;n++){
      if(st[v][n] <= 0.) continue;
      if(n > pNRow(v)) continue;
      if(Resp > st[v][n]) Resp = st[v][n];
    }
  }
  return log10(Resp > 0. ? Resp : 1.);  
}
double VarDatFile::pMax(int CoordY,int NVisMin,int NVisMax){
  double Resp = st[CoordY][NVisMin];
  for(int n=NVisMin;n<NVisMax;n++){
    if(Resp < st[CoordY][n])
      Resp = st[CoordY][n];
  }
  return Resp;
}
double VarDatFile::pMin(int CoordY,int NVisMin,int NVisMax){
  double Resp = st[CoordY][NVisMin];
  for(int n=NVisMin;n<NVisMax;n++){
    if(Resp > st[CoordY][n])
      Resp = st[CoordY][n];
  }
  return Resp;
}
double VarDatFile::pMaxLog(int CoordY,int NVisMin,int NVisMax){
  double Resp = st[CoordY][NVisMin];
  for(int n=NVisMin;n<NVisMax;n++){
    if(st[CoordY][n] <= 0.) continue;
    if(Resp < st[CoordY][n])
      Resp = st[CoordY][n];
  }
  return log10(Resp > 0. ? Resp : 1.);
}
double VarDatFile::pMinLog(int CoordY,int NVisMin,int NVisMax){
  double Resp = 0.;
  for(int n=NVisMin;n<NVisMax;n++){
    if(st[CoordY][n] <= 0.) continue;
    Resp = st[CoordY][n];
    break;
  }
  for(int n=NVisMin;n<NVisMax;n++){
    if(st[CoordY][n] <= 0.) continue;
    if(Resp > st[CoordY][n])
      Resp = st[CoordY][n];
  }
 return log10(Resp > 0. ? Resp : 1.);  
}
void VarDatFile::setXFormula(char *str){
  strcpy(XFormula,str);
}
void VarDatFile::setYFormula(char *str){
  strcpy(YFormula,str);
}
void VarDatFile::Reverse(){
  double Temp = 0.;
  for(int n=0;n<NMax/2;n++){
    for(int v=0;v<NVar;v++){
      Temp = st[v][n];
      st[v][n] = st[v][NMax-n-1];
      st[v][NMax-n-1] = Temp;
    }
  }
}
int VarDatFile::Smooth(double Fact,int CoordY,int NVisMin,int NVisMax){
  int NOrder = 3+1;
  //still not working for Fact != 1.
  Fact = 1.;
  int NOut = (int)(NMax*Fact);
  NMaxPunti = NOut;
  int NSmooth = 2;
  double *sw = (double *)calloc(NOut,sizeof(double));
  pMinMaxGlob(NVisMin,NVisMax);
  PuntiAlloc();
  double Max=st[CoordY][0];
  double Min=st[CoordY][0];
  for(int n=0;n<NMax;n++){
    if(Max < st[CoordY][n]) Max = st[CoordY][n];
    if(Min > st[CoordY][n]) Min = st[CoordY][n];
  }
  double DeltaIn=(Max-Min)/(double)(NMax-1);
  double DeltaOut=(Max-Min)/(double)(NOut-1);
  double *dArray = (double *)calloc(NMax+2*NOrder+2,sizeof(double));
  int CoordX = RefAbsc[CoordY];
  for(int vi=0;vi<NMax;vi++){
    Punti1[vi] = st[CoordX][vi];
  }
  for(int vi=0;vi<=NMax+2*NOrder;vi++){
    if(vi<NOrder){
      dArray[vi] = Min;
    }
    else if(vi > NMax){
      dArray[vi] = (vi-NOrder+2)*DeltaIn+Min;
    }
    else {
      dArray[vi] = (vi-NOrder+1)*DeltaIn+Min;
    }
  }
  Punti[0] = st[CoordY][0];
  Punti[NOut-1] = st[CoordY][NMax-1];
  for(int vo=1;vo<NOut-1;vo++){
    Punti[vo] = 0.;
    double x = DeltaOut*vo+Min;
    for(int vi=vo-1,vn=0;vi<vo+NOrder+1;vi++){
      vn = vi;
      if(vi < 0){
	vn = NMax + vi;
      }
      if(vi >= NMax){
	vn = vi - NMax;
      }
      double Blendx = Blend(dArray,x,vn,NOrder);
      Punti[vo] += Blendx*st[CoordY][vn];
    }
  }
  free(dArray);
  return NOut;
}
void VarDatFile::SmoothGauss(double Fact,int CoordY,int NVisMin,int NVisMax){
  //still not working for Fact != 1.
  Fact = 1.;
  int NOut = (int)(NMax*Fact);
  NMaxPunti = NOut;
  pMinMaxGlob(NVisMin,NVisMax);
  PuntiAlloc();
  int CoordX = RefAbsc[CoordY];
  for(int n=0;n<NMax;n++){
    Punti[n] = st[CoordY][n];
    Punti1[n] = st[CoordX][n];
  }
  Matrice Mask(5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  int NDim = 1;
  int IfMinImConv = 0;
  for(int v=0;v<NVar;v++){
    if(IsAbscissa(v)) continue;
    Mask.ConvoluteMatrix(st[v],NVisMax-NVisMin,NDim,IfMinImConv);
  }
}
void VarDatFile::DoubleDistFluct(){
  PuntiAlloc();
  SigErr(NVar < 8,"Two files are with the coordinate positions are expected\n");
  for(int n=0;n<NMax;n++){
    if(n > pNRow(0)) continue;
    if(n > pNRow(1)) continue;
    double Dist = 0.;
    for(int d=0;d<3;d++){
      Dist += SQR(st[1+d][n] - st[5+d][n]);
    }
    Dist = sqrt(Dist);
    Punti[n] = Dist;
    Punti1[n] = (double)n;
  }
}
