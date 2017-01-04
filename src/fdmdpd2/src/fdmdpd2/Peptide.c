#include "fdmdpd2.h"
static void PepVv1(struct beads *restrict b,int MemPos, int NPart);
static int BraketPos(const char *str,const char *mask,int *sPos,int *sLen){  
  if(!strncmp(str,mask,strlen(mask))) return 1;
  char *pInit = strpbrk(str,mask);
  int InitPos = pInit - str + 1;
  char *pOpen = strpbrk(str+InitPos,"([{");
  char *pClose = strpbrk(str+InitPos,"}])");
  if(pOpen == NULL || pClose == NULL) return 1;
  *sPos = pOpen - str +1;
  *sLen = pClose - str;
  //printf("%d %d %d %s\n",pInit,sPos,sLen,str+InitPos);
  return 0;
}
int Fetch(const char *str,const char *mask,const char *fmt, ... ){
  int sLen = 0;
  int sPos = 0;
  int IfContinue = BraketPos(str,mask,&sPos,&sLen);
  if(IfContinue) return 0;
  char *Field = (char *)calloc(sLen-sPos+1,sizeof(char));
  strncpy(Field,str+sPos,sLen-sPos);
  va_list args;
  va_start(args,fmt);
  vsscanf(Field,fmt,args);
  //  printf("Sotto %s trovato %s come %s\n",mask,Field,fmt);
  va_end(args);
  free(Field);
  return sPos;
}
/* int Fetch(const char *str,const char *mask,const char *fmt, ... ){ */
/*   int sLen = 0; */
/*   int sPos = 0; */
/*   if(BraketPos(str,mask,&sPos,&sLen)) return 0; */
/*   char *Field = strndupa(str + sPos, sLen); */
/*   printf("%d\n",sLen); */
/*   va_list args; */
/*   va_start(args,fmt); */
/*   vsscanf(Field,fmt,args); */
/*   printf("Sotto %s trovato %s come %s\n",mask,Field,fmt); */
/*   va_end(args); */
/*   return sPos; */
/* } */
/** Header for the outputNano.dat file */
void PepPrintHeader(FILE *FH, const struct beads *b){
  for(int n=0;n<b->NPep;n++){
    struct PEPTIDE *Pep = b->Pep + n;
    fprintf(FH,"# Pep g(%.2f %.2f %.2f) i[%d] d[%d %d] fn{%s}\n",Pep->Rad,Pep->Hei,Pep->kEl,n,Pep->NCircle,Pep->NSide,Pep->ArchFile);
  }   
}
static void PepDefine(FILE *FH,struct PEPTIDE *Pep,int nPep){
  double Geom[3];
  int NEdge[2];
  int NLink=0;
  char cLine[LINE_MAX];
  char Arch[60];
  fgets(cLine, LINE_MAX, FH);
  if( !Fetch(cLine,"i","%d",&Pep->Id));
  else Pep->Id = nPep;
  if( !Fetch(cLine,"g","%lf %lf %lf",Geom,Geom+1,Geom+2)){
    Pep->Rad = Geom[0];
    Pep->Hei = Geom[1];  
    Pep->kEl = Geom[2];
  }
  else{
    Pep->Rad = 1.;
    Pep->Hei = 4.;  
    Pep->kEl = 2000.;
  }
  if( !Fetch(cLine,"d","%d %d",NEdge,NEdge+1)){
    Pep->NSide = NEdge[0];
    Pep->NCircle = NEdge[1];
  }
  else{
    Pep->NSide = 10;
    Pep->NCircle = 12;
  }
  // Bulding list
  if( Fetch(cLine,"fn","%s",Arch) ){
    sprintf(Pep->ArchFile,"%s",Arch);
    FILE *Ciccia = fopen(Arch,"r");
    if(Ciccia == NULL) 
      fatal(EINVAL,"Architecture file %s not found",Arch);
    char LineArch[120];
    while(fgets(LineArch,120,Ciccia)){
      if(LineArch[0] != '#')
	NLink++;
    }
    Pep->NLink = NLink;
    Pep->Link = (int *)calloc(2*NLink,sizeof(*Pep->Link));
    if (Pep->Link == NULL) novm("Pep->Link");
    Pep->Elong = (double *)calloc(NLink,sizeof(*Pep->Elong));
    if (Pep->Elong == NULL) novm("Pep->Elong");
    Pep->kSpr = (double *)calloc(NLink,sizeof(*Pep->kSpr));
    if (Pep->kSpr == NULL) novm("Pep->kSpr");
    rewind(Ciccia);
    for(int l=0;fgets(LineArch,120,Ciccia);l++){
      //fgets(LineArch,120,Ciccia);
      if(LineArch[0] == '#'){
	l--;
	continue;
      }
      if(3==sscanf(LineArch,"%d %d %lf %lf",Pep->Link+l*2,Pep->Link+l*2+1,Pep->Elong+l,Pep->kSpr+l)){
	Pep->kSpr[l] = Pep->kEl;
      }
      //printf("--%d %d %lf %lf\n",Pep->Link[l*2],Pep->Link[l*2+1],Pep->Elong[l],Pep->kSpr[l]);
    }
    fclose(Ciccia);
  }
  else fatal(EINVAL,"Architecture file not specified");

}
void PepLoad(struct beads *b, FILE *FH){
  char buf[LINE_MAX];
  long fpos = ftell(FH);
  debug("Pep] Reading configurations");
  for (fgets(buf, sizeof(buf), FH); strcasestr(buf, "# Pep") == buf; b->NPep++)
    fgets(buf, sizeof(buf), FH);
  fseek(FH, fpos, SEEK_SET);
  assert(b->NPep > 0);
  b->Pep = calloc(b->NPep, sizeof(*b->Pep));
  if (b->Pep == NULL) novm("Pep");
  for(int n=0;n<b->NPep;n++){
    //fgetpos(FH,&PosTemp);
    PepDefine(FH,b->Pep+n,n);
  }
}
void PepFree(struct beads *b){
  if(b->NPep == 0) return;
  debug("Pep] Freeing");
  for (int n = 0; n < b->NPep; n++){
    free(b->Pep[n].Link);
    free(b->Pep[n].Elong);
  }
  free(b->Pep);
}
// Calculation of the elastic force bewteen the monomers composing the inclusion 
void PepCalcForces(struct beads *restrict b,int MemPos,const char *name,int NPart){
  int PepId = 0;
  sscanf(name+3,"%d",&PepId);
  struct PEPTIDE *Pep = b->Pep;
  for(int p=0;p<b->NPep;p++){
    Pep = b->Pep + p;
    if(Pep->Id == PepId) break;
  }
  for(int l=0;l<Pep->NLink;l++){
    b->e[1] += harmonic_spring(b,MemPos,Pep->Link[l*2+0],Pep->Link[l*2+1],Pep->kSpr[l],Pep->Elong[l],STRESS_CALC);
  }
  PepVv1(b,MemPos,NPart);
}
void CLinksCalcForces(struct beads *restrict b,int MemPos,const char *name,int NPart){
  int PepId = 0;
  sscanf(name+5,"%d",&PepId);
  struct PEPTIDE *Pep = b->Pep;
  for(int p=0;p<b->NPep;p++){
    Pep = b->Pep + p;
    if(Pep->Id == PepId) break;
  }
  for(int l=0;l<Pep->NLink;l++){
    b->e[1] += harmonic_spring(b,MemPos,Pep->Link[l*2+0],Pep->Link[l*2+1],Pep->kSpr[l],Pep->Elong[l],STRESS_CALC);
  }
  PepVv1(b,MemPos,NPart);
}
static void PepVv1(struct beads *restrict b,int MemPos, int NPart){
  struct TENS_PROF *Tens = &(b->Tens);
  double xMean[3] = {0.,0.,0.};
  VEC *xx = b->x_intra + MemPos;
  for(int p=0;p<NPart;p++){
    for(int d=0;d<3;d++){
      xMean[d] += xx[p][d];
    }
  }
  for(int d=0;d<3;d++){
    Tens->RefPos[d] = xMean[d]/(double)NPart;
  }
}
 
