#include "fdmdpd2.h"
static VEC2 *BPos = NULL;
static VEC2 *BFor = NULL;
// Serial
void ExtCalcForces(struct beads *restrict b){
  for(int e=0;e<b->NExt;e++){
    struct PEPTIDE *Ext = b->Ext;
    for(int l=0;l<Ext->NLink;l++){
      int p1 = Ext->Link[l*2+0];
      int p2 = Ext->Link[l*2+1];
      b->e[1] += harmonic_spring(b,0,Ext->Link[l*2+0],Ext->Link[l*2+1],Ext->kSpr[l],Ext->Elong[l],STRESS_OMIT);
    }
  }
}
// Parallel
/* void ExtCalcForces(struct beads *b){ */
/*   if(b->NExt < 1) return; */
/*   VEC2 *xv = b->xv + col_index * b->groupsize; */
/*   VEC  *f  = b->f  + col_index * b->groupsize; */
/*   VEC  *n  = b->nx + col_index * b->groupsize; */
/*   PASSPORT *passport = b->passport + col_index * b->groupsize; */
/*   if (BPos == NULL) novm("BPos"); */
/*   if (BFor == NULL) novm("BFor"); */
/*   for (int i = 0; i < b->groupsize; i++) { */
/*     int id = GET_ID(passport[i]); */
/*     if (!GET_EXISTS(passport[i])) continue; */
/*     for (int d = 0; d < 3; d++) { */
/*       BPos[id][d] = xv[i][d] + n[i][d] * b->l[d]; */
/*       BFor[id][d] =  f[i][d]; */
/*     } */
/*   } */
/*   MPI_Reduce(ismaster ? MPI_IN_PLACE : BPos, ismaster ? BPos : NULL, */
/*   	     ARRAY_SIZE(*xv) * b->nN, MPI_DOUBLE, MPI_SUM, 0, comm_grid); */
/*   MPI_Reduce(ismaster ? MPI_IN_PLACE : BFor, ismaster ? BFor : NULL, */
/*   	     ARRAY_SIZE(*xv) * b->nN, MPI_DOUBLE, MPI_SUM, 0, comm_grid); */
/*   double DistRel[3]; */
/*   double Dist = 0.; */
/*   for(int e=0;e<b->NExt;e++){ */
/*     struct PEPTIDE *Ext = b->Ext; */
/*     for(int l=0;l<Ext->NLink;l++){ */
/*       int p1 = Ext->Link[l*2+0]; */
/*       int p2 = Ext->Link[l*2+1]; */
/*       for(int d=0;d<3;d++){ */
/* 	DistRel[d] = BPos[p1][d] - BPos[p2][d]; */
/*       } */
/*       Dist = sqrt(SQR(DistRel[0])+SQR(DistRel[1])+SQR(DistRel[2])); */
/*       double Nrg = .5*Ext->kSpr[l]*SQR(Dist - Ext->Elong[l]); */
/*       double Force = -Ext->kSpr[l]*(Dist - Ext->Elong[l]); */
/*       printf("%d %d %lf %lf\n",p1,p2,Dist,Ext->Elong[l],Force); */
/*       b->e[1] += Nrg;continue; */
/*       for(int d=0;d<3;d++){ */
/* 	double Comp = Force*DistRel[d]/Dist; */
/* 	BFor[p1][d] += Comp; */
/* 	BFor[p2][d] -= Comp; */
/*       } */
/*     } */
/*   } */
/* for (int i = 0; i < b->groupsize; i++) { */
/*   int id = GET_ID(passport[i]); */
/*   if (!GET_EXISTS(passport[i])) continue; */
/*   for (int d = 0; d < 3; d++) { */
/*     f[i][d] = BFor[id][d]; */
/*   } */
/* } */
/* } */

static void ExtDefine(FILE *FH,struct PEPTIDE *Ext,int nExt);
void ExtPrintHeader(FILE *FH, const struct beads *b){
  for(int n=0;n<b->NExt;n++){
    struct PEPTIDE *Ext = b->Ext + n;
    fprintf(FH,"# Ext i[%d] fn{%s}\n",Ext[n].Id,Ext->ArchFile);
  }   
}
void ExtLoad(struct beads *b, FILE *FH){
  char buf[LINE_MAX];
  long fpos = ftell(FH);
  debug("Ext] Reading configurations");
  for (fgets(buf, sizeof(buf), FH); strcasestr(buf, "# Ext") == buf; b->NExt++)
    fgets(buf, sizeof(buf), FH);
  fseek(FH, fpos, SEEK_SET);
  assert(b->NExt > 0);
  b->Ext = calloc(b->NExt, sizeof(*b->Ext));
  if (b->Ext == NULL) novm("Ext");
  for(int n=0;n<b->NExt;n++){
    //fgetpos(FH,&PosTemp);
    ExtDefine(FH,b->Ext+n,n);
  }
  VEC2 *xv = b->xv + col_index * b->groupsize;
  VEC  *f  = b->f  + col_index * b->groupsize;
  VEC  *n  = b->nx + col_index * b->groupsize;
  PASSPORT *passport = b->passport + col_index * b->groupsize;
  BPos = calloc(b->nN, sizeof(*BPos));
  if (BPos == NULL) novm("BPos");
  BFor = calloc(b->nN, sizeof(*BFor));
  if (BFor == NULL) novm("BFor");
}
static void ExtDefine(FILE *FH,struct PEPTIDE *Ext,int nExt){
  double Geom[3];
  int NEdge[2];
  int NLink=0;
  char cLine[LINE_MAX];
  char Arch[60];
  fgets(cLine, LINE_MAX, FH);
  if( !Fetch(cLine,"i","%d",&Ext->Id));
  else Ext->Id = nExt;
  if( Fetch(cLine,"fn","%s",Arch) ){
    sprintf(Ext->ArchFile,"%s",Arch);
    FILE *Ciccia = fopen(Arch,"r");
    if(Ciccia == NULL) 
      fatal(EINVAL,"Architecture file %s not found",Arch);
    char LineArch[120];
    while(fgets(LineArch,120,Ciccia)){
      if(LineArch[0] != '#')
	NLink++;
    }
    Ext->NLink = NLink;
    Ext->Link = (int *)calloc(2*NLink,sizeof(*Ext->Link));
    if (Ext->Link == NULL) novm("Ext->Link");
    Ext->Elong = (double *)calloc(NLink,sizeof(*Ext->Elong));
    if (Ext->Elong == NULL) novm("Ext->Elong");
    Ext->kSpr = (double *)calloc(NLink,sizeof(*Ext->kSpr));
    if (Ext->kSpr == NULL) novm("Ext->kSpr");
    rewind(Ciccia);
    for(int l=0;fgets(LineArch,120,Ciccia);l++){
      //fgets(LineArch,120,Ciccia);
      if(LineArch[0] == '#'){
	l--;
	continue;
      }
      if(3==sscanf(LineArch,"%d %d %lf %lf",Ext->Link+l*2,Ext->Link+l*2+1,Ext->Elong+l,Ext->kSpr+l)){
	Ext->kSpr[l] = Ext->kEl;
      }
      //printf("--%d %d %lf %lf\n",Ext->Link[l*2],Ext->Link[l*2+1],Ext->Elong[l],Ext->kSpr[l]);
    }
    fclose(Ciccia);
  }
  else fatal(EINVAL,"Architecture file not specified");
}
void ExtFree(struct beads *b){
  if(b->NExt == 0) return;
  debug("Ext] Freeing");
  for (int n = 0; n < b->NExt; n++){
    free(b->Ext[n].Link);
    free(b->Ext[n].Elong);
  }
  free(b->Ext);
}
