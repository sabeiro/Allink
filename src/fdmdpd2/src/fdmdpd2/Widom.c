/* $Id: Tens.c 275 2010-11-18 16:44:51Z hoemberg $ */

#include "fdmdpd2.h"
#include "rand.h"
extern int save_beads(const char *fn, struct beads *restrict b);
extern double mdpd_EnergySingleCh(struct beads *restrict b,int Ch);
static void SwapChains(struct beads *restrict b,int c1,int c2);
static void SwapChains(struct beads *restrict b,int c1,int c2);
static void SwapLast(struct beads *restrict b,int c1,VEC2 *Lastx,VEC *BfLast);
static void InsertRest(struct beads *restrict b,VEC2 *PosTemp,int nChain,int StartPos);
static void InsertChainCorona(struct beads *restrict b,int nChain);
static void InsertChainBiased(struct beads *restrict b,int nChain);
static void InsertChain(struct beads *restrict b,int nChain);
static void LoadLastChain(char *filename,VEC2 *Lastx,VEC *BfLast,double *Edge,int Nppc);
static void WidomLoopRem(struct beads *restrict b);
static void WidomLoopAdd(struct beads *restrict b);
static void WidomNrg(struct beads *b);
int MCStepIn(struct beads *b);
int MCStepOut(struct beads *b);
static int WidomLoopMC(struct beads *b);
static void ReallocBeads(struct beads *restrict b);
static int IfMetropolis(double NewNrg,double ChemPot);

enum WidomCalcMode{
  WIDOM_IN = 1,
  WIDOM_OUT,
  WIDOM_NRG,
  WIDOM_NRG_CH,
  WIDOM_MC,
};

cfg_opt_t Widom_opts[] = {
  CFG_INT("NStep", 0, CFGF_NONE),
  CFG_STR("CalcMode", 0, CFGF_NONE),
  CFG_FLOAT("ChemPot", 0, CFGF_NONE),
  CFG_END()
};
//-------Loading-configuration---------------
void WidomLoad(struct beads *b,cfg_t *cfg){
  struct WIDOM *Widom = &(b->Widom);
  debug("Widom] Reading configurations");
  cfg_t *cWidom = cfg_getsec(cfg, "Widom");
  Widom->NStep = cfg_getint(cWidom,"NStep");
  Widom->CalcMode = 3;
  int NSect[3] = {3,3,3};
  const char *CalcMode = cfg_getstr(cWidom,"CalcMode");
  if(!strcmp(CalcMode,"out")){
    Widom->CalcMode = WIDOM_OUT;
  }
  else if(!strcmp(CalcMode,"in")){
    Widom->CalcMode = WIDOM_IN;
    //DomDec(Widom->Dc,b->l,NSect,b->xv,b->local_n*b->N[0],1.);
  }
  else if(!strcmp(CalcMode,"nrg")){
    Widom->CalcMode = WIDOM_NRG;
    //DomDec(Widom->Dc,b->l,NSect,b->xv,b->local_n*b->N[0],1.);
  }
  else if(!strcmp(CalcMode,"nrgCh")){
    Widom->CalcMode = WIDOM_NRG_CH;
    //DomDec(Widom->Dc,b->l,NSect,b->xv,b->local_n*b->N[0],1.);
  }
  else if(!strcmp(CalcMode,"mc")){
    Widom->CalcMode = WIDOM_MC;
    //DomDec(Widom->Dc,b->l,NSect,b->xv,b->local_n*b->N[0],1.);
  }
  else {
    Widom->CalcMode = 0;
  }
  Widom->ChemPot = cfg_getfloat(cWidom,"ChemPot");
  Widom->NAlloc = (groups_per_row+groups_per_col)*b->groupsize;
}
//--------------Deletion-routines------------------
static void SwapChains(struct beads *restrict b,int c1,int c2){
  double PosTemp[3];
  double BackFoldTemp[3];
  int p0 = col_index * b->groupsize;
  int Nppc = b->N[0];
  for(int p=p0+c1*Nppc,p2=p0+c2*Nppc;p<p0+(c1+1)*Nppc;p++,p2++){
    for(int d=0;d<3;d++){
      PosTemp[d] = b->xv[p][d];
      b->xv[p][d] = b->xv[p2][d];
      b->xv[p2][d] = PosTemp[d];
      BackFoldTemp[d] = b->nx[p][d];
      b->nx[p][d] = b->nx[p2][d];
      b->nx[p2][d] = BackFoldTemp[d];
    }
  }
  p0 = (groups_per_row + row_index) * b->groupsize;
  for(int p=p0+c1*Nppc,p2=p0+c2*Nppc;p<p0+(c1+1)*Nppc;p++,p2++){
    for(int d=0;d<3;d++){
      PosTemp[d] = b->xv[p][d];
      b->xv[p][d] = b->xv[p2][d];
      b->xv[p2][d] = PosTemp[d];
      BackFoldTemp[d] = b->nx[p][d];
      b->nx[p][d] = b->nx[p2][d];
      b->nx[p2][d] = BackFoldTemp[d];
    }
  }
}
static void SwapLast(struct beads *restrict b,int c1,VEC2 *Lastx,VEC *BfLast){
  int Nppc = b->N[0];
  double PosTemp;
  double BackFoldTemp;
  int pOne = col_index * b->groupsize;
  int pTwo = (groups_per_row + row_index) * b->groupsize;
  for(int p=c1*Nppc,p2=0;p<(c1+1)*Nppc;p++,p2++){
    for(int d=0;d<3;d++){
      PosTemp = Lastx[p2][d];
      Lastx[p2][d] = b->xv[p+pOne][d];
      b->xv[p+pOne][d] = PosTemp;
      b->xv[p+pTwo][d] = PosTemp;
      BackFoldTemp = BfLast[p2][d];
      BfLast[p2][d] = b->nx[p+pOne][d];
      b->nx[p+pOne][d] = BackFoldTemp;
      b->nx[p+pTwo][d] = BackFoldTemp;
    }
  }
}
//----------Insertion-routines-------------------
/*find the diblock limit of the chains and write 
  the distribution of the position of that point.
 */
static void StudySystem(struct beads *restrict b){
  struct WIDOM *Widom = &(b->Widom);
  int Nppc = b->N[0];
  int first = col_index * b->groupsize;
  int last = (col_index + 1) * b->groupsize;
  for(int p=first+1,t=0;p<Nppc;p++,t++){
    PASSPORT Pass1 = b->passport[p];
    int type = (int)GET_TYPE(Pass1);
    if(type != GET_TYPE(b->passport[p-1])){
      Widom->DiblockLim = t;
      break;
    }
  }
  double Norm = Nppc/(double)(last-first);
  for(int p=first+Widom->DiblockLim;p<last;p+=Nppc){
    int vz = (int)(b->xv[p][NORMAL]/b->l[NORMAL]*Widom->NBins);
    Widom->Distr[vz] += Norm;
  }
}
/* Write the position of the remaining beads from the position 
   of StartPos choosing the gaussian distribution of a spring chain */
static void InsertRest(struct beads *restrict b,VEC2 *PosTemp,int nChain,int StartPos){
  struct WIDOM *Widom = &(b->Widom);
  int Nppc = b->N[0];
  for(int p=StartPos+1;p<Nppc;p++){
    for(int d=0;d<3;d++){
      PosTemp[p][d] = PosTemp[p-1][d] + rng_gaussian(rng,0.,Widom->Var);
    }
  }
  for(int p=StartPos-1;p>=0;p--){
    for(int d=0;d<3;d++){
      PosTemp[p][d] = PosTemp[p+1][d] + rng_gaussian(rng,0.,Widom->Var);
    }
  }
  int p0 = col_index * b->groupsize;
  int pInit = (nChain)*Nppc+p0;
  for(int p=pInit,p1=0;p<pInit+Nppc;p++,p1++){
    if(p1 <= Widom->DiblockLim) SET_TYPE(b->passport[p],0);
    else SET_TYPE(b->passport[p],1);
    for(int d=0;d<3;d++){
      b->xv[p][d] = PosTemp[p1][d];
      b->nx[p][d] = floor(b->xv[p][d] / b->l[d]);
      b->xv[p][d] = b->xv[p][d] - b->nx[p][d] * b->l[d];      
    }
  }
  p0 = (groups_per_row + row_index) * b->groupsize;
  pInit = (nChain)*Nppc+p0;
  for(int p=pInit,p1=0;p<pInit+Nppc;p++,p1++){
    if(p1 <= Widom->DiblockLim) SET_TYPE(b->passport[p],0);
    else SET_TYPE(b->passport[p],1);
    for(int d=0;d<3;d++){
      b->xv[p][d] = PosTemp[p1][d];
      b->nx[p][d] = floor(b->xv[p][d] / b->l[d]);
      b->xv[p][d] = b->xv[p][d] - b->nx[p][d] * b->l[d];
    }
  }
}
/* The position of the first particle of the chain is distributed
   in a corona around the nanoparticle */
static void InsertChainCorona(struct beads *restrict b,int nChain){
  struct NANO *Nano = b->Nano;
  struct WIDOM *Widom = &(b->Widom);
  int Nppc = b->N[0];
  VEC2 *PosTemp = (VEC2 *)calloc(Nppc,sizeof(VEC2));
  double RadMin = Nano->Radius*1.0;
  double RadDelta = Nano->Radius*1.5-Nano->Radius*1.0;
  double theta = 2.*rng_uniform(rng)*M_PI;
  double ran = 2.*rng_uniform(rng) - 1.;
  double Thick = pow(1./3.,3.*rng_uniform(rng))*RadDelta+RadMin;
  PosTemp[0][TANG1]  = Thick*sqrt(1.-ran*ran)*cos(theta) 
    + Nano->Pos[TANG1];
  PosTemp[0][TANG2]  = Thick*sqrt(1.-ran*ran)*sin(theta) 
    + Nano->Pos[TANG2];
  PosTemp[0][NORMAL] = Thick*ran + Nano->Pos[NORMAL];
  InsertRest(b,PosTemp,nChain,0);
  free(PosTemp);
}
/* The position of the DiblockLim is the hydrophobic core */
static void InsertChainBiased(struct beads *restrict b,int nChain){
  VEC2 *PosTemp = (VEC2 *)calloc(b->N[0],sizeof(VEC2));
  double ZedDelta = b->l[NORMAL]*.7-b->l[NORMAL]*.3;
  double ZedMin = b->l[NORMAL]*.3;
  PosTemp[0][TANG1]  = rng_uniform(rng)*b->l[TANG1];
  PosTemp[0][TANG2]  = rng_uniform(rng)*b->l[TANG2];
  PosTemp[0][NORMAL] = rng_uniform(rng)*ZedDelta+ZedMin;
  //PosTemp[0][NORMAL] = .5*b->l[NORMAL] + Sz*rng_gaussian(rng,0.,GaussVar);
  InsertRest(b,PosTemp,nChain,0);
  free(PosTemp);
}
static void InsertChainDiblock(struct beads *restrict b,int nChain){
  struct WIDOM *Widom = &(b->Widom);
  int Nppc = b->N[0];
  VEC2 *PosTemp = (VEC2 *)calloc(Nppc,sizeof(VEC2));
  double ZedDelta = b->l[NORMAL]*.5-b->l[NORMAL]*.3;
  double ZedMin = b->l[NORMAL]*.3;
  PosTemp[Widom->DiblockLim][TANG1] = rng_uniform(rng)*b->l[TANG1];
  PosTemp[Widom->DiblockLim][TANG2] = rng_uniform(rng)*b->l[TANG2];
  double Random = rng_uniform(rng);
  if(Random < .5){
    ZedDelta = b->l[NORMAL]*.7-b->l[NORMAL]*.5;
    ZedMin = b->l[NORMAL]*.5;
  }
  PosTemp[Widom->DiblockLim][NORMAL] = Random*ZedDelta+ZedMin;
  InsertRest(b,PosTemp,nChain,Widom->DiblockLim);
  free(PosTemp);
}
static double CorrProbability(double Rand,double *Distr,int NBins){
  double Cum = 0.;
  for(int b=0;b<NBins;b++){
    Cum += Distr[b];
    if(Cum > Rand)
      return (double) b;
  }
  return 0.;
}
static void InsertChainProbDistr(struct beads *restrict b,int nChain){
  struct WIDOM *Widom = &(b->Widom);
  int Nppc = b->N[0];
  VEC2 *PosTemp = (VEC2 *)calloc(Nppc,sizeof(VEC2));
  double Zed = b->l[NORMAL]*CorrProbability(rng_uniform(rng),Widom->Distr,Widom->NBins);
  PosTemp[Widom->DiblockLim][TANG1] = rng_uniform(rng)*b->l[TANG1];
  PosTemp[Widom->DiblockLim][TANG2] = rng_uniform(rng)*b->l[TANG2];
  PosTemp[Widom->DiblockLim][NORMAL] = Zed;
  InsertRest(b,PosTemp,nChain,Widom->DiblockLim);
  free(PosTemp);
}
static void InsertChain(struct beads *restrict b,int nChain){
  int Nppc = b->N[0];
  VEC2 *PosTemp = (VEC2 *)calloc(Nppc,sizeof(VEC2));
  PosTemp[0][TANG1]  = rng_uniform(rng)*b->l[TANG1];
  PosTemp[0][TANG2]  = rng_uniform(rng)*b->l[TANG2];
  PosTemp[0][NORMAL] = rng_uniform(rng)*b->l[NORMAL];
  InsertRest(b,PosTemp,nChain,0);
  free(PosTemp);
}
static void LoadLastChain(char *filename,VEC2 *Lastx,VEC *BfLast,double *Edge,int Nppc){
  FILE *FileChain = fopen(filename,"r");
  if(FileChain == NULL) fatal(ENOENT, "Where is LastChain.dat?");
  int type;
  for(int p=0;p<Nppc;p++){
    if (fscanf(FileChain,"%lf %lf %lf %lf %lf %lf %d\n", 
	       &Lastx[p][0], &Lastx[p][1], &Lastx[p][2],
	       &Lastx[p][3], &Lastx[p][4], &Lastx[p][5], &type) != 7)
      fatal(EINVAL,"LastChain file to check");
    for(int d=0;d<3;d++){
      BfLast[p][d] = floor(Lastx[p][d] / Edge[d]);
      Lastx[p][d] -= BfLast[p][d] * Edge[d];
    }
  }
}
//-----------------Widom-routines-------------------
//just the self interacting energy
//how to touch the pairlist to include the others?
static double CalcNrgCh(struct beads *restrict b,int Ch){
  memset(b->f, 0, (groups_per_row + groups_per_col) *
	 b->groupsize * sizeof(*b->f));
  memset(b->e, 0, ARRAY_SIZE(b->e) * sizeof(*b->e));
  //RigidCalcForces(b);
  //mdpd_nonbonded_forces(b, b->mdpd, NBL_INVALID);
  //RigidChNrg(b,Ch);
  double Nrg = mdpd_EnergySingleCh(b,Ch);
  MPI_Reduce(ismaster ? MPI_IN_PLACE : b->e, ismaster ? b->e : NULL,ARRAY_SIZE(b->e), MPI_DOUBLE, MPI_SUM, 0, comm_grid);
  return b->e[0];
}
static double CalcNrgMon(struct beads *restrict b,int p1){
  memset(b->f, 0, (groups_per_row + groups_per_col) *
	 b->groupsize * sizeof(*b->f));
  memset(b->e, 0, ARRAY_SIZE(b->e) * sizeof(*b->e));
  //RigidCalcForces(b);
  //mdpd_nonbonded_forces(b, b->mdpd, NBL_INVALID);
  //RigidMonNrg(b,p1);
  mdpd_EnergySingleMon(b,p1);
  MPI_Reduce(ismaster ? MPI_IN_PLACE : b->e, ismaster ? b->e : NULL,ARRAY_SIZE(b->e), MPI_DOUBLE, MPI_SUM, 0, comm_grid);
  return b->e[0];
}
static void OneChainMore(struct beads *restrict b){
  b->groupsize += b->N[0];
  b->local_n++;
  b->n[0]++;
}
static void OneChainLess(struct beads *restrict b){
  b->groupsize -= b->N[0];
  b->local_n--;
  b->n[0]--;
}
//check on an ideal gas if a get the right density \mu = \ln\rho\Lambda^3
//is different for removal and insertion
static int IfMetropolisIn(double NrgDiff,double ChemPot,double Vol,double NPart){
  double Arg = ChemPot - NrgDiff;
  if((Vol*exp(Arg))/(NPart+1.) > rng_uniform(rng)) return 1;
  else return 0;
}
static int IfMetropolisOut(double NrgDiff,double ChemPot,double Vol,double NPart){
  if((NPart/exp(ChemPot))*Vol > rng_uniform(rng)) return 1;
  else return 0;
}
static int IfMetropolisMove(double NrgDiff){
  if(exp(NrgDiff) > rng_uniform(rng)) return 1;
  else return 0;
}
int MCStepIn(struct beads *b){
  struct WIDOM *Widom = &(b->Widom);
  double Vol = b->l[0]*b->l[1]*b->l[2];
  int nChain = b->local_n-1;
  double NPart = (double)nChain;
  OneChainMore(b);
  ReallocBeads(b);
  if(b->NNano == 1){
    InsertChainCorona(b,nChain);
  }
  else InsertChain(b,nChain);
  double NewNrg = CalcNrgCh(b,nChain);
  if( !IfMetropolisIn(NewNrg,Widom->ChemPot,Vol,NPart) ){
    OneChainLess(b);
    return 1;
  }
  fprintf(stderr,"Convenient chain %d %lf\r",nChain,NewNrg);
  return 0;
}
int MCStepOut(struct beads *b){
  struct WIDOM *Widom = &(b->Widom);
  double Vol = b->l[0]*b->l[1]*b->l[2];
  int nChain = b->local_n-1;
  double NPart = (double)nChain;
  int LastCh = b->local_n;
  if(b->local_n < 0) return 1;
  int Ch = (int)(rng_uniform(rng)*b->local_n);
  double NewNrg = CalcNrgCh(b,Ch);
  if( !IfMetropolisOut(NewNrg,Widom->ChemPot,Vol,NPart) ){
    return 1;
  }
  fprintf(stderr,"Chain removed %d %d %lf\r",Ch,b->local_n,NewNrg);
  //remove chain
  SwapChains(b,Ch,LastCh);
  OneChainLess(b);
  return 0;
}
int MCStepMove(struct beads *b){
  struct WIDOM *Widom = &(b->Widom);
  double Vol = b->l[0]*b->l[1]*b->l[2];
  int nChain = b->local_n-1;
  double NPart = (double)nChain;
  int LastCh = b->local_n;
  if(b->local_n < 0) return 1;
  int Ch = (int)(rng_uniform(rng)*b->local_n);
  double NewNrg = CalcNrgMon(b,Ch);
  if( !IfMetropolisOut(NewNrg,Widom->ChemPot,Vol,NPart) ){
    return 1;
  }
  fprintf(stderr,"Chain removed %d %d %lf\r",Ch,b->local_n,NewNrg);
  //remove chain
  SwapChains(b,Ch,LastCh);
  OneChainLess(b);
  return 0;
}
static void PrintCh(struct beads *restrict b,FILE *FConf,int nChain){
  int p0 = (groups_per_row + row_index) * b->groupsize;
  int pInit = (nChain)*b->N[0]+p0;
  for(int p=pInit;p<pInit+b->N[0];p++){
    double x =  b->xv[p][0] + b->nx[p][0] * b->l[0];
    double y =  b->xv[p][1] + b->nx[p][1] * b->l[1];
    double z =  b->xv[p][2] + b->nx[p][2] * b->l[2];
    fprintf(FConf,"%lf %lf %lf %lf %lf %lf %d\n",
	    x,y,z,0.,0.,0.,GET_TYPE(b->passport[p]));
  }
}
//---------------Widom-loops--------------------
static void WidomLoopRem(struct beads *restrict b){
  struct nblist nbl;
  int NStep = b->local_n;
  FILE *NrgFile = fopen("WidomRawOut.dat","w");
  if (NrgFile == NULL) fatal(EIO, "Couldn't open WidomRawOut.dat");
  forces_calc(b);
  double OldNrg = b->e[0];
  for(int s=0;s<NStep;s++){
    double Nrg = CalcNrgCh(b,s);
    fprintf(NrgFile,"%d %lf %lf\n",s,b->e[0],b->e[1]);
    printf("%d %lf %lf\n",s,OldNrg-b->e[0],b->e[1]);
    fflush(NrgFile);
    b->step = s;
  }
  fclose(NrgFile);
}
/* The snapshot file should have the number of 
   chain a unit smaller than the number of chains written.
   I just swap the last chain with the current one, 
   calculate and save the energy. */
// one block, one chain architecture
// highly inefficient
static void WidomLoopRemSwap(struct beads *restrict b){
  struct nblist nbl;
  int NStep = b->local_n;
  int Nppc = b->N[0];
  FILE *NrgFile = fopen("WidomRawOut.dat","w");
  if (NrgFile == NULL) fatal(EIO, "Couldn't open WidomRawOut.dat");
  forces_calc(b);
  VEC2 *Lastx;
  Lastx = calloc(Nppc, sizeof(*Lastx));
  if(Lastx == NULL) return novm("Last chain");
  VEC *BfLast;
  BfLast = (VEC *)calloc(Nppc, sizeof(*BfLast));
  if(BfLast == NULL) return novm("Last chain");
  LoadLastChain("LastChain.dat",Lastx,BfLast,Nppc,b->l);
  for(int s=0;s<NStep;s++){
    printf("%d\n",s);
    //swap the next chain with the one in the gray zone and recalculate the interactions
    SwapLast(b,s,Lastx,BfLast);
    b->l[0] += .0001;
    forces_calc(b);
    fprintf(NrgFile,"%d %lf %lf\n",s,b->e[0],b->e[1]);
    fflush(NrgFile);
    //printf("%d) %lf %lf %lf\n",s,b->e[0],b->e[1],b->e[2]);
    b->step = s;
  }
  fclose(NrgFile);
}
/* The snapshot file should have an additional chain whom 
   positions are changed every timestep */
static void WidomLoopAdd(struct beads *restrict b){
//  struct nblist nbl;
  struct WIDOM *Widom = &(b->Widom);
  int NStep = 10000;
  FILE *NrgFile = fopen("WidomRawIn.dat","w");
  Widom->NBins = 100;
  Widom->Distr = (double *)calloc(Widom->NBins,sizeof(double));
  StudySystem(b);
  forces_calc(b);
  OneChainMore(b);
  ReallocBeads(b);
  int nChain = b->local_n-1;
  for(int s=0;s<NStep;s++){
    if(b->NNano == 1){
      InsertChainCorona(b,nChain);
    }
    else
      //InsertChainDiblock(b,nChain);
      InsertChain(b,nChain);
    double Nrg = CalcNrgCh(b,nChain);
    fprintf(NrgFile,"%d %lf %lf\n",s,b->e[0],b->e[1]);
    fflush(NrgFile);
    b->step = s;
  }
  fclose(NrgFile);
}
static void WidomNrg(struct beads *b){
//  struct nblist nbl;
  b->e[0] = 0.;b->e[1] = 0.;b->e[2] = 0.;
  forces_calc(b);
  MPI_Reduce(ismaster ? MPI_IN_PLACE : b->e, ismaster ? b->e : NULL,ARRAY_SIZE(b->e), MPI_DOUBLE, MPI_SUM, 0, comm_grid);
  printf("%d) %lf %lf %lf\n",b->step,b->e[0],b->e[1],b->e[2]);
}
static int WidomLoopMCMove(struct beads *b){
  struct nblist nbl;
  struct WIDOM *Widom = &(b->Widom);
  int NStep = 1000000;
  Widom->NBins = 100;
  Widom->Distr = (double *)calloc(Widom->NBins,sizeof(double));
  StudySystem(b);
  forces_calc(b);
  int TrialIn = 0;
  int TrialOut = 0;
  for(int s=0;s<NStep;s++){
    int p = rng_uniform(rng)*b->local_n*b->N[0];
    double Pos[3] = {b->xv[p][0],b->xv[p][1],b->xv[p][2]};
    double Nrg0 = mdpd_EnergySingleMon(b,p);
    for(int d=0;d<3;d++)
      b->xv[p][d] += rng_gaussian(rng,0.,Widom->Var);
    double Nrg1 = mdpd_EnergySingleMon(b,p);
    double NrgDiff = Nrg1 - Nrg0;
    if(!IfMetropolisMove(NrgDiff)){
      for(int d=0;d<3;d++) b->xv[p][d] = Pos[d];
    }
    if( !(s%b->local_n) ) {
      b->step=s;
      b->time = b->step * b->dt;
    }
  }
  save_beads("result.dat",b);
  printf("Ratio In: %lf Out:%lf\n",(NStep-TrialIn)/(double)NStep,(NStep-TrialOut)/(double)NStep);
  return 0;
}
static int WidomLoopMC(struct beads *b){
//  struct nblist nbl;
  struct WIDOM *Widom = &(b->Widom);
  int NStep = 1000000;
  Widom->NBins = 100;
  Widom->Distr = (double *)calloc(Widom->NBins,sizeof(double));
  StudySystem(b);
  forces_calc(b);
  int TrialIn = 0;
  int TrialOut = 0;
  /* FILE *FConf = fopen("GaussianCoil.dat","w"); */
  /* fprintf(FConf, "# L=%lg %lg %lg t=%.8f blocks=%d\n", */
  /* 	  b->l[0], b->l[1], b->l[2], b->time, b->blocks); */
  /* mdpd_print_header(FConf, b->mdpd); */
  /* fprintf(FConf,"# n=%d N=%d name=%s\n",1000,b->N[0],"Ciccia"); */
  for(int s=0;s<NStep;s++){
    double Ran = (rng_uniform(rng));
    if(Ran > .5){
      if(MCStepIn(b)) TrialIn++;
      //else PrintCh(b,FConf,b->local_n-1);
    }
    else{
      if(MCStepOut(b)) TrialOut++;
    }
    if( !(s%b->local_n) ) {
      b->step = s;
      b->time = b->step * b->dt;
    }
  }
  save_beads("result.dat",b);
  printf("Ratio In: %lf Out:%lf\n",(NStep-TrialIn)/(double)NStep,(NStep-TrialOut)/(double)NStep);
  return 0;
}
int WidomLoop(struct beads *b){
  struct WIDOM *Widom = &(b->Widom);
  struct mdpd *mdpd = b->mdpd;
  Widom->Var = sqrt(SQR(mdpd_getRe(mdpd))/(double)(b->N[0]-1)/3.)/2.;
  if(Widom->CalcMode == WIDOM_IN)
    WidomLoopAdd(b);
  else if(Widom->CalcMode == WIDOM_OUT)
    WidomLoopRem(b);
  //WidomLoopRemSwap(b);
  else if(Widom->CalcMode == WIDOM_NRG)
    WidomNrg(b);
  else if(Widom->CalcMode == WIDOM_NRG_CH){
    FILE *FWidom = fopen("WidomChOut.dat","w");
    for(int c=0;c<b->local_n;c++){
      fprintf(FWidom,"%lf\n",CalcNrgCh(b,c));
    }
  }
  else if(Widom->CalcMode == WIDOM_MC)
    WidomLoopMCMove(b);
  //WidomLoopMC(b);
  else 
    return 1;
  return 0;
}
// memory 
void WidomFree(struct beads *b){
//  struct WIDOM *Widom = &(b->Widom);
}
static void ReallocBeads(struct beads *restrict b){
  struct WIDOM *Widom = &(b->Widom);
  for (int i = 0 ; i < b->blocks; i++) {
    b->nN += b->n[i] * b->N[i];
    /* for (int j = 0; j < b->n[i]; j++) { /\* molecules *\/ */
    /*   //if (is_chain_local(chain_count, i)) { */
    /*   b->local_n_b[b->local_n++] = i; */
    /* } */
  }
  size_t count = (groups_per_row+groups_per_col)*b->groupsize;
  if(Widom->NAlloc >= count) return;
  else Widom->NAlloc = count + 1000*b->N[0];
  b->f = realloc(b->f,Widom->NAlloc* sizeof(*b->f));
  if (b->f == NULL) novm("beads->f");
  b->xv = realloc(b->xv,Widom->NAlloc* sizeof(*b->xv));
  if (b->xv == NULL) novm("beads->xv");
  b->nx = realloc(b->nx,Widom->NAlloc* sizeof(*b->nx));
  if (b->nx == NULL) novm("beads->nx");
  b->passport = realloc(b->passport,Widom->NAlloc* sizeof(*b->passport));
  if (b->passport == NULL) novm("beads->passport");
  b->v2 = realloc(b->v2,Widom->NAlloc* sizeof(*b->v2));
  if (b->v2 == NULL) novm("beads->v2");
  b->x_intra = realloc(b->x_intra,b->groupsize* sizeof(*b->x_intra));
  if (b->x_intra == NULL) novm("beads->x_intra");
  b->f_intra = realloc(b->f_intra,b->groupsize* sizeof(*b->f_intra));
  if (b->f_intra == NULL) novm("beads->f_intra");
  count = 0;
  for (int i = 0; i < b->blocks; i++)
    count += b->n[i];
  b->local_n_b = realloc(b->local_n_b,count* sizeof(*b->local_n_b));
  if (b->local_n_b == NULL) novm("beads->local_n_b");
}
