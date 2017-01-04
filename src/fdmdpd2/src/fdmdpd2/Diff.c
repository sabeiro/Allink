#include "fdmdpd2.h"

cfg_opt_t Diff_opts[] = {
  CFG_INT("NBin", 0, CFGF_NONE),
  CFG_INT("NStep", 0, CFGF_NONE),
  CFG_END()
};
/*------------DEFINING--------------------------*/
void DiffLoad(struct beads *b,cfg_t *cfg){
  struct DIFF_PROF *Diff = &(b->Diff);
  debug("Diff] Reading configurations");
  //----------------Number-of-components-------------
  cfg_t *cTens = cfg_getsec(cfg, "Diff");
  Diff->NBin = cfg_getint(cTens,"NBin");
  Diff->NStep = cfg_getint(cTens,"NStep");
  Diff->OldPos = (double *)calloc(3*b->nN,sizeof(double));
  Diff->VolSlabInv = (double *)calloc(Diff->NBin,sizeof(double));
  Diff->Count = (double *)calloc(Diff->NBin,sizeof(double));
  Diff->PartOr = (int *)calloc(b->nN,sizeof(int));
  Diff->BoxRad = .5*MIN(b->l[TANG1],b->l[TANG2]);
  Diff->InitStep = b->step;
  double Volume1 = 0.;
  double HeightCyl = b->l[NORMAL];
  double LenMin = MIN(b->l[TANG1],b->l[TANG2]);
  double LenMax = MAX(b->l[TANG1],b->l[TANG2]);
  for(int s=0;s<Diff->NBin;s++){
    double Rad2 = (s+1.)/Diff->NBin*Diff->BoxRad;
    double Volume2 = M_PI * SQR(Rad2);
    if(Rad2 > LenMin*.5){
      double Pro2 = sqrt(SQR(Rad2)-SQR(.5*LenMin));
      double Angle2 = atan(.5*LenMin/Pro2);
      Volume2 = 2.*SQR(Rad2)*Angle2 + Pro2*LenMin;
    }
    double Volume = (Volume2 - Volume1)*HeightCyl;
    Diff->VolSlabInv[s] = 1. / Volume;
    Volume1 = Volume2;
  }
  int first = col_index*b->groupsize;
  int last = (col_index+1)*b->groupsize;
  VEC dr;
  struct NANO *Nano = b->Nano;
  for(int p = first;p<last;p++){
    double dr2 = 0.;
    for(int d=0;d<3;d++){
      Diff->OldPos[3*p+d] = b->xv[p][d];
      dr[d] = Nano->Pos[d] - b->xv[p][d];
      dr[d] -= floor( dr[d]/(b->l[d]) + .5)*b->l[d];
    }
    dr2 = sqrt(SQR(dr[TANG1]) + SQR(dr[TANG2]));
    int vr = (int)(dr2/Diff->BoxRad*Diff->NBin);
    Diff->PartOr[p] = vr;
    if(vr < 0 || vr >= Diff->NBin) continue;
    Diff->Count[vr] += 1.;
  }
  char FName[60];
  double Norm = 1./(double)b->nN;
  for(int v=0;v<Diff->NBin;v++){
    sprintf(FName,"MeanSqDispl%03d.dat",v);
    FILE *FWrite = fopen(FName,"w");
    fclose(FWrite);
  }
}
/*------------DEFINING--------------------------*/
void DiffFree(struct beads *b){
  struct DIFF_PROF *Diff = &(b->Diff);
  debug("Diff] Freeing");
  //----------------Number-of-components-------------
  free(Diff->OldPos);
  free(Diff->Count);
  free(Diff->PartOr);
}
int DiffFillProf(struct beads *b){
  struct DIFF_PROF *Diff = &(b->Diff);
  int s = b->step - Diff->InitStep;
  printf("step %d)/%d\n",s,Diff->NStep);
  if(s >= Diff->NStep){
    return 1;
  }
  struct NANO *Nano = b->Nano;
  int first = col_index*b->groupsize;
  int last = (col_index+1)*b->groupsize;
  double *Prof = (double *)calloc(Diff->NBin,sizeof(double));
  for(int p = first;p<last;p++){
    VEC dr;
    double dr2 = 0.;
    for(int d=0;d<3;d++){
      dr[d] = Nano->Pos[d] - b->xv[p][d];
      dr[d] -= floor( dr[d]/b->l[d] + .5)*b->l[d];
    }
    dr2 = sqrt(SQR(dr[TANG1]) + SQR(dr[TANG2]));
    int vr = (int)(dr2/Diff->BoxRad*Diff->NBin);
    if(vr < 0 || vr >= Diff->NBin) continue;
    int vr1 = Diff->PartOr[p];
    if(vr1 < 0 || vr1 >= Diff->NBin) continue;
    //if(vr != vr1) {printf("%d is out %d != %d\n",p,vr,vr1);}
    double Rad2 = 0.;
    for(int d=0;d<2;d++){
      Rad2 += SQR(b->xv[p][d] - Diff->OldPos[3*p+d]);
    }
    Prof[vr1] += Rad2;
  }
  char FName[60];
  for(int v=0;v<Diff->NBin;v++){
    sprintf(FName,"MeanSqDispl%03d.dat",v);
    FILE *FWrite = fopen(FName,"a");
    //fprintf(FWrite,"%d %lf\n",s,Prof[v]*Diff->VolSlabInv[v]);
    fprintf(FWrite,"%d %lf\n",s,Prof[v]/Diff->Count[v]);
    fclose(FWrite);
  }
  free(Prof);
  return 0;
}
