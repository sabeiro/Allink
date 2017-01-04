/* $Id: Tens.c 366 2013-08-08 14:19:08Z fuhrmans $ */

#include "fdmdpd2.h"

static void TensAddComponent(struct TENS_PROF *Tens,int v1,int v2,double *Pre,int *Comp,int NPre);
static void TensAddComponentRad(struct TENS_PROF *Tens,int v1,int v2,double *Pre,int *Comp,int NPre);
static void TensLineVolContrib(struct beads *b);
static void TensLineDensContrib(struct beads *b);
static void TensTiltVolContrib(struct beads *b);
static void TensTiltDensContrib(struct beads *b);
static void Tens3dVolContrib(struct beads *b);
static void Tens3dDensContrib(struct beads *b);
static void Tens3dAddCompPos(struct TENS_PROF *Tens,double *Pos1,double *Pos2,double *Pre,double *Edge);
static void Tens2dVolContrib(struct beads *b);
static void Tens2dDensContrib(struct beads *b);
static void Tens2dAddCompPos(struct TENS_PROF *Tens,double *Pos1,double *Pos2,double *Pre,double *Edge);
static void TensNoAdd(struct beads *b,double Force,PASSPORT prow,PASSPORT pcol) {}
static void TensLineAdd(struct beads *b,double Force,PASSPORT prow,PASSPORT pcol);
static void TensTiltAdd(struct beads *b,double Force,PASSPORT prow,PASSPORT pcol);
static void Tens3dAdd(struct beads *b,double Force,PASSPORT prow,PASSPORT pcol);
static void Tens2dAdd(struct beads *b,double Force,PASSPORT prow,PASSPORT pcol);
static void TensNoAddPos(struct beads *b, double ExtForce, VEC Pos1, VEC Pos2,int t1,int t2) {}
static void TensNoAddPre(struct beads *b, double *PreExt,PASSPORT prow,PASSPORT pcol){}
static void TensLineAddPos(struct beads *b, double ExtForce, VEC Pos1, VEC Pos2,int t1,int t2);
static void TensLineAddPre(struct beads *b, double *PreExt,PASSPORT prow,PASSPORT pcol);
static void TensTiltAddPos(struct beads *b, double ExtForce, VEC Pos1, VEC Pos2,int t1,int t2);
static void TensTiltAddPre(struct beads *b, double *PreExt,PASSPORT prow,PASSPORT pcol);
static void Tens3dAddPre(struct beads *b, double *PreExt,PASSPORT prow,PASSPORT pcol);
static void Tens3dAddPos(struct beads *b, double ExtForce, VEC Pos1, VEC Pos2,int t1,int t2);
static void Tens2dAddPre(struct beads *b, double *PreExt,PASSPORT prow,PASSPORT pcol);
static void Tens2dAddPos(struct beads *b, double ExtForce, VEC Pos1, VEC Pos2,int t1,int t2);
static int TypeIdx(int a, int b) __attribute__((const));
static double TiltDistance(double *Pos,double *NPos,double *Edge,double *Axis,double *dr,double *PosRel);
void (*TensStress)(struct beads *b,double Force,PASSPORT prow,PASSPORT pcol) = TensNoAdd;
void (*TensStressPos)(struct beads *b, double ExtForce, VEC Pos1, VEC Pos2,int t1,int t2) = TensNoAddPos;
void (*TensStressPre)(struct beads *b,double *PreExt,PASSPORT prow,PASSPORT pco) = TensNoAddPre;
static void TensLineAlloc(struct TENS_PROF *Tens);
static void Tens3dAlloc(struct TENS_PROF *Tens);
static void Tens2dAlloc(struct TENS_PROF *Tens);
static void TensLineWrite(struct beads *b);
static void Tens3dWrite(struct beads *b);
static void Tens3dWriteAv(struct beads *b);
static void Tens2dWrite(struct beads *b);
static void Tens2dWriteAv(struct beads *b);
static void TensLineFree(struct beads *b);
static void Tens3dFree(struct beads *b);
static void Tens2dFree(struct beads *b);

cfg_opt_t Tens_opts[] = {
  CFG_INT("NSlab", 0, CFGF_NONE),
  CFG_INT("NSlab3d", 0, CFGF_NONE),
  CFG_INT("NComp", 0, CFGF_NONE),
  CFG_FLOAT("NAverage", 0, CFGF_NONE),
  CFG_STR("CalcMode", 0, CFGF_NONE),
  CFG_END()
};
int TensLoop(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  if(!Tens->IfLoop) return 1;
  TensEnable(b);
  forces_calc(b);
  Tens->Step = b->step;
  TensWrite(b);
  return 0;
}
/*------------DEFINING--------------------------*/
void TensLoad(struct beads *b,cfg_t *cfg){
  struct TENS_PROF *Tens = &(b->Tens);
  debug("Tens] Reading configurations");
  //----------------Number-of-components-------------
  cfg_t *cTens = cfg_getsec(cfg, "Tens");
  Tens->NSlab = cfg_getint(cTens,"NSlab");  
  Tens->NComp = cfg_getint(cTens,"NComp");
  Tens->NAverage =  60.;
  Tens->NAverage = cfg_getfloat(cTens,"NAverage");
  const char *CalcMode = cfg_getstr(cTens,"CalcMode");
  Tens->NComp = 5;
  Tens->FNorma = 0.;
  for(int d=0;d<3;d++)
    Tens->RefPos[d] = .5*b->l[d];
  if(b->NNano > 0){
    struct NANO *Nano = b->Nano;
    //if(VAR_IF_TYPE(Nano->Shape,SHAPE_TIP))continue;
    for(int d=0;d<3;d++)
      Tens->RefPos[d] = Nano->Pos[d];
  }    
  Tens->RefAxis[TANG1] = 0.;
  Tens->RefAxis[TANG2] = 0.;
  Tens->RefAxis[NORMAL] = 1.;
  Tens->BoxRad = .5 * MIN(b->l[TANG1],b->l[TANG2]);
  Tens->NType = 1;//(TYPE_MAX * (TYPE_MAX + 1)) / 2;
  Tens->NLine = Tens->NComp*Tens->NType;
  Tens->IfLoop = 0;
  if(!strcmp(CalcMode,"line")){
    Tens->CalcMode = TENS_LINE;
    TensLineAlloc(Tens);
  }
  else if(!strcmp(CalcMode,"tilt")){
    Tens->CalcMode = TENS_TILT;
    TensLineAlloc(Tens);
  }
  else if(!strcmp(CalcMode,"no")){
    Tens->CalcMode = TENS_NO;
  }
  else if(!strcmp(CalcMode,"3d")){
    Tens->CalcMode = TENS_3D;
    Tens->NSlab = cfg_getint(cTens,"NSlab3d");
    Tens->NType = 1;
    Tens->NComp = 6;
    Tens->NLine = Tens->NComp*Tens->NType;
    Tens->Wrap[0] = 1;
    Tens->Wrap[1] = 1;
    Tens->Wrap[2] = 1;
    Tens3dAlloc(Tens);
  }
  else if(!strcmp(CalcMode,"2d")){
    Tens->CalcMode = TENS_2D;
    Tens->NSlab = cfg_getint(cTens,"NSlab3d");
    Tens->NType = 1;
    Tens->NComp = 6;
    Tens->NLine = Tens->NComp*Tens->NType;
    Tens->Wrap[0] = 0;
    Tens->Wrap[1] = 1;
    Tens->Wrap[2] = 1;
    Tens2dAlloc(Tens);
  }
  else if(!strcmp(CalcMode,"loop2d")){
    Tens->CalcMode = TENS_2D;
    Tens->NSlab = cfg_getint(cTens,"NSlab3d");
    Tens->NType = 1;
    Tens->NComp = 6;
    Tens->IfLoop = 1;
    Tens->NLine = Tens->NComp*Tens->NType;
    Tens2dAlloc(Tens);
  }
  else 
    fatal(EINVAL, "Calculation mode not valid %s |= Linear,Tilt,3d",CalcMode);
  debug("Tens] Allocating NSlab %d NComp %d NType %d",Tens->NSlab,Tens->NComp,Tens->NType);
}
//------------------Alloc--------------------------
void TensLineAlloc(struct TENS_PROF *Tens){
  Tens->Tension = calloc(Tens->NLine,sizeof(*Tens->Tension));
  if (Tens->Tension == NULL) novm("Tens->Tension");
  Tens->Count = calloc(Tens->NComp*TYPE_MAX,sizeof(*Tens->Count));
  if (Tens->Count == NULL) novm("Tens->Count");
  for(int l=0;l<Tens->NLine;l++){
    Tens->Tension[l] = calloc(Tens->NSlab,sizeof(**Tens->Tension));
    if (Tens->Tension[l] == NULL) novm("Tens->Tension");
  }
  for(int l=0;l<Tens->NComp*TYPE_MAX;l++){
    Tens->Count[l] = calloc(Tens->NSlab,sizeof(**Tens->Count));
    if (Tens->Count[l] == NULL) novm("Tens->Count");
  }
  Tens->VolSlabInv = calloc(Tens->NSlab,sizeof(*Tens->VolSlabInv));
  if (Tens->VolSlabInv == NULL) novm("Tens->VolSlabInv");
}
void Tens3dAlloc(struct TENS_PROF *Tens){
  Tens->Tension = calloc(Tens->NLine,sizeof(*Tens->Tension));
  if (Tens->Tension == NULL) novm("Tens->Tension");
  for(int l=0;l<Tens->NLine;l++){
    Tens->Tension[l] = calloc(CUBE(Tens->NSlab),sizeof(**Tens->Tension));
    if (Tens->Tension[l] == NULL) novm("Tens->Tension");
  }
  Tens->Count = calloc(TYPE_MAX,sizeof(*Tens->Count));
  if (Tens->Count == NULL) novm("Tens->Count");
  for(int l=0;l<TYPE_MAX;l++){
    Tens->Count[l] = calloc(CUBE(Tens->NSlab),sizeof(**Tens->Count));
    if (Tens->Count[l] == NULL) novm("Tens->Count");
  }
  Tens->VolSlabInv = calloc(Tens->NSlab,sizeof(*Tens->VolSlabInv));
  if (Tens->VolSlabInv == NULL) novm("Tens->VolSlabInv");
}
void Tens2dAlloc(struct TENS_PROF *Tens){
  Tens->Tension = calloc(Tens->NLine,sizeof(*Tens->Tension));
  if (Tens->Tension == NULL) novm("Tens->Tension");
  for(int l=0;l<Tens->NLine;l++){
    Tens->Tension[l] = calloc(SQR(Tens->NSlab),sizeof(**Tens->Tension));
    if (Tens->Tension[l] == NULL) novm("Tens->Tension");
  }
  Tens->Count = calloc(TYPE_MAX,sizeof(*Tens->Count));
  if (Tens->Count == NULL) novm("Tens->Count");
  for(int l=0;l<TYPE_MAX;l++){
    Tens->Count[l] = calloc(SQR(Tens->NSlab),sizeof(**Tens->Count));
    if (Tens->Count[l] == NULL) novm("Tens->Count");
  }
  Tens->VolSlabInv = calloc(Tens->NSlab,sizeof(*Tens->VolSlabInv));
  if (Tens->VolSlabInv == NULL) novm("Tens->VolSlabInv");
}
/* Cleanup */
void TensFree(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  if(Tens->CalcMode==TENS_NO) return;
  /////////////////////////////////////////////
  for(int l=0;l<Tens->NLine;l++)
    free(Tens->Tension[l]);
  for(int l=0;l<Tens->NComp*TYPE_MAX;l++)
    free(Tens->Count[l]);
  free(Tens->Tension);
  free(Tens->Count);
  free(Tens->VolSlabInv);
}
/** Disable calculation */
void TensDisable(){
  TensStress = TensNoAdd;
  TensStressPos = TensNoAddPos;
  TensStressPre = TensNoAddPre;
}
/** Enable calculation temporarily */
void TensEnable(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  if(Tens->CalcMode == TENS_LINE){
    TensStress = TensLineAdd;
    TensStressPos = TensLineAddPos;
    TensStressPre = TensLineAddPre;
    TensLineVolContrib(b);
    TensLineDensContrib(b);
  }
  else if(Tens->CalcMode == TENS_TILT){
    TensStress = TensTiltAdd;
    TensStressPos = TensTiltAddPos;
    TensStressPre = TensTiltAddPre;
    TensTiltVolContrib(b);
    TensTiltDensContrib(b);
  }
  else if(Tens->CalcMode == TENS_3D){
    TensStress = Tens3dAdd;
    TensStressPos = Tens3dAddPos;
    TensStressPre = Tens3dAddPre;
    Tens3dVolContrib(b);
    Tens3dDensContrib(b);
  }
  else if(Tens->CalcMode == TENS_2D){
    TensStress = Tens2dAdd;
    TensStressPos = Tens2dAddPos;
    TensStressPre = Tens2dAddPre;
    Tens2dVolContrib(b);
    Tens2dDensContrib(b);
  }
  else if(Tens->CalcMode == TENS_NO){
    TensDisable();
    return;
  }
  Tens->FNorma += 1.;
}
//############################VOLUME########################
/** Calculate the components of a vector projected on an axis , Axis should be already normalized */
static double ProjectOn(double *Axis,double *Pos){
  double Scalar = 0.;
  //double NormA = 0.;
  for(int d=0;d<3;d++){
    Scalar += Axis[d]*Pos[d];
    //NormA += Axis[d];
  }
  return Scalar; // /NormA
}
/** Matrix vector multiplication */
static void MatrVect(double *M,double *v,double *Resp){
  for(int c=0;c<3;c++){
    double Temp = 0.;
    for(int r=0;r<3;r++)
      Temp += M[3*c + r] * v[r];
    Resp[c] = Temp;
  }
}
/** Calculates the radial and normal volume for the s Slab*/
static void TensLineVolContrib(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  //debug("Tens] Volume Contribution");
  double Volume1 = 0.;
  double HeightCyl = b->l[NORMAL];
  double LenMin = MIN(b->l[TANG1], b->l[TANG2]);
  double LenMax = MAX(b->l[TANG1], b->l[TANG2]);
  for(int s=0;s<Tens->NSlab;s++){
    double Rad2 = (s+1.)/Tens->NSlab*Tens->BoxRad;
    double Volume2 = M_PI * SQR(Rad2);
    if(Rad2 > LenMin*.5){
      double Pro2 = sqrt(SQR(Rad2)-SQR(.5*LenMin));
      double Angle2 = atan(.5*LenMin/Pro2);
      Volume2 = 2.*SQR(Rad2)*Angle2 + Pro2*LenMin;
    }
    double Volume = (Volume2 - Volume1)*HeightCyl;
    Tens->VolSlabInv[s] = 1. / Volume;
    Volume1 = Volume2;
  }
  Tens->VolSlabInvPerp = 1/(b->l[0]*b->l[1]*b->l[2]);
}
static void TensTiltVolContrib(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  //debug("Tens] Volume Contribution");
  double Volume1 = 0.;
  double HeightCyl = b->l[NORMAL];
  double LenMin = .5*MIN(b->l[TANG1], b->l[TANG2]);
  double LenMax = .5*MAX(b->l[TANG1], b->l[TANG2]);
  for(int s=0;s<Tens->NSlab;s++){
    double Rad2 = (s+1.)/Tens->NSlab*Tens->BoxRad;
    double Volume2 = M_PI * SQR(Rad2);
    if(Rad2 > LenMin*.5){
      double Pro2 = sqrt(SQR(Rad2)-SQR(LenMin));
      double Angle2 = atan(LenMin/Pro2);
      Volume2 = 2.*SQR(Rad2)*Angle2 + 2.*Pro2*LenMin;
    }
    double Volume = (Volume2 - Volume1)*HeightCyl;
    Tens->VolSlabInv[s] = Volume > 0 ? 1. / Volume : 1.;
    Volume1 = Volume2;
  }
  Tens->VolSlabInvPerp = 1/(b->l[0]*b->l[1]*b->l[2]);
  double Zed[3] = {0.,0.,1.};
  double RotAxis[3] = {0.,0.,0.};
  // Projected on the moving axis
  double NanoAxis[3] = {Tens->RefAxis[0],Tens->RefAxis[1],Tens->RefAxis[2]};
  double RotAngle = AngleVect(NanoAxis,Zed);
  VectProd(NanoAxis,Zed,RotAxis);
  RotMatr9(Tens->Rot,RotAxis,RotAngle);
}
/** Calculates the radial and normal volume for the s Slab*/
static void Tens3dVolContrib(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  //debug("Tens] Volume Contribution");
  double Volume1 = 0.;
  double HeightCyl = b->l[NORMAL];
  double LenMin = MIN(b->l[TANG1], b->l[TANG2]);
  double LenMax = MAX(b->l[TANG1], b->l[TANG2]);
  for(int s=0;s<Tens->NSlab;s++){
    double Rad2 = (s+1.)/Tens->NSlab*LenMax*.5;
    double Volume2 = M_PI * SQR(Rad2);
    if(Rad2 > LenMin*.5){
      double Pro2 = sqrt(SQR(Rad2)-SQR(.5*LenMin));
      double Angle2 = atan(.5*LenMin/Pro2);
      Volume2 = 2.*SQR(Rad2)*Angle2 + Pro2*LenMin;
    }
    double Volume = (Volume2 - Volume1)*HeightCyl;
    Tens->VolSlabInv[s] = 1. / Volume;
    Volume1 = Volume2;
  }
  Tens->VolSlabInvPerp = 1/(b->l[0]*b->l[1]*b->l[2]);
  double Zed[3] = {0.,0.,1.};
  double RotAxis[3] = {0.,0.,0.};
  // Projected on the moving axis
  //double NanoAxis[3] = {Nano->Axis[0],Nano->Axis[1],Nano->Axis[2]};
  // Projected on the fixed normal
  double NanoAxis[3] = {0.,0.,1.};
  double RotAngle = AngleVect(NanoAxis,Zed);
  VectProd(NanoAxis,Zed,RotAxis);
  RotMatr9(Tens->Rot,RotAxis,RotAngle);
}
static void Tens2dVolContrib(struct beads *b){
  TensLineVolContrib(b);
}
//######################DENSITY################################
/** Sums over all particles and adds a contribution of density in the s Slab */
static void TensLineDensContrib(struct beads *b){
  struct TENS_PROF *Tens = &b->Tens;
  //debug("Tens] Density Contribution");
  double LengthInv = Tens->NSlab/b->l[NORMAL];
  double BoxRadInv = Tens->NSlab/Tens->BoxRad;
  double VolInv = Tens->NSlab/(b->l[0]*b->l[1]*b->l[2]);
  int pFirst = col_index * b->groupsize; 
  int pLast = (col_index + 1) * b->groupsize;
  // calculate the density
  for(int p=pFirst;p<pLast;p++){
    int type = GET_TYPE(b->passport[p]);
    double r1x = b->xv[p][TANG1] - Tens->RefPos[TANG1];
    r1x -= floor( r1x/(b->l[TANG1]) + .5)*b->l[TANG1];
    double r1y = b->xv[p][TANG2] - Tens->RefPos[TANG2];
    r1y -= floor( r1y/(b->l[TANG2]) + .5)*b->l[TANG2];
    double Dist = hypot(r1x, r1y);
    double Height = b->xv[p][NORMAL];
    int sn = (int) ( Height * LengthInv);
    //assert(sn >= 0 && sn < Tens->NSlab);
    Tens->Count[type+0*TYPE_MAX][sn] += VolInv;
    int sr = (int)(Dist*BoxRadInv);
    //assert(sr >= 0 && sr < Tens->NSlab);
    if(sr < 0 || sr >= Tens->NSlab) continue;
    Tens->Count[type+1*TYPE_MAX][sr] += Tens->VolSlabInv[sr];
  }
}
/** Sums over all particles and adds a contribution of density in the s Slab */
static void TensTiltDensContrib(struct beads *b){
  struct TENS_PROF *Tens = &b->Tens;
  //debug("Tens] Density Contribution");
  double LengthInv = Tens->NSlab/b->l[NORMAL];
  Tens->BoxRad = .5 * sqrt( SQR(b->l[0]) + SQR(b->l[1]) + SQR(b->l[2]) );
  double BoxRadInv = Tens->NSlab/Tens->BoxRad;
  double VolInv = Tens->NSlab/(b->l[0]*b->l[1]*b->l[2]);
  int pFirst = col_index * b->groupsize; 
  int pLast = (col_index + 1) * b->groupsize;
  // calculate the density
  for(int p=pFirst;p<pLast;p++){
    int type = GET_TYPE(b->passport[p]);
    double dr[3];
    double PosRel[3] = {b->xv[p][0],b->xv[p][1],b->xv[p][2]};
    double Rotated[3] = {b->xv[p][0],b->xv[p][1],b->xv[p][2]};
    double Dist = TiltDistance(b->xv[p],Tens->RefPos,b->l,Tens->RefAxis,dr,PosRel);
    ProjectOn(Tens->RefAxis,PosRel);
    MatrVect(Tens->Rot,PosRel,Rotated);
    double Height = Rotated[2] + Tens->RefPos[2];
    int sn = (int) ( Height * LengthInv);
    //assert(sn >= 0 && sn < Tens->NSlab);
    Tens->Count[type+0*TYPE_MAX][sn] += VolInv;
    int sr = (int)(Dist*BoxRadInv);
    if(sr < 0 || sr >= Tens->NSlab) continue;
    Tens->Count[type+1*TYPE_MAX][sr] += Tens->VolSlabInv[sr];
  }
}
static void Tens3dDensContrib(struct beads *b){
  struct TENS_PROF *Tens = &b->Tens;
  double Cm[3];
  double PosP1[3];
  int pFirst = col_index * b->groupsize; 
  int pLast = (col_index + 1) * b->groupsize;
  // calculate the density
  for(int p=pFirst;p<pLast;p++){
    int type = GET_TYPE(b->passport[p]);
    for(int d=0;d<3;d++){
      PosP1[d] = b->xv[p][d] - (b->l[d]*.5 - Tens->RefPos[d]);
      PosP1[d] -= floor(PosP1[d]/b->l[d])*b->l[d];
    }
    int vx = (int)(PosP1[TANG1]/b->l[TANG1]*Tens->NSlab);
    int vy = (int)(PosP1[TANG2]/b->l[TANG2]*Tens->NSlab);
    int vz = (int)(PosP1[NORMAL]/b->l[NORMAL]*Tens->NSlab);
    assert(vx >= 0 && vx < Tens->NSlab);
    assert(vy >= 0 && vy < Tens->NSlab);
    assert(vz >= 0 && vz < Tens->NSlab);
    /* if(vx < 0 || vx >= Tens->NSlab) continue; */
    /* if(vy < 0 || vy >= Tens->NSlab) continue; */
    /* if(vz < 0 || vz >= Tens->NSlab) continue; */
    int v = (vx*Tens->NSlab+vy)*Tens->NSlab+vz;
    Tens->Count[type][v] += 1.;
  }
}
static void Tens2dDensContrib(struct beads *b){
  struct TENS_PROF *Tens = &b->Tens;
  //debug("Tens] Density Contribution");
  double LengthInv = Tens->NSlab/b->l[NORMAL];
  double BoxRadInv = Tens->NSlab/Tens->BoxRad;
  double VolInv = Tens->NSlab/(Tens->BoxRad*b->l[NORMAL]);
  int pFirst = col_index * b->groupsize; 
  int pLast = (col_index + 1) * b->groupsize;
  // calculate the density
  for(int p=pFirst;p<pLast;p++){
    int type = GET_TYPE(b->passport[p]);
    double r1x = b->xv[p][TANG1] - Tens->RefPos[TANG1];
    r1x -= floor( r1x/(b->l[TANG1]) + .5)*b->l[TANG1];
    double r1y = b->xv[p][TANG2] - Tens->RefPos[TANG2];
    r1y -= floor( r1y/(b->l[TANG2]) + .5)*b->l[TANG2];
    double Dist = hypot(r1x, r1y);
    double Height = (b->xv[p][NORMAL] - Tens->RefPos[NORMAL]) + b->l[NORMAL]*.5 ;
    Height -= floor(Height/b->l[NORMAL])*b->l[NORMAL];
    int vr = (int)(Dist/Tens->BoxRad*Tens->NSlab);
    int vz = (int)(Height/b->l[NORMAL]*Tens->NSlab);
    if(vr < 0 || vr >= Tens->NSlab) continue;
    assert(vz >= 0 && vz < Tens->NSlab);
    int v = vr*Tens->NSlab+vz;
    Tens->Count[type][v] += VolInv*Tens->VolSlabInv[vr];
  }
}
//#################PRESSURE#############################

/** Adds the contribution of a line tension in every slab it crosses with respect to the normal and radial distance.
 *  See: Frenkel/Smit chapter 17, page 472 */
static int TypeIdx(int a, int b){
  return 0;
  assert(a >= 0 && a < TYPE_MAX);
  assert(b >= 0 && b < TYPE_MAX);
  return a > b ?
    b * TYPE_MAX - ((b > 1) ? (b * (b - 1) / 2) : 0) + a - b:
    a * TYPE_MAX - ((a > 1) ? (a * (a - 1) / 2) : 0) + b - a;
}
static void TensLineAdd(struct beads *b, double ExtForce,PASSPORT prow,PASSPORT pcol){
  int Id1 = GET_POS(prow);
  int Id2 = GET_POS(pcol);
  int t1  = GET_TYPE(prow);
  int t2  = GET_TYPE(pcol);
  TensLineAddPos(b,ExtForce,b->xv[Id1],b->xv[Id2],t1,t2);
}
static void TensTiltAdd(struct beads *b, double ExtForce,PASSPORT prow,PASSPORT pcol){
  int Id1 = GET_POS(prow);
  int Id2 = GET_POS(pcol);
  int t1  = GET_TYPE(prow);
  int t2  = GET_TYPE(pcol);
  TensTiltAddPos(b,ExtForce,b->xv[Id1],b->xv[Id2],t1,t2);
}
void TensLineAddPos(struct beads *b, double ExtForce,double *Pos1,double *Pos2,int t1,int t2)
{
  struct TENS_PROF *Tens = &(b->Tens);
  double LengthInv = Tens->NSlab/b->l[NORMAL];
  double BoxRadInv = Tens->NSlab/Tens->BoxRad;
  int TypeSum = TypeIdx(t1,t2);
  double PosRel1[3];
  double PosRel2[3];
  double DistRel[4] = {0.,0.,0.,0.};
  for(int d=0;d<3;d++){
    DistRel[d] = Pos2[d] - Pos1[d];
    DistRel[d] -= floor( DistRel[d]/(b->l[d]) + .5)*b->l[d];
    PosRel1[d] = Pos1[d] - Tens->RefPos[d];
    PosRel1[d] -= floor( PosRel1[d]/(b->l[d]) + .5)*b->l[d];
    PosRel2[d] = Pos2[d] - Tens->RefPos[d];
    PosRel2[d] -= floor( PosRel2[d]/(b->l[d]) + .5)*b->l[d];
    DistRel[3] += SQR(DistRel[d]);
  }
  double x1d = sqrt( SQR(PosRel1[TANG1]) + SQR(PosRel1[TANG2] ));
  double x2d = sqrt( SQR(PosRel2[TANG1]) + SQR(PosRel2[TANG2]) );
  DistRel[3] = sqrt(DistRel[3]);
  double Pre[3];
  int Comp[3];
  int vn1 = (int) ( (Pos1[NORMAL]) * LengthInv);
  int vn2 = (int) ( (Pos2[NORMAL]) * LengthInv);
  int vr1 = (int) ( (x1d) * BoxRadInv);
  int vr2 = (int) ( (x2d) * BoxRadInv);
  //-----------------(n n)--------------------------
  Pre[0] = ExtForce*SQR(DistRel[NORMAL])/DistRel[3];
  Comp[0] = TypeSum + 0 * Tens->NType;
  //-----------------(n r)--------------------------
  Pre[1] = .5*ExtForce*( SQR(DistRel[TANG1]) + SQR(DistRel[TANG2]) )/DistRel[3];
  Comp[1] = TypeSum + 1 * Tens->NType;
  if(vn1*vn2 > 0 && vn1 <= Tens->NSlab && vn2 <= Tens->NSlab)
    TensAddComponent(Tens,vn1,vn2,Pre,Comp,2);
  //-----------------(r l1)-------------------------
  Pre[0] = ExtForce * SQR(DistRel[TANG1]) / DistRel[3];
  Comp[0] = TypeSum + 2 * Tens->NType;
  //-----------------(r l2)-------------------------
  Pre[1] = ExtForce * SQR(DistRel[TANG2]) / DistRel[3];
  Comp[1] = TypeSum + 3 * Tens->NType;
  //-----------------(r n)--------------------------
  Pre[2] = ExtForce * SQR(DistRel[NORMAL]) / DistRel[3];
  Comp[2] = TypeSum + 4 * Tens->NType;
  if(vr1*vr2 > 0 && vr1 <= Tens->NSlab && vr2 <= Tens->NSlab)
    TensAddComponentRad(Tens,vr1,vr2,Pre,Comp,3);
}
void TensLineAddPre(struct beads *b, double *PreExt,PASSPORT prow,PASSPORT pcol){
  struct TENS_PROF *Tens = &(b->Tens);
  int Id1 = GET_POS(prow);
  int Id2 = GET_POS(pcol);
  int t1  = GET_TYPE(prow);
  int t2  = GET_TYPE(pcol);
  double *Pos1 = b->xv[Id1];
  double *Pos2 = b->xv[Id2];
  double LengthInv = Tens->NSlab/b->l[NORMAL];
  double BoxRadInv = Tens->NSlab/Tens->BoxRad;
  int TypeSum = TypeIdx(t1,t2);
  double PosRel1[3];
  double PosRel2[3];
  for(int d=0;d<3;d++){
    PosRel1[d] = Pos1[d] - Tens->RefPos[d];
    PosRel1[d] -= floor( PosRel1[d]/(b->l[d]) + .5)*b->l[d];
    PosRel2[d] = Pos2[d] - Tens->RefPos[d];
    PosRel2[d] -= floor( PosRel2[d]/(b->l[d]) + .5)*b->l[d];
  }
  double x1d = sqrt( SQR(PosRel1[TANG1]) + SQR(PosRel1[TANG2] ));
  double x2d = sqrt( SQR(PosRel2[TANG1]) + SQR(PosRel2[TANG2]) );
  double Pre[3];
  int Comp[3];
  int vn1 = (int) ( (Pos1[NORMAL]) * LengthInv);
  int vn2 = (int) ( (Pos2[NORMAL]) * LengthInv);
  int vr1 = (int) ( (x1d) * BoxRadInv);
  int vr2 = (int) ( (x2d) * BoxRadInv);
  //-----------------(n n)--------------------------
  Pre[0] = PreExt[NORMAL];
  Comp[0] = TypeSum + 0 * Tens->NType;
  //-----------------(n r)--------------------------
  Pre[1] = .5*(PreExt[TANG1]+PreExt[TANG2]);
  Comp[1] = TypeSum + 1 * Tens->NType;
  if(vn1*vn2 > 0 && vn1 <= Tens->NSlab && vn2 <= Tens->NSlab)
    TensAddComponent(Tens,vn1,vn2,Pre,Comp,2);
  //-----------------(r l1)-------------------------
  Pre[0] = PreExt[TANG1];
  Comp[0] = TypeSum + 2 * Tens->NType;
  //-----------------(r l2)-------------------------
  Pre[1] = PreExt[TANG2];
  Comp[1] = TypeSum + 3 * Tens->NType;
  //-----------------(r n)--------------------------
  Pre[2] = PreExt[NORMAL];
  Comp[2] = TypeSum + 4 * Tens->NType;
  if(vr1*vr2 > 0 && vr1 <= Tens->NSlab && vr2 <= Tens->NSlab)
    TensAddComponentRad(Tens,vr1,vr2,Pre,Comp,3);
}
/* The radial and normal distances are expressed in the reference system
   of the inclusion, translated wrt RefPos and rotated wrt RefAxis */
static double TiltDistance(double *Pos,double *NPos,double *Edge,double *Axis,double *dr,double *PosRel){
  //Calculate the joining PosRel Axis
  for(int d=0;d<3;d++){
    PosRel[d] = NPos[d] - Pos[d];
    PosRel[d] -= floor( PosRel[d]/(Edge[d]) + .5)*Edge[d];
  }
  double Dist2 = Dist2AxisPerp(PosRel,Axis,dr);
  //Project PosRel on Axis
  double Scalar = 0.;
  double NormP = 0.;
  for(int d=0;d<3;d++){
    Scalar += Axis[d]*Pos[d];
    NormP += SQR(Pos[d]);
  }
  double Cos = Scalar/sqrt(NormP);
  for(int d=0;d<3;d++){
    PosRel[d] = Cos*Axis[d];
  }
  return sqrt(Dist2);
}
void TensTiltAddPos(struct beads *b, double ExtForce,double *Pos1,double *Pos2,int t1,int t2)
{
  struct TENS_PROF *Tens = &(b->Tens);
  double LengthInv = Tens->NSlab/b->l[NORMAL];
  double BoxRadInv = Tens->NSlab/Tens->BoxRad;
  int TypeSum = TypeIdx(t1,t2);
  double Dist1[3], PosRel1[3];
  double Dist2[3], PosRel2[3];
  double DistRel[4] = {0.,0.,0.,0.};
  for(int d=0;d<3;d++){
    DistRel[d] = Pos2[d] - Pos1[d];
    DistRel[d] -= floor( DistRel[d]/(b->l[d]) + .5)*b->l[d];
    DistRel[3] += SQR(DistRel[d]);
  }
  double x1d = TiltDistance(Pos1,Tens->RefPos,b->l,Tens->RefAxis,Dist1,PosRel1);
  double x2d = TiltDistance(Pos2,Tens->RefPos,b->l,Tens->RefAxis,Dist2,PosRel2);
  DistRel[3] = sqrt(DistRel[3]);
  double Pre[3];
  int Comp[3];
  int vn1 = (int) ( (Pos1[NORMAL]) * LengthInv);
  int vn2 = (int) ( (Pos2[NORMAL]) * LengthInv);
  int vr1 = (int) ( (x1d) * BoxRadInv);
  int vr2 = (int) ( (x2d) * BoxRadInv);
  //-----------------(n n)--------------------------
  Pre[0] = ExtForce*SQR(DistRel[NORMAL])/DistRel[3];
  Comp[0] = TypeSum + 0 * Tens->NType;
  //-----------------(n r)--------------------------
  Pre[1] = .5*ExtForce*( SQR(DistRel[TANG1]) + SQR(DistRel[TANG2]) )/DistRel[3];
  Comp[1] = TypeSum + 1 * Tens->NType;
  if(vn1*vr2 > 0 && vn1 <= Tens->NSlab && vn2 <= Tens->NSlab)
    TensAddComponent(Tens,vn1,vn2,Pre,Comp,2);
  //-----------------(r l1)-------------------------
  Pre[0] = ExtForce * SQR(DistRel[TANG1]) / DistRel[3];
  Comp[0] = TypeSum + 2 * Tens->NType;
  //-----------------(r l2)-------------------------
  Pre[1] = ExtForce * SQR(DistRel[TANG2]) / DistRel[3];
  Comp[1] = TypeSum + 3 * Tens->NType;
  //-----------------(r n)--------------------------
  Pre[2] = ExtForce * SQR(DistRel[NORMAL]) / DistRel[3];
  Comp[2] = TypeSum + 4 * Tens->NType;
  if(vr1*vr2 > 0 && vr1 <= Tens->NSlab && vr2 <= Tens->NSlab)
    TensAddComponentRad(Tens,vr1,vr2,Pre,Comp,3);
}
void TensTiltAddPre(struct beads *b, double *PreExt,PASSPORT prow,PASSPORT pcol){
  struct TENS_PROF *Tens = &(b->Tens);
  int Id1 = GET_POS(prow);
  int Id2 = GET_POS(pcol);
  int t1  = GET_TYPE(prow);
  int t2  = GET_TYPE(pcol);
  double *Pos1 = b->xv[Id1];
  double *Pos2 = b->xv[Id2];
  double LengthInv = Tens->NSlab/b->l[NORMAL];
  double BoxRadInv = Tens->NSlab/Tens->BoxRad;
  int TypeSum = TypeIdx(t1,t2);
  double Dist1[3], PosRel1[3];
  double Dist2[3], PosRel2[3];
  for(int d=0;d<3;d++){
    PosRel1[d] = Pos1[d] - Tens->RefPos[d];
    PosRel1[d] -= floor( PosRel1[d]/(b->l[d]) + .5)*b->l[d];
    PosRel2[d] = Pos2[d] - Tens->RefPos[d];
    PosRel2[d] -= floor( PosRel2[d]/(b->l[d]) + .5)*b->l[d];
  }
  double x1d = TiltDistance(Pos1,Tens->RefPos,b->l,Tens->RefAxis,Dist1,PosRel1);
  double x2d = TiltDistance(Pos2,Tens->RefPos,b->l,Tens->RefAxis,Dist2,PosRel2);
  double Pre[3];
  int Comp[3];
  int vn1 = (int) ( (Pos1[NORMAL]) * LengthInv);
  int vn2 = (int) ( (Pos2[NORMAL]) * LengthInv);
  int vr1 = (int) ( (x1d) * BoxRadInv);
  int vr2 = (int) ( (x2d) * BoxRadInv);
  //-----------------(n n)--------------------------
  Pre[0] = PreExt[NORMAL];
  Comp[0] = TypeSum + 0 * Tens->NType;
  //-----------------(n r)--------------------------
  // conventional, there is no definition.
  Pre[1] = .5*(PreExt[TANG1]+PreExt[TANG2]);
  Comp[1] = TypeSum + 1 * Tens->NType;
  if(vn1*vn2 > 0 && vn1 <= Tens->NSlab && vn2 <= Tens->NSlab)
    TensAddComponent(Tens,vn1,vn2,Pre,Comp,2);
  //-----------------(r l1)-------------------------
  Pre[0] = PreExt[TANG1];
  Comp[0] = TypeSum + 2 * Tens->NType;
  //-----------------(r l2)-------------------------
  Pre[1] = PreExt[TANG2];
  Comp[1] = TypeSum + 3 * Tens->NType;
  //-----------------(r n)--------------------------
  Pre[2] = PreExt[NORMAL];
  Comp[2] = TypeSum + 4 * Tens->NType;
  if(vr1*vr2 > 0 && vr1 <= Tens->NSlab && vr2 <= Tens->NSlab)
    TensAddComponentRad(Tens,vr1,vr2,Pre,Comp,3);
}
static void TensAddComponent(struct TENS_PROF *Tens,int v1,int v2,double *Pre,int *Comp,int NPre){
  int NSlab = Tens->NSlab;
  int NSlabHalf = (int)(Tens->NSlab*.5);
  if(v1 == v2){
    for(int p=0;p<NPre;p++)
      Tens->Tension[Comp[p]][v1] += Pre[p];
    return;
  }
  int vMax = v1 > v2 ? v1 : v2;
  int vMin = v1 < v2 ? v1 : v2;
  if( vMax - vMin > NSlabHalf)
    {
    double NSlabInv = 1./(double)(NSlab - vMax + vMin + 1);
    for(int p=0;p<NPre;p++){
      for(int v = 0;v <= vMin;v++)
    	Tens->Tension[Comp[p]][v]  += Pre[p]*NSlabInv;
      for(int v = vMax;v < NSlab;v++)
      	Tens->Tension[Comp[p]][v]  += Pre[p]*NSlabInv;
    }
  }
  else 
    {
    double NSlabInv = 1./(double)(vMax + 1 - vMin);
    for(int p=0;p<NPre;p++)
      for(int v = vMin;v <= vMax;v++)
	Tens->Tension[Comp[p]][v]  += Pre[p]*NSlabInv;
  }
}
static void TensAddComponentRad(struct TENS_PROF *Tens,int v1,int v2,double *Pre,int *Comp,int NPre){
  int NSlab = Tens->NSlab;
  int NSlabHalf = (int)(Tens->NSlab*.5);
  if(v1 == v2){
    for(int p=0;p<NPre;p++)
      Tens->Tension[Comp[p]][v1] += Pre[p];
    return;
  }
  int vMax = MAX(v1, v2);
  int vMin = MIN(v1, v2);
  if( vMax - vMin > NSlabHalf){
    double NSlabInv = 1./(double)(NSlab - vMax + vMin + 1);
    for(int p=0;p<NPre;p++) {
      for(int v = vMax;v < NSlab;v++)
        Tens->Tension[Comp[p]][v]  += Pre[p]*Tens->VolSlabInv[v]*NSlabInv;
      for(int v = 0;v <= vMin;v++)
        Tens->Tension[Comp[p]][v]  += Pre[p]*Tens->VolSlabInv[v]*NSlabInv;
    }
  }
  else{
    double NSlabInv = 1./(double)(vMax + 1 - vMin);
    for(int v = vMin;v <= vMax;v++)
      for(int p=0;p<NPre;p++)
	Tens->Tension[Comp[p]][v]  += Pre[p]*Tens->VolSlabInv[v]*NSlabInv;
  }
}
//-------------3d-----------------
static void Tens3dAdd(struct beads *b, double ExtForce,PASSPORT prow,PASSPORT pcol){
  int Id1 = GET_POS(prow);
  int Id2 = GET_POS(pcol);
  int t1  = GET_TYPE(prow);
  int t2  = GET_TYPE(pcol);
  Tens3dAddPos(b,ExtForce,b->xv[Id1],b->xv[Id2],t1,t2);
}
void Tens3dAddPos(struct beads *b, double Force,double *Pos1,double *Pos2,int t1,int t2)
{
  struct TENS_PROF *Tens = &(b->Tens);
  double DistRel[4] = {0.,0.,0.,0.};
  double PosP1[3];
  double PosP2[3];
  for(int d=0;d<3;d++){
    PosP1[d] = Pos1[d] - (b->l[d]*.5 - Tens->RefPos[d]);
    PosP2[d] = Pos2[d] - (b->l[d]*.5 - Tens->RefPos[d]);
    PosP1[d] -= floor(PosP1[d]/b->l[d])*b->l[d];
    PosP2[d] -= floor(PosP2[d]/b->l[d])*b->l[d];
    DistRel[d] =  Pos2[d] - Pos1[d];
    DistRel[d] -= floor( DistRel[d]/(b->l[d]) + .5)*b->l[d];
    DistRel[3] += SQR(DistRel[d]);
  }
  DistRel[3] = sqrt(DistRel[3]);
  double Pre[6];
  Pre[0] = Force*DistRel[0]*DistRel[0]/DistRel[3];
  Pre[1] = Force*DistRel[1]*DistRel[1]/DistRel[3];
  Pre[2] = Force*DistRel[2]*DistRel[2]/DistRel[3];
  Pre[3] = Force*DistRel[0]*DistRel[1]/DistRel[3];
  Pre[4] = Force*DistRel[0]*DistRel[2]/DistRel[3];
  Pre[5] = Force*DistRel[1]*DistRel[2]/DistRel[3];
  Tens3dAddCompPos(Tens,PosP1,PosP2,Pre,b->l);
}
void Tens3dAddPre(struct beads *b, double *Pre,PASSPORT prow,PASSPORT pcol){
  struct TENS_PROF *Tens = &(b->Tens);
  int Id1 = GET_POS(prow);
  int Id2 = GET_POS(pcol);
  int t1  = GET_TYPE(prow);
  int t2  = GET_TYPE(pcol);
  double *Pos1 = b->xv[Id1];
  double *Pos2 = b->xv[Id2];
  Tens3dAddCompPos(Tens,Pos1,Pos2,Pre,b->l);
}
/** Divide the tension line in NPoint and add to every small cube crosse by the line the correspondent pressure */
static void Tens3dAddCompPos(struct TENS_PROF *Tens,double *Pos1,double *Pos2,double *Pre,double *Edge){
  int NPoint = 100;
  double PointInv = 1./(double)NPoint;
  double Deltav[3];
  double Pos[3];
  int vCurr[3];
  for(int d=0;d<3;d++){
    // starts from the smallest position
    Pos[d] = Pos1[d] < Pos2[d] ? Pos1[d] : Pos2[d];
    Deltav[d] = fabs((Pos2[d]-Pos1[d])/(double)NPoint);
    // forward/backward?
    if( fabs(Pos2[d] - Pos1[d]) > Edge[d]*.5){
      Deltav[d] = -fabs((MAX(Pos1[d],Pos2[d])-Edge[d]-Pos[d])/(double)NPoint);
    }
  }
  for(int p=0;p<NPoint;p++){
    for(int d=0;d<3;d++){
      vCurr[d] = (int)(Pos[d] * Tens->NSlab/Edge[d]);
      //assert(vCurr[d] >= 0 && vCurr[d] < Tens->NSlab);
      if(vCurr[d] < 0 || vCurr[d] >= Tens->NSlab) continue;
      // update the current position forward/backward
      Pos[d] += Deltav[d];
      // wrap the position to the other end
      if(Pos[d] < 0.){
	if(!Tens->Wrap[d]) return;
	Pos[d] = Edge[d] + Pos[d];
      }
    }
    int vTot = (vCurr[0]*Tens->NSlab+vCurr[1])*Tens->NSlab+vCurr[2];
    for(int c=0;c<Tens->NComp;c++){
      Tens->Tension[c][vTot] += Pre[c]*PointInv;
    }
  }
}
//-------------2d-----------------
static void Tens2dAdd(struct beads *b, double ExtForce,PASSPORT prow,PASSPORT pcol){
  int Id1 = GET_POS(prow);
  int Id2 = GET_POS(pcol);
  int t1  = GET_TYPE(prow);
  int t2  = GET_TYPE(pcol);
  Tens2dAddPos(b,ExtForce,b->xv[Id1],b->xv[Id2],t1,t2);
}
void Tens2dAddPos(struct beads *b, double Force,double *Pos1,double *Pos2,int t1,int t2)
{
  struct TENS_PROF *Tens = &(b->Tens);
  double DistRel[4] = {0.,0.,0.,0.};
  double PosP1[3];
  double PosP2[3];
  for(int d=0;d<2;d++){
    PosP1[d] = Pos1[d] - Tens->RefPos[d];
    PosP2[d] = Pos2[d] - Tens->RefPos[d];
    PosP1[d] -= floor(PosP1[d]/b->l[d])*b->l[d];
    PosP2[d] -= floor(PosP2[d]/b->l[d])*b->l[d];
    DistRel[d] = Pos2[d] - Pos1[d];
    DistRel[d] -= floor( DistRel[d]/(b->l[d]) + .5)*b->l[d];
    DistRel[3] += SQR(DistRel[d]);
  }
  PosP1[NORMAL] = Pos1[NORMAL] - Tens->RefPos[NORMAL] + b->l[NORMAL]*.5 ;
  PosP2[NORMAL] = Pos2[NORMAL] - Tens->RefPos[NORMAL] + b->l[NORMAL]*.5 ;
  PosP1[NORMAL] -= floor(PosP1[NORMAL]/b->l[NORMAL])*b->l[NORMAL];
  PosP2[NORMAL] -= floor(PosP2[NORMAL]/b->l[NORMAL])*b->l[NORMAL];
  DistRel[NORMAL] = Pos2[NORMAL] - Pos1[NORMAL];
  DistRel[NORMAL] -= floor( DistRel[NORMAL]/(b->l[NORMAL]) + .5)*b->l[NORMAL];
  DistRel[3] += SQR(DistRel[NORMAL]);
  DistRel[3] = sqrt(DistRel[3]);
  double Pre[6];
  Pre[0] = Force*DistRel[0]*DistRel[0]/DistRel[3];
  Pre[1] = Force*DistRel[1]*DistRel[1]/DistRel[3];
  Pre[2] = Force*DistRel[2]*DistRel[2]/DistRel[3];
  Pre[3] = Force*DistRel[0]*DistRel[1]/DistRel[3];
  Pre[4] = Force*DistRel[0]*DistRel[2]/DistRel[3];
  Pre[5] = Force*DistRel[1]*DistRel[2]/DistRel[3];
  Tens2dAddCompPos(Tens,PosP1,PosP2,Pre,b->l);
}
void Tens2dAddPre(struct beads *b, double *Pre,PASSPORT prow,PASSPORT pcol){
  struct TENS_PROF *Tens = &(b->Tens);
  int Id1 = GET_POS(prow);
  int Id2 = GET_POS(pcol);
  int t1  = GET_TYPE(prow);
  int t2  = GET_TYPE(pcol);
  double *Pos1 = b->xv[Id1];
  double *Pos2 = b->xv[Id2];
  double PosP1[3];
  double PosP2[3];
  for(int d=0;d<2;d++){
    PosP1[d] = Pos1[d] - Tens->RefPos[d];
    PosP2[d] = Pos2[d] - Tens->RefPos[d];
    PosP1[d] -= floor(PosP1[d]/b->l[d])*b->l[d];
    PosP2[d] -= floor(PosP2[d]/b->l[d])*b->l[d];
  }
  PosP1[NORMAL] = Pos1[NORMAL] - Tens->RefPos[NORMAL] + b->l[NORMAL]*.5;
  PosP2[NORMAL] = Pos2[NORMAL] - Tens->RefPos[NORMAL] + b->l[NORMAL]*.5;
  PosP1[NORMAL] -= floor(PosP1[NORMAL]/b->l[NORMAL])*b->l[NORMAL];
  PosP2[NORMAL] -= floor(PosP2[NORMAL]/b->l[NORMAL])*b->l[NORMAL];
  Tens2dAddCompPos(Tens,PosP1,PosP2,Pre,b->l);
}
/** Divide the tension line in NPoint and add to every small cube crosse by the line the correspondent pressure */
static void Tens2dAddCompPos(struct TENS_PROF *Tens,double *ExtPos1,double *ExtPos2,double *Pre,double *ExtEdge){
  int NPoint = 100;
  double PointInv = 1./(double)NPoint;
  double Deltav[2];
  double Pos[2];
  double Pos1[2];
  double Pos2[2];
  int vCurr[2];
  double Edge[2];
  Pos1[0] = hypot(ExtPos1[TANG1],ExtPos1[TANG2]);
  Pos2[0] = hypot(ExtPos2[TANG1],ExtPos2[TANG2]);
  Pos1[1] = ExtPos1[NORMAL];
  Pos2[1] = ExtPos2[NORMAL];
  Edge[0] = Tens->BoxRad;
  if(Pos1[0] >= Edge[0] || Pos2[0] >= Edge[0]) return;
  Edge[1] = ExtEdge[NORMAL];
  for(int d=0;d<2;d++){
    // starts from the smallest position
    Pos[d] = Pos1[d] < Pos2[d] ? Pos1[d] : Pos2[d];
    Deltav[d] = fabs((Pos2[d]-Pos1[d])/(double)NPoint);
    // forward/backward?
    if( fabs(Pos2[d] - Pos1[d]) > Edge[d]*.5){
      Deltav[d] = -fabs((MAX(Pos1[d],Pos2[d])-Edge[d]-Pos[d])/(double)NPoint);
    }
  }
  for(int p=0;p<NPoint;p++){
    for(int d=0;d<2;d++){
      vCurr[d] = (int)(Pos[d] * Tens->NSlab/Edge[d]);
      assert(vCurr[d] >= 0 && vCurr[d] < Tens->NSlab);
      // update the current position forward/backward
      Pos[d] += Deltav[d];
      // wrap the position to the other end
      if(Pos[d] < 0.){
	if(!Tens->Wrap[d]) return;
	Pos[d] = Edge[d] + Pos[d];
      }      //if(Deltav[d]<0.) printf("%lf %d %lf\n",Pos[d],d,Deltav[d]);
    }
    int vTot = (vCurr[0]*Tens->NSlab+vCurr[1]);
    for(int c=0;c<Tens->NComp;c++){
      Tens->Tension[c][vTot] += Pre[c]*PointInv*Tens->VolSlabInv[vCurr[0]];
    }
  }
}
//###############-WRITE-SUM-################################
void TensWrite(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  Tens->Step = (int)(b->step/(double)Tens->NAverage);
  if(Tens->CalcMode == TENS_3D)
    Tens3dWrite(b);
  else if(Tens->CalcMode == TENS_2D)
    Tens2dWrite(b);
  else if(Tens->CalcMode == TENS_NO) ;
  else
    TensLineWrite(b);
}
/** Overwrites the Tension*.dat files */
void TensLineWrite(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  if(Tens->FNorma <= 0.) return;
  if(Tens->NSlab == 0) return;
  double Vol = (b->l[0] * b->l[1] * b->l[2]);
  double InvVolume = -Tens->NSlab/Vol;
  double InvValues = 1. / Tens->NSlab;
  double InvNorma = 1. / Tens->FNorma;
  char stringa[60];
  FILE *FileToWrite;
  char Legend[256];
  debug("Tens] Summing the values of every processor");
  for(int l=0;l<Tens->NLine;l++){
    MPI_Reduce(ismaster ? MPI_IN_PLACE : Tens->Tension[l], ismaster ? Tens->Tension[l] : NULL, Tens->NSlab, MPI_DOUBLE, MPI_SUM, 0, comm_grid);
    if (!ismaster) {
      memset(Tens->Tension[l], 0, Tens->NSlab * sizeof(**Tens->Tension));
    }
  }
  for(int l=0;l<Tens->NComp*TYPE_MAX;l++){
    MPI_Reduce(ismaster ? MPI_IN_PLACE : Tens->Count[l], ismaster ? Tens->Count[l] : NULL, Tens->NSlab, MPI_DOUBLE, MPI_SUM, 0, comm_grid);
    if (!ismaster) {
      memset(Tens->Count[l], 0, Tens->NSlab * sizeof(**Tens->Count));
    }
  }
  if (!ismaster) return;
  //header
  debug("Tens] Writing the Tension.dat files");
  sprintf(Legend,"# 1=x ");
  //double Pressure[6];for(int c=0;c<6;c++)Pressure[c] = 0.;
  for(int t=0;t<TYPE_MAX;t++)
    for(int tt=t;tt<TYPE_MAX;tt++)
      sprintf(Legend+strlen(Legend),"%d=%d-%d ",TypeIdx(t,tt)+2,t,tt);
  for(int t=0;t<TYPE_MAX;t++)
    sprintf(Legend+strlen(Legend),"%d=%d ",Tens->NType+t+2,t);
  for(int c=0;c<2;c++){
    sprintf(stringa,"TensionS%09dZ%d.dat",Tens->Step,c);
    if( (FileToWrite = fopen(stringa,"w")) == NULL){
      fatal(EIO, "Could not open the Tension file");
    }
    mdpd_print_header(FileToWrite,b->mdpd);
    RigidPrintHeader(FileToWrite,b);
    fprintf(FileToWrite,"# %s\n",Legend);
   //body
    for(int v=0;v<Tens->NSlab;v++){
      fprintf(FileToWrite,"%lf ",v*InvValues*b->l[NORMAL]);
      for(int t=0;t<Tens->NType;t++){
	fprintf(FileToWrite," %g ",Tens->Tension[t+Tens->NType*c][v]*InvNorma*InvVolume);
	//Pressure[c] += Tens->Tension[t+Tens->NType*c][v]*InvNorma;
      }
      for(int t=0;t<TYPE_MAX;t++) 
	fprintf(FileToWrite," %g ",Tens->Count[t+TYPE_MAX*0][v]*InvNorma);
      fprintf(FileToWrite,"\n");
    }
    fclose(FileToWrite);
  }
  //for(int c=0;c<6;c++) printf("%d) %g %g\n",Tens->Step,Pressure[c]/Vol,b->virial[c]/Vol);
  for(int c=2;c<Tens->NComp;c++){
    sprintf(stringa,"TensionS%09dR%d.dat",Tens->Step,c-2);
    if( (FileToWrite = fopen(stringa,"w")) == NULL){
      fatal(EIO, "Could not open the Tension file");
    }
    mdpd_print_header(FileToWrite,b->mdpd);
    RigidPrintHeader(FileToWrite,b);
    fprintf(FileToWrite,"# %s\n",Legend);
    for(int v=0;v<Tens->NSlab;v++){
      fprintf(FileToWrite,"%lf ",v*InvValues*Tens->BoxRad);
      for(int t=0;t<Tens->NType;t++)
	fprintf(FileToWrite," %g ",-Tens->Tension[t+Tens->NType*c][v]*InvNorma);
      for(int t=0;t<TYPE_MAX;t++) fprintf(FileToWrite," %g ",Tens->Count[t+TYPE_MAX*1][v]*InvNorma);
      fprintf(FileToWrite,"\n");
    }
    fclose(FileToWrite);
  }
  //preparing for a new file
  for(int l=0;l<Tens->NComp;l++){
    memset(Tens->Tension[l],0.,Tens->NSlab*sizeof(double));
  }
  for(int l=0;l<TYPE_MAX;l++){
    memset(Tens->Count[l],0.,Tens->NSlab*sizeof(double));
  }
  Tens->FNorma = 0.;
}
/** Overwrites the Tension*.dat files */
void Tens3dWrite(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  debug("Tens] Summing the values of every processor");
  for(int l=0;l<Tens->NLine;l++){
    //MPI_Reduce(ismaster ? MPI_IN_PLACE : Tens->Tension[l], ismaster ? Tens->Tension[l] : NULL,CUBE(Tens->NSlab), MPI_DOUBLE, MPI_SUM, 0, comm_grid);
    MPI_Allreduce(MPI_IN_PLACE,Tens->Tension[l],CUBE(Tens->NSlab),MPI_DOUBLE,MPI_SUM,comm_grid);  
  }
  for(int l=0;l<TYPE_MAX;l++){
    //MPI_Reduce(ismaster ? MPI_IN_PLACE : Tens->Count[l], ismaster ? Tens->Count[l] : NULL, CUBE(Tens->NSlab), MPI_DOUBLE, MPI_SUM, 0, comm_grid);
    MPI_Allreduce(MPI_IN_PLACE,Tens->Count[l],CUBE(Tens->NSlab),MPI_DOUBLE,MPI_SUM,comm_grid);
  }
  Tens3dWriteAv(b);
  //preparing for a new file
  for(int l=0;l<Tens->NComp;l++){
    memset(Tens->Tension[l],0.,CUBE(Tens->NSlab)*sizeof(double));
  }
  for(int l=0;l<TYPE_MAX;l++){
    memset(Tens->Count[l],0.,CUBE(Tens->NSlab)*sizeof(double));
  }
  Tens->FNorma = 0.;
}
void Tens3dWriteAv(struct beads *b){
  if (!ismaster) return;
  struct TENS_PROF *Tens = &(b->Tens);
  if(Tens->FNorma <= 0.) return;
  if(Tens->NSlab == 0) return;
  debug("Tens] Writing the Tension.dat files");
  double Vol = b->l[0] * b->l[1] * b->l[2];
  double rho = mdpd_getrhocoex(b->mdpd) / 10.;
  double InvVolume = rho*CUBE(Tens->NSlab)/Vol;
  double InvValues = 1. / Tens->NSlab;
  double InvNorma = 1. / Tens->FNorma;
  int NFile = 3;//6;
  char stringa[60];
  FILE *FileToWrite;
  // print header
  char Legend[256];
  sprintf(Legend,"# 1=x 2=y 3=z ");
  for(int t=0;t<TYPE_MAX;t++)
    sprintf(Legend+strlen(Legend),"%d=%d ",Tens->NType+t+2,t);
  double Pressure[6];for(int c=0;c<6;c++)Pressure[c] = 0.;
  for(int c=0;c<NFile;c++){
    sprintf(stringa,"TensionS%09dC%d.xvl",Tens->Step,c);
    if( (FileToWrite = fopen(stringa,"w")) == NULL){
      fatal(EIO, "Could not open the Tension file");
    }
    fprintf(FileToWrite,"# l(%.1f %.1f %.1f) r(%.2f %.2f %.2f) v[%d] d[color]\n",b->l[0],b->l[1],b->l[2],Tens->RefPos[0],Tens->RefPos[1],Tens->RefPos[2],Tens->NSlab);
    mdpd_print_header(FileToWrite,b->mdpd);
    RigidPrintHeader(FileToWrite,b);
    fprintf(FileToWrite,"# %s\n",Legend);
    // print entries
    for(int vx=0;vx<Tens->NSlab;vx++){
      double x = vx*InvValues*b->l[TANG1];
      for(int vy=0;vy<Tens->NSlab;vy++){
	double y = vy*InvValues*b->l[TANG2];
	for(int vz=0;vz<Tens->NSlab;vz++){
	  double z = vz*InvValues*b->l[NORMAL];
	  int v = (vx*Tens->NSlab+vy)*Tens->NSlab+vz;
	  Pressure[c] += Tens->Tension[c][v]*InvNorma;
	  if(Tens->Count[0][v] <= 0 && Tens->Count[1][v] <= 0 && ABS(Tens->Tension[c][v]) <= 0)
	    continue;
	  fprintf(FileToWrite,"{x(%.3f %.3f %.3f)",x,y,z);
	  fprintf(FileToWrite," v( %lf %.2f %.2f)}\n",
		  -Tens->Tension[c][v]*InvNorma*InvVolume,
		  Tens->Count[0][v]*InvNorma*InvVolume,
		  Tens->Count[1][v]*InvNorma*InvVolume);
	}
      }
    }
    fclose(FileToWrite);
  }
  //Check
  for(int c=0;c<6;c++) printf("%d) %g %g\n",Tens->Step,Pressure[c]*InvNorma/Vol,b->virial[c]/Vol);
}
/** Overwrites the Tension*.dat files */
void Tens2dWrite(struct beads *b){
  struct TENS_PROF *Tens = &(b->Tens);
  debug("Tens] Summing the values of every processor");
  for(int l=0;l<Tens->NLine;l++){
    MPI_Allreduce(MPI_IN_PLACE,Tens->Tension[l],SQR(Tens->NSlab),MPI_DOUBLE,MPI_SUM,comm_grid);  
  }
  for(int l=0;l<TYPE_MAX;l++){
    MPI_Allreduce(MPI_IN_PLACE,Tens->Count[l],SQR(Tens->NSlab),MPI_DOUBLE,MPI_SUM,comm_grid);
  }
  Tens2dWriteAv(b);
  //preparing for a new file
  for(int l=0;l<Tens->NComp;l++){
    memset(Tens->Tension[l],0.,SQR(Tens->NSlab)*sizeof(double));
  }
  for(int l=0;l<TYPE_MAX;l++){
    memset(Tens->Count[l],0.,SQR(Tens->NSlab)*sizeof(double));
  }
  Tens->FNorma = 0.;
}
void Tens2dWriteAv(struct beads *b){
  if (!ismaster) return;
  struct TENS_PROF *Tens = &(b->Tens);
  if(Tens->FNorma <= 0.) return;
  if(Tens->NSlab == 0) return;
  debug("Tens] Writing the Tension.dat files");
  double Vol = b->l[0] * b->l[1] * b->l[2];
  double rho = mdpd_getrhocoex(b->mdpd) / 10.;
  double InvVolume = rho*CUBE(Tens->NSlab)/Vol;
  double InvValues = 1. / (double)Tens->NSlab;
  double InvNorma = 1. / (double)Tens->FNorma;
  int NFile = 3;//6;
  char stringa[60];
  FILE *FileToWrite;
  // print header
  char Legend[256];
  sprintf(Legend,"# 1=x 2=y 3=z ");
  for(int t=0;t<TYPE_MAX;t++)
    sprintf(Legend+strlen(Legend),"%d=%d ",Tens->NType+t+2,t);
  for(int c=0;c<NFile;c++){
    sprintf(stringa,"TensionS%09dL%d.xvl",Tens->Step,c);
    if( (FileToWrite = fopen(stringa,"w")) == NULL){
      fatal(EIO, "Could not open the Tension file");
    }
    fprintf(FileToWrite,"# l(%.1f %.1f %.1f) r(%.2f %.2f %.2f) v[%d] d[color]\n",Tens->BoxRad,b->l[NORMAL],1.,Tens->RefPos[0],Tens->RefPos[1],Tens->RefPos[2],Tens->NSlab);
    mdpd_print_header(FileToWrite,b->mdpd);
    RigidPrintHeader(FileToWrite,b);
    fprintf(FileToWrite,"# %s\n",Legend);
    // print entries
    for(int vx=0;vx<Tens->NSlab;vx++){
      double x = vx*InvValues*Tens->BoxRad;
      for(int vy=0;vy<Tens->NSlab;vy++){
	double y = vy*InvValues*b->l[NORMAL];
	int v = (vx*Tens->NSlab+vy);
	if(Tens->Count[0][v] <= 0 && Tens->Count[1][v] <= 0 && ABS(Tens->Tension[c][v]) <= 0)
	  continue;
	fprintf(FileToWrite,"{x(%.4f %.4f %.4f)",x,y,0.);
	fprintf(FileToWrite," v( %lf %.4f %.4f)}\n",
		Tens->Tension[c][v]*InvNorma*InvVolume,
		Tens->Count[0][v]*InvNorma*InvVolume,
		Tens->Count[1][v]*InvNorma*InvVolume);
      }
    }
    fclose(FileToWrite);
  }
}

