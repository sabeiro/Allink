/* $Id: Rigid.c 366 2013-08-08 14:19:08Z fuhrmans $ */
#include "fdmdpd2.h"
#include "rand.h"
/// Quaternion
typedef double QUADRI[4];
//-----------------distances----------------------
static double (*RigidDistance)(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidCylDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidPillDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidTiltDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidSphericalDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidElipsDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidParabDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidTorusDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidNoDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){return 0.;}
static double RigidWallDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidTiltWallDist(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidPoreDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static double RigidJanusDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano);
static void RigidCalcMinimal(struct beads *restrict b,struct NANO *Nano);
//----------------forces--------------------------
static void DefForceParam(struct NANO *Nano);
static double RigidForceLJ39(struct NANO *Nano,double Dr,double *Pot,double Sign);
static double RigidForceOld(struct NANO *Nano,double Dr,double *Pot,double Sign);
static double RigidForceCutOff(struct NANO *Nano,double Dr,double *Pot,double Sign);
static double RigidForceHamaker(struct NANO *Nano,double Dr,double *Pot,double Sign);
static double RigidForceCyl(struct NANO *Nano,double Dr,double *Pot,double Sign);
static double RigidForceLinear(struct NANO *Nano,double Dr,double *Pot,double Sign);
static double RigidForceHarmonic(struct NANO *Nano,double Dr,double *Pot,double Sign);
static double RigidForceNone(struct NANO *Nano,double Dr,double *Pot,double Sign){return 0.;}
static void UpdateAVel(struct NANO *Nano,double dt);
static double (*RigidForce)(struct NANO *Nano,double Dr,double *Pot,double Sign);
static void PrintForce(struct beads *b);
static void PrintForce3d(struct beads *b);
//------------------weighting-functions----------------
static double WeightFuncSquare(const double r);
static double WeightFuncCubic(const double r);
static double WeightFuncSpline(const double r);
static double WeightFuncNo(const double r){return 1.;}
static double (*WeightFunc)(double r) = WeightFuncNo;
static double DerWeightFuncSquare(const double r);
static double DerWeightFuncCubic(const double r);
static double DerWeightFuncSpline(const double r);
static double DerWeightFuncNo(const double r){return 0.;}
static double (*DerWeightFunc)(double r) = DerWeightFuncNo;
static double SignFuncNo(const double r,int type){return 1.;};
static double SignFuncType(const double r,int type){return type == 0 ? -1. : 1.;};
static double SignFuncJanus(const double r,int type);
static double SignFuncCyl(const double r,int type);
static double SignFuncWall(const double r,int type){return type == 0 ? 1. : -1.;};
static double SignFuncTorus(const double r,int type){return type == 0 ? 1. : 0.;};
static double SignFuncBoth(const double r,int type){return 1.;};
static double (*SignFunc)(double r,int type) = SignFuncNo;
static void PorePos(struct beads *restrict b);
static void StalkPos(struct beads *restrict b);
//---------------Transformation-------------------------
//void VectProd(double *v,double *u,double *Resp);
static void MatrVect(double *M,double *v,double *Resp);
static void MatrMatrix(double *M1,double *M2,double *Resp);
static void CopyMatrix(double *Resp,const double *M2);
static void Swap(double *a,double *b);
static void Transpose(double *M1);
static void RotMatrix(double *M,QUAT q);
static void RotMatrix16(double *M,QUAT q);
static double Dist2Axis(double *Pos,double *Axis);
//double Dist2AxisPerp(double *Pos,double *Axis, double *Perp);
static double ProjectOn(double *Axis,double *Pos);
static double ProjectOnResp(double *Axis,double *Pos,double *Resp);
static void PrintMatrix(double *M1);
FILE *Statistics;
//void RotMatr9(double *M,double *Axis,double Angle);
//#########################DEFINING/ALLOCATING###########################3
/** Point to the corrispondence distance and velocity verlet 
    functions */
static void RigidPointShape(int Shape){
  SignFunc = SignFuncType;
  WeightFunc = WeightFuncNo;
  DerWeightFunc = DerWeightFuncNo;
  RigidForce = RigidForceCutOff;
  if(VAR_IF_TYPE(Shape,SHAPE_NONE)){
    RigidDistance = RigidNoDistance;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_SPH)){
    RigidDistance = RigidSphericalDistance;
    RigidForce = RigidForceHamaker;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_TIP)){
    //RigidDistance = RigidParabDistance;
    RigidDistance = RigidElipsDistance;
    RigidForce = RigidForceLJ39;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_CYL)){
    RigidDistance = RigidCylDistance;
    //RigidForce = RigidForceLJ39;
    //RigidForce = RigidForceCyl;
    //WeightFunc = WeightFuncSquare;
    //DerWeightFunc = DerWeightFuncSquare;
    //SignFunc = SignFuncCyl;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_PILL)){
    RigidDistance = RigidPillDistance;
    WeightFunc = WeightFuncSquare;
    DerWeightFunc = DerWeightFuncSquare;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_TILT)){
    RigidDistance = RigidTiltDistance;
    //WeightFunc = WeightFuncSquare;
    //DerWeightFunc = DerWeightFuncSquare;
    SignFunc = SignFuncCyl;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_WALL)){
    RigidDistance = RigidTiltWallDist; 
    SignFunc = SignFuncWall;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_PORE)){
    RigidDistance = RigidTiltDistance;
    RigidForce = RigidForceNone;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_STALK)){
    RigidDistance = RigidTorusDistance;
    RigidForce = RigidForceNone;
    //SignFunc = SignFuncTorus;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_EXT)){
    RigidDistance = RigidTiltDistance;
    RigidForce = RigidForceHarmonic;
    //WeightFunc = WeightFuncSquare;
    //DerWeightFunc = DerWeightFuncSquare;
    SignFunc = SignFuncType;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_JANUS)){
    RigidDistance = RigidJanusDistance;
    //WeightFunc = WeightFuncCubic;
    //SignFunc = SignFuncJanus;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_TORUS)){
    RigidDistance = RigidTorusDistance;
    RigidForce = RigidForceNone;
    //SignFunc = SignFuncTorus;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_HARM)){
    RigidDistance = RigidTiltWallDist;
    RigidForce = RigidForceHarmonic;
    SignFunc = SignFuncBoth;
  }
  else if(VAR_IF_TYPE(Shape,SHAPE_UMBR)){
    RigidDistance = RigidWallDistance;
    RigidForce = RigidForceHarmonic;
    SignFunc = SignFuncBoth;
  }
  else
    fatal(EINVAL, "Invalid shape of the inclusion");
}
void DefForceParam(struct NANO *Nano){
  double Pot[2];
  Nano->CutOff = Nano->Radius*3.;
  RigidPointShape(Nano->Shape);
  if(RigidForce == RigidForceCutOff)
    Nano->CutOff = Nano->Radius+Nano->Coating*2.;
  Nano->Baseline = 0.;
  RigidForce(Nano,Nano->CutOff,Pot,0);
  Nano->Baseline = -Pot[0];
  Nano->DistThr = 0.;
  for(int i=10000;i>0;i--){
    double x = (Nano->CutOff)/10000.*i;
    double For = RigidForce(Nano,x,Pot,0);
    Pot[0] *= WeightFunc(0.);
    if(For >= Nano->ForThr){
      Nano->DistThr = x;
      Nano->PotThr = Pot[0];
      break;
    }
  }
  printf("DefForce %lf %lf %lf\n",Nano->DistThr,Nano->PotThr,Nano->Baseline);
}
static void RigidDefine(FILE *FH,struct NANO *Nano,double dt,int n){
  double Pos[3];
  double Axis[3];
  double Char[5];
  char Shape[10];
  int iShape=0;
  char cLine[256];
  fgets(cLine, LINE_MAX, FH);
  if( !Fetch(cLine,"x","%lf %lf %lf",Pos,Pos+1,Pos+2)){
    Nano->Pos[0] = 0.;
    Nano->Pos[1] = 0.;
    Nano->Pos[2] = 0.;
  }
  else{
    Nano->Pos[0] = Pos[0];
    Nano->Pos[1] = Pos[1];
    Nano->Pos[2] = Pos[2];
    Nano->OldPos[0] = Pos[0];
    Nano->OldPos[1] = Pos[1];
    Nano->OldPos[2] = Pos[2];
  }    
  if( !Fetch(cLine,"a","%lf %lf %lf",Axis,Axis+1,Axis+2)){
    Nano->Axis[TANG1]  = 0.;
    Nano->Axis[TANG2]  = 0.;
    Nano->Axis[NORMAL] = 1.;
  }
  else{
    Nano->Axis[0] = Axis[0];
    Nano->Axis[1] = Axis[1];
    Nano->Axis[2] = Axis[2];
    Normalize(Nano->Axis);
  }
  if( !Fetch(cLine,"c","%lf %lf %lf %lf",Char,Char+1,Char+2,Char+3))
    fatal(EINVAL,"Rigid characteristic are not specified");
  if( !Fetch(cLine,"s","%s",Shape))
    fatal(EINVAL,"Rigid shape is not specified");
  iShape = 0;
  if(!strcmp(Shape,"no")){
    VAR_ADD_TYPE(iShape,SHAPE_NONE);
  }
  else if(!strcmp(Shape,"sph"))
    VAR_ADD_TYPE(iShape,SHAPE_SPH);
  else if(!strcmp(Shape,"tip"))
    VAR_ADD_TYPE(iShape,SHAPE_TIP);
  else if(!strcmp(Shape,"cyl")){
    VAR_ADD_TYPE(iShape,SHAPE_CYL);
    VAR_ADD_TYPE(iShape,SHAPE_HEI);
    VAR_ADD_TYPE(iShape,SHAPE_SMOOTH);
  }
  else if(!strcmp(Shape,"tilt")){
    VAR_ADD_TYPE(iShape,SHAPE_TILT);
    VAR_ADD_TYPE(iShape,SHAPE_HEI);
    VAR_ADD_TYPE(iShape,SHAPE_SMOOTH);
  }
  else if(!strcmp(Shape,"pill")){
    VAR_ADD_TYPE(iShape,SHAPE_PILL);
    VAR_ADD_TYPE(iShape,SHAPE_HEI);
    VAR_ADD_TYPE(iShape,SHAPE_SMOOTH);
  }
  else if(!strcmp(Shape,"wall"))
    VAR_ADD_TYPE(iShape,SHAPE_WALL);
  else if(!strcmp(Shape,"pore"))
    VAR_ADD_TYPE(iShape,SHAPE_PORE);
  else if(!strcmp(Shape,"stalk"))
    VAR_ADD_TYPE(iShape,SHAPE_STALK);
  else if(!strcmp(Shape,"torus"))
    VAR_ADD_TYPE(iShape,SHAPE_TORUS);
  else if(!strcmp(Shape,"ext"))
    VAR_ADD_TYPE(iShape,SHAPE_EXT);
  else if(!strcmp(Shape,"janus"))
    VAR_ADD_TYPE(iShape,SHAPE_JANUS);
  else if(!strcmp(Shape,"harm"))
    VAR_ADD_TYPE(iShape,SHAPE_HARM);
  else if(!strcmp(Shape,"umbr"))
    VAR_ADD_TYPE(iShape,SHAPE_UMBR);
  else
    fatal(EINVAL,"Rigid shape %s not valid",Shape);
  RigidPointShape(iShape);
  Nano->Mass = 10.;
  Nano->Radius = Char[0];
  Nano->CutOff = Nano->Radius*3.;
  Nano->Viscosity = 5.;
  Nano->Hamaker = Char[1];
  Nano->Height = Char[2];
  Nano->Coating = Char[3];
  Nano->Shape = iShape;
  Nano->Pot[0] = 0.;
  Nano->Pot[1] = 0.;
  Nano->Inter = 0;
  Nano->ForThr = 50.;
  Nano->SwitchOff = 1.;
  DefForceParam(Nano);
  for(int d=0;d<3;d++){
    Nano->Vel[d] = 0.;
    Nano->AVel[d] = 0.;
  }
  Nano->Gamma = 3.*3.14*Nano->Radius*Nano->Viscosity;
  Nano->Zeta = sqrt(12.*2.*1.*Nano->Gamma/Nano->Mass/dt);
  char FileName[60];
  sprintf(FileName,"outputNano%d.dat",n);
  Nano->Save = fopen(FileName,"w");
  if(Nano->Save == NULL) novm("Nano.Save");
  if(ismaster){
    fprintf(Nano->Save,"# r %lf h %lf m %lf z %lf g %lf dt %lf\n",Nano->Radius,Nano->Hamaker,Nano->Mass,Nano->Zeta,Nano->Gamma,dt);
    fprintf(Nano->Save,"#1Time 2Rad 3Inter 4Pot 5Chem 6Cons 7AMom 8VelLat 9VelNorm 10VelAng 11Pre0 12Pre1 13Pre2\n");
    fflush(Nano->Save);
  }
}
/** Header for the outputNano.dat file */
void RigidPrintHeader(FILE *FH, const struct beads *b)
{
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    if(Nano->Shape != 0){
      char Shape[10];
      sprintf(Shape,"no");
      if(VAR_IF_TYPE(Nano->Shape,SHAPE_SPH)) sprintf(Shape,"sph");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_CYL)) sprintf(Shape,"cyl");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_PILL)) sprintf(Shape,"pill");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_TILT)) sprintf(Shape,"tilt");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_WALL)) sprintf(Shape,"wall");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_PORE)) sprintf(Shape,"pore");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_TIP)) sprintf(Shape,"tip");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_EXT)) sprintf(Shape,"ext");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_JANUS)) sprintf(Shape,"janus");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_TORUS)) sprintf(Shape,"torus");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_STALK)) sprintf(Shape,"stalk");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_HARM)) sprintf(Shape,"harm");
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_UMBR)) sprintf(Shape,"umbr");
      fprintf(FH,"# Rigid x(%.2f %.2f %.2f) a(%.3g %.3g %.3g) c(%.2f %.2f %.2f %lf) s{%s} \n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],Nano->Axis[0],Nano->Axis[1],Nano->Axis[2],Nano->Radius,Nano->Hamaker,Nano->Height,Nano->Coating,Shape);
    }
  }   
}
void PrintForce(struct beads *b){
  printf("Printing the potential and exit\n");
  struct NANO *Nano = b->Nano;
  double Pot[2] = {0.,0.};
  double RnpMin = 0.5;
  double RnpMax = 6.0;
  double Delta = (RnpMax-RnpMin)/11.;
  int Nppc = 32;//mdpd_getN(b)
  double Rho = mdpd_getrhocoex(b->mdpd)/(double)Nppc*CUBE(mdpd_getRe(b->mdpd));
  double HamakerVal[11] = {1.,1.5,1.7,2.,2.5,3.0,4.,6.,10.,14.,16.};
  double RadVal[10] = {0.5,1.,1.5,1.7,2.0,2.5,3.0,4.0,5.0,6.0};
  FILE *NrgFile = fopen("Eps2InterfNrg.dat","w");
  fprintf(NrgFile,"#Ref tension gamma = kBT/b^2 sqrt(chi/6) = 2.88 (30) = 4.08 (60)\n");
  for(int h=0;h<11;h++){
    Nano->Hamaker = HamakerVal[h];
    for(int r = 0;r<10;r++){
      double Rnp = RadVal[r];
      double InterfNrg = 0.;
      double InterfCount=0.;
      /* char File2Save[60]; */
      //sprintf(File2Save,"PotentialR%03.1fH%.1f.dat",Rnp,Nano->Hamaker);
      //FILE *Ciccia = fopen(File2Save,"w");
      //fprintf(Ciccia,"#x PotCut ForceCut PotHam ForceHam\n");
      Nano->Radius = Rnp;
      int NStep = 1000;
      double OldPot = 0.;
      double StartPos = 0.;
      for(int i=0;i<NStep;i++){
	double x = (Nano->CutOff)/(double)NStep*i;
	double For = RigidForce(Nano,x,Pot,0);
	//fprintf(Ciccia,"%lf %lf %lf\n",x,Pot[0],For);
	if(Pot[0] < 0.){
	  InterfNrg += Pot[0]*x*x*Rho*4.*M_PI;
	  InterfCount += 1.;
	  if(OldPot > 0.)
	    StartPos = x;
	}
	OldPot = Pot[0];
      }
      double Norm = InterfCount/(Nano->CutOff - StartPos);
      fprintf(NrgFile,"%lf %lf %lf %lf\n",Nano->Radius,Nano->Hamaker,InterfNrg/Norm,InterfNrg*Rho/(SQR(Nano->Radius)*4.*M_PI)/InterfCount);
      //fclose(Ciccia);
    }
  }
  fclose(NrgFile);
}
void PrintForce3d(struct beads *b){
  //printf("Printing the potential and exit\n");
  struct NANO *Nano = b->Nano;
  double Pot[2] = {0.,0.};
  FILE *CicciaPot = fopen("Potential3d.dat","w");
  FILE *CicciaFor = fopen("Force3d.dat","w");
  //FILE *Cicciax = fopen("Pot3dx.dat","w");
  //FILE *Cicciaz = fopen("Pot3dz.dat","w");
  for(double y = 0.;y<b->l[1];y+=b->l[1]*.005){
    for(double z = 0.;z<b->l[2];z+=b->l[2]*.005){
      double dr[3] = {0.,0.,0.};
      double Pos[3] = {Nano->Pos[0],y,z};
      double NPos[3];
      for(int d=0;d<3;d++)
	NPos[d] = Nano->Pos[d]-floor(Nano->Pos[d]/b->l[d])*b->l[d];
      double dr2 = RigidDistance(Pos,NPos,b->l,dr,Nano);
      double Sign = SignFunc(Pos[NORMAL],0);
      //if(dr2>SQR(Nano->CutOff))continue;
      double Dr = sqrt(dr2);
      double DrInv = 1./Dr;
      //Cons = Sign;
      double Cons = RigidForce(Nano,Dr,Pot,Sign);
      Cons = WeightFunc(Pos[NORMAL])*Cons;// - DerWeightFunc(Pos[NORMAL])*Pot[0];
      Pot[0] *= WeightFunc(Pos[NORMAL]);
      Cons *= Nano->SwitchOff;
      if(Cons > 10.) Cons = 10.;
      if(Pot[0] > 10.) Pot[0] = 10.;
      fprintf(CicciaPot,"%lf %lf %lf\n",y,z,Pot[0]);
      fprintf(CicciaFor,"%lf %lf %lf\n",y,z,Cons*dr[2]*DrInv);
      //fprintf(Cicciax,"%lf %lf %lf\n",y,z,Cons*dr[1]*DrInv);
      //fprintf(Cicciaz,"%lf %lf %lf\n",y,z,Cons*dr[2]*DrInv);
    }
  }
  fclose(CicciaPot);
  fclose(CicciaFor);
  //fclose(Cicciax);
  //fclose(Cicciaz);
  //exit(0);
}
void RigidLoad(struct beads *b, FILE *FH){
  printf("NewRigidNewMain\n");
  char buf[256];
  struct NANO *Nano;
  debug("Rigid] Reading configurations");
  double EdgeMin = MIN(MIN(b->l[0],b->l[1]),b->l[2]);
  fpos_t PosTemp;
  fgetpos(FH,&PosTemp);
  b->NNano = 0;
  do {
    fgets(buf, sizeof(buf), FH);
    if(strstr(buf, "# Rigid") == buf)
      b->NNano++;
    else 
      break;
  } while (1==1);
  fsetpos(FH,&PosTemp);
  assert(b->NNano >= 0);
  b->Nano = calloc(b->NNano,sizeof(*b->Nano));
  if (b->Nano == NULL) novm("Nano");
  for(int n=0;n<b->NNano;n++){
    Nano = b->Nano + n;
    fgetpos(FH,&PosTemp);
    RigidDefine(FH,Nano,b->dt,n);
    Nano->CutOff = MIN(EdgeMin,Nano->CutOff);
  }
  RigidPointShape(Nano->Shape);
  if(ismaster) PrintForce3d(b);
  Statistics = fopen("RigidStatistics.dat","w");
}
/**
   Frees the nanoparticle
*/
void RigidFree(struct beads *b)
{
  if(b->NNano == 0) return;
  debug("Rigid] Freeing");
  for (int i = 0; i < b->NNano; i++)
    fclose(b->Nano[i].Save);
  free(b->Nano);
  fclose(Statistics);
}
//####################THE#FORCE#LOOP###########################
static void Mutual(struct beads *restrict b){
  return;
  //----------Mutual-interaction------------
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_TIP)) continue;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_WALL)) continue;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_EXT)) continue;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_STALK)) continue;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_PORE)) continue;
    for(int nn=n+1;nn<b->NNano;nn++){
      struct NANO *Nano1 = b->Nano + nn;
      if(VAR_IF_TYPE(Nano1->Shape,SHAPE_TIP)) continue;
      if(VAR_IF_TYPE(Nano1->Shape,SHAPE_WALL)) continue;
      if(VAR_IF_TYPE(Nano1->Shape,SHAPE_EXT)) continue;
      if(VAR_IF_TYPE(Nano1->Shape,SHAPE_STALK)) continue;
      if(VAR_IF_TYPE(Nano1->Shape,SHAPE_PORE)) continue;
      double dr[3] = {0.,0.,0.};
      double dr2 = 0.;
      // FIXME: Back-fold
      for(int d=0;d<3;d++){
	dr[d] = Nano->Pos[d] - Nano1->Pos[d];
	dr[d] -= floor(dr[d]/b->l[d] + .5)*b->l[d];
      }
      dr2 = SQR(dr[TANG1]) + SQR(dr[TANG2]);
      double CutOff = 1.2*(Nano1->Radius + Nano->Radius);
      if(dr2 > SQR(CutOff)) continue;
      double Dr = sqrt(dr2);
      double Pot[2] = {0.,0.};
      //double Cons = RigidForce(Nano1,Dr,Pot,1);
      double Cons = -500.*(Dr - CutOff);
      Nano1->Pot[0] += Pot[0];
      Nano1->Pot[1] += Pot[1];
      Nano->Pot[0] += Pot[0];
      Nano->Pot[1] += Pot[1];
      for (int d = 0; d < 2; d++) {
	//if(Nano->Shape == SHAPE_TILT && d == NORMAL) continue;
	Nano1->Force[d] -= Cons*dr[d] / Dr;
	Nano->Force[d]  += Cons*dr[d] / Dr;
	Nano1->Pre[d]  += Nano1->Force[d] * dr[d];
	Nano->Pre[d]   += Nano->Force[d]  * dr[d];
      }
      TensStressPos(b,Cons,Nano1->Pos,Nano->Pos,2,2);
    }
  }
}
void RigidCalcForces1(struct beads *restrict b){
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    int NNano = b->NNano;
    int NPart = b->N[0]*b->n[0];
    int Rank = 0;
    int Size = 0;
    double PosBf[3];
    for(int d = 0;d<3;d++){
      PosBf[d] = Nano->Pos[d] - floor( Nano->Pos[d]/b->l[d] )*b->l[d];
    }
    double CutOff = Nano->Radius*3.;
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Size);
    for(int d=0;d<4;d++)
      Nano->ForceRMS[d] = 0.;
    for(int d=0;d<3;d++){
      Nano->Force[d] = 0.;
      Nano->Pot[d] = 0.;
      Nano->Pre[d] = 0.;
    }
    Nano->Inter = 0;
    const VEC lh = {.5 * b->l[0], .5 * b->l[1], .5 * b->l[2]};
    size_t count = (groups_per_row + groups_per_col) * b->groupsize;
    count = 0;
    for (int i = 0; i < b->blocks; i++)
      count += b->n[i];
    int first = col_index*b->groupsize; 
    int last = (col_index+1)*b->groupsize;
    for (int p = first;p<last;p++){
      double Dis = 0.;
      double Cons = 0.;
      double Pot = 0.;
      int type = GET_TYPE(b->passport[p]);
      int pId = GET_ID(b->passport[p]);
      VEC dr;
      double dr2=0.;
      VEC2 *xv = b->xv + pId;
      double xv_bf[3];
      for(int d=0;d<3;d++){
	dr[d] = Nano->Pos[d] - xv[0][d];
	dr[d] -= floor( dr[d]/(b->l[d]) + .5)*b->l[d];
	xv_bf[d] = xv[0][d];
	xv_bf[d] -= floor( xv[0][d]/b->l[d])*b->l[d];
      }
      dr2 = SQR(dr[TANG1]) + SQR(dr[TANG2]);// + SQR(dr[2]);
      if(dr2 > SQR(CutOff) ) continue;
      if(xv_bf[NORMAL] > Nano->Pos[NORMAL] + Nano->Height*.5){
	dr2 += SQR(xv_bf[NORMAL] - Nano->Pos[NORMAL] + Nano->Height*.5);
      }
      if( xv_bf[NORMAL] < Nano->Pos[NORMAL] - Nano->Height*.5){
	dr2 += SQR(xv_bf[NORMAL] - Nano->Pos[NORMAL] - Nano->Height*.5);
      }
      dr[NORMAL] = 0.;
      double Dr = sqrt(dr2);
      if(dr2 > SQR(CutOff) ) continue;
      double idr = fabs(1./(Dr-Nano->Radius+1.));
      if(dr2 < SQR(.6*Nano->Radius)){//LJ cut off 
	Cons = 100.;
	Pot = -100. / idr + 1000.*(1.+.6*Nano->Radius) + 0.0369862318104964;
	Nano->Pot[1] += 1000.*+.6*Nano->Radius;
      }
      else {
	double idr2 = SQR(idr);
	double idr4 = SQR(idr2);
	double idr6 = CUBE(idr2);
	if(type == 1){
	  Cons = Nano->Hamaker * idr4 * (9. * idr4 * idr2 + 3.);
	  Pot = Nano->Hamaker * idr * idr2 * (idr6 + 1.) + 0.0369862318104964;
	  Nano->Pot[1] += Nano->Hamaker*3/(Nano->Radius)*idr*idr2*(3*idr6+1.);
	}
	else{
	  Cons = Nano->Hamaker * idr4 * (9. * idr4 * idr2 - 3.);
	  Pot = Nano->Hamaker * idr * idr2 * (idr6 - 1.) + 0.0369862318104964;
	  Nano->Pot[1] += Nano->Hamaker*3/(Nano->Radius)*idr*idr2*(3*idr6-1.);
	}
	//b->e[1] += Pot;
	Nano->Pot[0] += Pot;
      }
      //printf("%d/%d) id %d t %d (%lf %lf) %lf\n",p,b->step,pId,type,Cons,Pot,dr2);
      /* #ifndef WALL */
      /*       for(d=0;d<3;d++)  */
      /* #endif// */
      for(int d=0;d<3;d++) {
	b->f[pId][d] -= (Cons)*dr[d]*idr/Nano->Radius;
	Nano->Force[d] += (Cons/Nano->Mass)*dr[d]*idr/Nano->Radius;
	// -----------------pressure--------------
	//Nano->Pre[d] += Nano->Force[d] * dr[d];
	b->virial[d] -= Cons * SQR(dr[d])/idr;
	/* b->f[pId][d] -= (Cons)*dr[d]/Dr; */
	/* Nano->Force[d] += (Cons/Nano->Mass)*dr[d]/Dr; */
	/* b->virial[d] -= Cons * SQR(dr[d])/Dr; */
      }
      Nano->Inter++;
      //printf("%d\n",Nano->Inter);
    }
    if(ismaster){
      for(int c=0;c<4;c++)
	Nano->ForceRMS[c] = 0.;
      double Ran = 0.;
      double Dis = 0.;
      for(int d=0;d<3;d++){
	Ran = Nano->Zeta*(2.*rng_uniform(rng) - 1.);
	Dis = - Nano->Gamma*Nano->Vel[d];
	Nano->ForceRMS[0] += SQR(Nano->Force[d]);
	Nano->ForceRMS[1] += SQR(Dis);
	Nano->ForceRMS[2] += SQR(Ran);
	Nano->Force[d] += (Ran + Dis);
	Nano->ForceRMS[3] += SQR(Nano->Force[d]);
      }
    }
    double VelLat = sqrt(SQR(Nano->Vel[TANG1]) + SQR(Nano->Vel[TANG2]));
    double VelNorm = fabs(Nano->Vel[NORMAL]);
    fprintf(Nano->Save,"%d %lf %lf %d %lf %lf %lf %lf %lf %lf %lf\n",b->step,Nano->Radius,Nano->Pot[0],Nano->Inter,Nano->ForceRMS[0],Nano->Pot[1],VelLat,VelNorm,Nano->Pre[0],Nano->Pre[1],Nano->Pre[2]);
    //if (ismaster)
    {
      /* Nano particle - Nano particle interaction */
      /*     for(int i=0;i<b->NNano;i++) { */
      /*       for(int j=i+1;j<b->NNano;j++) { */
      /* 	double dr[3], dr2=0., f=0., idr, tmp; */	
      /* 	for(int d=0;d<3;d++) { */
      /* 	  dr[d] = Nano->Pos[3*j+d] - Nano->Pos[3*i+d]; */
      /* 	  while(dr[d] >  lh[d]) dr[d] -= b->l[d]; */
      /* 	  while(dr[d] < -lh[d]) dr[d] += b->l[d]; */
      /* 	  dr2 += SQR(dr[d]); */
      /* 	} */
      /* 	if(dr2 > 3.*Nano->Radius) continue; */
      /* 	idr = 1./sqrt(dr2); */
      /* 	//	f = force_nano_nano(idr,Nano); */
      /* 	//	tmp = .5*potential_nano_nano(idr,Nano); */
      /* 	Nano->u[i] += tmp; */
      /* 	Nano->u[j] += tmp; */	
      /* 	for(int d=0;d<3;d++) { */
      /* 	  Nano->f[3*j+d] += f * dr[d]; */
      /* 	  Nano->f[3*i+d] -= f * dr[d]; */
      /* 	} */
      /*       } */
      /*     } */
    }
    MPI_Allreduce(MPI_IN_PLACE, Nano->Force, 3*NNano, MPI_DOUBLE, MPI_SUM,comm_grid);
    MPI_Allreduce(MPI_IN_PLACE, Nano->Pot, 3*NNano, MPI_DOUBLE, MPI_SUM,comm_grid);
  }
  return 0;
}
void RigidCalcForces(struct beads *restrict b){
  struct TENS_PROF *Tens = &(b->Tens);
  if(Tens->CalcMode!=TENS_NO) IfTens = 1;
  if(b->NNano == 0 ) return;
  int first = col_index * b->groupsize; 
  int last = (col_index + 1) * b->groupsize;
  double Pot[2] = {0.,0.};
  double TotNrg = 0.;
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    //------------Inizialize-Nano-------------
    RigidPointShape(Nano->Shape);
    for(int d=0;d<3;d++){
      Nano->AMom[d] = 0.;
      Nano->AMomTemp[d] = 0.;
      Nano->Force[d] = 0.;
      Nano->Pre[d] = 0.;
    }
    Nano->Pot[0] = 0.;
    Nano->Pot[1] = 0.;
    Nano->Inter = 0;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_SPH)){
      if(1==0){
	Nano->SwitchOff = Nano->Coating;
      }
    }
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_NONE))  continue;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_STALK)) continue;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_TORUS)) continue;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_PORE))  continue;
    /* if(VAR_IF_TYPE(Nano->Shape,SHAPE_CYL)){ */
    /*   RigidCalcMinimal(b,Nano); */
    /*   continue; */
    /* } */
    //------------Calculate-interactions------
    for(int p = first; p <last ; p++){
      PASSPORT Pass1 = b->passport[p];
      int type = (int)GET_TYPE(Pass1);
      double dr[3] = {0.,0.,0.};
      double Pos[3] = {b->xv[p][0],b->xv[p][1],b->xv[p][2]};
      double NPos[3]= {Nano->Pos[0],Nano->Pos[1],Nano->Pos[2]};
      for(int d=0;d<3;d++)
	NPos[d] = Nano->Pos[d]-floor(Nano->Pos[d]/b->l[d])*b->l[d];
      double dr2 = RigidDistance(Pos,NPos,b->l,dr,Nano);
      if(dr2>SQR(Nano->CutOff)) continue;
      double Sign = SignFunc(Pos[NORMAL],type);
      double Dr = sqrt(dr2);
      double DrInv = 1./Dr;
      double Cons = RigidForce(Nano,Dr,Pot,Sign);
      Cons = WeightFunc(Pos[NORMAL])*Cons;// - DerWeightFunc(Pos[NORMAL])*Pot[0];
      //printf("%d %lf\n",p,Pot[0]);
      Pot[0] *= WeightFunc(Pos[NORMAL]);
      Cons *= Nano->SwitchOff;
      b->e[4] += Pot[0]*Nano->SwitchOff;
      b->e[0] += Pot[0]*Nano->SwitchOff;
      Nano->Pot[0] += Pot[0]*Nano->SwitchOff;
      Nano->Pot[1] += Pot[1]*Nano->SwitchOff;
      Nano->Inter++;
      TotNrg += Pot[0];
      double idr = fabs(1./(Dr-Nano->Radius+1.));
      for(int d=0;d<3;d++) {
	/* b->f[p][d] -=  Cons*dr[d]*idr; */
	/* Nano->Force[d] += Cons*dr[d]*idr; */
	/* b->virial[d] -= Cons*SQR(dr[d])/idr; */
	b->f[p][d]     -= Cons * dr[d] * DrInv;
	Nano->Force[d] += Cons * dr[d] * DrInv;
	b->virial[d]   -= Cons * SQR(dr[d])*DrInv;
      	Nano->Pre[d]   += Cons * SQR(dr[d]) * DrInv;
      	Nano->AMom[d]  += Nano->AMomTemp[d]*Cons;
      }
      b->virial[3] += dr[0] * dr[1] * Cons * DrInv;
      b->virial[4] += dr[0] * dr[2] * Cons * DrInv;
      b->virial[5] += dr[1] * dr[2] * Cons * DrInv;
      TensStressPos(b,Cons,b->xv[p],NPos,type,2);
    }
    MPI_Allreduce(MPI_IN_PLACE, Nano->Force, 3, MPI_DOUBLE, MPI_SUM, comm_grid);
    MPI_Allreduce(MPI_IN_PLACE, Nano->Pot, 2, MPI_DOUBLE, MPI_SUM, comm_grid);
  }
  if(ismaster){
    Mutual(b);
    for(int n=0;n<b->NNano;n++){
      struct NANO *Nano = b->Nano + n;
      //----------thermostat---------------------
      for(int c=0;c<4;c++) Nano->ForceRMS[c] = 0.;
      for(int d=0;d<3;d++){
	double Ran = Nano->Zeta  * (2.*rng_uniform(rng) - 1.);
	double Dis = - Nano->Gamma*Nano->Vel[d];
	Nano->ForceRMS[0] += sqrt(SQR(Nano->Force[d]));
	Nano->ForceRMS[1] += sqrt(SQR(Dis));
	Nano->ForceRMS[2] += sqrt(SQR(Ran));
	Nano->ForceRMS[3] += SQR(Nano->AMom[d]);
	Nano->Force[d] += (Ran + Dis);
      }
      //------------update-the-file--------------------
      double VelLat = sqrt(SQR(Nano->Vel[TANG1]) + SQR(Nano->Vel[TANG2]));
      double VelNorm = Nano->Vel[NORMAL];
      double VelAng = sqrt( SQR(Nano->AVel[0]) + SQR(Nano->AVel[1]) + SQR(Nano->AVel[2]) );
      fprintf(Nano->Save,"%.3f %d %.5g %.5g ",b->time,Nano->Inter,Nano->Pot[0],Nano->Pot[1]);
      fprintf(Nano->Save,"%.5g %.5g ",Nano->Force[NORMAL],Nano->ForceRMS[0]);
      //fprintf(Nano->Save,"%.3g %.3g %.3g ",VelLat,VelNorm,VelAng);
      //fprintf(Nano->Save,"%.3g %.3g %.3g ",b->virial[0],b->virial[1],b->virial[2]);
      fprintf(Nano->Save,"%lf %lf %lf ",Nano->Pos[0],Nano->Pos[1],b->Nano->Pos[2]);
      fprintf(Nano->Save,"%.1f %.1f %.1f\n",Nano->Pre[0],Nano->Pre[1],Nano->Pre[2]);
    }
  }
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    MPI_Allreduce(MPI_IN_PLACE, Nano->Force, 3, MPI_DOUBLE, MPI_SUM,comm_grid);
    MPI_Allreduce(MPI_IN_PLACE, Nano->Pot, 2, MPI_DOUBLE, MPI_SUM,comm_grid);
  }
}
void RigidChNrg(struct beads *restrict b,int Ch){
  int pFirst = col_index * b->groupsize + b->N[b->local_n_b[0]]*Ch; 
  int pLast = pFirst + b->N[b->local_n_b[0]];
  double Pot[2] = {0.,0.};
  int Type = 0;
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    if(Nano->Shape == SHAPE_NONE){break;}
    RigidPointShape(Nano->Shape);
    double NPos[3];
    Nano->SwitchOff = Nano->Coating;
    for(int d=0;d<3;d++)
      NPos[d] = Nano->Pos[d] - floor(Nano->Pos[d]/b->l[d])*b->l[d];
    for(int p = pFirst; p <pLast ; p++){
      PASSPORT Pass1 = b->passport[p];
      int type = (int)GET_TYPE(Pass1);
      double dr[3] = {0.,0.,0.};
      double Pos[3] = {b->xv[p][0],b->xv[p][1],b->xv[p][2]};
      double dr2 = RigidDistance(Pos,NPos,b->l,dr,Nano);
      if(dr2 > SQR(Nano->CutOff)) continue;
      double Dr = sqrt(dr2);
      double Cons = RigidForce(Nano,Dr,Pot,type+Type)*Nano->SwitchOff;
      b->e[0] += Pot[0]*Nano->SwitchOff;
    }
  }
}
void RigidMonNrg(struct beads *restrict b,int p){
  double Cons = 0.;
  double Pot[2] = {0.,0.};
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    if(Nano->Shape == SHAPE_NONE){break;}
    RigidPointShape(Nano->Shape);
    double NPos[3];
    Nano->SwitchOff = Nano->Coating;
    for(int d=0;d<3;d++)
      NPos[d] = Nano->Pos[d] - floor(Nano->Pos[d]/b->l[d])*b->l[d];
    PASSPORT Pass1 = b->passport[p];
    int type = (int)GET_TYPE(Pass1);
    double dr[3] = {0.,0.,0.};
    double Pos[3] = {b->xv[p][0],b->xv[p][1],b->xv[p][2]};
    double dr2 = RigidDistance(Pos,NPos,b->l,dr,Nano);
    if(dr2 > SQR(Nano->CutOff)) continue;
    double Dr = sqrt(dr2);
    Cons = RigidForce(Nano,Dr,Pot,type)*Nano->SwitchOff;
    b->e[0] += Pot[0]*Nano->SwitchOff;
  }
}
//#######################DISTANCES#############################
static double (*RigidDistance)(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano) = NULL;
/** Distance from a cylinder */
static double RigidPillDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  double PosRel[3] = {0.,0.,0.};
  for(int d=0;d<3;d++){
    dr[d] = NPos[d] - Pos[d];   
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
    PosRel[d] = dr[d];
  }
  double dr2 = SQR(dr[TANG1]) + SQR(dr[TANG2]);
  /* The cylinder ends with a cupola, if the particle is above the
     height of the cylinder, the cylinder it is seen as
     two spheres in \pm Nano->Height*.5 */
  if(dr[NORMAL] > Nano->Height*.5){
    dr[NORMAL] = dr[NORMAL] - Nano->Height*.5;
    dr2 += SQR(dr[NORMAL]);
  }
  else if(dr[NORMAL] < - Nano->Height*.5){
    dr[NORMAL] = dr[NORMAL] + Nano->Height*.5;
    dr2 += SQR(dr[NORMAL]);
  }
  else
    dr[NORMAL] = 0.;
  Pos[NORMAL] = fabs(dr[NORMAL]/Nano->Radius);
  /*Position of the surface point intersected by the normal 
    connecting the axis with the monomer position */
  double HeiOnAxis = ProjectOn(Nano->Axis,PosRel);
  double AxPos[4] = {0.,0.,0.,0.};
  for(int d=0;d<3;d++){
    AxPos[d] = -Nano->Axis[d]*HeiOnAxis + dr[d];
    AxPos[3] += SQR(AxPos[d]);
  }
  AxPos[3] = sqrt(AxPos[3]);
  for(int d=0;d<3;d++){
    NPos[d] = Nano->Pos[d] + Nano->Axis[d]*HeiOnAxis + Nano->Radius*AxPos[d]/AxPos[3];
  }
  return dr2;
}
static double RigidCylDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  double dr2 = 0.;
  for(int d=0;d<3;d++){
    dr[d] = NPos[d] - Pos[d];
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
    Pos[d] -= floor(Pos[d]/Edge[d])*Edge[d];
  }
  dr2 = SQR(dr[TANG1]) + SQR(dr[TANG2]);// + SQR(dr[2]);
  //if(dr2 > SQR(Nano->Radius*3) ) return 100.;
  if(Pos[NORMAL] > NPos[NORMAL] + Nano->Height*.5){
    dr2 += SQR(Pos[NORMAL] - NPos[NORMAL] + Nano->Height*.5);
  }
  if( Pos[NORMAL] < Nano->Pos[NORMAL] - Nano->Height*.5){
    dr2 += SQR(Pos[NORMAL] - NPos[NORMAL] - Nano->Height*.5);
  }
  dr[NORMAL] = 0.;
  return dr2;
}
static double RigidCylDistance1(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  double PosRel[3] = {0.,0.,0.};
  for(int d=0;d<3;d++){
    dr[d] = NPos[d] - Pos[d];
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
    PosRel[d] = dr[d];
  }
  // if and how much above the cylinder (for the smoothing)
  if(dr[NORMAL] > Nano->Height*.5){
    /* Pos[NORMAL] = fabs((dr[NORMAL]-Nano->Height*.5)); */
    dr[NORMAL] = 100.;
    /* dr[NORMAL] -= Nano->Height*.5; */
    /* Pos[NORMAL] = fabs((dr[NORMAL]-Nano->Height*.5)/(1.5*Nano->Radius)); */
    /* dr[NORMAL] = (dr[NORMAL] - (Nano->Height*.5+Nano->Radius)); */
    /* for(int d=0;d<3;d++) dr[d]  *= .65; */
    /* if(dr[NORMAL] > 0.) dr[NORMAL] = 100.; */
  }
  else if(dr[NORMAL] < -Nano->Height*.5){
    /* Pos[NORMAL] = fabs((dr[NORMAL]+Nano->Height*.5)); */
    dr[NORMAL] = 100.;
    /* dr[NORMAL] += Nano->Height*.5; */
    /* Pos[NORMAL] = fabs((dr[NORMAL]+Nano->Height*.5)/(1.5*Nano->Radius)); */
    /* dr[NORMAL] = (dr[NORMAL] + (Nano->Height*.5+Nano->Radius)); */
    /* for(int d=0;d<3;d++) dr[d]  *= .65; */
    /* if(dr[NORMAL] < 0.) dr[NORMAL] = 100.; */
  }
  else{
    Pos[NORMAL] = 0.;
    dr[NORMAL]  = 0.;
  }
  if(IfTens){
    /*Position of the surface point intersected by the normal 
      connecting the axis with the monomer position */
    double HeiOnAxis = ProjectOn(Nano->Axis,PosRel);
    double AxPos[4] = {0.,0.,0.,0.};
    for(int d=0;d<3;d++){
      AxPos[d] = -Nano->Axis[d]*HeiOnAxis + dr[d];
      AxPos[3] += SQR(AxPos[d]);
    }
    AxPos[3] = sqrt(AxPos[3]);
    for(int d=0;d<3;d++){
      NPos[d] = Nano->Pos[d] + Nano->Axis[d]*HeiOnAxis + Nano->Radius*AxPos[d]/AxPos[3];
    }
  }
  double dr2 = SQR(dr[TANG1]) + SQR(dr[TANG2]) + SQR(dr[NORMAL]);
  return dr2;
}
/** Distance form a sphere */
static double RigidSphericalDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  double dr2 = 0.;
  double PosBf[3];
  for(int d=0;d<3;d++){
    PosBf[d] = Pos[d] - floor( Pos[d]/Edge[d] )*Edge[d];
    dr[d] = NPos[d] - PosBf[d];
    if(dr[d] >  .5*Edge[d]) dr[d] -= Edge[d];
    if(dr[d] < -.5*Edge[d]) dr[d] += Edge[d];
  }
  if(IfTens){
    /* position on the sphere of the conjunction between the nanoparticle centere and the bead */
    double Norm = SQR(Nano->Radius)/(SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]));
    Norm = sqrt(Norm);
    for(int d=0;d<3;d++){
      NPos[d] = NPos[d] + dr[d]*Norm;
    }
  }
  dr2 = SQR(dr[TANG1]) + SQR(dr[TANG2]) + SQR(dr[NORMAL]);
  return dr2;
}
/** Distance form a elipsoid */
static double RigidElipsDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  double dr2 = 0.;
  double Ratio = Nano->Height > 0. ? Nano->Radius/Nano->Height : 1.;
  for(int d=0;d<3;d++){
    double PosBf = Pos[d];
    dr[d] = NPos[d] - PosBf;
    if(d==NORMAL) continue;
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
  }
  dr[NORMAL] = Ratio*dr[NORMAL];
  dr2 = SQR(dr[TANG1]) + SQR(dr[TANG2]) + SQR(dr[NORMAL]);
  return dr2;
}
/** Distance form a parabola */
static double RigidParabDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  for(int d=0;d<3;d++){
    double PosBf = Pos[d] - floor( Pos[d]/Edge[d] )*Edge[d];
    dr[d] = NPos[d] - PosBf;
    if(d==2) continue;
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
  }
  double r2 = SQR(dr[TANG1])+SQR(dr[TANG2]);
  return SQR(-4.*Nano->Height*dr[NORMAL]+r2+3.*Nano->Height+4.*Nano->Height*Nano->Pos[NORMAL]);
  double r = sqrt(r2);
  double z = Pos[NORMAL];//Nano[n].Pos[NORMAL] - Pos[NORMAL];
  double a = -4.*Nano->Height;
  double b = 1.;
  double c = -2.*0.;
  double d = 3.*Nano->Height+4.*Nano->Height*Nano->Pos[NORMAL];
  Pos[NORMAL] = 0.;
  return SQR(a*z + b*r2 + c*r + d);
}
/** Distance from a planar wall */
static double RigidTiltWallDist(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  for(int d=0;d<3;d++){
    dr[d] = NPos[d] - Pos[d];
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
  }
  double a = Nano->Axis[0];
  double b = Nano->Axis[1];
  double c = Nano->Axis[2];
  double d = -Nano->Axis[0]*NPos[0] - Nano->Axis[1]*NPos[1] - Nano->Axis[2]*NPos[2];
  double Dist2 = a*dr[0]+b*dr[1]+c*dr[2]+d;
  Dist2 = SQR(Dist2)/(SQR(a)+SQR(b)+SQR(c));
  return Dist2;
}
/** Distance from a planar wall */
static double RigidWallDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  double dr2=0.;
  dr[NORMAL] = NPos[NORMAL] - Pos[NORMAL];
  dr[NORMAL] -= floor( dr[NORMAL]/(Edge[NORMAL]) + .5)*Edge[NORMAL];
  dr2 = SQR(dr[NORMAL]);
  NPos[TANG1] = Pos[TANG1];
  NPos[TANG2] = Pos[TANG2];
  return dr2;
}
/** Distance from cylindrical external potential */
static double RigidPoreDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  for(int d=0;d<3;d++){
    dr[d] = Pos[d] - NPos[d];
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
  }
  Pos[NORMAL] = dr[NORMAL]/Nano->Height;
  dr[NORMAL] = 0.;
  return  SQR(dr[TANG1]) + SQR(dr[TANG2]);
}
static double RigidTorusDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  for(int d=0;d<3;d++){
    dr[d] = Pos[d] - NPos[d];
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
  }
  double Radxy = sqrt(SQR(dr[TANG1]) + SQR(dr[TANG2]));
  double Temp = SQR(Nano->Height - Radxy) - SQR(Nano->Radius) + SQR(dr[NORMAL]);
  return SQR(Temp);
}
static double RigidJanusDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  double PosRel[3] = {0.,0.,0.};
  for(int d=0;d<3;d++){
    dr[d] = NPos[d] - Pos[d];
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
    PosRel[d] = dr[d];
  }
  PosRel[TANG2] = 2.*dr[TANG2];
  // writes on dr the distance from the axis
  double r2 = Dist2AxisPerp(PosRel,Nano->Axis,dr);
  Pos[NORMAL] = dr[NORMAL] - 4.*Nano->Radius;
  return r2;
}
/** Calculate the distance to the surface of a tilted cylinder */
static double RigidTiltDistance(double *Pos,double *NPos,double *Edge,double *dr,struct NANO *Nano){
  double PosRel[3] = {0.,0.,0.};
  for(int d=0;d<3;d++){
    dr[d] = NPos[d] - Pos[d];
    dr[d] -= floor( dr[d]/(Edge[d]) + .5)*Edge[d];
    PosRel[d] = dr[d];
  }
  // writes on dr the distance from the axis
  double r2 = Dist2AxisPerp(PosRel,Nano->Axis,dr);
  /* update the temporal angolar momentum (to be multiplicated
     by the inertial tensor */
  VectProd(dr,PosRel,Nano->AMomTemp);
  /* The cylinder ends with a cupola, if the particle is above the 
     height of the cylinder, the cylinder is seen as
     two spheres in \pm Nano->Height*.5 */
  // if and how much above the cylinder (for the smoothing)
  Pos[NORMAL] = 0.;
  double HeiOnAxis = ProjectOn(Nano->Axis,PosRel);
  if(fabs(HeiOnAxis) > Nano->Height*.5){
    double Sign = HeiOnAxis > 0. ? 1. : -1.;
    r2 = 0.;
    for(int d=0;d<3;d++){
      dr[d] = - (Sign*(Nano->Height*.5+Nano->Radius)*Nano->Axis[d]) + PosRel[d];
      dr[d]  *= .65;
    }
    for(int d=0;d<3;d++){
      r2 += SQR(dr[d]);
    }
   Pos[NORMAL] = fabs(HeiOnAxis - Sign*.5*Nano->Height);
  }
  /* if( fabs(HeiOnAxis) > .5*Nano->Height){ */
  /*   double Sign = HeiOnAxis > 0. ? 1. : -1.; */
  /*   r2 = 0.; */
  /*   for(int d=0;d<3;d++){ */
  /*     dr[d] = - (Sign*Nano->Height*.5*Nano->Axis[d]) + PosRel[d]; */
  /*     r2 += SQR(dr[d]); */
  /*   } */
  /*   //Smoothing if Pos[NORMAL] > 0. */
  /*   Pos[NORMAL] = fabs(HeiOnAxis - Sign*.5*Nano->Height); */
  /* } */
  return r2;
}
//############################FORCES###############################
static double WeightFuncSpline(const double x){
  if(fabs(x) < 1.)
    return 2.*fabs(x)*x*x - 3.*x*x + 1.;
  else
    return 0.;
}
static double DerWeightFuncSpline(const double x){
  if(fabs(x) < 1.)
    return 6.*x*x - 6.*fabs(x);
  else
    return 0.;
}
static double WeightFuncSquare(const double x){
  if(fabs(x) < 1.)
    return x*x - 2.*fabs(x) + 1.;
  else
    return 0.;
}
static double DerWeightFuncSquare(const double x){
  if(fabs(x) < 1.)
    return 2.*x - 2.;
  else
    return 0.;
}
static double WeightFuncCubic(const double x){
  if(ABS(x) < 1.)
    return x*(-SQR(x)+1.);
  else
    return 0.;
}
static double DerWeightFuncCubic(const double x){
  if(ABS(x) < 1.)
    return -3.*x*x+1;
  else
    return 0.;
}
static double SignFuncJanus(const double x,int type){
  double c = 1.5;
  if(type > 0) return 1.;
  if(x <= -1.)
    return -1.;
  else if(x < 1.)
    //return (1.-c)*CUBE(2.*x-1.) + c*(2.*x-1.);
    return (1.-c)*CUBE(x) + c*(x);
  return 1.;
}
static double SignFuncCyl(const double x,int type){
  double c = 1.5;
  if(type > 0) return 1.;
  if(x <= 0.)
    return -1.;
  else if(x < 1.)
    return (1.-c)*CUBE(2.*x-1.) + c*(2.*x-1.);
  return 1.;
}
/** Definition of a standard Lennard-Jones 39 Potential */
static double RigidForceOld(struct NANO *Nano,double Dr,double *Pot,double Sign){
  double Cons = 0.;
  double OffSet = 1.;
  if( Dr <  Nano->Radius - 0.33*OffSet ){
    Pot[0] = 33.*Nano->Hamaker - 478.*(Dr - Nano->Radius + 0.33*OffSet)*Nano->Hamaker;
    Pot[1] = 478.*Nano->Hamaker;
    Cons = 478.*Nano->Hamaker;
    return Cons;
  }
  double idr = Nano->Radius/Dr;
  double idr2 = SQR(idr);
  double idr4 = SQR(idr2);
  double idr6 = CUBE(idr2);
  double Hamaker = Nano->Hamaker;
  Pot[0] = idr*idr2*Hamaker*(idr6    + Sign*1.) + Nano->Baseline;
  Cons   = idr4    *Hamaker*(9.*idr6 + Sign*3.)/Nano->Radius;
  Pot[1] = idr4    *Hamaker*(9.*idr6 + Sign*3.);
  return Cons;
}
static double RigidForceCyl(struct NANO *Nano,double Dr,double *Pot,double Sign){
  double idr = fabs(1./(Dr-Nano->Radius+1.));
  if(SQR(Dr) < SQR(.6*Nano->Radius)){//LJ cut off 
    Pot[0] = -100. / idr + 1000.*(1.+.6*Nano->Radius) + 0.0369862318104964;
    Pot[1] = 1000.*+.6*Nano->Radius;
    return 100.;
  }
  double idr2 = SQR(idr);
  double idr4 = SQR(idr2);
  double idr6 = CUBE(idr2);
  Pot[0] = Nano->Hamaker * idr * idr2 * (idr6 + Sign*1.) + 0.0369862318104964;
  Pot[1] = Nano->Hamaker*3/(Nano->Radius)*idr*idr2*(3*idr6 + Sign*1.);
  return Nano->Hamaker * idr4 * (9. * idr4 * idr2 + Sign*3.);
}
/** Linear attractive/repulsive potential */
static double RigidForceLinear(struct NANO *Nano,double Dr,double *Pot,double Sign){
  if( Dr < Nano->Radius)
    return 1.*(Nano->Hamaker*(- Dr + Nano->Radius)/(Nano->Radius) - Nano->Hamaker*.6666);
  else 
    return Nano->Hamaker*(  Dr - 3.*Nano->Radius)/(3.*Nano->Radius);
}
/** Hamker Physica IV 10 p 1058 (1937) */
static double RigidForceHamaker(struct NANO *Nano,double Dr,double *Pot,double Sign){
  double SurfTens = 3.*Nano->Hamaker;
  double Cons = 0.;
  double Sigma  = 1.0;//0.04;
  double Slope = 1.00203;
  double Intercept = 0.31739;
  double Rnp = Nano->Radius/Slope-Intercept;
  if( Dr <  Nano->DistThr){
    Pot[0] = Nano->PotThr - Nano->ForThr*(Dr - Nano->DistThr);
    Pot[1] = Nano->PotThr;
    return Nano->ForThr;
  }
  double DrInv  = 1./Dr;
  double DrP1 = 1./(Dr+Rnp);
  double DrP3 = CUBE(DrP1);
  double DrP6 = SQR(DrP3);
  double DrP9 = DrP6*DrP3;
  double DrM1 = 1./(-Dr+Rnp);
  double DrM3 = CUBE(DrM1);
  double DrM6 = SQR(DrM3);
  double DrM9 = DrM6*DrM3;
  double PreRep = Sigma*1./360.*DrInv*SurfTens;
  double PreAttr = 1./12.*DrInv*SurfTens;
  double Rep = (9.0*Rnp+Dr)*DrP9 - (9.0*Rnp-Dr)*DrM9;
  double Attr= (3.0*Rnp+Dr)*DrP3 - (3.0*Rnp-Dr)*DrM3;
  double RepP  = -DrP9+9.*(9.*Rnp+Dr)*DrP9*DrP1
    -(9.*Rnp+Dr)*DrP9*DrInv;
  double RepM  = -DrM9+9.*(9.*Rnp-Dr)*DrM9*DrM1
    -(9.*Rnp-Dr)*DrM9*DrInv;
  double AttrP = -DrP3+3.*(3.*Rnp+Dr)*DrP3*DrP1
    -(3.*Rnp+Dr)*DrP3*DrInv;
  double AttrM = -DrM3+3.*(3.*Rnp-Dr)*DrM3*DrM1
    -(3.*Rnp-Dr)*DrM3*DrInv;
  double RepChem = 9.*DrP9 - (9.*Rnp+Dr)*9.*DrP9*DrP1 
    - 9.*DrM9 + 9.*(9.*Rnp-Dr)*DrM9*DrM1;
  double AttrChem = 3.*DrP3 - (3.*Rnp+Dr)*3.*DrP3*DrP1 
    - 3.*DrM9 + 3.*(3.*Rnp-Dr)*DrM3*DrM1;
  Pot[0] = PreRep*Rep + Sign*PreAttr*Attr + Nano->Baseline;
  Pot[1] = PreRep*RepChem + Sign*PreAttr*AttrChem;
  Cons = PreRep*(RepP+RepM) + Sign*PreAttr*(AttrP+AttrM);
  return Cons;
}
/** 39 Lennard-Jones potential translated by the radius of the inclusion and a fixed coating ruled by a */
static double RigidForceCutOff(struct NANO *Nano,double Dr,double *Pot,double Sign){
  double Cons = 0.;
  double Strength = pow(Nano->Coating,9.)*Nano->Hamaker;
  double Hamaker  = pow(Nano->Coating,3.)*Nano->Hamaker;
  if( Dr <  Nano->DistThr){
    Pot[0] = Nano->PotThr - Nano->ForThr*(Dr - Nano->DistThr);
    Pot[1] = Nano->ForThr;
    Cons = Nano->ForThr;
    return Cons;
  }
  double idr = 1./(Dr-Nano->Radius+Nano->Coating);
  double idr2 = idr*idr;
  double idr4 = idr2*idr2;
  double idr6 = idr2*idr4;
  Pot[0] = idr*idr2*(Strength*idr6   +Sign*Hamaker)+Nano->Baseline;
  Cons   = idr4    *(Strength*9.*idr6+Sign*3.*Hamaker);
  Pot[1] = idr4    *(Strength*9.*idr6+Sign*3.*Hamaker);
  return Cons;
}
// Harmonic potential
static double RigidForceHarmonic(struct NANO *Nano,double Dr,double *Pot,double Sign){
  if(Dr > Nano->Radius){Pot[0] = 0.; return 0.;}
  if( Dr <  Nano->DistThr ){
    Pot[0] = Nano->PotThr - Nano->ForThr*(Dr - Nano->DistThr);
    Pot[1] = Nano->ForThr;
    return Nano->ForThr;
  }
  double a = -Nano->Hamaker/SQR(Nano->Radius);
  double b = 0.;
  double c = Nano->Hamaker;
  double Cons = - 2.*a*Dr - b;
  Pot[0] = a*Dr*Dr + b*Dr + c;
  return Sign*Cons;
}
static double RigidForceLJ39(struct NANO *Nano,double Dr,double *Pot,double Sign){
  double Cons = 0.;
  double NanoRadInv = 1./Nano->Radius;
  if( Dr <  Nano->DistThr ){
    Pot[0] = Nano->PotThr - Nano->ForThr*(Dr - Nano->DistThr);
    Pot[1] = Nano->ForThr;
    return Nano->ForThr;
  }
  double idr = Nano->Radius/(Dr);
  double idr2 = idr*idr;
  double idr4 = idr2*idr2;
  double idr6 = idr2*idr4;
  Pot[0] = Nano->Hamaker*idr*idr2*(idr6+Sign*1.)+Nano->Baseline;
  Cons   = Nano->Hamaker*idr4*(9.*idr6 +Sign*3.)*NanoRadInv;
  Pot[1] = Nano->Hamaker*idr4*(9.*idr6 +Sign*3.)*NanoRadInv;
  return Cons;
}
//##########################VELOCITY#VERLET#######################
/** 
    First step of the velocity Verlet integration
*/
/** 
    First step of the velocity Verlet integration
*/
void RigidVv1Umbr(struct beads *restrict b){
  double dt = b->dt;
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    if(n==0) continue;
    else if(n==1){
      Nano->Force[1] -= 50.*(Nano->Pos[1] - Nano->OldPos[1]);
      double tmp = .5*Nano->Force[1]*dt;
      Nano->Pos[1] += (Nano->Vel[1]+tmp)*dt;
      Nano->Vel[1] += tmp;
    }
  }
}
void RigidVv1(struct beads *restrict b){
  if(b->NNano == 0 )return;
  for(int n=0;n<b->NNano;n++){
    struct NANO *Nano = b->Nano + n;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_NONE)){
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_SPH)){
      for(int d=0;d<3;d++){
	double tmp = .5*Nano->Force[d]*b->dt/Nano->Mass;
	Nano->Vel[d] += tmp;
	Nano->Pos[d] += (Nano->Vel[d])*b->dt;
      }
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_TIP)){
      Nano->Pos[NORMAL] -= Nano->Coating;
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_CYL)){
      for(int d=0;d<3;d++){
      /* for(int d=0;d<2;d++){ */
      	double tmp = .5*Nano->Force[d]*b->dt/Nano->Mass;
      	Nano->Vel[d] += tmp;
      	Nano->Pos[d] += (Nano->Vel[d])*b->dt;
      }
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_PILL)){
      for(int d=0;d<3;d++){
	double tmp = .5*Nano->Force[d]*b->dt/Nano->Mass;
	Nano->Vel[d] += tmp;
	Nano->Pos[d] += (Nano->Vel[d])*b->dt;
      }
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_TILT)){
      double PosTop[3]; // Position of the top of the cylinder
      double VelTop[3]; // Velocity at the top
      for(int d=0;d<3;d++){
	PosTop[d] = Nano->Axis[d]*Nano->Height*.5;
	double Ran = Nano->Zeta  * (2.*rng_uniform(rng) - 1.);
	double Dis = -Nano->Gamma*Nano->AVel[d];
	Nano->AMom[d] += Ran + Dis;
      }
      //---------Equation-of-motion--------------------
      for(int d=0;d<3;d++){
	double tmp = .5*Nano->Force[d]*b->dt/Nano->Mass;
	Nano->Vel[d] += tmp;
	Nano->Pos[d] += (Nano->Vel[d])*b->dt;
      }
      //---------Angular---------------
      UpdateAVel(Nano,b->dt);
      VectProd(Nano->AVel,PosTop,VelTop);
      for(int d=0;d<3;d++){
	PosTop[d] += VelTop[d]*b->dt;
	Nano->Axis[d] = PosTop[d];
      }
      Normalize(Nano->Axis);
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_WALL)){
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_JANUS)){
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_UMBR)){
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_TORUS)){
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_EXT)){
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_PORE)){
      //PorePos(b);
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_STALK)){
    }
  }
}
void RigidVv2(struct beads *restrict b){
  if(b->NNano == 0 ) return ;
  double VolumeInv = 1./(b->l[0]*b->l[1]*b->l[2]);
  for(int n=0;n<b->NNano;n++) {
    struct NANO *Nano = b->Nano + n;
    for(int d=0;d<3;d++){
      Nano->Vel[d] += .5*Nano->Force[d]*b->dt/Nano->Mass;
      Nano->Pre[d] += SQR(Nano->Vel[d])*VolumeInv;
    }
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_TILT))
      UpdateAVel(Nano,b->dt);
    Nano->v2 = SQR(Nano->Vel[0]) + SQR(Nano->Vel[1]) + SQR(Nano->Vel[2]);
  }
  struct TENS_PROF *Tens = &(b->Tens);
  struct NANO *Nano = b->Nano;
  for(int d=0;d<3;d++){
    Tens->RefPos[d] = Nano->Pos[d];
    Tens->RefAxis[d] = Nano->Axis[d];
  }
}
static void PorePos(struct beads *restrict b){
  const int NGrid = 20;
  double *Plot = calloc(NGrid*NGrid,sizeof(double));
  int first = col_index * b->groupsize; 
  int last = (col_index + 1) * b->groupsize;
  double NPart = (double)(last-first);
  double Cm[3] = {0.,0.,0.};
  //FILE *Ciccia = fopen("Ciccia.dat","w");
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      Plot[gx*NGrid+gy] = 1.;
    }
  }
  for(int p=first;p<last;p++){
    int gx = (int)(b->xv[p][TANG1]/b->l[TANG1]*NGrid);
    int gy = (int)(b->xv[p][TANG2]/b->l[TANG2]*NGrid);
    Plot[gx*NGrid+gy] += 1000.;
    Cm[NORMAL] += b->xv[p][NORMAL];
  }
  double Dx = .5*b->l[TANG1]/(double)NGrid;
  double Dy = .5*b->l[TANG2]/(double)NGrid;
  double Count = 0.;
  for(int gx=0;gx<NGrid;gx++){
    double x = gx/(double)NGrid*b->l[TANG1] + Dx;
    for(int gy=0;gy<NGrid;gy++){
      double y = gy/(double)NGrid*b->l[TANG2] + Dy;
      //Plot[gx*NGrid+gy] /= 100.;
      //fprintf(Ciccia,"%lf %lf %lf\n",x,y,10./Plot[gx*NGrid+gy]);
      Cm[TANG1] += x/(Plot[gx*NGrid+gy]);
      Cm[TANG2] += y/(Plot[gx*NGrid+gy]);
      Count     +=1./(Plot[gx*NGrid+gy]);
    }
  }
  Cm[0] /= Count;
  Cm[1] /= Count;
  Cm[2] /= NPart;
  struct NANO *Nano = b->Nano;
  for(int d=0;d<3;d++){
    if(Cm[d] < b->l[d])
      Nano->Pos[d] = Cm[d];
  }
  //fprintf(Ciccia,"%l %lf %lf\n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2]);
  //printf("%lf %lf %lf\n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2]);  
  free(Plot);
  //fclose(Ciccia);exit(0);
}
static void StalkPos(struct beads *restrict b){
  int first = col_index * b->groupsize;
  int last = (col_index + 1) * b->groupsize;
  double NPart = (double)(last-first);
  double Cm[3] = {0.,0.,0.};
  const int NGrid = 20;
  double *Plot = calloc(NGrid*NGrid,sizeof(double));
  //FILE *Ciccia = fopen("Ciccia.dat","w");
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      Plot[gx*NGrid+gy] = 1.;
    }
  }
  for(int p=first;p<last;p++){
    int gx = (int)(b->xv[p][TANG1]/b->l[TANG1]*NGrid);
    int gy = (int)(b->xv[p][TANG2]/b->l[TANG2]*NGrid);
    Plot[gx*NGrid+gy] += 1000.;
    Cm[NORMAL] += b->xv[p][NORMAL];
  }
  double Dx = .5*b->l[TANG1]/(double)NGrid;
  double Dy = .5*b->l[TANG2]/(double)NGrid;
  double Count = 0.;
  for(int gx=0;gx<NGrid;gx++){
    double x = gx/(double)NGrid*b->l[TANG1] + Dx;
    for(int gy=0;gy<NGrid;gy++){
      double y = gy/(double)NGrid*b->l[TANG2] + Dy;
      //Plot[gx*NGrid+gy] /= 100.;
      //fprintf(Ciccia,"%lf %lf %lf\n",x,y,10./Plot[gx*NGrid+gy]);
      Cm[TANG1] += x*(Plot[gx*NGrid+gy]);
      Cm[TANG2] += y*(Plot[gx*NGrid+gy]);
      Count     +=   (Plot[gx*NGrid+gy]);
    }
  }
  Cm[0] /= Count;
  Cm[1] /= Count;
  Cm[2] /= NPart;
  struct NANO *Nano = b->Nano;
  for(int d=0;d<3;d++){
    if(Cm[d] < b->l[d])
      Nano->Pos[d] = Cm[d];
  }
  //fclose(Ciccia);exit(0);
  //fprintf(Ciccia,"%lf %lf %lf\n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2]);
  //printf("%lf %lf %lf\n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2]);  
  free(Plot);
}
static void UpdateAVel(struct NANO *Nano,double dt){
  //-----------Trasformation-of-the-inertia-tensor---
  double InTensor[9]; //Inverse inertia tensor
  double Rot[9]; // Rotation matrix
  double RotT[9]; // Transpose of the rotation matrix
  double Resp1[9]; // Response
  double Resp2[9]; // Response
  double AVel[3]; // Angular velocity
  double Axis[3]; // New rotation axis
  double Normal[3]; // Normal to the membrane
  Normal[TANG1] = 0.; 
  Normal[TANG2] = 0.;
  Normal[NORMAL] = 1.;
  //VectProd(Nano->Axis,Normal,Axis);
  VectProd(Normal,Nano->Axis,Axis);
  double Angle = VectAngle(Nano->Axis,Normal);
  QUAT Quat = Quaternion(Axis,Angle);
  //RotMatrix(Rot,Quat);
  RotMatr9(Rot,Axis,Angle);
  CopyMatrix(RotT,Rot);
  Transpose(RotT);
  for(int i=0;i<9;i++)
    InTensor[i] = 0.;
  InTensor[0] = 1./(Nano->Mass*(.25*SQR(Nano->Radius)+1./12.*SQR(Nano->Height) ) );
  InTensor[4] = InTensor[0];
  InTensor[8] = 1./(.5*Nano->Mass*SQR(Nano->Radius));
  // Calculating r^t I^-1 r L = \dot w
  MatrMatrix(RotT,InTensor,Resp1);
  MatrMatrix(Resp1,Rot,Resp2);
  MatrVect(Resp2,Nano->AMom,AVel);
  /* PrintMatrix(Rot); */
  /* PrintMatrix(InTensor); */
  /* PrintMatrix(Resp2); */
  /* printf("%lf %lf %lf\n",AVel[0],AVel[1],AVel[2]); */
  for(int d=0;d<3;d++){
    Nano->AVel[d] += AVel[d]*.5*dt;
  }
}
/* //####################MISCELLANEOUS########################## */
/* static double *VectProd(double *u,double *v){ */
/*   double Res[3]; */
/*   for(int d=0;d<33;d++){ */
/*     int NUno = (d+1)%3; */
/*     int NDue = (d+2)%3; */
/*     Resp[d] = u[NUno] * v[NDue] - u[NDue]*v[NUno]; */
/*   } */
/*   return Resp; */
/* } */

double NormVect(double *Axis){
  double Norma = 0.;
  for(int d=0;d<3;d++){
    Norma += SQR(Axis[d]);
  }
  if(Norma > 0. ) return sqrt(Norma);
  return 1.;
}
/** Normalize a vector */
void Normalize(double *Axis){
  double NormaInv = 1./NormVect(Axis);
  for(int d=0;d<3;d++)
    Axis[d] *= NormaInv;
}
/// Vectorial product
void VectProd(double *v,double *u,double *Resp){
  Resp[0] = v[1]*u[2] - v[2]*u[1];
  Resp[1] = v[2]*u[0] - v[0]*u[2];
  Resp[2] = v[0]*u[1] - v[1]*u[0];
}
/// Define a quaternion as a rotation of an angle Angle around a Axis
double VectAngle(double *u,double *v){
  double Resp = 0.;
  for(int d=0;d<3;d++)
    Resp += u[d]*v[d];
  Resp /= NormVect(u)*NormVect(v);
  return acos(Resp);
}
QUAT Quaternion(double *Axis,double Angle){
  QUAT q;
  q.w = cos(Angle*.5);
  double NormaInv = 1./NormVect(Axis);
  double Sin = sin(Angle*.5);
  q.x = Sin*Axis[0]*NormaInv;
  q.y = Sin*Axis[1]*NormaInv;
  q.z = Sin*Axis[2]*NormaInv;
  return q;
}
/** Matrix vector multiplication */
void MatrVect(double *M,double *v,double *Resp){
  for(int c=0;c<3;c++){
    double Temp = 0.;
    for(int r=0;r<3;r++)
      Temp += M[3*c + r] * v[r];
    Resp[c] = Temp;
  }
}
/** Matrix Matrix multiplication */
void MatrMatrix(double *M1,double *M2,double *Resp){
  for(int c=0;c<3;c++){
    for(int r=0;r<3;r++){
      double Temp = 0.;
      for(int i=0;i<3;i++)
	Temp += M1[3*i+r] * M2[3*c+i];
      Resp[3*c+r] = Temp;
    }
  }
}
/** Copy a matrix */
void CopyMatrix(double *Resp,const double *M2){
  for(int c=0;c<3;c++)
    for(int r=0;r<3;r++)
      Resp[3*c+r] = M2[3*c+r];
}
/** print a matrix */
void PrintMatrix(double *M1){
  for(int c=0;c<3;c++){
    printf("|");
    for(int r=0;r<3;r++){
      printf("%lf ",M1[3*c+r]);
    }
    printf("|\n");
  }
}
void Swap(double *a,double *b){
  double Temp = *a;
  *a = *b;
  *b = Temp;
}
/** Transpone */
void Transpose(double *M1){
  for(int c=0;c<3;c++)
    for(int r=c+1;r<3;r++){
      double Temp = M1[3*c+r];
      M1[3*c+r] = M1[3*r+c];
      M1[3*r+c] = Temp;
    }
}
/** Create a rotation matrix from an axis and angle */
void RotMatr9(double *M,double *Axis,double Angle){
  double c = cos(Angle);
  double s = sin(Angle);
  double t = 1. - c;
  double Norm = 0.;
  for(int d=0;d<3;d++)
    Norm += SQR(Axis[d]);
  Norm = Norm > 0. ? sqrt(Norm) : 1.;
  double x = Axis[0]/Norm;
  double y = Axis[1]/Norm;
  double z = Axis[2]/Norm;
  int NRow = 3;
  M[NRow*0+0]  = t*x*x + c;
  M[NRow*0+1]  = t*x*y + z*s;
  M[NRow*0+2]  = t*x*z - y*s;
  
  M[NRow*1+0]  = t*x*y - z*s;
  M[NRow*1+1]  = t*y*y + c;
  M[NRow*1+2]  = t*y*z + x*s;
  
  M[NRow*2+0]  = t*x*z + y*s;
  M[NRow*2+1]  = t*y*z - x*s;
  M[NRow*2+2]  = t*z*z + c;
}
void RotMatr16(double *M,double *Axis,double Angle){
  double c = cos(Angle);
  double s = sin(Angle);
  double t = 1. - c;
  double Norm = 0.;
  for(int d=0;d<3;d++)
    Norm += SQR(Axis[d]);
  Norm = Norm > 0. ? sqrt(Norm) : 1.;
  double x = Axis[0]/Norm;
  double y = Axis[1]/Norm;
  double z = Axis[2]/Norm;
  int NRow = 4;
  M[NRow*0+0]  = t*x*x + c;
  M[NRow*0+1]  = t*x*y + z*s;
  M[NRow*0+2]  = t*x*z - y*s;
  M[NRow*0+3]  = 0.;
    
  M[NRow*1+0]  = t*x*y - z*s;
  M[NRow*1+1]  = t*y*y + c;
  M[NRow*1+2]  = t*y*z + x*s;
  M[NRow*1+3]  = 0.;
    
  M[NRow*2+0]  = t*x*z + y*s;
  M[NRow*2+1]  = t*y*z - x*s;
  M[NRow*2+2]  = t*z*z + c;
  M[NRow*2+3]  = 0.;
    
  M[NRow*3+0]  = 0.;
  M[NRow*3+1]  = 0.;
  M[NRow*3+2]  = 0.;
  M[NRow*3+3]  = 1.;
}
/** Create a rotation matrix from a quaternion */
void RotMatrix(double *M,QUAT q){
  // first row
  M[0] = SQR(q.w) + SQR(q.x) - SQR(q.y) - SQR(q.z);
  M[1] = 2.*q.x*q.y + 2.*q.w*q.z;
  M[2] = 2.*q.x*q.z - 2.*q.w*q.y;
  //
  M[3] = 2.*q.x*q.y - 2.*q.w*q.z;
  M[4] = SQR(q.w) - SQR(q.x) + SQR(q.y) - SQR(q.z);
  M[5] = 2.*q.z*q.y + 2.*q.w*q.x;
  //
  M[6] = 2.*q.x*q.z + 2.*q.w*q.y;
  M[7] = 2.*q.y*q.z - 2.*q.w*q.x;
  M[8] = SQR(q.w) - SQR(q.x) - SQR(q.y) + SQR(q.z);
}
void RotMatrix16(double *M,QUAT q){
  // first row
  M[0] = 1.-2.*q.x*q.x-2.*q.z*q.z;
  M[1] = 2.*q.x*q.y + 2.*q.w*q.z;
  M[2] = 2.*q.x*q.z - 2.*q.w*q.y;
  M[3] = 0.;
  //
  M[4] = 2.*q.x*q.y - 2.*q.w*q.z;
  M[5] = 1.-2.*q.x*q.x-2.*q.z*q.z;
  M[6] = 2.*q.x*q.y + 2.*q.w*q.x;
  M[7] = 0.;
  //
  M[8] = 2.*q.x*q.z + 2.*q.w*q.y;
  M[9] = 2.*q.y*q.z - 2.*q.w*q.x;
  M[10] = 1.-2.*q.x*q.x-2.*q.y*q.y;
  M[11] = 0.;
  //
  M[12] = 0.;
  M[13] = 0.;
  M[14] = 0.;
  M[15] = 1.;
}
/** Calculate the distance perpendicular to the axis of the nanoparticle using sin th = |a^b|/(|a||b|), Axis should be already normalized */
double Dist2Axis(double *Pos,double *Axis){
  double Perp1[3];
  double Norm1 = 0.;
  //double NormA = 0.;
  VectProd(Pos,Axis,Perp1);
  for(int d=0;d<3;d++){
    Norm1 += SQR(Perp1[d]);
    //NormA += SQR(Axis[d]);
  }
  return Norm1;///(NormA);
}
/** Calculate the distance and the vector componentes perpendicular to the axis of the nanoparticle using sin th = |a^b|/(|a||b|) , Axis should be already normalized */
double Dist2AxisPerp(double *Pos,double *Axis, double *Perp){
  double Perp1[3];
  //double NormA = 0.;
  double Norm1 = 0.;
  //double Norm2 = 0.;
  double Dist = 0.;
  VectProd(Pos,Axis,Perp1);
  VectProd(Axis,Perp1,Perp);
  for(int d=0;d<3;d++){
    Norm1 += SQR(Perp1[d]);
    //NormPos += SQR(Pos[d]);
    //NormA += SQR(Axis[d]);
    //Norm2 += SQR(Perp[d]);
  }
  Dist = Norm1;// /(NormA);
  //double Fact = sqrt( Norm1/(NormPos*Norm2) );// /sqrt(NormA);
  /* for(int d=0;d<3;d++){ */
  /*   Perp[d] = Perp[d]; */
  /* } */
  return Dist;
}
/** Calculate the components of a vector projected on an axis , Axis should be already normalized */
double ProjectOn(double *Axis,double *Pos){
  double Scalar = 0.;
  //double NormA = 0.;
  for(int d=0;d<3;d++){
    Scalar += Axis[d]*Pos[d];
    //NormA += Axis[d];
  }
  return Scalar; // /NormA;
}
/** Calculate the components of a vector projected on an axis , Axis should be already normalized */
double ProjectOnResp(double *Axis,double *Pos,double *Resp){
  double Scalar = 0.;
  //double NormA = 0.;
  double NormP = 0.;
  for(int d=0;d<3;d++){
    Scalar += Axis[d]*Pos[d];
    //NormA += Axis[d];
    NormP += SQR(Pos[d]);
  }
  double Cos = Scalar/sqrt(NormP);// /sqrt(NormA);
  for(int d=0;d<3;d++){
    Resp[d] = Cos*Axis[d];// /NormA
  }
  return Scalar;
}
double AngleVect(double *u,double *v){
  double Resp = 0.;
  for(int d=0;d<3;d++){
    Resp += u[d]*v[d];
  }
  Resp /= NormVect(u);
  Resp /= NormVect(v);
  //half circle is missing
  return cos(Resp);
}
