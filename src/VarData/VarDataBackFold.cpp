/***********************************************************************
VarData: This  Program reads and writes a specific file format 
storing all the information relative to a set of equal structure
polymers to the CHAIN, PART and GENERAL structures. It provides 
two different ways to backfold the coordinates, a fanction that 
creates an initial system with different option and some function
for the data analisys. The first calculate the distribution of the
monomer in the box, the second the distribution of the bonds.
Copyright (C) 2008 by Giovanni Marelli <sabeiro@virgilio.it>


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
***********************************************************************/
#include "../include/VarData.h"
;
int VarData::BfDefChain(){
  if(VAR_IF_TYPE(SysType,VAR_CHAIN_DEF)){
    for(int c=0,b=0;c<Gen->NChain;c++){
      for(int d=0;d<3;d++){
	Ch[c].Pos[d] -= floor(Ch[c].Pos[d]*pInvEdge(d))*Gen->Edge[d];
      }
    }
    return 0;
  }
  VarMessage("BackFold.DefChains");
  Vettore Ax0(1.,0.,0.);
  Vettore Ax2(0.,0.,1.);
  for(int b=0,NCh=0;b<pNBlock();NCh+=Block[b++].NChain){
    double *xc = (double *)calloc(pNPCh(b),sizeof(double));
    double *yc = (double *)calloc(pNPCh(b),sizeof(double));
    double *zc = (double *)calloc(pNPCh(b),sizeof(double));
    //if(strcasestr(Block[b].Name, "PEP") == Block[b].Name) continue;
    int ChType = CHAIN_POLY;
    if(Block[b].Asym == 0) ChType = CHAIN_ADDED;
    //printf("\n%s type %d asym %d ppc %d #chain %d offset %d %d\n",Block[b].Name,ChType,Block[b].Asym,pNPCh(b),Block[b].NChain,NCh,Block[b].NChain+NCh);
    for(int c=NCh;c<NCh+Block[b].NChain;c++){
      Ch[c].Type = 0;
      if(Block[b].Asym == 0){
	VAR_ADD_TYPE(Ch[c].Type,CHAIN_ADDED);
      }
      for(int d=0;d<3;d++){
	Ch[c].Pos[d] = 0.;
	Ch[c].Vel[d] = 0.;
      }
      int p1 = Block[b].InitIdx + (c-NCh)*pNPCh(b);
      for(int p=p1;p<p1+pNPCh(b);p++){
	for(int d=0;d<3;d++){
	  Ch[c].Pos[d] += Pm[p].Pos[d];
	  Ch[c].Vel[d] += Pm[p].Vel[d];
	}
	}
      for(int d=0;d<3;d++){
	Ch[c].Pos[d] /= (double)pNPCh(b);
	Ch[c].Vel[d] /= (double)pNPCh(b);
      }
      for(int ppc=0;ppc<pNPCh(b);ppc++){
	xc[ppc] = Pm[p1+ppc].Pos[0];
	yc[ppc] = Pm[p1+ppc].Pos[1];
	zc[ppc] = Pm[p1+ppc].Pos[2];
      }
      RETTA rzx = Mat->InterRett(zc,xc,pNPCh(b));
      RETTA rzy = Mat->InterRett(zc,yc,pNPCh(b));
      double x1 = Pm[p1].Pos[2]*rzx.m + rzx.q;
      double x2 = Pm[p1+pNPCh(b)-1].Pos[2]*rzx.m + rzx.q;
      double y1 = Pm[p1].Pos[2]*rzy.m + rzy.q;
      double y2 = Pm[p1+pNPCh(b)-1].Pos[2]*rzy.m + rzy.q;
      Ch[c].Dir[2] = Pm[p1].Pos[2] - Pm[p1+pNPCh(b)-1].Pos[2];
      Ch[c].Dir[1] = y1 - y2;
      Ch[c].Dir[0] = x1 - x2;
      Vettore ChDir(Ch[c].Dir[0],Ch[c].Dir[1],Ch[c].Dir[2]);
      Ch[c].Angle = Ax0.Angle(&Ax2,&ChDir);
      if(Ch[c].Angle > .5*M_PI)
	VAR_ADD_TYPE(Ch[c].Type,CHAIN_UP);
      else
	VAR_ADD_TYPE(Ch[c].Type,CHAIN_DOWN);
      Vettore Origin(pCm(0)-Ch[c].Pos[0],pCm(1)-Ch[c].Pos[1],pCm(2)-Ch[c].Pos[2]);
      double Angle = Origin.Angle(&ChDir);
      if(Angle < .5*M_PI)
	VAR_ADD_TYPE(Ch[c].Type,CHAIN_OUTER);
      else 
	VAR_ADD_TYPE(Ch[c].Type,CHAIN_INNER);
    }
    free(xc);
    free(yc);
    free(zc);
  }
  VAR_ADD_TYPE(SysType,VAR_CHAIN_DEF);
  return 0;
}
void VarData::BfPep(){
  int nBlock = 0;
  for(int n=0;n<pNNano();n++){
    if(Nano[n].Shape != SHAPE_CLUSTER) continue;
    int nChain = 0;
    for(int b=0;b<Gen->NBlock;b++){
      if(!strncmp(Block[b].Name,"PEP",3)){
	if(b > nBlock){
	  nBlock = b;
	  break;
	}
      }
      nChain += Block[b].NChain;
    }
    int p1 = Block[nBlock].InitIdx;
    double Cm[3] = {0.,0.,0.};
    double NCm = 0.;
    double *xc = (double *)calloc(pNPCh(nBlock),sizeof(double));
    double *yc = (double *)calloc(pNPCh(nBlock),sizeof(double));
    double *zc = (double *)calloc(pNPCh(nBlock),sizeof(double));
    for(int ppc=0;ppc<pNPCh(nBlock);ppc++){
      xc[ppc] = Pm[p1+ppc].Pos[0];
      yc[ppc] = Pm[p1+ppc].Pos[1];
      zc[ppc] = Pm[p1+ppc].Pos[2];
    }
    RETTA rzx = Mat->InterRett(zc,xc,pNPCh(nBlock));
    RETTA rzy = Mat->InterRett(zc,yc,pNPCh(nBlock));
    double x1 = Pm[p1].Pos[2]*rzx.m + rzx.q;
    double x2 = Pm[p1+pNPCh(nBlock)-1].Pos[2]*rzx.m + rzx.q;
    double y1 = Pm[p1].Pos[2]*rzy.m + rzy.q;
    double y2 = Pm[p1+pNPCh(nBlock)-1].Pos[2]*rzy.m + rzy.q;
    Nano[n].Axis[2] = Pm[p1].Pos[2] - Pm[p1+pNPCh(nBlock)-1].Pos[2];
    Nano[n].Axis[1] = y1 - y2;
    Nano[n].Axis[0] = x1 - x2;      
    for(int p=Block[nBlock].InitIdx;p<Block[nBlock].EndIdx;p++){
      for(int d=0;d<3;d++) Cm[d] += Pm[p].Pos[d];
      NCm += 1.;
    }
    for(int d=0;d<3;d++) Nano[n].Pos[d] = Cm[d]/NCm;
    SetNanoBkf(n);
    double Norm = 0.;
    for(int d=0;d<3;d++){
      //SigErr(isnan(Nano[n].Axis[d]),"Wrong nano axis %d %lf ",d,Nano[n].Axis[d]);
      Norm += SQR(Nano[n].Axis[d]);
    }
    Norm = sqrt(Norm);
    for(int d=0;d<3;d++){
      Nano[n].Axis[d] /= -Norm;
    }
    free(xc);
    free(yc);
    free(zc);
  }
}
int VarData::BfEdge(){
  double Edge[3][2];
  if( VAR_IF_TYPE(SysType,VAR_EDGE) ) return 0;
  Edge[0][0] = Pm[0].Pos[0];Edge[0][1] = Pm[0].Pos[0];
  Edge[1][0] = Pm[0].Pos[1];Edge[1][1] = Pm[0].Pos[1];
  Edge[2][0] = Pm[0].Pos[2];Edge[2][1] = Pm[0].Pos[2];
  for(int p=0;p<Gen->NPart;p++){
    for(int d=0;d<3;d++){
      if(Edge[d][1] < Pm[p].Pos[d] )
	Edge[d][1] = Pm[p].Pos[d];
      if(Edge[d][0] >  Pm[p].Pos[d])
	Edge[d][0] = Pm[p].Pos[d];
    }
  }
  for(int d=0;d<3;d++)
    SetEdge(Edge[d][1] - Edge[d][0],d);
    // SetEdge(Edge[d][1],d);
  VAR_ADD_TYPE(SysType,VAR_EDGE);
  return 0;
}
void VarData::ShiftRef(int BackFold){
  VarMessage("BackFold");
  if(BackFold == BF_SKIP) return ;
  BfEdge();
  BfPep();
  double Shift[3];
  for(int d=0;d<3;d++){
    Shift[d] = -ShiftPos[d]*Gen->Edge[d];
  }
  for(int c=0;c<Gen->NChain;c++){
    for(int d=0;d<3;d++){
      Ch[c].Pos[d] -= Shift[d];
      //Ch[c].Pos[d] -= floor(Ch[c].Pos[d]/Gen->Edge[d])*Gen->Edge[d];
    }
  }
  if(BackFold == BF_CHAIN){
    BfDefChain();
    for(int p=0;p<Gen->NPart;p++){
      for(int d=0;d<3;d++){
	Pm[p].Pos[d] -= Shift[d];
	Pm[p].Bkf[d] = -floor((Ch[Pm[p].CId].Pos[d]-Shift[d])*pInvEdge(d))*Gen->Edge[d];
      }
    }
    for(int c=0;c<Gen->NChain;c++){
      for(int d=0;d<3;d++){
	Ch[c].Pos[d] -= floor(Ch[c].Pos[d]*pInvEdge(d))*pEdge(d);
      }
    }
  }
  else if(BackFold == BF_PART){
    for(int p=0;p<Gen->NPart;p++){
      for(int d=0;d<3;d++){
	Pm[p].Pos[d] -= Shift[d];
	Pm[p].Bkf[d] = -floor(Pm[p].Pos[d]*pInvEdge(d))*Gen->Edge[d];
      }
    }
  }
  else if(BackFold == BF_NANO){
    for(int d=0;d<3;d++){
      ShiftPos[d] = pNanoPos(0,d)*pInvEdge(d);
      Shift[d] = (ShiftPos[d]-.5)*Gen->Edge[d];
    }
    for(int p=0;p<Gen->NPart;p++){
      for(int d=0;d<3;d++){
  	Pm[p].Pos[d] -= Shift[d];
  	Pm[p].Bkf[d] = -floor(Pm[p].Pos[d]*pInvEdge(d))*Gen->Edge[d];
      }
    }
  }
  else if(BackFold == BF_TILT){
    int n=0;
    for(int d=0;d<3;d++){
      ShiftPos[d] = pNanoPos(n,d);
      Shift[d] = ShiftPos[d]-.5*Gen->Edge[d];
    }
    Vettore NanoAx(Nano[n].Axis[0],Nano[n].Axis[1],Nano[n].Axis[2]);
    NanoAx.Set(0.,CNorm);
    NanoAx.Normalize();
    Vettore Ex(0.,0.,0.);
    Vettore Zed(0.,0.,0.);
    Zed.Set(1.,CNorm);
    Ex.Set(1.,CLat1);
    //double Angle = atan(Nano[0].Axis[CLat2]/Nano[0].Axis[CLat1]);//Zed.Angle(&NanoAx,&Ex);
    double Angle = Ex.Angle(&NanoAx);
    if(isnan(Angle)) Angle = 0.;
    int NRow = 4;
    Matrice M(Zed.x,.5*Angle,NRow);
    // printf("angle %lf\n",Angle*360./DUE_PI);
    // M.Print();
    Vettore PosOld(NRow);
    Vettore PosNew(NRow);
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
	PosOld.Set(Pm[p].Pos[d]-pNanoPos(0,d),d);
      }
      for(int r=0;r<NRow;r++){
	double Temp=0.;
	for(int c=0;c<NRow;c++)
	  Temp += M.Val(r,c)*PosOld.Val(c);
	PosNew.Set(Temp,r);
      }
      //PosNew = M.Mult(PosOld);
      for(int d=0;d<3;d++){
	Pm[p].Pos[d] = PosNew.Val(d) +.5*pEdge(d);
	Pm[p].Bkf[d] = -floor(Pm[p].Pos[d]*pInvEdge(d))*Gen->Edge[d];
      }
    }
    Nano[0].Axis[CLat2] = -sqrt(1.-SQR(Nano[0].Axis[CNorm]));
    Nano[0].Axis[CLat1] = 0.;
  }
  for(int n=0;n<MAX(1,Gen->NNano);n++){
    for(int d=0;d<3;d++){
      Nano[n].Pos[d] -= Shift[d];
    }
    SetNanoBkf(n);
  }
}
//obsolete?
void VarData::DefBlock(int *NChStep,int How){
  if(How == VAR_OPPOSED){
    Gen->NBlock = 4;
    Block = (BLOCK *)realloc(Block,Gen->NBlock*sizeof(*Block));
    sprintf(Block[0].Name,"Lower1");
    sprintf(Block[1].Name,"Upper1");
    sprintf(Block[2].Name,"Lower2");
    sprintf(Block[3].Name,"Upper2");
  }
  else if(How == VAR_VESICLE || How == VAR_TUBE){
    Gen->NBlock = 2;
    Block = (BLOCK *)realloc(Block,Gen->NBlock*sizeof(*Block));
    sprintf(Block[0].Name,"Inner");
    sprintf(Block[1].Name,"Outer");
  }
  Block[0].InitIdx= 0;
  Block[0].EndIdx = NChStep[0]*Gen->NPCh;
  Block[0].NChain = NChStep[0];
  Block[0].NPCh = Gen->NPCh;
  for(int b=1;b<Gen->NBlock;b++){
    Block[b].InitIdx= Block[b-1].EndIdx;
    Block[b].EndIdx = NChStep[b]*Gen->NPCh+Block[b-1].EndIdx;
    Block[b].NChain = NChStep[b];
    Block[b].NPCh = Gen->NPCh;
  }
  for(int b=0;b<Gen->NBlock;b++)
    printf("%d %d %d %d\n",Block[b].NChain,Block[b].InitIdx,Block[b].EndIdx,Block[b].NChain);
}
//obsolete?
bool VarData::BackFold(int How){
  VarMessage("BackFold");
  if(How == BF_SKIP) return 0;
  BfEdge();
  double ShiftNano[3];
  SetNanoBkf(0);  
  if(How == BF_PART)
    for(int p=0;p<Gen->NPart;p++)
      for(int d=0;d<3;d++)
	Pm[p].Bkf[d] = -floor(Pm[p].Pos[d]/Gen->Edge[d])*Gen->Edge[d];
  else if(How == BF_CHAIN)
    for(int p=0;p<Gen->NPart;p++){
      int c = Pm[p].CId;
      if(c < 0 || c >= Gen->NChain) continue;
      for(int d=0;d<3;d++)
	Pm[p].Bkf[d] = -floor(Ch[c].Pos[d]/Gen->Edge[d])*Gen->Edge[d];
    }
  else if(How == BF_NANO)
    for(int p=0;p<Gen->NPart;p++)
      for(int d=0;d<3;d++)
	Pm[p].Pos[d] += ShiftNano[d];
  BfDefChain();
  return 0;
}
//Obsolete
bool VarData::ShiftSys(int How){
  if(How == SHIFT_NO) return 0;
  double Shift[3] = {0.,0.,0.};
  SetNanoBkf(0);
  if(How == SHIFT_CM){//Cm
    for(int d=0;d<3;d++)
      Gen->Cm[d] = 0.;
    for(int p=0;p<Gen->NPart;p++){
      Gen->Cm[0] += Pm[p].Pos[0];
      Gen->Cm[1] += Pm[p].Pos[1];
      Gen->Cm[2] += Pm[p].Pos[2];
    }
    for(int d=0;d<3;d++){
      Gen->Cm[d] /= (double)Gen->NPart;
      Shift[d] = (Gen->Edge[d]*.5-Gen->Cm[d]);
      Nano->Pos[d] += Shift[d]; 
    }
    SetNanoBkf(0);
  }
  else if(How == SHIFT_NANO){//Nano Particle
    for(int d=0;d<3;d++){
      Shift[d] = (Gen->Edge[d]*.5 - pNanoPos(0,d));
      Nano->Pos[d] += Shift[d]; 
    }
    SetNanoBkf(0);
  }
  else if(How == SHIFT_CM_NANO){//x,y Nano; z Cm
    for(int d=0;d<3;d++)
      Gen->Cm[d] = 0.;
    for(int p=0;p<Gen->NPart;p++){
      Gen->Cm[2] += Pm[p].Pos[2];
    }
    Gen->Cm[0] = pNanoPos(0,0);
    Gen->Cm[1] = pNanoPos(0,1);
    Gen->Cm[2] /= (double)Gen->NPart;
    for(int d=0;d<3;d++){
      Shift[d] = (Gen->Edge[d]*.5 - Gen->Cm[d]);
      if(d < 2) Nano->Pos[d] += Shift[d];
    }
    SetNanoBkf(0);
  }
  else return 1;
  for(int p=0;p<Gen->NPart;p++){
    // if(Pm[p].Typ == 2){
    //   Pm[p].Pos[0] = Nano->PosBf[0];
    //   Pm[p].Pos[1] = Nano->PosBf[1];
    //   Pm[p].Pos[2] = Nano->PosBf[2];
    //   continue;
    // }
    Pm[p].Pos[0] += Shift[0];
    Pm[p].Pos[1] += Shift[1];
    Pm[p].Pos[2] += Shift[2];
  }
  //printf("NanoPos %lf %lf %lf %lf %lf %lf\n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],Nano->PosBf[0],Nano->PosBf[1],Nano->PosBf[2]);
  return 0;
}
void VarData::BackBone(double *Line,int NBin){
  double *Count = new double[NBin];
  double *LineS = new double[NBin];
  for(int b=0;b<NBin;b++){
    Line[b] = 0.;
    Count[b] = 0.;
  }
  CLat1 = 1;
  CLat2 = 0;
  //average
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Pos[2] < .02*pEdge(2))continue;
    int b = (int)(Pm[p].Pos[CLat1]*pInvEdge(CLat1)*NBin);
    if(b<0 || b >= NBin)continue;
    double Weight = Pm[p].Pos[2];
    Line[b] += Pm[p].Pos[CLat2]*Weight;
    Count[b] += Weight;
  }
  for(int b=0;b<NBin;b++){
    Line[b] /= Count[b] > 0. ? Count[b] : 1.;
  }
  //smooth
  for(int v=0;v<NBin;v++) LineS[v] = Line[v];
  InterBSpline1D(LineS,Line,NBin,NBin);
  for(int v=0;v<NBin;v++) LineS[v] = Line[v];
  InterBSpline1D(LineS,Line,NBin,NBin);
  // //reweighting
  // for(int b=0;b<NBin;b++){
  //   Line[b] = 0.;
  //   Count[b] = 0.;
  // }
  // double NormDist = 1./(.5*pEdge(CLat2));
  // for(int p=0;p<pNPart();p++){
  //   if(Pm[p].Pos[2] < .02*pEdge(2))continue;
  //   int b = (Pm[p].Pos[CLat1]*pInvEdge(CLat1)*NBin);
  //   if(b < 0 || b >= NBin)continue;
  //   double Dist = SQR(Pm[p].Pos[CLat2] - LineS[b]);
  //   if(Dist > SQR(.1)) continue;
  //   int bm1 = b-3;
  //   if(bm1 < 0) bm1 = b;
  //   double Distm1 = SQR(Pm[p].Pos[CLat2] - LineS[bm1]);
  //   if(Distm1 > SQR(.3)) continue;
  //   int bp1 = b+3;
  //   if(bp1 > NBin-1) bp1 = b;
  //   double Distp1 = SQR(Pm[p].Pos[CLat2] - LineS[bp1]);
  //   if(Distp1 > SQR(.3)) continue;
  //   double Weight = 100.*(1. - Dist*Distm1*Distp1*CUBE(NormDist));
  //   Line[b] += Pm[p].Pos[CLat2]*Weight;
  //   Count[b] += Weight;
  // }
  // for(int b=0;b<NBin;b++){
  //   Line[b] /= Count[b] > 0. ? Count[b] : 1.;
  // }
  // //smooth
  // for(int v=0;v<NBin;v++) LineS[v] = Line[v];
  // InterBSpline1D(LineS,Line,NBin,NBin);
  // for(int v=0;v<NBin;v++) LineS[v] = Line[v];
  // InterBSpline1D(LineS,Line,NBin,NBin);
  //fourier
#ifdef USE_FFTW
  fftw_complex *FouOut = (fftw_complex *)fftw_malloc(NBin*sizeof(fftw_complex));
  fftw_complex *FouIn  = (fftw_complex *)fftw_malloc(NBin*sizeof(fftw_complex));
  fftw_plan direct = fftw_plan_dft_1d(NBin,
  				     FouIn,FouOut,FFTW_FORWARD,FFTW_PATIENT);
  fftw_plan reverse = fftw_plan_dft_1d(NBin,
  				     FouOut,FouIn,FFTW_BACKWARD,FFTW_PATIENT);
  for(int b=0;b<NBin;b++) FouIn[b][0] = Line[b];
  fftw_execute(direct);
  int NComp = NBin;// - 2;//NBin/2;
  for(int b=NBin-1;b>=NComp;b--){
    FouOut[b][0] = 0.;
    FouOut[b][1] = 0.;
  }
  fftw_execute(reverse);
  for(int b=0;b<NBin;b++) Line[b] = FouIn[b][0]/(double)NBin;  
#endif
  //write
  if(1==0){
    FILE *FOut = fopen("BackBone.dat","w");
    fprintf(FOut,"#l(%lf %lf %lf) v[%d] d[part]\n",pEdge(CLat1),pEdge(CLat2),pEdge(CNorm),NBin);
    int NPart = 0;
    for(int p=0;p<pNPart();p++){
      if(Pm[p].Pos[2] < .02*pEdge(2)) continue;
      fprintf(FOut,"{t[0 0 0] x(%lf %lf %lf)}\n",Pm[p].Pos[CLat1],Pm[p].Pos[CLat2],Pm[p].Pos[CNorm]);
      NPart++;
    }
    for(int b=0;b<NBin;b++){
      double x = (b+.5) * pEdge(CLat1)/(double)NBin;
      double y = b * pEdge(CLat2)/(double)NBin;
      fprintf(FOut,"{t[0 0 1] x(%lf %lf 0.) l[%d]} \n",x,Line[b],NPart+1);
      NPart++;
    }
    fclose(FOut);
  }
  FILE *FOut1 = fopen("LineProf.dat","w");
  for(int b=0;b<NBin;b++){
    double x = (b+.5) * pEdge(CLat1)/(double)NBin;
    double y = b * pEdge(CLat2)/(double)NBin;
    fprintf(FOut1,"%lf %lf\n",x,Line[b]);
  }
  fclose(FOut1);  
}
void VarData::StalkLineProf(double *Line,int NBin){
  //-------------------alloc
  Vettore Ax0(1.,0.,0.);
  Vettore Ax2(0.,0.,1.);
  double Dist[4];
  BfDefChain();
  //double CmStalk[3] = {0.,0.,0.};
  double *CountStalk =  (double *)calloc(NBin,sizeof(double));
  int NSample = 8;
  int NIter = 6;
  Matrice Mask(5);
  Mask.FillGaussian(.5,3.);
  int NDim = 1;
  int IfMinImConv = 1;
  int IfSpline = 1;
  int IfGaussian = 1;
  double InvNSample = 1./(double)NSample;
  double **PlotUp    = (double **)calloc(NSample,sizeof(double));
  double **PlotUpS   = (double **)calloc(NSample,sizeof(double));
  double **CountUp   = (double **)calloc(NSample,sizeof(double));
  double **PlotDown  = (double **)calloc(NSample,sizeof(double));
  double **PlotDownS = (double **)calloc(NSample,sizeof(double));
  double **CountDown = (double **)calloc(NSample,sizeof(double));
  double *LineS = (double *)calloc(NBin,sizeof(double));
  for(int s=0;s<NSample;s++){
    PlotUp[s]    = (double *)calloc(NSample,sizeof(double));
    PlotUpS[s]   = (double *)calloc(NSample,sizeof(double));
    CountUp[s]   = (double *)calloc(NSample,sizeof(double));
    PlotDown[s]  = (double *)calloc(NSample,sizeof(double));
    PlotDownS[s] = (double *)calloc(NSample,sizeof(double));
    CountDown[s] = (double *)calloc(NSample,sizeof(double));
  }
  for(int vx=0;vx<NBin;vx++){
    Line[vx] = 0.;
  }
  //---------------defines the two midplanes
  //center of mass of the upper layer of the lower membrane
  for(int b=0,cOff=0;b<pNBlock();b++,cOff+=Block[b].NChain){
    if(strcasestr(Block[b].Name, "TT0") != Block[b].Name) continue;
    for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
      if(!VAR_IF_TYPE(Ch[Pm[p].CId].Type,CHAIN_UP))continue;
      if(Pm[p].Typ == 0) continue;
      double Posx = Pm[p].Pos[CLat1] - floor(Pm[p].Pos[CLat1]*pInvEdge(CLat1))*pEdge(CLat1);
      double Posy = Pm[p].Pos[CLat2] - floor(Pm[p].Pos[CLat2]*pInvEdge(CLat2))*pEdge(CLat2);
      int sx = (int)(Posx*pInvEdge(CLat1)*NSample);
      int sy = (int)(Posy*pInvEdge(CLat2)*NSample);
      PlotDown[sx][sy] += Pm[p].Pos[CNorm];
      CountDown[sx][sy] += 1.;
    }
  }
  //center of mass of the lower layer of the upper membrane
  for(int b=0,cOff=0;b<pNBlock();b++,cOff+=Block[b].NChain){
    if(strcasestr(Block[b].Name, "TT1") != Block[b].Name) continue;
    for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
      if(!VAR_IF_TYPE(Ch[Pm[p].CId].Type,CHAIN_DOWN))continue;
      if(Pm[p].Typ == 0) continue;
      //if(Pm[p].Pos[2] < .5*pEdge(2))continue;
      double Posx = Pm[p].Pos[CLat1] - floor(Pm[p].Pos[CLat1]*pInvEdge(CLat1))*pEdge(CLat1);
      double Posy = Pm[p].Pos[CLat2] - floor(Pm[p].Pos[CLat2]*pInvEdge(CLat2))*pEdge(CLat2);
      int sx = (int)(Posx*pInvEdge(CLat1)*NSample);
      int sy = (int)(Posy*pInvEdge(CLat2)*NSample);
      PlotUp[sx][sy] += Pm[p].Pos[CNorm];
      CountUp[sx][sy] += 1.;
    }
  }
  for(int sx=0;sx<NSample;sx++){
    for(int sy=0;sy<NSample;sy++){
      PlotUp[sx][sy] /= CountUp[sx][sy] > 0. ? CountUp[sx][sy] : 1.;
      PlotDown[sx][sy] /= CountDown[sx][sy] > 0. ? CountDown[sx][sy] : 1.;
    }
  }
  //InterBSpline2D(PlotUp,PlotUpS,NSample,NSample);
  //InterBSpline2D(PlotDown,PlotDownS,NSample,NSample);
  for(int sx=0;sx<NSample;sx++){
    for(int sy=0;sy<NSample;sy++){
      PlotUpS[sx][sy] = PlotUp[sx][sy];
      PlotDownS[sx][sy] = PlotDown[sx][sy];
    }
  }
  //----------------count the chains between the two midplanes
  double AngleMax = 1./(.5*M_PI);
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ == 1) continue;
    int sx = (int)(Pm[p].Pos[CLat1]*pInvEdge(CLat1)*NSample);
    if(sx < 0 || sx >= NSample) continue;
    int sy = (int)(Pm[p].Pos[CLat2]*pInvEdge(CLat2)*NSample);
    if(sy < 0 || sy >= NSample) continue;
    if(Pm[p].Pos[CNorm] > PlotUpS[sx][sy] || Pm[p].Pos[CNorm] < PlotDownS[sx][sy]) continue;
    int vx = (int)(Pm[p].Pos[CLat2]*pInvEdge(CLat2)*NBin);
    if(vx < 0 || vx >= NBin) continue;
    int c = Pm[p].CId;
    Vettore ChDir(Ch[c].Dir[0],Ch[c].Dir[1],Ch[c].Dir[2]);
    double Angle = Ax0.Angle(&Ax2,&ChDir);
    if(Angle < .7 || Angle > 2.3) continue;
    double Weight = 1.;//1. - SQR((Angle-.5*M_PI)*AngleMax);
    Line[vx] += Pm[p].Pos[CLat1]*Weight;
    CountStalk[vx] += Weight;
  }
  for(int vx=0;vx<NBin;vx++){
    Line[vx] /= CountStalk[vx] > 0. ? CountStalk[vx] : 1.;
  }
  //-------------------------smoothing
  if(IfSpline){
    for(int v=0;v<NBin;v++) LineS[v] = Line[v];
    InterBSpline1D(LineS,Line,NBin,NBin);
    for(int v=0;v<NBin;v++) LineS[v] = Line[v];
    InterBSpline1D(LineS,Line,NBin,NBin);
  }
  if(IfGaussian){
    Mask.ConvoluteMatrix(Line,NBin,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Line,NBin,NDim,IfMinImConv);
  }
  // // //------------------------Squared-average------------
  for(int i=0;i<NIter;i++){
    for(int v=0;v<NBin;v++){
      LineS[v] = Line[v];
      Line[v] = 0.;
      CountStalk[v] = 0.;
    }
    double NormDist = 1./(.5*pEdge(CLat2));
    for(int p=0;p<pNPart();p++){
      if(Pm[p].Typ == 1) continue;
      int sx = (int)(Pm[p].Pos[CLat1]*pInvEdge(CLat1)*NSample);
      if(sx < 0 || sx >= NSample) continue;
      int sy = (int)(Pm[p].Pos[CLat2]*pInvEdge(CLat2)*NSample);
      if(sy < 0 || sy >= NSample) continue;
      if(Pm[p].Pos[CNorm] > PlotUpS[sx][sy] || Pm[p].Pos[CNorm] < PlotDownS[sx][sy]) continue;
      int vx = (int)(Pm[p].Pos[CLat2]*pInvEdge(CLat2)*NBin);
      if(vx < 0 || vx >= NBin) continue;
      double Dist = SQR(Pm[p].Pos[CLat1] - LineS[vx]);
      if(Dist > SQR(3.)) continue;
      double Weight = 100.*(1. - SQR((Pm[p].Pos[CLat1]-LineS[vx])*NormDist));
      Line[vx] += Pm[p].Pos[CLat1]*Weight;
      CountStalk[vx] += Weight;
    }
    double Average = 0.;
    for(int vx=0;vx<NBin;vx++){
      Line[vx] /= CountStalk[vx] > 0. ? CountStalk[vx] : 1.;
      Average += Line[vx];
    }
    // //-------------------------smoothing
    if(IfSpline){
      for(int v=0;v<NBin;v++) LineS[v] = Line[v];
      InterBSpline1D(LineS,Line,NBin,NBin);
      for(int v=0;v<NBin;v++) LineS[v] = Line[v];
      InterBSpline1D(LineS,Line,NBin,NBin);
    }
    if(IfGaussian){
      Mask.ConvoluteMatrix(Line,NBin,NDim,IfMinImConv);
      Mask.ConvoluteMatrix(Line,NBin,NDim,IfMinImConv);
    }
  }
  //-------------------------print-temp-files
  if(1==0){
    FILE *F2Write = fopen("TwoMid.dat","w");
    FILE *FWrite = fopen("StalkLine.dat","w");
    FILE *FWrite1 = fopen("StalkLine1.dat","w");
    fprintf(F2Write,"#l(%lf %lf %lf) v[%d] d[part]\n",pEdge(0),pEdge(1),pEdge(2),NSample);
    WriteSurf(F2Write,PlotUpS,NSample,0);
    WriteSurf(F2Write,PlotDownS,NSample,SQR(NSample));
    for(int p=0;p<pNPart();p++){
      if(Pm[p].Typ == 1) continue;
      int sx = (int)(Pm[p].Pos[CLat1]*pInvEdge(CLat1)*NSample);
      if(sx < 0 || sx >= NSample) continue;
      int sy = (int)(Pm[p].Pos[CLat2]*pInvEdge(CLat2)*NSample);
      if(sy < 0 || sy >= NSample) continue;
      if(Pm[p].Pos[CNorm] > PlotUpS[sx][sy] || Pm[p].Pos[CNorm] < PlotDownS[sx][sy]) continue;
      int c = Pm[p].CId;
      Vettore ChDir(Ch[c].Dir[0],Ch[c].Dir[1],Ch[c].Dir[2]);
      double Angle = Ax0.Angle(&Ax2,&ChDir);
      if(Angle < .7 || Angle > 2.3) continue;
      fprintf(F2Write,"{t[0 0 1] x(%lf %lf %lf)}\n",Pm[p].Pos[0],Pm[p].Pos[1],Pm[p].Pos[2]);
      fprintf(FWrite,"%lf %lf \n",Pm[p].Pos[1],Pm[p].Pos[0]);
    }
    for(int vx=0;vx<NBin;vx++){
      fprintf(FWrite1,"%lf %lf \n",vx*pEdge(CLat2)/(double)NBin,Line[vx]);
    }
    fclose(FWrite);
    fclose(F2Write);
    //exit(0);
  }
  for(int s=0;s<NSample;s++){
    free(PlotUp[s]);
    free(PlotUpS[s]);
    free(CountUp[s]);
    free(PlotDown[s]);
    free(PlotDownS[s]);
    free(CountDown[s]);
  }
  free(PlotUp);
  free(PlotUpS);
  free(CountUp);
  free(CountDown);
  free(PlotDown);
  free(PlotDownS);
  free(CountStalk);
}
void VarData::StalkPos2(double *OldPos,double *CmStalk){
  Vettore Ax0(1.,0.,0.);
  Vettore Ax2(0.,0.,1.);
  double Dist[4];
  BfDefChain();
  double CountStalk = 1.;
  for(int b=0,cOff=0;b<pNBlock();b++,cOff+=Block[b].NChain){
    if(strncmp(Block[b].Name,"TT",2)) continue;
    for(int c=cOff;c<cOff+Block[b].NChain;c++){
      Vettore ChDir(Ch[c].Dir[0],Ch[c].Dir[1],Ch[c].Dir[2]);
      double Angle = Ax0.Angle(&Ax2,&ChDir);
      if(Ch[c].Pos[2] < OldPos[2] - 1.5 || Ch[c].Pos[2] > OldPos[2] + 1.5) continue;
      if(Angle < 1. || Angle > 2.) continue;
      for(int p=0;p<Block[b].NPCh;p++){
      	int pCurr = Block[b].InitIdx + (c-cOff)*Block[b].NPCh + p;
      	if(Pm[p].Typ == 1) continue;
  	for(int d=0;d<3;d++){
  	  Dist[d] = Pm[p].Pos[d] - OldPos[d];
  	  Dist[d] -= floor(Dist[d]*pInvEdge(d))*pEdge(d);
  	}
  	Dist[3] = (SQR(Dist[0])+SQR(Dist[1])+SQR(Dist[2]));
  	double Weight = .00001/(Dist[3]*SQR(Angle-.5*M_PI));
  	//double Weight = .1/(Dist[3]);
  	if(Weight > 1.) continue;
  	//printf("%lf %lf %lf\n",Angle,Dist[3],Weight);
  	for(int d=0;d<3;d++){
  	  CmStalk[d] += Pm[p].Pos[d]*Weight;//(Pm[p].Pos[d]-OldPos[d])*Weight;
  	}
  	CountStalk += Weight;
      }
    }
  }
  if(CountStalk <= 0.) CountStalk = 1.;
  for(int d=0;d<3;d++){
    CmStalk[d] /= CountStalk;
  }
}
void VarData::StalkPos3(double *OldPos,double *CmStalk){
  int NSample = 32;
  int NGrid = NSample-1;
  double VolEl = pVol()/(double)CUB(NSample);
  double *Plot3d = (double *)calloc(CUBE(NSample),sizeof(double));
  double *Plot2d = (double *)calloc(SQR(NGrid),sizeof(double));
  double Min = 0.;
  double Max = 0.;
  double CountStalk = 0.;
  VAR_TRIANGLE *Triang = NULL;
  double Dist[4];
  for(int d=0;d<3;d++) CmStalk[d] = 0.;
  for(int p=0;p<pNPart();p++){
    //if(Pm[p].Typ != 0) continue;
    double Posx = Pm[p].Pos[0] - floor(Pm[p].Pos[0]*pInvEdge(0))*pEdge(0);
    double Posy = Pm[p].Pos[1] - floor(Pm[p].Pos[1]*pInvEdge(1))*pEdge(1);
    double Posz = Pm[p].Pos[2] - floor(Pm[p].Pos[2]*pInvEdge(2))*pEdge(2);
    int sx = (int)(Posx*pInvEdge(0)*NSample);
    int sy = (int)(Posy*pInvEdge(1)*NSample);
    int sz = (int)(Posz*pInvEdge(2)*NSample);
    int sTot = (sx*NSample+sy)*NSample+sz;
    Plot3d[sTot] += VolEl;
    CmStalk[CNorm] += Pm[p].Pos[CNorm];
    if(Max < Plot3d[sTot]) Max = Plot3d[sTot];
    if(Min > Plot3d[sTot]) Min = Plot3d[sTot];
  }
  double IsoLevel = .1*Max;
  int NTri = 0;
  Triang = MarchingCubes(Plot3d,NSample,IsoLevel,&NTri);  
  free(Plot3d);
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      Plot2d[gx*NGrid+gy] = 1.;
    }
  }
  for(int t=0;t<NTri;t++){
    for(int v=0;v<3;v++){
      for(int d=0;d<3;d++){
	Dist[d] = OldPos[d] - Triang[t].p[v].x[d];
      }
      Dist[3] = SQR(Dist[0])+SQR(Dist[1])+SQR(Dist[2]);
      int gx = (int)(Triang[t].p[v].x[CLat1]*pInvEdge(CLat1)*NGrid);
      int gy = (int)(Triang[t].p[v].x[CLat2]*pInvEdge(CLat2)*NGrid);
      if(gx < 0 || gx >= NGrid) continue;
      if(gy < 0 || gy >= NGrid) continue;
      Plot2d[gx*NGrid+gy] += sqrt(Dist[3]);
    }
  }
  double Dx = .5*pEdge(CLat1)/(double)NGrid;
  double Dy = .5*pEdge(CLat2)/(double)NGrid;
  double Count = 0.;
  if(1==0){
    FILE *Ciccia = fopen("Ciccia.dat","w");
    for(int gx=0;gx<NGrid;gx++){
      double x = gx/(double)NGrid*pEdge(CLat1) + Dx;
      for(int gy=0;gy<NGrid;gy++){
	double y = gy/(double)NGrid*pEdge(CLat1) + Dy;
	double Weight = 1./(Plot2d[gx*NGrid+gy]);
	fprintf(Ciccia,"%lf %lf %lf\n",x,y,Weight*10000.);
      }
    }
    fclose(Ciccia);
    exit(0);
  }
  for(int gx=0;gx<NGrid;gx++){
    double x = gx/(double)NGrid*pEdge(CLat1) + Dx;
    for(int gy=0;gy<NGrid;gy++){
      double y = gy/(double)NGrid*pEdge(CLat1) + Dy;
      double Weight = 1./(Plot2d[gx*NGrid+gy]*SQR(x-OldPos[0])*SQR(y-OldPos[1]));
      CmStalk[CLat1] += x*Weight;
      CmStalk[CLat2] += y*Weight;
      CountStalk     +=   Weight;
      //printf("  %lf %lf %lf\n",CmStalk[CLat1],CmStalk[CLat2],Weight);
    }
  }
  CmStalk[0] /= CountStalk;
  CmStalk[1] /= CountStalk;
  CmStalk[2] /= pNPart();
  free(Plot2d);
}
int VarData::StalkPos4(double *OldPos,double *CmStalk){
  Nano->Shape = SHAPE_TORUS;
  Point2Shape(Nano->Shape);
  int NInt = 30;
  int NBin = 12;
  double kElPhob = 40000.;
  double kElPhil = -1.;//-2.;
  double MinRad = .2;//.4;
  double MinHei = .1;//2.;
  double RadMax = .3;
  double RadMin = .3;//.4;
  double HeiMax = 6.;
  double HeiMin = 1.0;
  double GainRad = 100.;
  double GainHei = 2000.;
  double MoveStep = .05;
  double RadStep = .01;
  double HeiStep = .1;
  //Yuliya
  double *Count = (double *)calloc(NBin*NBin,sizeof(double));
  // first configuration
  for(int d=0;d<3;d++) Nano->Pos[d] = OldPos[d];
  Nano->Rad = MAX(MIN(RadMax,OldPos[3] + 1.0 ),RadMin);
  Nano->Height = MAX(MIN(HeiMax,OldPos[4] + 0.),HeiMin);
  double OldNrg = 0.;
  for(int p=0;p<pNPart();p++){
    double Dist2 = NanoDist2(Pm[p].Pos,0);
    if(Dist2 > SQR(2.*Nano->Rad)) continue;
    if(Pm[p].Typ == 0) OldNrg += kElPhob*Dist2;
    if(Pm[p].Typ == 1) OldNrg += kElPhil*Dist2;
  }
  OldNrg += (GainRad*SQR(Nano->Rad-MinRad) + GainHei*SQR(Nano->Height-MinHei));
  // fake Monte Carlo
  for(int i=0;i<NInt;i++){
    //change position
    double OldNPos[5] = {Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],Nano->Rad,Nano->Height};
    for(int m=0;m<3;m++){
      //change the torus
      if(m==0){
      	for(int d=0;d<3;d++){
      	  Nano->Pos[d] += MoveStep*(2.*Mat->Casuale()-1.);
      	}
	SetNanoBkf(0);
      }
      if(m==1){
      	Nano->Rad += RadStep*(2.*Mat->Casuale()-1.);
      	if(Nano->Rad < RadMin || Nano->Rad > RadMax){
      	  Nano->Rad = OldNPos[3];
      	  continue;
      	}
      }
      if(m==2){
      	Nano->Height += HeiStep*(2.*Mat->Casuale()-1.);
      	if(Nano->Height < HeiMin || Nano->Height > HeiMax){
      	  Nano->Height = OldNPos[4];
      	  continue;
      	}
      }
      // recalc the energy
      double Nrg = 0.;
      for(int p=0;p<pNPart();p++){
	double Dist2 = NanoDist2(Pm[p].Pos,0);
	if(Dist2 > SQR(2.*Nano->Rad)) continue;
	if(Pm[p].Typ == 0) Nrg += kElPhob*Dist2;
	if(Pm[p].Typ == 1) Nrg += kElPhil*Dist2;
      }
      Nrg += (GainRad*SQR(Nano->Rad-MinRad) + GainHei*SQR(Nano->Height-MinHei));
      // accept/remove
      if(exp(OldNrg-Nrg) > Mat->Casuale()){
	if(m==0){
	  printf("NewNrg \t%lf Rad ___ Hei ___ Pos %.2f %.2f %.2f\n",Nrg-OldNrg,Nano->Pos[0]-OldNPos[0],Nano->Pos[1]-OldNPos[1],Nano->Pos[2]-OldNPos[2]);
	}
	if(m==1){
	  printf("NewNrg \t%lf Rad %.3f Hei ___ Pos ___ ___ ___\n",Nrg-OldNrg,Nano->Rad);
	}
	if(m==2){
	  printf("NewNrg \t%lf Rad ___ Hei %.3f Pos ___ ___ ___\n",Nrg-OldNrg,Nano->Height);
	}
	OldNrg = Nrg;
      }
      else{
	if(m==0){
	  for(int d=0;d<3;d++){
	    Nano->Pos[d] = OldNPos[d];
	  }
	  SetNanoBkf(0);
	}
 	if(m==1){
	  Nano->Rad = OldNPos[3];
	}
 	if(m==2){
	  Nano->Height = OldNPos[4];
	}
      }
      //printf("%lf %lf %lf -> %lf %lf %lf\n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],OldNPos[0],OldNPos[1],OldNPos[2]);
      //printf("%lf-> %lf %lf-> %lf\n",SRad,Nano->Rad,LRad,Nano->Height);  
    }
  }
  double Cm[3];
  double CountCm = 0;
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ == 1) continue;
    double Pos[3];
    for(int d=0;d<3;d++){
      Pos[d] = Pm[p].Pos[d] - Nano->Pos[d];
      Pos[d] -= floor(Pos[d]*pInvEdge(d))*pEdge(d);
    }
    double Rad = sqrt( SQR(Pos[CLat1]) + SQR(Pos[CLat2]) );
    if(Rad > Nano->Height) continue;
    for(int d=0;d<3;d++){
      Cm[d] += Pm[p].Pos[d];
    }
    CountCm += 1.;
    if(Pm[p].Pos[CNorm] < Nano->Pos[CNorm] - Nano->Rad || Pm[p].Pos[CNorm] > Nano->Pos[CNorm] + Nano->Rad) continue;
    int vx = (int)(Pos[CLat1]/Nano->Height*NBin);
    vx += NBin/2;
    if(vx < 0 || vx >= NBin) continue;
    int vy = (int)(Pos[CLat2]/Nano->Height*NBin);
    vy += NBin/2;
    if(vy < 0 || vy >= NBin) continue;
    Count[vx*NBin+vy] += 1.;
  }
  double Area = 0.;
  for(int vx=0;vx<NBin;vx++){
    for(int vy=0;vy<NBin;vy++){
      if(Count[vx*NBin+vy] < 1.) continue;
      Area += 1.;
    }
  }
  if(CountCm <= 0.){
    printf("No particles in the torus\n");
    return 1;
  }
  Nano->Area = SQR(Nano->Height)*Area/(double)(SQR(NBin));
  for(int d=0;d<3;d++){
    Cm[d] /= CountCm;
    //Nano->Pos[d] = Cm[d];
    CmStalk[d] = Nano->Pos[d];
  }
  SetNanoBkf(0);
  printf("Pos %lf %lf %lf Area %lf Count %lf\n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],Nano->Area,CountCm);
  free(Count);
  return 0;
}
#include <Cubo.h>
double VarData::NormalWeight(VAR_TRIANGLE *Triang,double *WeightL,int NGrid,int NTri){
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  NeiVertex VList(NTri,3,NGrid,Edge);
  double CountStalk = 0.;
  Vettore Ax(1.,0.,0.);
  Vettore n1(3);
  Vettore n2(3);
  double Max = 0.;
  for(int t=0;t<NTri;t++){
    for(int v=0;v<3;v++){
      int vCurr = Triang[t].v[v];
      VList.Add(vCurr,t,Triang[t].p[v].x);
    }
  }
  VList.Reorder();
  VList.SetCounters();
  for(int t=0;t<NTri;t++){
    for(int v=0;v<3;v++){
      int vCurr = Triang[t].v[v];
      double Weight = 0.;
      for(VList.SetCounters(vCurr);VList.IfItCell(vCurr);VList.IncrCurr(vCurr)){
	int tt = VList.VertCurr(vCurr);
	int p = tt*pNPCh();
	for(int d=0;d<3;d++){
	  n1.Set(Pm[p].Vel[d],d);
	}
	if(n1.Norm() <= 0.) continue;
	Weight += n2.Angle(&Ax,&n1);
      }
      WeightL[vCurr] = Weight;
      if(Max < Weight) Max = Weight;
    }
  }
  return Max;
}
/** Algorithm to connect al the vertices in a single chain, many weird cases are not covered */
void VarData::ConnectLineChain(VAR_LINE *Triang,int NGrid,int NTri){
  SetNPart(2*NTri);
  int *Exist = (int *)calloc(NGrid*NGrid,sizeof(int));
  int NPart = 0;
  DdLinkedList *Pc = new DdLinkedList(Gen->Edge,pNPart(),1.5);
  double InvNGrid = 1./(double)NGrid;
  for(int t=0;t<NTri;t++){
    for(int v=0;v<2;v++){
      int vx = (int)(Triang[t].p[v].x[0]*pEdge(0)*InvNGrid);
      int vy = (int)(Triang[t].p[v].x[1]*pEdge(1)*InvNGrid);
      int vv = vx*NGrid+vy;
      Exist[vv] += 1;
      if(Exist[vv] > 1) continue;
      Pc->AddPart(NPart,Triang[t].p[v].x);
      for(int d=0;d<3;d++) Pm[NPart].Pos[d] = Triang[t].p[v].x[d];
      NPart++;
    }
  }
  SetNPart(NPart);
  for(int p=0;p<pNPart();p++){
    Ln[p].Link[0] = 0;
  }
  double DistRel[4];
  char FName[60];
  int link = 0;
  double Pos[5] = {0.,.5*pCm(CLat2),pCm(CLat2),1.5*pCm(CLat2),pEdge(CLat2)};
  double pList[4];
  for(int p=0,c=0;p<pNPart();p++){
    if(Pm[p].Pos[CLat1] <= 0.9 && Pm[p].Pos[CLat2] > Pos[c] + 3.){
      Pos[c+1] = Pm[p].Pos[CLat2] + .1;
      pList[c] = p;
      c++;
      if(c==4) break;
    }
  }
  for(int p=pList[0],c=0;c<4;c++,p=pList[c]){
    sprintf(FName,"Chain%d.dat",c);
    FILE *FChain = fopen(FName,"w");
    fprintf(FChain,"%lf %lf\n",Pm[p].Pos[0],Pm[p].Pos[1]);
    link = p;
    for(int p1=0;p1<pNPart();p1++){
      for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
	if(p1 == Pc->p2Curr) continue;
	//printf("%d %d %d %d\n",p,p1,Pc->p2Curr,link);
	if(link == Pc->p2Curr){
	  fprintf(FChain,"%lf %lf\n",Pm[p1].Pos[0],Pm[p1].Pos[1]);
	  Ln[p1].Link[0] = link;
	  link = p1;
	  break;
	}
      }
    }
  }
}
/** Algorithm to connect al the vertices in a single chain*/
void VarData::ConnectLineChain3(VAR_LINE *Triang,int NGrid,int NTri){
  SetNPart(2*NTri);
  int *Called = (int *)calloc(pNPart(),sizeof(int));
  int *Exist = (int *)calloc(NGrid*NGrid,sizeof(int));
  double *DirPrev = (double *)calloc(3*pNPart(),sizeof(int));
  int NChain = 0;
  int NPart = 0;
  double InvNGrid = 1./(double)NGrid;
  DdLinkedList *Pc = new DdLinkedList(Gen->Edge,pNPart(),2.5);
  for(int t=0;t<NTri;t++){
    for(int v=0;v<2;v++){
      int vx = (int)(Triang[t].p[v].x[0]*pEdge(0)*InvNGrid);
      int vy = (int)(Triang[t].p[v].x[1]*pEdge(1)*InvNGrid);
      int vv = vx*NGrid+vy;
      Exist[vv] += 1;
      if(Exist[vv] > 1) continue;
      Pc->AddPart(NPart,Triang[t].p[v].x);
      for(int d=0;d<3;d++) Pm[NPart].Pos[d] = Triang[t].p[v].x[d];
      NPart++;
    }
  }
  SetNPart(NPart);
  for(int p=0;p<pNPart();p++){
    Ln[p].NLink = 0;
    Ln[p].Link[0] = -1;
    Exist[p] = 0;
    Pm[p].CId = p;
  }
  // for(int p=0;p<pNPart();p++){
  //   Ln[p].NLink = 1;
  //   int link = Pc->FindClosest(p);
  //   if(link == -1) Ln[p].NLink = 0;
  //   Ln[p].Link[0] = link;
  // }
  double DistRel[4];
  for(int p=0;p<pNPart();p++){
    Ln[p].NLink = 1;
    double Closest = 1000000.;
    int pClosest = -1;
    for(Pc->SetCurr(p);Pc->IfCurr();Pc->NextCurr()){
      int link = Pc->p2Curr;
      if(link == p) continue;
      if(Ln[link].Link[0] == p) continue;
      Pc->Dist2Curr(DistRel);
      double Weight = DirPrev[link*3+0]*DistRel[0] + DirPrev[link*3+1]*DistRel[1];
      if(Weight >= 0.) Weight = 0.9;
      else if(Weight < 0.) Weight = 1.1;
      if(DistRel[3]*Weight < Closest){
  	Closest = DistRel[3];
  	pClosest = link;
  	for(int d=0;d<3;d++){
	  DirPrev[3*p+d] = DistRel[d];
  	}
      }
    }
    if(pClosest == -1){
      Ln[p].NLink = 0;
      continue;
    }
    Ln[p].Link[0] = pClosest;
    Exist[pClosest] += 1;
    if(Exist[pClosest] > 1 && Ln[pClosest].NLink == 1) Ln[p].NLink = 0;
  }
  NChain = 0;
  for(int p=0;p<pNPart();p++){
    int link = Ln[p].Link[0];
    if(Pm[p].CId > NChain){
      NChain++;
      Pm[p].CId = NChain;
    }
    Pm[link].CId = Pm[p].CId;
  }
  NChain = 0;
  for(int p=pNPart();p>=0;p--){
    int link = Ln[p].Link[0];
    if(Pm[p].CId > NChain){
      NChain++;
      Pm[p].CId = NChain;
    }
    Pm[link].CId = Pm[p].CId;
  }
  NChain = 0;
  for(int p=0;p<pNPart();p++){
    int link = Ln[p].Link[0];
    if(Pm[p].CId > NChain){
      NChain++;
      Pm[p].CId = NChain;
    }
    Pm[link].CId = Pm[p].CId;
  }
  SetNChain(NChain);
  return;
  //isolate multiple connected points
  // for(int p=0;p<pNPart();p++){
  //   int link = Ln[p].Link[0];
  //   Called[link]++;
  //   if(Called[link] > 1){
  //     Ln[p].NLink = 0;
  //     NPart--;
  //   }
  // }
  for(int p=0;p<pNPart()-1;p++){
    if(Ln[p].NLink == 0) continue;
    int link = Ln[p].Link[0];
    int MemPos = -1;
    for(int pp=p+1;pp<pNPart();pp++){
      if(Pm[pp].Idx == link){
	MemPos = pp;
	break;
      }
    }
    if(MemPos == -1){
      Ln[p].NLink = 0;
      continue;
    }
    printf("%d %d %d %d\n",p,link,MemPos,Pm[p+1].Idx);
    //if(p < pNPart() - 2) SwapPart(p+1,p+2);
    SwapPart(p+1,MemPos);
    printf("%d %d\n",Pm[p+1].Idx,Pm[MemPos].Idx);
  }
  for(int p=0;p<pNPart()-1;p++){
    Ln[p].Link[0] = p+1;
  }
  // for(int p=0;p<pNPart();p++){
  //   int link = Ln[p].Link[0];
  //   for(int pp=0;pp<pNPart();pp++){
  //     if(Pm[pp].Idx == link){
  // 	Ln[p].Link[0] = pp;
  // 	break;
  //     }
  //   }
  // }
  SetNPart(NPart);
  NChain = 0;
  for(int p=0;p<pNPart();p++){
    Pm[p].CId = NChain;
    if(Ln[p].NLink == 0) NChain++;
  }
  SetNChain(NChain);
  free(Exist);
  free(Called);
  free(DirPrev);
}
/** Algorithm to connect al the vertices in a single chain, many weird cases are not coverd*/
void VarData::ConnectLineChain2(VAR_LINE *Triang,int NGrid,int NTri){
  //can't follow degenerate triangles..
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  NeiVertex VList(NTri,2,NGrid,Edge);
  double CountStalk = 0.;
  Vettore Ax(1.,0.,0.);
  Vettore n1(3);
  Vettore n2(3);
  double Max = 0.;
  int NTria = 0;
  for(int t=0;t<NTri;t++){
    if(Triang[t].p[0].x[0] == Triang[t].p[1].x[0] && Triang[t].p[0].x[1] == Triang[t].p[1].x[1]) continue;
    for(int v=0;v<2;v++){
      int vCurr = Triang[t].v[v];
      VList.Add(vCurr,NTria,Triang[t].p[v].x);
      for(int d=0;d<3;d++) Pm[NTria*2+v].Pos[d] = Triang[t].p[v].x[d];
      NTria++;
    }
  }
  VList.Reorder();
  SetNPart(2*NTri);
  //VList.Print();
  VList.SetCounters();
  for(int p=0;p<pNPart();p++){
    Ln[p].NLink = 0;
  }
  for(int t=0;t<NTri;t++){
    for(int v=0;v<2;v++){
      int vCurr = Triang[t].v[v];
      for(VList.SetCounters(vCurr);VList.IfItCell(vCurr);VList.IncrCurr(vCurr)){
	int tt = VList.TriaCurr(vCurr);
	int p = tt*2+v;
	if(vCurr > tt){
	  Ln[p].NLink = 1;
	  Ln[p].Link[0] = vCurr;
	}
      }
    }
  }
  int nChain=0;
  for(int c=0;c<pNChain();c++){
    Ch[c].NPCh = 0;
  }
  for(int p=0;p<pNPart();p++){
    if(Pm[p].CId <= nChain) continue;
    Ch[nChain].InitBead = p;
    int pp1 = Ln[p].Link[0];
    int pp = 0;
    printf("------ %d %d\n",p,pp1);
    for(pp=pp1;Ln[pp].NLink > 0;pp=Ln[pp].Link[0]){
      printf(" %d %d\n",pp,Ch[nChain].NPCh);
      Pm[pp1].CId = nChain;
      Ch[nChain].NPCh++;
      pp = Ln[pp].Link[0];
      if(pp == pp1) break;
    }
    Ch[nChain].EndBead = pp;
    nChain++;
  }
  SetNChain(nChain);
  for(int c=0;c<pNChain();c++){
    printf("%d %d %d %d\n",c,Ch[c].InitBead,Ch[c].EndBead,Ch[c].NPCh);
    for(int p=Ch[c].InitBead,ppc=0;ppc<Ch[c].NPCh;p=Ln[p].Link[0]){
      printf("%d %d \n",c,p);
    }
  }
}
int VarData::StalkPos(double *OldPos){
  double CmStalk[3] = {OldPos[0],OldPos[1],OldPos[2]};
  //pPos(OldPos);
  //if(StalkPos4(OldPos,CmStalk)) return 1;
  int NTrials = 10;
  int n=0;
  while(StalkPos4(OldPos,CmStalk)){
    n++;
    if(n > NTrials) continue;
  }
  //pPos(CmStalk);
  for(int d=0;d<3;d++){
    //Nano->Pos[d]   = CmStalk[d];// + OldPos[d];
    //Nano->PosBf[d] = Nano->Pos[d];// - floor(Nano->Pos[d]*pInvEdge(d))*pEdge(d);
    OldPos[d] = pNanoPos(0,d);;
  }
  OldPos[3] = Nano->Rad;
  OldPos[4] = Nano->Height;
  Nano->Axis[CLat1] = 0.;
  Nano->Axis[CLat2] = 0.;
  Nano->Axis[CNorm] = 1.;
  // Nano->Rad     = .5;
  // Nano->Height  = 4.;
  // Nano->Hamaker = 1.;
  // Nano->OffSet  = 1.;
  // Nano->Shape = SHAPE_STALK;
  return 0;
}
double VarData::PorePos(){
  const int NGrid = 26;
  double *Plot = (double *)calloc(NGrid*NGrid,sizeof(double));
  double *Plot2 = (double *)calloc(NGrid*NGrid,sizeof(double));
  double Cm[3]  = {0.,0.,0.};
  double Cm1[3] = {0.,0.,0.};
  double Incr = 500.;
  double Threshold = SQR(SQR(1./Incr));
  int NReweight = 2;
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      Plot[gx*NGrid+gy] = 1.;
    }
  }
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ != 0 ) continue;
    int gx = (int)(Pm[p].Pos[CLat1]*pInvEdge(CLat1)*NGrid);
    int gy = (int)(Pm[p].Pos[CLat2]*pInvEdge(CLat2)*NGrid);
    if(gx < 0 || gx >= NGrid) continue;
    if(gy < 0 || gy >= NGrid) continue;
    Plot[gx*NGrid+gy] += Incr;
  }
  double Dx = .5*pEdge(CLat1)/(double)NGrid;
  double Dy = .5*pEdge(CLat2)/(double)NGrid;
  double Count = 0.;
  //FILE *Ciccia = fopen("PoreGrid.dat","w");
  double Area = 0.;
  //weigthing by neighbours
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      Plot2[gx*NGrid+gy] = Plot[gx*NGrid+gy];
    }
  }
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      if(gx == 0){ Plot[gx*NGrid+gy] = 1./SQR(SQR(Plot2[gx*NGrid+gy]));continue;}
      if(gy == 0){ Plot[gx*NGrid+gy] = 1./SQR(SQR(Plot2[gx*NGrid+gy]));continue;}
      if(gx == NGrid-1){ Plot[gx*NGrid+gy] = 1./SQR(SQR(Plot2[gx*NGrid+gy]));continue;}
      if(gy == NGrid-1){ Plot[gx*NGrid+gy] = 1./SQR(SQR(Plot2[gx*NGrid+gy]));continue;}
      Plot[gx*NGrid+gy] = 1./(Plot2[(gx-1)*NGrid+gy]*Plot2[(gx+1)*NGrid+gy]*Plot2[gx*NGrid+(gy-1)]*Plot2[gx*NGrid+(gy+1)]);
    }
  }
  //above the threshold
  for(int gx=0;gx<NGrid;gx++){
    double x = gx/(double)NGrid*pEdge(CLat1) + Dx;
    for(int gy=0;gy<NGrid;gy++){
      double y = gy/(double)NGrid*pEdge(CLat2) + Dy;
      //Plot[gx*NGrid+gy] /= 100.;
      //fprintf(Ciccia,"%lf %lf %lf 0\n",x,y,10.*Plot[gx*NGrid+gy]);
      Cm[CLat1] +=  x*Plot[gx*NGrid+gy];
      Cm[CLat2] +=  y*Plot[gx*NGrid+gy];
      Count     += Plot[gx*NGrid+gy];
      if(Plot[gx*NGrid+gy] > Threshold ){
	Area += 1.;//Plot[gx*NGrid+gy];
	//fprintf(Ciccia,"%lf %lf %lf %d\n",x,y,1.,2);
      }
      //else fprintf(Ciccia,"%lf %lf %lf %d\n",x,y,.5,0);
    }
  }
  Area = Area/(double)SQR(NGrid)*pEdge(CLat1)*pEdge(CLat2);
  Area /= 1.;//Count;
  Nano->Rad = sqrt(Area/M_PI);
  Cm[0] /= Count;
  Cm[1] /= Count;
  //reweighting wrt to the former position
  for(int r=0;r<NReweight;r++){
    Count = 0.;
    Cm1[0] = Cm[0];Cm1[1] = Cm[1];
    Cm[0] = 0.;Cm[1] = 0.;
    Area = 0.;
    //printf("%lf %lf\n",Cm1[0],Cm1[1]);
    for(int gx=0;gx<NGrid;gx++){
      double x = gx/(double)NGrid*pEdge(CLat1) + Dx;
      for(int gy=0;gy<NGrid;gy++){
	double y = gy/(double)NGrid*pEdge(CLat2) + Dy;
	double Dist = (SQR(x-Cm1[CLat1]) + SQR(y-Cm1[CLat2]));
	if(Dist > SQR(1.1*Nano->Rad)) continue;
	double Weight = 1./pow(Dist,.2);
	//fprintf(Ciccia,"%lf %lf %lf 1\n",x,y,10.*Plot[gx*NGrid+gy]*Weight);
	Cm[CLat1] +=  x*Plot[gx*NGrid+gy]*Weight;
	Cm[CLat2] +=  y*Plot[gx*NGrid+gy]*Weight;
	Count     +=  Plot[gx*NGrid+gy]*Weight;
	if(Plot[gx*NGrid+gy] > Threshold ) Area += 1.;
     }
    }
    Cm[0] /= Count;
    Cm[1] /= Count;
    Area = Area/(double)SQR(NGrid)*pEdge(CLat1)*pEdge(CLat2);
    //Nano->Rad = sqrt(Area/M_PI);
  }
  //assigning
  Cm[2] = pCm(2);
  for(int d=0;d<3;d++){
    Nano->Pos[d] = Cm[d];
  }
  SetNanoBkf(0);
  //asymmetry
  double AreaOut = 0.;
  double AreaIn = 0.;
  for(int gx=0;gx<NGrid;gx++){
    double x = gx/(double)NGrid*pEdge(CLat1) + Dx;
    for(int gy=0;gy<NGrid;gy++){
      double y = gy/(double)NGrid*pEdge(CLat2) + Dy;
      double Dist = SQR(x-Cm[0]) + SQR(y-Cm[1]);
      if(Dist <= SQR(Nano->Rad)){
	AreaIn += 1.;
	if(Plot[gx*NGrid+gy] < Threshold ) AreaOut += 1.;
	//fprintf(Ciccia,"%lf %lf %lf %d\n",x,y,0.,1);
      }
      else{
	if(Plot[gx*NGrid+gy] > Threshold ) AreaOut += 1.;
      }
    }
  }
  //fclose(Ciccia);
  free(Plot);
  return AreaOut/AreaIn;
}
