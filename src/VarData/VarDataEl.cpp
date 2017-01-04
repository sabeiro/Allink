/***********************************************************************
VarDataEl: Elaboration functions for the VarData class. This functions
provides a simple manipulation of the data read by [Open]. The
options are provided to elaborate different system's shapes.
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

int VarData::RadDistr(int Values,double *Plot,double Border[2],int How){
  return 1;
}
char *VarData::SysState(){
  char *Info;
  Info = (char *)malloc(256*sizeof(char));
  Gen->Temp = 0.;
  double VelSquare[3];
  for(int d=0;d<3;d++){
    Gen->Pre[d] = 0.;
    VelSquare[d] = 0.;
  }
  for(int p=0;p<Gen->NPart;p++){
    for(int d=0;d<3;d++){
      Gen->Pre[d] += SQUARE(Pm[p].Vel[d]);
    }
  }
  Gen->Temp = (Gen->Pre[0] + Gen->Pre[1] + Gen->Pre[2])/(3.*Gen->NPart);
  for(int d=0;d<3;d++){
    Gen->Pre[d] /= (Gen->Edge[0]*Gen->Edge[1]*Gen->Edge[2]);
  }
  Gen->SurfTens = Gen->Edge[2]*(Gen->Pre[2] - .5*(Gen->Pre[0] + Gen->Pre[1]) );
  sprintf(Info,"Pre0\tPre1\tPre2\tSurfTens\tTemp\n%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n",Gen->Pre[0],Gen->Pre[1],Gen->Pre[2],Gen->SurfTens,Gen->Temp);
  return Info;
}
MOMENTI VarData::SampleSurfaceMem(int NSample){
  if(IfPlotMem){
    MOMENTI m1;
    m1.Num = 0;
  }
  PlotMem = (double *) calloc(SQR(NSample),sizeof(double));
  double *PlotR = (double *) calloc(SQR(NEdge),sizeof(double));
  //LoadDensFile(PlotR,NEdge);
  MOMENTI m1 = SampleSurfacePart(PlotR,NEdge,0);
  InterBSpline2D(PlotR,PlotMem,NEdge,NSample);
  if(IfPlotMem){
    free(PlotR);
  }
  IfPlotMem = 1;
  return m1;
}
void VarData::SampleSurface(double *Plot,int NSample,int Type){
  double Average=0.;
  double NAverage=0;
  double **Norma = (double **)calloc(NSample,sizeof(double));
  for(int v=0;v<NSample;v++){
    Norma[v] = (double *)calloc(NSample,sizeof(double));
    for(int vv=0;vv<NSample;vv++){
      Plot[v*NSample+vv] = 0.;
    }
  }
  if(Type == 0){
    for(int c=0;c<Gen->NChain;c++){
      //if(!CHAIN_IF_TYPE(Ch[c].Type,NChType))continue;
      Average += Ch[c].Pos[CNorm];NAverage += 1.;
	int v = (int)(Ch[c].Pos[CLat1]*pInvEdge(CLat1)*NSample);
	if( v < 0 || v >= NSample) continue;
	int vv = (int)(Ch[c].Pos[CLat2]*pInvEdge(CLat2)*NSample);
	if( vv < 0 || vv >= NSample) continue;
	Plot[v*NSample+vv] += Ch[c].Pos[CNorm];
	Norma[v][vv] += 1.;
	//printf("%d %d %lf %lf\n",v,vv,Plot[v][vv],Norma[v][vv]);
    }
  }
  else {
    for(int p=0;p<Gen->NPart;p++){
      Average += Pm[p].Pos[CNorm];NAverage += 1.;
	int v = (int)(Pm[p].Pos[CLat1]*pInvEdge(CLat1)*NSample);
	if( v < 0 || v >= NSample) continue;
	int vv = (int)(Pm[p].Pos[CLat2]*pInvEdge(CLat2)*NSample);
	if( vv < 0 || vv >= NSample) continue;
	Plot[v*NSample+vv] += Pm[p].Pos[CNorm];
	Norma[v][vv] += 1.;
	//printf("%d %d %lf %lf\n",v,vv,Plot[v][vv],Norma[v][vv]);
    }
  }
  Average /= (double)NAverage;
  for(int v=0;v<NSample;v++){
    for(int vv=0;vv<NSample;vv++){
      if(Norma[v][vv] > 0.){
	Plot[v*NSample+vv] /= Norma[v][vv];
      }
      else Plot[v*NSample+vv] = Average;
    }
  }
  for(int v=0;v<NSample;v++)
    free(Norma[v]);
  free(Norma);
}
MOMENTI VarData::SampleSurface(Matrice *Plot,int NSample,int Type){
  MOMENTI m1;
  double Average=0.;
  double NAverage=0;
  double **Norma = (double **)calloc(NSample,sizeof(double));
  for(int v=0;v<NSample;v++){
    Norma[v] = (double *)calloc(NSample,sizeof(double));
    for(int vv=0;vv<NSample;vv++){
      Plot->Set(v,vv,0);
    }
  }
  for(int c=0;c<Gen->NChain;c++){
    if(!CHAIN_IF_TYPE(Ch[c].Type,NChType))continue;
    Average += Ch[c].Pos[CNorm];NAverage += 1.;
      int v = (int)(Ch[c].Pos[CLat1]*pInvEdge(CLat1)*NSample);
      if( v < 0 || v >= NSample) continue;
      int vv = (int)(Ch[c].Pos[CLat2]*pInvEdge(CLat2)*NSample);
      if( vv < 0 || vv >= NSample) continue;
      Plot->Add(v,vv,Ch[c].Pos[CNorm]);
      Norma[v][vv] += 1.;
      //printf("%d %d %'lf\n",v,vv,Plot[v][vv]);
  }
  Average /= (double)NAverage;
  for(int v=0;v<NSample;v++){
    for(int vv=0;vv<NSample;vv++){
      if(Norma[v][vv] > 0.)
	Plot->Set(v,vv,Plot->Val(v,vv)/Norma[v][vv]);
      //pasticcio
      else Plot->Set(v,vv,Average);
      //printf("%d %d %'lf\n",v,vv,Plot[v][vv]);
    }
  }
  free(Norma);
  m1.Uno = Average;
  m1.Num = SQR(NSample);
  return m1;
}
void VarData::LoadDensFile(double **Plot,int NSample){
  double Round = 0.001;
  double *Count = (double *)calloc(3*NSample*NSample,sizeof(double));
  for(int p=0;p<pNPart();p++){
    int sr = (int)((pPos(p,0)+Round)*pInvEdge(0)*NSample);
    if(sr < 0 || sr >= NSample) continue;
    int sz = (int)((pPos(p,1)+Round)*pInvEdge(1)*NSample);
    if(sz < 0 || sz >= NSample) continue;
    int t = pType(p);
    Plot[t][sr*NSample+sz] += pPos(p,2);
    Count[(sr*NSample+sz)*3+t] += 1.;
  }  
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  for(int t=0;t<3;t++){
    for(int s=0;s<NSample*NSample;s++){
      if(Count[s*3+t] > 0.)
	Plot[t][s] /= Count[s*3+t];
    }
  }
  int NDim = 2;
  int IfMinImConv = 1;
  for(int t=0;t<3;t++){
    Mask.ConvoluteMatrix(Plot[t],NSample,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot[t],NSample,NDim,IfMinImConv);
  }
  free(Count);
}
MOMENTI VarData::SampleSurfacePart(double *Plot,int NSample,int Type){
  double Round = 0.001;
  double Average=0.;
  double NAverage=0;
  double Min = 10000000.;
  double Max =-10000000.;
  double **Norma = (double **)calloc(NSample,sizeof(double));
  MOMENTI m1;
  for(int v=0;v<NSample;v++){
    Norma[v] = (double *)calloc(NSample,sizeof(double));
    for(int vv=0;vv<NSample;vv++){
      Plot[v*NSample+vv] = 0;
    }
  }
  for(int p=0;p<Gen->NPart;p++){
    if(Pm[p].Typ != Type) continue;
    int v = (int)((Pm[p].Pos[CLat1]+Round)*pInvEdge(CLat1)*NSample);
    if( v < 0 || v >= NSample) continue;
    int vv = (int)((Pm[p].Pos[CLat2]+Round)*pInvEdge(CLat2)*NSample);
    if( vv < 0 || vv >= NSample) continue;
    Plot[v*NSample+vv] += Pm[p].Pos[CNorm];
    Norma[v][vv] += 1.;
    Average += Pm[p].Pos[CNorm];NAverage += 1.;
  }
  Average /= (double)NAverage;
  for(int v=0;v<NSample;v++){
    for(int vv=0;vv<NSample;vv++){
      if(Norma[v][vv] > 0.){
	Plot[v*NSample+vv] /= Norma[v][vv];
      }
      else Plot[v*NSample+vv] = Average;
      if(Plot[v*NSample+vv] < Min) Min = Plot[v*NSample+vv];
      if(Plot[v*NSample+vv] > Max) Max = Plot[v*NSample+vv];
    }
  }
  for(int v=0;v<NSample;v++)
    free(Norma[v]);
  free(Norma);
  m1.Min = Min;
  m1.Max = Max;
  m1.Uno = Average;
  m1.Num = SQR(NSample);
  return m1;
}
int VarData::SpatialDerivative(Matrice *Surface,Matrice *Resp,SPLINE Weight,int NSample){
  SampleSurface(Surface,NSample,NChType);
  Matrice *Mask = new Matrice(Weight);
  for(int h=0;h<Resp->Size();h++)
    for(int w=0;w<Resp->Size();w++)
      Resp->Set(h,w,0.);
  Mat->ApplyFilter(Surface,Resp,Mask);
  Mask->Transpose();
  Mat->ApplyFilter(Surface,Resp,Mask);
  delete Mask;
  return 0;
}
Properties VarData::SysProperties(){
  Properties *Pr = new Properties();
  double StructPhil = 0.;
  double StructPhob = 0.;
  for(int c=0;c<Gen->NChain;c++){
    for(int cc=0;cc<Gen->NChain;cc++){
      Pr->ChDiff += QUAD((Ch[c].Pos[CLat1]-Ch[cc].Pos[CLat1]));
      Pr->ChDiff += QUAD((Ch[c].Pos[CLat2]-Ch[cc].Pos[CLat2]));
    }
    double CmPhil[3] = {0.,0.,0.};
    double CmPhob[3] = {0.,0.,0.};
    double RadPhil=0.;
    double RadPhob=0.;
    double NPhil=0.;
    double NPhob=0.;
    double q=1./Gen->Edge[3];
    double Dist[3] = {0.,0.,0.}; 
    for(int p=c*pNPCh();p<pNPCh()*c+Block[0].Asym;p++){//Phob
      for(int d=0;d<3;d++){
	CmPhob[d] += Pm[p].Pos[d];
	for(int pp=c*Gen->NPCh;pp<Gen->NPCh*c+Block[0].Asym-1;pp++){//Phob
	  Dist[d] = Pm[p].Pos[d] - Pm[pp].Pos[d];
	  if(p == pp) StructPhob += 1.;
	  else 
	    //StructPhob += Dist[d]*sin(q*Dist[d])/(q);
	    StructPhob += exp(-QUAD(( q*Dist[d] ))/6. );
	}
	if(p == Gen->NPCh*c+Block[0].Asym-1) continue;
	Pr->RePhob += QUAD(( Pm[p].Pos[d] - Pm[p+1].Pos[d]));
      }
    }
    for(int p=c*Gen->NPCh+Block[0].Asym;p<Gen->NPCh*(c+1)-1;p++){//Phil
      for(int d=0;d<3;d++){
	CmPhil[d] += Pm[p].Pos[d];
	for(int pp=c*Gen->NPCh+Block[0].Asym;pp<Gen->NPCh*(c+1)-1;pp++){//Phil
	  Dist[d] = Pm[p].Pos[d] - Pm[pp].Pos[d];
	  if(p == pp) StructPhil += 1.;
	  else 
	    //	    StructPhil += sin(q*Dist[d])/(q*Dist[d]);
	    StructPhil += exp(-QUAD(( q*Dist[d] ))/6. );
	}
	if( p == Gen->NPCh*(c+1)-1) continue;
	Pr->RePhil += QUAD(( Pm[p].Pos[d] - Pm[p+1].Pos[d]));
      }
    }
    for(int d=0;d<3;d++){
      CmPhob[d] /= (double) Block[0].Asym;
      CmPhil[d] /= (double) Gen->NPCh-Block[0].Asym;
    }
    for(int p=c*Gen->NPCh;p<Gen->NPCh*c+Block[0].Asym;p++){//Phob
      for(int d=0;d<3;d++)
	RadPhob += QUAD(( CmPhob[d] - Pm[p].Pos[d]));
    }
    for(int p=c*Gen->NPCh+Block[0].Asym;p<Gen->NPCh*(c+1);p++){//Phil
      for(int d=0;d<3;d++)
	RadPhil += QUAD(( CmPhil[d] - Pm[p].Pos[d]));
    }
    Pr->GyrPhob += RadPhob / (double)(Block[0].Asym*Gen->NChain);
    Pr->GyrPhil += RadPhil / (double)((Gen->NPCh - Block[0].Asym)*Gen->NChain);
  }
  Pr->RePhob /= Gen->NChain;
  Pr->RePhil /= Gen->NChain;
  Pr->ChDiff /= (double)QUAD((Gen->NChain));
  Pr->FactPhob += StructPhob/(double)(Gen->NChain);//*QUAD((Asym)));
  Pr->FactPhil += StructPhil/(double)(Gen->NChain);//*QUAD((Gen->NPCh-Asym)));
  Pr->Print();
  return *Pr;
}
int VarData::Folding(){
  double *Interp1 = (double *)calloc(Gen->NPCh,sizeof(double));
  double *Interp2 = (double *)calloc(Gen->NPCh,sizeof(double));
  double *InterpN = (double *)calloc(Gen->NPCh,sizeof(double));
  RETTA r1;
  RETTA r2;
  for(int c=0;c<Gen->NChain;c++){
    double Comp[3]={0.,0.,0.};
    for(int p= c*Gen->NPCh,pp=0;p<(c+1)*Gen->NPCh;p++,pp++){
      Interp1[pp] = Pm[p].Pos[CLat1];
      Interp2[pp] = Pm[p].Pos[CLat2];
      InterpN[pp] = Pm[p].Pos[CNorm];
      if(p == (c+1)*Gen->NPCh - 1) break;
      double Segm[3]={0.,0.,0.};
      double Rad=0.;
      for(int d=0;d<3;d++){
	Segm[d] = Pm[p].Pos[d] - Pm[p+1].Pos[d];
	Rad += QUAD((Segm[d]));
      }
      Rad = sqrt(Rad);
      for(int d=0;d<3;d++){
	Comp[d] += Segm[d] / (Rad*Gen->NPCh);
      }
    }
    if( ASS((Comp[CNorm])) < .7)
      Ch[c].Type |= CHAIN_FLABBY;
    else 
      Ch[c].Type |= CHAIN_STRETCH;
    r1 = Mat->InterRett(Interp1,InterpN,Gen->NPCh);
    r2 = Mat->InterRett(Interp2,InterpN,Gen->NPCh);
    if( POS(r1.m) < .5 || POS(r2.m) < .5)
      Ch[c].Type |= CHAIN_TILTED;
  }
  free(Interp1);
  free(Interp2);
  free(InterpN);
  return 0;
}
void VarData::ChangeNChain(int NChain,int nBlock){
  int NChainTemp = 0;
  Block[nBlock].NChain = NChain;
  Block[nBlock].EndIdx = Block[nBlock].InitIdx + Block[nBlock].NChain*Block[nBlock].NPCh;
  NChainTemp = Block[0].NChain;
  for(int b=1;b<Gen->NBlock;b++){
    Block[b].InitIdx = Block[b-1].EndIdx;
    Block[b].EndIdx = Block[b].InitIdx + Block[b].NChain*Block[b].NPCh;
    NChainTemp += Block[b].NChain;
  }
  Gen->NChain = NChainTemp;
}
void VarData::SwapChain(int c1,int c2){
  SwapChain(c1,c2,0);
}
void VarData::SwapChain(int c1,int c2,int b){
  int p1 = Block[b].InitIdx+c1*Block[b].NPCh;
  int p2 = Block[b].InitIdx+c2*Block[b].NPCh;
  double Temp[3];
  for(int p=0;p<Block[b].NPCh;p++){
    for(int d=0;d<3;d++){
      Temp[d] = Pm[p+p1].Pos[d];
      Pm[p+p1].Pos[d] = Pm[p+p2].Pos[d];
      Pm[p+p2].Pos[d] = Temp[d];
      Temp[d] = Pm[p+p1].Vel[d];
      Pm[p+p1].Vel[d] = Pm[p+p2].Vel[d];
      Pm[p+p2].Vel[d] = Temp[d];
    }
  }
  for(int d=0;d<3;d++){
    Temp[d] = Ch[c1].Pos[d];
    Ch[c1].Pos[d] = Ch[c2].Pos[d];
    Ch[c2].Pos[d] = Temp[d];
    Temp[d] = Ch[c1].Vel[d];
    Ch[c1].Vel[d] = Ch[c2].Vel[d];
    Ch[c2].Vel[d] = Temp[d];
  }
}
void VarData::SwapPart(int p1,int p2){
  double Temp[3];
  for(int d=0;d<3;d++){
    Temp[d] = Pm[p1].Pos[d];
    Pm[p1].Pos[d] = Pm[p2].Pos[d];
    Pm[p2].Pos[d] = Temp[d];
    Temp[d] = Pm[p1].Vel[d];
    Pm[p1].Vel[d] = Pm[p2].Vel[d];
    Pm[p2].Vel[d] = Temp[d];
  }
  int tmp = Ln[p1].NLink;
  Ln[p1].NLink = Ln[p2].NLink;
  Ln[p2].NLink = tmp;
  for(int l=0;l<MAX(Ln[p1].NLink,Ln[p2].NLink);l++){
    tmp = Ln[p2].Link[l];
    Ln[p2].Link[l] = Ln[p1].Link[l];
    Ln[p1].Link[l] = tmp;
  }
  tmp = Pm[p1].CId;
  Pm[p1].CId = Pm[p2].CId;
  Pm[p2].CId = tmp;
  tmp = Pm[p1].Typ;
  Pm[p1].Typ = Pm[p2].Typ;
  Pm[p2].Typ = tmp;
  tmp = Pm[p1].Idx;
  Pm[p1].Idx = Pm[p2].Idx;
  Pm[p2].Idx = tmp;
}
double VarData::TwoPartDist2(int p1,int p2,double *DistRel){
  return TwoPartDist2(Pm[p1].Pos,p2,DistRel);
}
double VarData::TwoPartDist2(double *Pos,int p2,double *DistRel){
  for(int d=0;d<3;d++){
    DistRel[d] = Pos[d] - Pm[p2].Pos[d];
    DistRel[d] -= floor(DistRel[d]*pInvEdge(d) + .5)*pEdge(d);
  }
  return DistRel[3] = (SQR(DistRel[0])+SQR(DistRel[1])+SQR(DistRel[2]));
}
double VarData::TwoPartDist(double *Pos,int p2,double *DistRel){
  return DistRel[3] = sqrt(TwoPartDist2(Pos,p2,DistRel));
}
double VarData::TwoPartDist(int p1,int p2,double *DistRel){
  return TwoPartDist(Pm[p1].Pos,p2,DistRel);
}
int VarData::TwoPartDist(int p1,int p2,double *DistRel,double CutOff){
  return TwoPartDist(Pm[p1].Pos,p2,DistRel,CutOff);
}
int VarData::TwoPartDist(double *Pos,int p2,double *DistRel,double CutOff){
  for(int d=0;d<3;d++){
    DistRel[d] = Pos[d] - Pm[p2].Pos[d];
    DistRel[d] -= floor(DistRel[d]*pInvEdge(d) + .5)*pEdge(d);
  }
  DistRel[3] = (SQR(DistRel[0])+SQR(DistRel[1])+SQR(DistRel[2]));
  if(DistRel[3] > SQR(CutOff)) return 0;
  DistRel[3] = sqrt(DistRel[3]);
  return 1;
}
void VarData::ShiftBlock(Vettore *Shift,int b){
  for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
    if(!(Pm[p].CId%2)) continue;
    for(int d=0;d<3;d++){
      Pm[p].Pos[d] += Shift->Val(d);
      Pm[p].Pos[d] -= floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);
    }
  }
}
void VarData::RotateBlock(Vettore *Axis,Vettore *Origin,int b){
  double Norm = Axis->Normalize();
  if(Norm <= 0.){printf("Cannot rotate block, the axis is null\n"); return;};
  Vettore Proj(3);
  Vettore Pos(3);
  for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
    for(int d=0;d<3;d++){
      Pos.Set(Pm[p].Pos[d]-Origin->Val(d),d);
    }
    Proj.Copy(&Pos);
    //Proj.ProjOnAxis(Axis);
    Proj.PerpTo(&Pos,Axis);
    for(int d=0;d<3;d++){
      Pm[p].Pos[d] = Origin->Val(d) + Pos[d] + 2.*Proj[d];
      Pm[p].Pos[d] -= floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);
    }
  }
}
void VarData::MirrorBlock(Vettore *S1,Vettore *S2,Vettore *S3,int b){
  Vettore PS(3);
  Vettore Pos(3);
  for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
    for(int d=0;d<3;d++){
      Pos.Set(Pm[p].Pos[d],0);
    }
    PS.ProjOnSurf(S1,S2,S3,&Pos);
    for(int d=0;d<3;d++){
      Pm[p].Pos[d] = PS.Val(d);
    }
  }  
}
void VarData::Transform(int b){
  if(b>=pNBlock())return;
  Vettore Axis(1.,0.,1.);
  Axis.Normalize();
  Vettore Shift(0.,0.,-.6*pEdge(2));
  Vettore Px1(.5*pEdge(0),.5*pEdge(1),.1*pEdge(2));
  Vettore Px2(.5*pEdge(0),.5*pEdge(1),.3*pEdge(2));
  Vettore Px3(.5*pEdge(0),.5*pEdge(1),.7*pEdge(2));
  int nNano = b-1;
  Vettore Origin(pNanoPos(nNano,0),pNanoPos(nNano,1),pNanoPos(nNano,2));
  Origin.Print();
  //ShiftBlock(&Shift,b);
  RotateBlock(&Axis,&Origin,b);
  //MirrorBlock(&Px1,&Px2,&Px3,b);
}
void VarData::Point2Shape(int iShape){
  if(VAR_IF_TYPE(iShape,SHAPE_NONE))
    Nano_Dist = &VarData::FieldNo;
  else if(VAR_IF_TYPE(iShape,SHAPE_SPH))
    Nano_Dist = &VarData::FieldSphere;
  else if(VAR_IF_TYPE(iShape,SHAPE_CYL))
    //Nano_Dist = &VarData::FieldCyl;
    Nano_Dist = &VarData::FieldTransMem;
  else if(VAR_IF_TYPE(iShape,SHAPE_CLUSTER))
    Nano_Dist = &VarData::FieldCyl;
  else if(VAR_IF_TYPE(iShape,SHAPE_PILL))
    Nano_Dist = &VarData::FieldCyl;
  else if(VAR_IF_TYPE(iShape,SHAPE_TILT))
    //Nano_Dist = &VarData::FieldSphere;
    Nano_Dist = &VarData::FieldTransMem;
  else if(VAR_IF_TYPE(iShape,SHAPE_WALL))
    Nano_Dist = &VarData::FieldTiltWall;
  else if(VAR_IF_TYPE(iShape,SHAPE_PORE))
    Nano_Dist = &VarData::FieldSphere;
  else if(VAR_IF_TYPE(iShape,SHAPE_EXT)) 
    Nano_Dist = &VarData::FieldSphere;
  else if(VAR_IF_TYPE(iShape,SHAPE_JANUS))
    Nano_Dist = &VarData::FieldJanus;
  else if(VAR_IF_TYPE(iShape,SHAPE_STALK))
    Nano_Dist = &VarData::FieldTorus;
  else if(VAR_IF_TYPE(iShape,SHAPE_TIP))
    //Nano_Dist = &VarData::FieldParab;
    Nano_Dist = &VarData::FieldElips;
  else if(VAR_IF_TYPE(iShape,SHAPE_TORUS))
    Nano_Dist = &VarData::FieldTorus;
  else if(VAR_IF_TYPE(iShape,SHAPE_HARM))
    Nano_Dist = &VarData::FieldTiltWall;
  else if(VAR_IF_TYPE(iShape,SHAPE_UMBR))
    Nano_Dist = &VarData::FieldTiltWall;
  else if(VAR_IF_TYPE(iShape,SHAPE_BOUND))
    Nano_Dist = &VarData::FieldBound;
  else{
    Nano_Dist = &VarData::FieldNo;
    //printf("No distance fuction assigned to the shape of the nano shape %d\n",iShape);
  }
}
double VarData::NanoDist2(double x,double y,double z,int n){
  double Pos[3] = {x,y,z};
  return NanoDist2(Pos,n);
}
double VarData::FieldNo(double *Pos,int n){
  return 100000.;
}
double VarData::FieldSphere(double *Pos,int n){
  return (SQR(pNanoPos(n,0)-Pos[0])+SQR(pNanoPos(n,1)-Pos[1])+SQR(pNanoPos(n,2)-Pos[2]));
}
double VarData::FieldElips(double *Pos,int n){
  double PosRel[3];
  for(int d=0;d<3;d++){
    PosRel[d] = pNanoPos(n,d) - Pos[d];
    PosRel[d] -= floor(PosRel[d]*pInvEdge(d) + .5)*pEdge(d);
  }
  double Ratio = Nano[n].Height>0.?Nano[n].Rad/Nano[n].Height:.5;
  double Lat = SQR(PosRel[0])+SQR(PosRel[1]);
  double Norm = SQR(Ratio*(PosRel[2]));
  return (Lat+Norm);
}
double VarData::FieldParab(double *Pos,int n){
  double PosRel[3];
  for(int d=0;d<3;d++){
    PosRel[d] = pNanoPos(n,d) - Pos[d];
    //PosRel[d] -= floor(PosRel[d]*pInvEdge(d) + .5)*pEdge(d);
  }
  double r2 = SQR(PosRel[0])+SQR(PosRel[1]);
  return SQR(-4.*Nano[n].Height*Pos[CNorm]+r2+3.*Nano[n].Height+4.*Nano[n].Height*Nano[n].Pos[CNorm]);
  double r = sqrt(r2);
  double z = Pos[CNorm];//Nano[n].Pos[CNorm] - Pos[CNorm];
  double a = -4.*Nano[n].Height;
  double b = 1.;
  double c = -2.*0.;
  double d = 3.*Nano[n].Height+4.*Nano[n].Height*Nano[n].Pos[CNorm];
  return SQR(a*z + b*r2 + c*r + d);
}
double VarData::FieldCyl(double *Pos,int n){
  double Dist2 = SQR(pNanoPos(n,0)-Pos[0])+SQR(pNanoPos(n,1)-Pos[1]);
  if(Pos[2] > Nano[n].Height*.5+pNanoPos(n,2))
    Dist2 += SQR(Pos[2]-Nano[n].Height*.5-pNanoPos(n,2));
  else if(Pos[2] < -Nano[n].Height*.5+Nano[n].Pos[2])
    Dist2 += SQR(Pos[2]+Nano[n].Height*.5+pNanoPos(n,2));
  return (Dist2);
}
double VarData::FieldTransMem(double *Pos,int n){
  double PosRel[3];
  for(int d=0;d<3;d++){
    PosRel[d] = pNanoPos(n,d) - Pos[d];
    PosRel[d] -= floor(PosRel[d]*pInvEdge(d) + .5)*pEdge(d);
  }
  // double RadLat = SQR(PosRel[CLat1]) + SQR(PosRel[CLat2]);
  // if(PosRel[CNorm] > .5*Nano[n].Height + 1.) return 500.;
  // if(PosRel[CNorm] < -.5*Nano[n].Height - 1.) return 500.;
  // if(PosRel[CNorm] > .5*Nano[n].Height){
  //   PosRel[CNorm] -= .5*Nano[n].Height - .5*Nano[n].Rad;
  //   if(SQR(PosRel[CNorm]) > RadLat )
  //     return - RadLat + SQR(PosRel[CNorm]);
  //   else
  //     return RadLat - SQR(PosRel[CNorm]);
  // }
  // if(PosRel[CNorm] < -.5*Nano[n].Height){
  //   PosRel[CNorm] += .5*Nano[n].Height - .5*Nano[n].Rad;
  //   if(SQR(PosRel[CNorm]) > RadLat )
  //     return - RadLat + SQR(PosRel[CNorm]);
  //   else
  //     return RadLat - SQR(PosRel[CNorm]);
  // }
  // PosRel[CNorm] = 0.;
  // return SQR(PosRel[CLat1]) + SQR(PosRel[CLat2]) + SQR(PosRel[CNorm]);
  // if and how much above the cylinder (for the smoothing)
  if(PosRel[CNorm] > Nano[n].Height*.5){
    PosRel[CNorm] = (PosRel[CNorm] - (Nano[n].Height*.5+Nano[n].Rad));
    for(int d=0;d<3;d++)
      PosRel[d]  *= .65;//*Nano[n].Rad;
    if(PosRel[CNorm] > 0.) PosRel[CNorm] = 100.;
  }
  else if(PosRel[CNorm] < -Nano[n].Height*.5){
    PosRel[CNorm] = (PosRel[CNorm] + (Nano[n].Height*.5+Nano[n].Rad));
    for(int d=0;d<3;d++)
      PosRel[d]  *= .65;//*Nano[n].Rad;
    if(PosRel[CNorm] < 0.) PosRel[CNorm] = 100.;
  }
  else{
    PosRel[CNorm]  = 0.;
  }
  double Dist2 = SQR(PosRel[0]) + SQR(PosRel[1]) + SQR(PosRel[2]);
  return (Dist2);
}
double VarData::FieldTorus(double *Pos,int n){
  double PosRel[3];
  for(int d=0;d<3;d++){
    PosRel[d] = pNanoPos(n,d) - Pos[d];
    PosRel[d] -= floor(PosRel[d]*pInvEdge(d) + .5)*pEdge(d);
  }
  double Radxy = sqrt(SQR(PosRel[CLat1]) + SQR(PosRel[CLat2]));
  double Temp = SQR(Nano[n].Height - Radxy) - SQR(Nano[n].Rad) + SQR(PosRel[CNorm]);
  return SQR(Temp);
}
double VarData::FieldTilt(double *Pos,int n){
  Vettore NanoAxis(Nano[n].Axis[0],Nano[n].Axis[1],Nano[n].Axis[2]);
  Vettore PosRel(3);
  Vettore Dist(3);
  double dr[3];
  for(int d=0;d<3;d++){
    PosRel.Set(pNanoPos(n,d) - Pos[d],d);
  }
  double r2 = fabs(Dist.PerpTo3(&PosRel,&NanoAxis));
  double HeiOnAxis = PosRel.ProjOnAxis(&NanoAxis);
  if(fabs(HeiOnAxis) > .5*Nano[n].Height){
    double Sign = HeiOnAxis > 0. ? 1. : -1.;
    r2 = 0.;
    for(int d=0;d<3;d++){
      dr[d] = - (Sign*Nano[n].Height*.5*Nano[n].Axis[d]) + PosRel[d];
      r2 += SQR(dr[d]);
    }
  }
  return r2;
}
double VarData::FieldTiltWall(double *Pos,int n){
  double a = Nano[n].Axis[0];
  double b = Nano[n].Axis[1];
  double c = Nano[n].Axis[2];
  double d = -Nano[n].Axis[0]*Nano[n].Pos[0] - Nano[n].Axis[1]*Nano[n].Pos[1] - Nano[n].Axis[2]*Nano[n].Pos[2];
  double Dist2 = a*Pos[0]+b*Pos[1]+c*Pos[2]+d;
  Dist2 = SQR(Dist2)/(SQR(a)+SQR(b)+SQR(c));
  return Dist2;
}
double VarData::FieldBound(double *Pos,int n){
  double dr[3] = {Pos[0],Pos[1],Pos[2]};
  dr[CNorm] = 0.;
  if(Pos[CLat1] > .5*pEdge(CLat1)) dr[CLat1] = pEdge(CLat1) - Pos[CLat1];
  if(Pos[CLat2] > .5*pEdge(CLat2)) dr[CLat2] = pEdge(CLat2) - Pos[CLat2];
  if(dr[CLat1] < dr[CLat2] ) dr[CLat2] = 0.;
  else dr[CLat1] = 0.;
  return SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]);
  // return MIN(SQR(dr[CLat1]),SQR(dr[CLat2]));
}
double VarData::FieldJanus(double *Pos,int n){
  Vettore NanoAxis(Nano[n].Axis[0],Nano[n].Axis[1],Nano[n].Axis[2]);
  Vettore PosRel(3);
  Vettore Dist(3);
  for(int d=0;d<3;d++){
    PosRel.Set(pNanoPos(n,d) - Pos[d],d);
  }
  PosRel.Set(PosRel.Val(2)*2.,2);
  double r2 = fabs(Dist.PerpTo3(&PosRel,&NanoAxis));
  return SQR(r2);
}

// # include "../Matematica/cvt.H"

// int VarData::Voronoi(){
//   int NDim = 2;
//   int NCell = Gen->NChain;
//   int Batch = 1000;
//   int DontCreate = 4;
//   int SampleType = 1;//Halton
//   int NSample = 10000;
//   int NItMax = 40;
//   int NItFix = 1;
//   int Seme = 45322643;
//   int NIt = 0;
//   double ItDiff = 0.;
//   double Energy = 0.;
//   double *Pos = (double *)calloc(NCell*NDim,sizeof(double));
//   for(int c=0;c<NCell;c++){
//     for(int d=0;d<NDim;d++){
//       Pos[c*NDim+d] = Ch[c].Pos[d];
//     }
//   }
//   cvt(NDim,NCell,Batch,DontCreate,SampleType,NSample,NItMax,NItFix,&Seme,Pos,&NIt,&ItDiff,&Energy);
// }

// #ifndef VOROPP
// #include "../share/voro++/src/voro++.cc"
// #include "../../share/include/container.hh"

// int VarData::Voronoi(){
//   double xMin = 0.;
//   double xMax = Gen->Edge[CLat1];
//   double yMin = 0.;
//   double yMax = Gen->Edge[CLat2];
//   double zMin = .4*Gen->Edge[CNorm];
//   double zMax = .6*Gen->Edge[CNorm];
//   double LengthScale = 1./(xMax*2.6);
//   int nxf = (int)( (xMax-xMin)*LengthScale ) + 1;
//   int nyf = (int)( (yMax-yMin)*LengthScale ) + 1;
//   int nzf = (int)( (zMax-zMin)*LengthScale ) + 1;
//   bool xperiodic = true;
//   bool yperiodic = true;
//   bool zperiodic = false;
//   const int memory=8;
//   container con(xMin,xMax,yMin,yMax,zMin,zMax,nxf,nyf,nzf,xperiodic,yperiodic,zperiodic,memory);
//   for(int c=0;c<Gen->NChain;c++){
//     if(!CHAIN_IF_TYPE(Ch[c].Type,NChType)) continue;
//     con.put(0,Ch[c].Pos[CLat1],Ch[c].Pos[CLat2],Ch[c].Pos[CNorm]);
//   }
//   //con.draw_cells_gnuplot("random_points_v.gnu");
//   //con.draw_cells_pov("pack_six_cube_v.pov");
//   //con.draw_particles_pov("pack_six_cube_p.pov");
//   return 0;
// }
// #else
// int VarData::Voronoi(){
//   printf("Voronoi library not supplied\n");
//   return 1;
// }
// #endif
