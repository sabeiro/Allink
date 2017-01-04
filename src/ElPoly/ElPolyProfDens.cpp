/***********************************************************************
ElPolyProfDens: 
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
#include "ElPoly.h"
/**
a 
*/
int ElPoly::DensProf(int NBin,int NSample,int Coord){
  int NType = 3;
  double *DensityAv = (double *)calloc(NType*NBin,sizeof(double));
  double FreeE = 0.;
  double ChDensMean = 0.;
  double ChDensSqr  = 0.;
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  double *VolEl = (double *)calloc(NBin,sizeof(double));
  VolumeCircSlab(VolEl,NBin);
  if(Coord > 3 || Coord < 0) return 1;
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    //printf("Opening: %s\n",cFile[f]);
    if(Coord < 3)
      if(OpenRisk(cFile[f],BF_PART))return 0;
      else if(Coord == 3){
	if(OpenRisk(cFile[f],BF_NANO))return 0;
      }
    double RefZed = .5*pCm(CNorm);
    if(pNNano() > 0) RefZed = pNanoPos(0,CNorm);
    // if(Coord == 2 && VAR_IF_TYPE(SysType,VAR_MEMBRANE) )
    //   ShiftSys(SHIFT_CM);
    /// sum on small patches and shift wrt the weighted average
    if(Coord < 3){
      //DensityProfile(Coord,NBin,NType,DensityAv);
      BfDefChain();
      int NPartition = NSample;
      double *Density = (double *)calloc(NPartition*NPartition*NType*NBin,sizeof(double));
      double Pos[4];
      char cLine[STRSIZE];
      for(int b=0;b<pNBlock();b++){
	if(!strncmp(Block[b].Name,"PEP",3) ) continue;
	if(!strcmp(Block[b].Name,"WATER") )  continue;
	for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
	  int t = pType(p);
	  for(int d=0;d<3;d++) Pos[d] = pPos(p,d);
	  Pos[CNorm] = pPos(p,CNorm) - RefZed + .5*pEdge(CNorm);
	  int vx = (int)(NPartition*Pos[CLat1]*pInvEdge(CLat2));
	  if(vx < 0 || vx >= NBin) continue;
	  int vy = (int)(NPartition*Pos[CLat2]*pInvEdge(CLat2));
	  if(vy < 0 || vy >= NBin) continue;
	  int vz = (int)(NBin*Pos[CNorm]*pInvEdge(CNorm));
	  if(vz < 0 || vz >= NBin) continue;
	  if(CHAIN_IF_TYPE(Ch[pChain(p)].Type,CHAIN_ADDED)) t = 2;
	  if(!strncmp(Block[b].Name,"ADDED",5)) t = 2;
	  if(!strncmp(Block[b].Name,"CHOL",4))  t = 2;
	  Density[((vz*NType+t)*NPartition+vy)*NPartition+vx] += 1.;
	}
      }
      for(int vx=0;vx<NPartition;vx++){
	for(int vy=0;vy<NPartition;vy++){
	  double Media = 0.;
	  double Weight = 0.;
	  for(int v=0;v<NBin;v++){
	    Weight += Density[((v*NType+0)*NPartition+vy)*NPartition+vx];
	    Media += v*Density[((v*NType+0)*NPartition+vy)*NPartition+vx];
	  }
	  int vMedia = Weight > 0 ? (int)(Media/Weight) : NBin/2;
	  int vOffset = vMedia - NBin/2;
	  for(int t=0;t<NType;t++){
	    for(int v=0;v<NBin;v++){
	      int vOther = v + vOffset;
	      if(vOther < 0)
		vOther = NBin - vOther;
	      if(vOther >= NBin)
		vOther = vOther - NBin;
	      DensityAv[v*NType+t] += Density[((vOther*NType+t)*NPartition+vy)*NPartition+vx];
	    }
	  }
	}
      }
      free(Density);
    }
    else{
      Vettore NanoAxis(Nano[0].Axis[0],Nano[0].Axis[1],Nano[0].Axis[2]);
      Vettore PosRel(3);
      Vettore Dist(3);
      for(int b=0;b<pNBlock();b++){
	if(!strncmp(Block[b].Name,"PEP",3) ) continue;
	for(int p=0;p<pNPart();p++){
	  for(int d=0;d<3;d++){
	    PosRel.Set(remainder(pNanoPos(0,d)-pPos(p,d),pEdge(d)),d);
	  }
	  double RadDist = ABS(Dist.PerpTo3(&PosRel,&NanoAxis));
	  int v = (int)(RadDist/pEdge(3)*NBin);
	  if(v < 0 || v >= NBin) continue;
	  int t = pType(p);
	  if(CHAIN_IF_TYPE(Ch[pChain(p)].Type,CHAIN_ADDED)) t = 2;
	  DensityAv[v*NType+t] += 1.;
	}
      }
    }
    ChDensMean += pNChain()/(pEdge(CLat1)*pEdge(CLat2));
    ChDensSqr  += SQR(pNChain()/(pEdge(CLat1)*pEdge(CLat2)));
  }
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,DensityAv,NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    printf("\n");
    double NormaF = (NFile[1] - NFile[0]);
    double VolumeSlab = (pEdge(0)*pEdge(1)*pEdge(2))/(double)NBin;
    for(int v=0;v<NBin;v++){
      for(int t=0;t<NType;t++){
	if(Coord < 3) 
	  DensityAv[v*NType+t] /= VolumeSlab*NormaF;
	else 
	  DensityAv[v*NType+t] /= NormaF*VolEl[v];
      }
    }
    char FileName[256];
    sprintf(FileName,"DensityChi%.0fKappa%.0fRho%.0f.dat",pchiN(),pkappaN(),prho());
    FILE *FileToWrite = fopen(FileName,"w");
    char cSystem[STRSIZE];
    SysDef(cSystem);
    double ChDensErr = sqrt(ChDensSqr - SQR(ChDensMean)/NormaF )/NormaF;
    fprintf(FileToWrite,"# %s\n",cSystem);
    fprintf(FileToWrite,"# ChainDens %lf %lf\n",ChDensMean/NormaF,ChDensErr);
    for(int v=0;v<NBin;v++){
      fprintf(FileToWrite,"%lf %lf %lf %lf\n",(double)v/(double)NBin*pEdge(Coord),DensityAv[v*NType],DensityAv[v*NType+1],DensityAv[v*NType+2]);
    }
    fclose(FileToWrite);
#ifdef OMPI_MPI_H
  }
#endif
  free(DensityAv);
  return 0;
}
/**
a  to be completed
*/
void ElPoly::DensProfNormalSlab(int NBin,int NSample,int Coord){
  int NType = 3;
  double *DensityAv = (double *)calloc(NType*NBin,sizeof(double));
  double FreeE = 0.;
  double ChDensMean = 0.;
  double ChDensSqr  = 0.;
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  double *VolEl = (double *)calloc(NBin,sizeof(double));
  VolumeCircSlab(VolEl,NBin);
  if(Coord > 3 || Coord < 0) return;
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    //printf("Opening: %s\n",cFile[f]);
    if(Coord < 3)
      if(OpenRisk(cFile[f],BF_PART))return;
      else if(Coord == 3){
	if(OpenRisk(cFile[f],BF_NANO))return;
      }
    double RefZed = .5*pCm(CNorm);
    if(pNNano() > 0) RefZed = pNanoPos(0,CNorm);
    // if(Coord == 2 && VAR_IF_TYPE(SysType,VAR_MEMBRANE) )
    //   ShiftSys(SHIFT_CM);
    /// sum on small patches and shift wrt the weighted average
    if(Coord < 3){
      //DensityProfile(Coord,NBin,NType,DensityAv);
      BfDefChain();
      int NPartition = NSample;
      double *Density = (double *)calloc(NPartition*NPartition*NType*NBin,sizeof(double));
      double Pos[4];
      char cLine[STRSIZE];
      for(int b=0;b<pNBlock();b++){
	if(!strncmp(Block[b].Name,"PEP",3) ) continue;
	if(!strcmp(Block[b].Name,"WATER") )  continue;
	for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
	  int t = pType(p);
	  for(int d=0;d<3;d++)
	    Pos[d] = pPos(p,d);
	  Pos[CNorm] = pPos(p,CNorm) - RefZed + .5*pEdge(CNorm);
	  int vx = (int)(NPartition*Pos[CLat1]*pInvEdge(CLat2));
	  if(vx < 0 || vx >= NBin) continue;
	  int vy = (int)(NPartition*Pos[CLat2]*pInvEdge(CLat2));
	  if(vy < 0 || vy >= NBin) continue;
	  int vz = (int)(NBin*Pos[CNorm]*pInvEdge(CNorm));
	  if(vz < 0 || vz >= NBin) continue;
	  if(CHAIN_IF_TYPE(Ch[pChain(p)].Type,CHAIN_ADDED)) t = 2;
	  if(!strncmp(Block[b].Name,"ADDED",5)) t = 2;
	  if(!strncmp(Block[b].Name,"CHOL",4))  t = 2;
	  Density[((vz*NType+t)*NPartition+vy)*NPartition+vx] += 1.;
	}
      }
      for(int vx=0;vx<NPartition;vx++){
	for(int vy=0;vy<NPartition;vy++){
	  double Media = 0.;
	  double Weight = 0.;
	  for(int v=0;v<NBin;v++){
	    Weight += Density[((v*NType+0)*NPartition+vy)*NPartition+vx];
	    Media += v*Density[((v*NType+0)*NPartition+vy)*NPartition+vx];
	  }
	  int vMedia = Weight > 0 ? (int)(Media/Weight) : NBin/2;
	  int vOffset = vMedia - NBin/2;
	  for(int t=0;t<NType;t++){
	    for(int v=0;v<NBin;v++){
	      int vOther = v + vOffset;
	      if(vOther < 0)
		vOther = NBin - vOther;
	      if(vOther >= NBin)
		vOther = vOther - NBin;
	      DensityAv[v*NType+t] += Density[((vOther*NType+t)*NPartition+vy)*NPartition+vx];
	    }
	  }
	}
      }
      free(Density);
    }
    else{
      Vettore NanoAxis(Nano[0].Axis[0],Nano[0].Axis[1],Nano[0].Axis[2]);
      Vettore PosRel(3);
      Vettore Dist(3);
      for(int b=0;b<pNBlock();b++){
	if(!strncmp(Block[b].Name,"PEP",3)) continue;
	for(int p=0;p<pNPart();p++){
	  for(int d=0;d<3;d++){
	    PosRel.Set(remainder(pNanoPos(0,d)-pPos(p,d),pEdge(d)),d);
	  }
	  double RadDist = ABS(Dist.PerpTo3(&PosRel,&NanoAxis));
	  int v = (int)(RadDist/pEdge(3)*NBin);
	  if(v < 0 || v >= NBin) continue;
	  int t = pType(p);
	  if(CHAIN_IF_TYPE(Ch[pChain(p)].Type,CHAIN_ADDED))
	    t = 2;
	  DensityAv[v*NType+t] += 1.;
	}
      }
    }
    ChDensMean += pNChain()/(pEdge(CLat1)*pEdge(CLat2));
    ChDensSqr  += SQR(pNChain()/(pEdge(CLat1)*pEdge(CLat2)));
  }
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,DensityAv,NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    printf("\n");
    double NormaF = (NFile[1] - NFile[0]);
    double VolumeSlab = (pEdge(0)*pEdge(1)*pEdge(2))/(double)NBin;
    for(int v=0;v<NBin;v++){
      for(int t=0;t<NType;t++){
	if(Coord < 3) 
	  DensityAv[v*NType+t] /= VolumeSlab*NormaF;
	else 
	  DensityAv[v*NType+t] /= NormaF*VolEl[v];
      }
    }
    char FileName[256];
    sprintf(FileName,"DensityChi%.0fKappa%.0fRho%.0f.dat",pchiN(),pkappaN(),prho());
    FILE *FileToWrite = fopen(FileName,"w");
    char cSystem[STRSIZE];
    SysDef(cSystem);
    double ChDensErr = sqrt(ChDensSqr - SQR(ChDensMean)/NormaF )/NormaF;
    fprintf(FileToWrite,"# %s\n",cSystem);
    fprintf(FileToWrite,"# ChainDens %lf %lf\n",ChDensMean/NormaF,ChDensErr);
    for(int v=0;v<NBin;v++){
      fprintf(FileToWrite,"%lf %lf %lf %lf\n",(double)v/(double)NBin*pEdge(Coord),DensityAv[v*NType],DensityAv[v*NType+1],DensityAv[v*NType+2]);
    }
    fclose(FileToWrite);
#ifdef OMPI_MPI_H
  }
#endif
  free(DensityAv);
  return;
}
/**
a 
*/
void ElPoly::SlabProf(int NBin,int nNano,int Coord1){
  if(Coord1 != CLat1 && Coord1 != CLat2){
    printf("%d is not a lateral coordinate\n",Coord1,CLat1,CLat2);
  }
  //SigErr(pNNano() < 2,"The number of nano should be at least two");
  int Coord2 = Coord1 == CLat1 ? CLat2 : CLat1;
  double *Dens = (double *) calloc(NBin,sizeof(double));
  double *Temp = (double *) calloc(NBin,sizeof(double));
  double *Prof = (double *) calloc(6*NBin,sizeof(double));
  double **Plot = (double **)calloc(3,sizeof(double));
  for(int t=0;t<3;t++){
    Plot[t] = (double *)calloc(NBin*NBin,sizeof(double));
  }
  double *Count = (double *) calloc(6*NBin,sizeof(double));
  Vettore Pos(3);
  Vettore PosRel(3);
  //FILE *Ciccia = fopen("Ciccia.dat","w");
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    //--------------open,-back-fold----------------------
    if(OpenRisk(cFile[f],BF_NO)) return;
    if(pNNano() == 1){
      SetNNano(2);
      Nano[1].Pos[0] = 0.;
      Nano[1].Pos[1] = 0.;
      Nano[1].Pos[2] = .5*pEdge(CNorm);
    }
    Vettore Axis(3);
    for(int d=0;d<3;d++){
      Axis.Set(Nano[1].Pos[d]-Nano[0].Pos[d],d);
    }
    Axis.Set(0.,CNorm);
    Axis.Normalize();
    double InvEdge = (pEdge(CLat1)*Axis.Val(0) + pEdge(CLat2)*Axis.Val(1))/(Axis.Val(0) + Axis.Val(1));
    InvEdge = 1./InvEdge;
    //Nano[nNano].Pos[CLat1] = 0.;
    //Nano[nNano].Pos[CLat2] = 0.;
    for(int b=0;b<pNBlock();b++){
      int TypePlus = 0;
      if(!strncmp(Block[b].Name,"PEP",3) ) TypePlus = 1;
      if(!strcmp(Block[b].Name,"WATER") )  continue;
      for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
	int t = 2*pType(p);
	for(int d=0;d<3;d++){
	  double Dist = Pm[p].Pos[d] - Nano[nNano].Pos[d];
	  Dist -= floor(Dist*pInvEdge(d)+.5)*pEdge(d);
	  Pos.Set(Dist,d);
	}
	Pos.Set(0.,CNorm);
	double LatDist = PosRel.PerpTo(&Pos,&Axis);
	if(LatDist > Nano[nNano].Rad) continue;
	double Dist = Pos.ProjOnAxis(&Axis);
	Pos.Set(pPos(p,CNorm)-Nano[nNano].Pos[CNorm],CNorm);
	int t1 = TypePlus ? 2 : pType(p);
	//fprintf(Ciccia,"%lf %lf %lf %d\n",pPos(p,0),pPos(p,1),pPos(p,2),t1);
	//fprintf(Ciccia,"%lf %lf %lf %d\n",Pos.Val(0)+Nano[nNano].Pos[0]-5.,Pos.Val(1)+Nano[nNano].Pos[1]+5.,Pos.Val(2)+Nano[nNano].Pos[2],t1);
	if(pPos(p,CNorm) > pNanoPos(nNano,CNorm)) t += 1;
	int v = (int)(Dist*InvEdge*NBin);
	if(v < 0 || v >= NBin) continue;
	int z = (int)(pPos(p,CNorm)*pInvEdge(CNorm)*NBin);
	if(z < 0 || z >= NBin) continue;
	Prof[v*6+t] += pPos(p,CNorm) - pCm(CNorm);
	Count[v*6+t] += 1.;
	Plot[t1][v*NBin+z] += 1.;
	if(!TypePlus) Dens[v] += 1.;
      }
    }
  }
  //fclose(Ciccia);
  printf("\n");
  //normalize
  for(int b=0;b<NBin;b++){
    for(int t=0;t<6;t++){
      Prof[b*6+t] /= Count[b*6+t] > 0. ? Count[b*6+t] : 1.;
    }
  }
  double Norm = 1./(NFile[1]-NFile[0]);
  for(int v=0;v<SQR(NBin);v++){
    for(int t=0;t<3;t++){
      Plot[t][v] *= Norm;
    }
  }
  Matrice Mask1(5);
  Mask1.FillGaussian(.5,3.);
  Mask1.Print();
  int NDim = 1;
  int IfMinImConv = 1;
  Mask1.ConvoluteMatrix(Dens,NBin,NDim,IfMinImConv);
  //write profile
  FILE *FSlab = fopen("SlabProfile.dat","w");
  FILE *FDens = fopen("SlabDensProfile.dat","w");
  double InvNFile = 1./(double)(NFile[1]-NFile[0]);
  for(int b=0;b<NBin;b++){
    double x = b/(double)NBin*pEdge(Coord1);
    //fprintf(FSlab,"%lf %lf %lf %lf %lf %lf %lf\n",x,Prof[b*6  ],Prof[b*6+1],Prof[b*6+2],Prof[b*6+3],Prof[b*6+4],Prof[b*6+5],Prof[b*6+6]);
    fprintf(FSlab,"%lf %lf %lf %lf %lf \n",x,Prof[b*6  ],Prof[b*6+1],Prof[b*6+2],Prof[b*6+3],Prof[b*6+4]);
    fprintf(FDens,"%lf %lf\n",x,Dens[b]*InvNFile);
  }
  fclose(FSlab);
  fclose(FDens);
  //smooth and write dens
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  NDim = 2;
  IfMinImConv = 1;
  for(int t=0;t<3;t++){
    Mask.ConvoluteMatrix(Plot[t],NBin,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot[t],NBin,NDim,IfMinImConv);
  }
  FILE *FSlabProf = fopen("SlabDens.xzd","w");
  double LatDim[3] = {pEdge(Coord1),pEdge(CNorm),30.};
  PrintDens(FSlabProf,Plot,LatDim,NBin);
  fclose(FSlabProf);
  free(Prof);
  free(Count);
  free(Temp);
  free(Dens);
}
/**
a 
*/
void ElPoly::PrintDens(FILE *FileToWrite,double **Plot,double *LatDim,int NBin){
  int link[4] = {0,0,0,0};
  double InvNBin = 1./(double)NBin;
  fprintf(FileToWrite,"#l(%lf %lf %lf) v[%d] d[%s]\n",LatDim[0],LatDim[1],LatDim[2],NBin,ChooseDraw(EL_DENS));
  double Dr = InvNBin*LatDim[0];
  double Dz = InvNBin*LatDim[1];
  for(int t=0,p=0,c=0;t<3;t++){
    for(int vv=0;vv<NBin-1;vv++){
      for(int v=0;v<NBin-1;v++){
	link[0] = (v+0)*NBin + (vv+0);
	if(Plot[t][link[0]] < .02) continue;
	//fprintf(FileToWrite,"{t[%d %d %d] x(%lf %lf %lf)}\n",p,c,t, v*InvNBin*pEdge(3),(vv)*InvNBin*(Border[1]-Border[0]),(Plot[(v*NBin+vv)*NType+t]));
	//continue;
	//------------defines-the-squares---------------------
	link[1] = (v+0)*NBin + (vv+1);
	link[2] = (v+1)*NBin + (vv+0);
	link[3] = (v+1)*NBin + (vv+1);
	int pRef = p;
	for(int lx=0;lx<2;lx++){
	  for(int ly=0;ly<2;ly++){
	    int l = 2*lx+ly;
	    double NanoAdded = Plot[2][(link[l])];//+Plot[3][(link[l])];
	    double Phob = t == 0 ? Plot[0][(link[l])] : 0.;
	    double Phil = t == 1 ? Plot[1][(link[l])] : 0.;
	    double r = (v +lx)*Dr + .5*Dr;
	    double z = (vv+ly)*Dz + .5*Dz;
	    double dens = Plot[t][(link[l])];//+Plot[(v*NBin+vv)*NType+1]+Plot[(v*NBin+vv)*NType+2]);
	    int l1 = pRef + (p+1)%4;
	    int l2 = pRef + (p+2)%4;
	    int l3 = pRef + (p+3)%4;
	    fprintf(FileToWrite,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf) l[%d] l[%d] l[%d] }\n",p,c,t,r,z,dens,NanoAdded,Phob,Phil,l1,l2,l3);
	    p++;
	  }
	}
	c++;
      }
    }
  }
}
/**
a 
*/
void ElPoly::CartDens(int NBin,int nNano){
  int NType = 5;
  double **Plot = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<NType;t++){
    Plot[t] = (double *)calloc(NBin*NBin,sizeof(double));
  }
  double NBinInv = 1./(double)NBin;
  double InvNFile = 1./(double)(NFile[1]-NFile[0]);
  double LatDim[3] = {pEdge(CLat1),pEdge(CNorm),30.};
  char FName[20];
  double Width = 2.;
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    //--------------open,-back-fold----------------------
    if(OpenRisk(cFile[f],NBackFold)) return;
    for(int p=0;p<pNPart();p++){
      double Dx = fabs(pCm(CLat1) - pPos(p,CLat1));
      double Dy = fabs(pCm(CLat2) - pPos(p,CLat2));
      if(Dy > Width) continue;
      int vx = (int)(2.*Dx*pInvEdge(CLat1)*NBin);
      if(vx<0 || vx>=NBin) continue;
      int vz = (int)(pPos(p,CNorm)*pInvEdge(CNorm)*NBin);
      if(vz<0 || vz>=NBin) continue;
      int t = pType(p);
      //if((p%Block[0].Asym)!=0)continue;
      // if(VAR_IF_TYPE(Ch[Pm[p].CId].Type,CHAIN_INNER)) t = 0;
      // else if(VAR_IF_TYPE(Ch[Pm[p].CId].Type,CHAIN_OUTER)) t = 1;
      // else continue;
      Plot[t][vx*NBin+vz] += 1.;
    }
  }
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  for(int v=0;v<NBin*NBin;v++) for(int t=0;t<NType;t++) Plot[t][v] *= InvNFile;
  int NDim = 2;
  int IfMinImConv = 1;
  for(int t=0;t<NType;t++){
    Mask.ConvoluteMatrix(Plot[t],NBin,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot[t],NBin,NDim,IfMinImConv);
  }
  sprintf(FName,"CartDens.dat");
  FILE *FileToWrite = fopen(FName,"w");
  PrintDens(FileToWrite,Plot,LatDim,NBin);
  int NAvx = 4;
  int NAvz = 4;
  double NAv = 1./(double)(NAvx*NAvz);
  double Avx = 0.;
  double Avz = 0.;
  double Av = 0.;
  for(int vx=0;vx<NBin*NBin;vx++) Plot[3][vx] = 1.;
  for(int t=0;t<2;t++){
    int p = 0;
    sprintf(FName,"LineLayer%d.dat",t);
    FILE *FLine = fopen(FName,"w");
    if(1==1){
      fprintf(FLine,"#l(%lf %lf %lf) v[%d] d[part]\n",pEdge(0),pEdge(1),pEdge(2),NBin);
      for(int vx=0;vx<NBin;vx++){
	for(int vz=0;vz<NBin;vz++){
	  fprintf(FLine,"{t[%d %d %d] x(%lf %lf %lf)}\n",p++,0,t,(double)(vx*pEdge(CLat1)/(double)NBin),(double)(vz*pEdge(CNorm)/(double)NBin),Plot[t][vx*NBin+vz],t);
	}
      }
    }
    for(int vx=0;vx<NBin;vx+=NAvx){
      for(int vz=0;vz<NBin;vz+=NAvz){
	for(int vvx=0;vvx<NAvx;vvx++){
	  for(int vvz=0;vvz<NAvz;vvz++){
	    int vTot = (vx+vvx)*NBin + vz+vvz;
	    if(vTot >= SQR(NBin)) break;
	    Plot[3][vTot] *= Plot[t][vTot];
	    double Weight = Plot[t][vTot];
	    if(Weight < 1.0) continue;
	    Avx += (double)(vx*pEdge(CLat1)/(double)NBin)*Weight;
	    Avz += (double)(vz*pEdge(CNorm)/(double)NBin)*Weight;
	    Av  += Weight;
	  }
	}
	//fprintf(FLine,"%lf %lf %lf %d\n",Avx/Av,Avz/Av,Av/(double)NAverage,2);
	fprintf(FLine,"{t[%d %d %d] x(%lf %lf %lf) l[%d]}\n",p++,0,2,Avx/Av,Avz/Av,Av*NAv,p+1);
	Avx = 0.;
	Avz = 0.;
	Av  = 0.;
      }
    }
    fclose(FLine);
  }
  sprintf(FName,"LineLayer2.dat");
  FILE *FLine = fopen(FName,"w");
  for(int vx=0;vx<NBin;vx++){
    for(int vz=0;vz<NBin;vz++){
      double Weight = Plot[3][vx*NBin+vz];
      if(Weight < 1.1) continue;
      fprintf(FLine,"%lf %lf %lf\n",(double)(vx*pEdge(CLat1)/(double)NBin),(double)(vz*pEdge(CNorm)/(double)NBin),Weight);      
    }
  }
  for(int t=0;t<NType;t++) free(Plot[t]);
  free(Plot);
}
/**
a 
*/
void ElPoly::RadDistrF(int NBin,int How,int nNano){
  if(nNano > pNNano()){
    printf("The specified nanoparticle doesn't exist\n");
    return ;
  }
  int NType = 5;
  double *Dens = (double *)calloc(NBin,sizeof(double));
  double **Plot = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<NType;t++){
    Plot[t] = (double *)calloc(NBin*NBin,sizeof(double));
  }
  //--------The-first-file-defines-the-reference-distance--------
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  double RefPos[3] = {pCm(CLat1),pCm(CLat2),pCm(CNorm)};
  double LatDim[3] = {pEdge(3),pEdge(CNorm),30.};
  if(How >= 6){
    LatDim[0] = pEdge(CLat1);
    LatDim[1] = pEdge(CLat2);
  }
  double InvNBin  = 1./(double)NBin;
  double InvNFile = 1./(double)(NFile[1]-NFile[0]);
  double OldPos[3] = {pNanoPos(nNano,0),pNanoPos(nNano,1),pNanoPos(nNano,2)};
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    //--------------open,-back-fold----------------------
    if(OpenRisk(cFile[f],NBackFold)) return;
    BfDefChain();
    if(How == 5) PorePos();
    RefPos[CLat1] = pNanoPos(nNano,CLat1);
    RefPos[CLat2] = pNanoPos(nNano,CLat2);
    RefPos[CNorm] = pNanoPos(nNano,CNorm);
    if(How == 2){
      RefPos[CLat1] = pCm(CLat1);
      RefPos[CLat2] = pCm(CLat2);
      RefPos[CNorm] = pCm(CNorm);
    }
    if(How == 0){
      for(int d=0;d<3;d++){
	RefPos[d] = pCm(d);
      }
    }
    if(How < 6){
      //SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
      //----------define-the-distance,-fill-Plot---------------
      double TangInv = ABS(Nano[nNano].Axis[CNorm]) > 0. ? sqrt(QUAD(Nano[nNano].Axis[CLat1])+QUAD(Nano[nNano].Axis[CLat2]))/Nano[nNano].Axis[CNorm] : 0.;
      if(How == 0) TangInv = 0.;
      Vettore Rotated(3);
      Vettore Dist(3);
      Vettore PosRel(3);
      Vettore Zed(3);
      Zed.Set(0.,CLat1);Zed.Set(0.,CLat2);Zed.Set(1.,CNorm);
      //Vettore NanoAxis(1.,0.,0.);
      Vettore NanoAxis(Nano[nNano].Axis[0],Nano[nNano].Axis[1],Nano[nNano].Axis[2]);
      if(How == 0){
	NanoAxis.Set(0.,CLat1);
	NanoAxis.Set(0.,CLat2);
	NanoAxis.Set(1.,CNorm);
      }
      Vettore RotAxis(3);
      RotAxis.VetV(&NanoAxis,&Zed);
      double Angle = Zed.Angle(&NanoAxis,&Zed);
      Matrice M(RotAxis.x,Angle,3);
      int vzNano = (int)(pNanoPos(nNano,CNorm)*InvNBin*pInvEdge(CNorm));
      for(int b=0;b<pNBlock();b++){
	int IfAdded = 0;
	if(!strncmp(Block[b].Name,"PEP",3) ) IfAdded = 2;
	if(!strcmp(Block[b].Name,"STUFFING") ) continue;
	if(!strcmp(Block[b].Name,"WATER") ) continue;
	if(!strcmp(Block[b].Name,"CHOL") ) IfAdded = 1;
	if(!strcmp(Block[b].Name,"OIL") )  continue;//IfAdded = 1;
	if(!strcmp(Block[b].Name,"ADDED") )  continue;//IfAdded = 1;
	for(int p=Block[b].InitIdx,link=0;p<Block[b].EndIdx;p++){
	  // double RadCh = SQR(Ch[Pm[pPos(.CId,0)-Nano[nNano].PosBf[0])+SQR(Ch[Pm[pPos(.CId,1)-Nano[nNano].PosBf[1]);
	  // if( RadCh < SQR(Nano[nNano].Rad)){continue;}
	  //finding the position on the axis of the nanoparticle

	  int c = pChain(p);
	  int p1 = c*pNPCh(b);
	  //if(fabs(Pm[p1].Pos[CNorm] - Nano[0].Pos[CNorm]) < 2.5) continue;
	  if(fabs(Pm[p1].Pos[CNorm] - Nano[0].Pos[CNorm]) < 3.5) continue;

	  // if(!strcmp(Block[b].Name,"TT0"))
	  //   if(VAR_IF_TYPE(Ch[c].Type,CHAIN_INNER)) continue;
	  // if(!strcmp(Block[b].Name,"TT1"))
	  //   if(VAR_IF_TYPE(Ch[c].Type,CHAIN_OUTER)) continue;
	  // if(VAR_IF_TYPE(Ch[c].Type,CHAIN_OUTER)) continue;
	

	  for(int d=0;d<3;d++){
	    double Dist = pPos(p,d) - RefPos[d];
	    Dist -= floor( Dist*pInvEdge(d) + .5)*pEdge(d);
	    PosRel.Set(Dist,d);
	  }
	  double RadDist = fabs(Dist.PerpTo3(&PosRel,&NanoAxis));
	  // double RadDist = hypot(PosRel.Val(CLat1),PosRel.Val(CLat2));
	  //correspondent bin on the radial distance
	  int vr = (int)(NBin*RadDist*pInvEdge(3));
	  //Projecting the position on the axis and rotating to the normal
	  PosRel.ProjOnAxis(&NanoAxis);
	  M.Mult(PosRel,Rotated);
	  //double PosNorm = Rotated[CNorm] + .5*pEdge(CNorm) - floor(PosNorm*pInvEdge(CNorm))*pEdge(CNorm);
	  double PosNorm = PosRel.Val(CNorm) + .5*pEdge(CNorm);
	  int vz = (int)(NBin*(PosNorm)*pInvEdge(CNorm));
	  if( vz < 0 || vz >= NBin){
	    printf("Axis %lf %lf %lf \n",Nano[nNano].Axis[0],Nano[nNano].Axis[1],Nano[nNano].Axis[2]);
	    printf("Normal %lf > %lf %d > %d\n",PosNorm,pEdge(CNorm),vz,NBin);
	    continue;
	  }
	  if( vr < 0 || vr >= NBin){
	    continue;
	  }
	  //shifting wrt the inclination of the peptide
	  // v += (int)( (vv-NBin/2)*TangInv);
	  // if( v < 0) v += NBin;
	  // if( v >= NBin) v -= NBin;
	  int t = pType(p);
	  if(IfAdded==1){
	    t = 2;
	    Plot[0][vr*NBin+vz] += 1.;
	  }
	  if(IfAdded==2) t = 3;
	  else Dens[vr] += 1.;
	  Plot[t][vr*NBin+vz] += 1.;
	}
      }
    }
    else if(How >= 6){
      for(int b=0;b<pNBlock();b++){
	int IfAdded = 0;
	if(!strncmp(Block[b].Name,"PEP",3) ) continue;//IfAdded = 1;
	//if(!strncmp(Block[b].Name,"ADDED",3) ) IfAdded = 1;
	if(!strcmp(Block[b].Name,"STUFFING") ) continue;
	if(!strcmp(Block[b].Name,"WATER") ) continue;
	if(!strcmp(Block[b].Name,"CHOL") ) IfAdded = 1;
	if(!strcmp(Block[b].Name,"OIL") )  IfAdded = 1;
	if(!strcmp(Block[b].Name,"ADDED") )  IfAdded = 1;
	for(int p=Block[b].InitIdx,link=0;p<Block[b].EndIdx;p++){
	  double Dist = pPos(p,CLat1) - RefPos[CLat1];
	  Dist -= floor(Dist*pInvEdge(CLat1)+.5)*pEdge(CLat1) - .5*pEdge(CLat1);
	  int vx  = (int)(NBin*Dist*pInvEdge(CLat1));
	  Dist = pPos(p,CLat2) - RefPos[CLat2];
	  Dist -= floor(Dist*pInvEdge(CLat2)+.5)*pEdge(CLat2) - .5*pEdge(CLat2);
	  int vy = (int)(NBin*Dist*pInvEdge(CLat2));
	  int t = pType(p);
	  if( IfAdded ) Plot[2][vx*NBin+vy] += 1.;
	  else Plot[t][(vx*NBin+vy)] += 1.;
	}
      }
    }
  }
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,Plot, NType*NBin*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    printf("\n");
    //------------Normalize----------------------
    if(How < 6){
      double *VolContr = (double *)calloc(NBin,sizeof(double));
      VolumeCircSlab(VolContr,NBin);
      for(int v=0;v<NBin;v++){
	Dens[v] /= VolContr[v]*(NFile[1] - NFile[0]);
	for(int vv=0;vv<NBin;vv++){
	  for(int t=0,n=0;t<NType;t++){
	    Plot[t][(v*NBin+vv)] /= VolContr[v]*(NFile[1]-NFile[0]);
	  }
	}
      }
      free(VolContr);
    }
    if(How >= 6){
      for(int v=0;v<NBin;v++){
	Dens[v] /= (NFile[1] - NFile[0]);
	for(int vv=0;vv<NBin;vv++){
	  for(int t=0,n=0;t<NType;t++){
	    Plot[t][(v*NBin+vv)] *= InvNFile;
	  }
	}
      }
    }
    Matrice Mask(5,5);
    Mask.FillGaussian(.5,3.);
    Mask.Print();
    int NDim = 2;
    int IfMinImConv = 0;
    if(How == 6) IfMinImConv = 1;
    for(int t=0;t<NType;t++){
      Mask.ConvoluteMatrix(Plot[t],NBin,NDim,IfMinImConv);
      Mask.ConvoluteMatrix(Plot[t],NBin,NDim,IfMinImConv);
    }
    //--------------------writes-the-files---------------------
    char FileName[256];
    char FileName2[256];
    if(How == 0){//Cm Reference
      sprintf(FileName,"RadDistrCmR%.1fH%.1fh%.1f",Nano[nNano].Rad,Nano[nNano].Hamaker,Nano[nNano].Height);
    }
    else if(How==1 || How == 4 || How == 5){//Nano reference
      sprintf(FileName,"RadDistrNanoR%.1fH%.1fh%.1f",Nano[nNano].Rad,Nano[nNano].Hamaker,Nano[nNano].Height);
    }      
    else if(How == 2){//Cm Nano Reference
      sprintf(FileName,"RadDistrCmNR%.1fH%.1fh%.1f",Nano[nNano].Rad,Nano[nNano].Hamaker,Nano[nNano].Height);
    }
    else if(How == 6){//Cm Nano Reference
      sprintf(FileName,"CartDistrNanoR%.1fH%.1fh%.1f",Nano[nNano].Rad,Nano[nNano].Hamaker,Nano[nNano].Height);
    }
    sprintf(FileName2,"%s.rzd",FileName);
    FILE *FileToWrite = fopen(FileName2,"w");
    PrintDens(FileToWrite,Plot,LatDim,NBin);
    fclose(FileToWrite);
    sprintf(FileName,"DensProfR%.1fH%.1fh%.1f.dat",Nano[nNano].Rad,Nano[nNano].Hamaker,Nano[nNano].Height);
    FILE *FToWrite = fopen(FileName,"w");
    for(int v=0;v<NBin;v++){
      double r = v*InvNBin*pEdge(3);
      fprintf(FToWrite,"%lf %lf\n",r,Dens[v]);
    }
#ifdef OMPI_MPI_H
  }
#endif
  for(int t=0;t<NType;t++)
    free(Plot[t]);
  free(Plot);
}
/** 
    Calculate the (r,z) profile of the bond length distribution. 
*/
void ElPoly::BondDistr(int NSample){
  int NType = 3;
  double **Plot = (double **)calloc(NType,sizeof(double));
  double **Count = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<3;t++){
    Plot[t] = (double *)calloc(SQR(NSample),sizeof(double));
    Count[t] = (double *)calloc(SQR(NSample),sizeof(double));
  }
  double InvNSample = 1./(double)NSample;
  double Cm[4] = {0.,0.,0.,0.};
  double Dist[4] = {0.,0.,0.,0.};
  double Bond[4] = {0.,0.,0.,0.};
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    OpenRisk(cFile[f],BF_CHAIN);
    for(int b=0;b<pNBlock();b++){
      int IfAdded = 0;
      if(!strncmp(Block[b].Name,"PEP",3) ) continue;// IfAdded = 1;
      for(int p=Block[b].InitIdx,link=0;p<Block[b].EndIdx-1;p++){
	if(pChain(p) != pChain(p+1)) continue;
	Dist[3] = 0.;
	Bond[3] = 0.;
	for(int d=0;d<3;d++){
	  Cm[d]   = .5*(pPos(p,d) + pPos(p+1,d));
	  Bond[d] = Pm[p].Pos[d] - Pm[p+1].Pos[d];
	  Bond[3] += SQR(Bond[d]);
	  Dist[d] = Cm[d] - Nano[0].Pos[d];
	  Dist[d] -= floor( Dist[d]*pInvEdge(d) + .5)*pEdge(d);
	}
	double DistR = sqrt(SQR(Dist[CLat1])+SQR(Dist[CLat2]));
	int t = pType(p);
	if(t != pType(p+1)) continue;
	if(IfAdded) t = 2;
	double z = .5*pEdge(CNorm) + Dist[CNorm];
	double r = DistR;
	int vz = (int)(z*pInvEdge(CNorm)*NSample);
	if(vz < 0 || vz >= NSample) continue;
	int vr = (int)(r*pInvEdge(3)*NSample);
	if(vr < 0 || vr >= NSample) continue;
	Plot[t][vr*NSample+vz] += sqrt(Bond[3]);
	Count[t][vr*NSample+vz] += 1.;
      }
    }
  }
  // for(int t=0;t<NType;t++){
  //   for(int v=0;v<SQR(NSample);v++){
  //     Plot[t][v] /= Count[t][v] > 0. ? 1./Count[t][v] : 1.;
  //   }
  // }
  double *VolContr = (double *)calloc(NSample,sizeof(double));
  VolumeCircSlab(VolContr,NSample);
  for(int v=0;v<NSample;v++){
    for(int vv=0;vv<NSample;vv++){
      for(int t=0,n=0;t<NType;t++){
	Plot[t][(v*NSample+vv)] /= VolContr[v]*(NFile[1]-NFile[0]);
	//double Norm = Count[t][v*NSample+vv] > 1. ? 1./Count[t][v*NSample+vv] : 1.;
	//Plot[t][(v*NSample+vv)] *= Norm;
      }
    }
  }
  free(VolContr);
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  int NDim = 2;
  int IfMinImConv = 0;
  for(int t=0;t<NType;t++){
    Mask.ConvoluteMatrix(Plot[t],NSample,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot[t],NSample,NDim,IfMinImConv);
  }
  char FileName2[60];
  sprintf(FileName2,"BondDistr.rzd");
  FILE *FileToWrite = fopen(FileName2,"w");
  double LatDim[3] = {pEdge(3),pEdge(CNorm),15.};
  PrintDens(FileToWrite,Plot,LatDim,NSample);
  fclose(FileToWrite);
  for(int t=0;t<NType;t++){
    free(Plot[t]);
    free(Count[t]);
  }
  free(Plot);
  free(Count);
}
/** 
    Calculate the (r,z) profile of the splay angle distribution. 
*/
void ElPoly::SplayDistr(int NSample){
  int NType = 3;
  double **Plot = (double **)calloc(NType,sizeof(double));
  double **Count = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<3;t++){
    Plot[t] = (double *)calloc(SQR(NSample),sizeof(double));
    Count[t] = (double *)calloc(SQR(NSample),sizeof(double));
  }
  double InvNSample = 1./(double)NSample;
  double Cm[4] = {0.,0.,0.,0.};
  double CmP[6] = {0.,0.,0.,0.,0.,0.};
  double E2EP[6] = {0.,0.,0.,0.,0.,0.};
  double GyrP[2] = {0.,0.};
  double DistBA[4] = {0.,0.,0.,0.};
  double DistCB[4] = {0.,0.,0.,0.};
  int IfMartini = 0;
  int HalfLim = (int)(Block[0].Asym*.5);
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  double Direct[6] = {0.,0.,0.,0.,0.,0.};
  //int pLimit[4] = {4,8,8,13};//Martini
  //double TailLength[2] = {1./4.,1./5.};//Martini     
  int pLimit[4] = {0,4,4,9};//DFT
  double TailLength[2] = {1./4.,1./5.};//DFT
  Vettore Ax0(1.,0.,0.);
  Vettore Ax2(0.,0.,1.);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    OpenRisk(cFile[f],BF_CHAIN);
    for(int b=0;b<pNBlock();b++){
      double *xc = (double *)calloc(pNPCh(b),sizeof(double));
      double *yc = (double *)calloc(pNPCh(b),sizeof(double));
      double *zc = (double *)calloc(pNPCh(b),sizeof(double));
      if(!strncmp(Block[b].Name,"PEP",3) ) continue;
      if(!strncmp(Block[b].Name,"ADDED",3) ) continue;
      for(int c=0;c<Block[b].NChain;c++){
	for(int i=0;i<2;i++){
	  GyrP[i] = 0.;
	  int p1 = Block[b].InitIdx + c*Block[b].NPCh + pLimit[i*2];
	  int p2 = Block[b].InitIdx + c*Block[b].NPCh + pLimit[i*2+1];
	  int NBead = 0;
	  for(int d=0;d<3;d++){
	    CmP[i*3+d] = 0.;
	    E2EP[i*3+d] = Pm[p1].Pos[d] - Pm[p2].Pos[d];
	  }
	  for(int p=p1;p<p2;p++){
	    int ppc = p-p1;
	    xc[ppc] = Pm[p].Pos[0];
	    yc[ppc] = Pm[p].Pos[1];
	    zc[ppc] = Pm[p].Pos[2];
	    NBead += 1;
	    for(int d=0;d<3;d++){
	      CmP[i*3+d] += Pm[p].Pos[d];
	    }
	  }
	  for(int d=0;d<3;d++){
	    CmP[i*3+d] /= (double) NBead;
	  }
	  for(int p=p1;p<p2;p++){
	    for(int d=0;d<3;d++){
	      GyrP[i] += SQR(Pm[p].Pos[d]-CmP[i*3+d]);
	    }
	  }
	  RETTA rzx = Mat->InterRett(zc,xc,NBead);
	  RETTA rzy = Mat->InterRett(zc,yc,NBead);
	  double x1 = Pm[p1].Pos[2]*rzx.m + rzx.q;
	  double x2 = Pm[p2].Pos[2]*rzx.m + rzx.q;
	  double y1 = Pm[p1].Pos[2]*rzy.m + rzy.q;
	  double y2 = Pm[p2].Pos[2]*rzy.m + rzy.q;
	  Direct[2*i+2] = Pm[p1].Pos[2] - Pm[p2].Pos[2];
	  Direct[2*i+1] = y1 - y2;
	  Direct[2*i+0] = x1 - x2;
	}
	Vettore Dir1(Direct[0],Direct[1],Direct[2]);
	Vettore Dir2(Direct[3],Direct[4],Direct[5]);
	double Angle = Ax0.Angle(&Dir1,&Dir2);
	if(isnan(Angle)) continue;
	//int pA = Block[b].InitIdx + c*Block[b].NPCh + 3;//Martini
	int pA = Block[b].InitIdx + c*Block[b].NPCh + 3;//DFT
	for(int d=0;d<3;d++){
	  Cm[d] = Nano[0].Pos[d] - Pm[pA].Pos[d];
	  CmP[0+d] = CmP[0+d] - Nano[0].Pos[d];
	  CmP[0+d] -= floor( CmP[0+d]*pInvEdge(d) + .5)*pEdge(d);
	  CmP[3+d] = CmP[3+d] - Nano[0].Pos[d];
	  CmP[3+d] -= floor( CmP[3+d]*pInvEdge(d) + .5)*pEdge(d);
	}
	double E2E1 = sqrt(SQR(E2EP[0+0])+SQR(E2EP[0+1])+SQR(E2EP[0+2]));
	double E2E2 = sqrt(SQR(E2EP[3+0])+SQR(E2EP[3+1])+SQR(E2EP[3+2]));
	double DistR = sqrt(SQR(Cm[CLat1])+SQR(Cm[CLat2]));
	double z = .5*pEdge(CNorm) + Cm[CNorm];
	double r = DistR;
	int vz = (int)(z*pInvEdge(CNorm)*NSample);
	if(vz < 0 || vz >= NSample) continue;
	int vr = (int)(r*pInvEdge(3)*NSample);
	if(vr < 0 || vr >= NSample) continue;
	Plot[0][vr*NSample+vz] += Angle;
	Count[0][vr*NSample+vz] += 1.;
	r = sqrt(SQR(CmP[0+CLat1])+SQR(CmP[0+CLat2]));
	z = .5*pEdge(CNorm) + CmP[0+CNorm];
	vz = (int)(z*pInvEdge(CNorm)*NSample);
	if(vz < 0 || vz >= NSample) continue;
	vr = (int)(r*pInvEdge(3)*NSample);
	if(vr < 0 || vr >= NSample) continue;
	// Plot[1][vr*NSample+vz] += Dir1.Norm()*TailLength[0];
	Plot[1][vr*NSample+vz] += E2E1*TailLength[0];
	Count[1][vr*NSample+vz] += 1.;
	Plot[2][vr*NSample+vz] += sqrt(GyrP[0]*TailLength[0]);
	Count[2][vr*NSample+vz] += 1.;
	r = sqrt(SQR(CmP[3+CLat1])+SQR(CmP[3+CLat2]));
	z = .5*pEdge(CNorm) + CmP[3+CNorm];
	vz = (int)(z*pInvEdge(CNorm)*NSample);
	if(vz < 0 || vz >= NSample) continue;
	vr = (int)(r*pInvEdge(3)*NSample);
	if(vr < 0 || vr >= NSample) continue;
	//Plot[1][vr*NSample+vz] += Dir2.Norm()*TailLength[1];
	Plot[1][vr*NSample+vz] += E2E2*TailLength[1];
	Count[1][vr*NSample+vz] += 1.;
	Plot[2][vr*NSample+vz] += sqrt(GyrP[1]*TailLength[1]);
	Count[2][vr*NSample+vz] += 1.;
      }
      free(xc);
      free(yc);
      free(zc);
    }
  }
  double *VolContr = (double *)calloc(NSample,sizeof(double));
  VolumeCircSlab(VolContr,NSample);
  for(int v=0;v<NSample;v++){
    for(int vv=0;vv<NSample;vv++){
      for(int t=0,n=0;t<NType;t++){
	Plot[t][(v*NSample+vv)] /= VolContr[v]*(NFile[1]-NFile[0]);
	//double Norm = Count[t][v*NSample+vv] > 1. ? 1./Count[t][v*NSample+vv] : 1.;
	//Plot[t][(v*NSample+vv)] *= Norm;
      }
    }
  }
  free(VolContr);
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  int NDim = 2;
  int IfMinImConv = 0;
  for(int t=0;t<NType;t++){
    Mask.ConvoluteMatrix(Plot[t],NSample,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot[t],NSample,NDim,IfMinImConv);
  }
  char FileName2[60];
  sprintf(FileName2,"SplayDistr.rzd");
  FILE *FileToWrite = fopen(FileName2,"w");
  double LatDim[3] = {pEdge(3),pEdge(CNorm),15.};
  PrintDens(FileToWrite,Plot,LatDim,NSample);
  fclose(FileToWrite);
  for(int t=0;t<NType;t++){
    free(Plot[t]);
    free(Count[t]);
  }
  free(Plot);
  free(Count);
}
/**
a 
*/
void ElPoly::RadDens2Thick(int NBin){
  double *ThickProf = (double *)calloc(3*NBin,sizeof(double));
  double *Count = (double *)calloc(3*NBin,sizeof(double));
  double **Plot = (double **)calloc(3,sizeof(double));
  for(int t=0;t<3;t++){
    Plot[t] = (double *)calloc(NBin*NBin,sizeof(double));
  }
  LoadDensFile(Plot,NBin);
  double Cm[2] = {0.,0.};
  double CmCount[2] = {0.,0.};
  double InvNBin = 1./(double)NBin;
  for(int sr=0;sr<NBin;sr++){
    for(int sz=0;sz<NBin/2;sz++){
      double Pos = (sz-NBin/2)*pEdge(1)*InvNBin;
      double Weight = Plot[0][sr*NBin+sz]*Plot[1][sr*NBin+sz];
      ThickProf[sr*3  ] += Pos*Weight;
      Count[sr*3  ] += Weight;
    }
    for(int sz=NBin/2;sz<NBin;sz++){
      double Pos = (sz-NBin/2)*pEdge(1)*InvNBin;
      double Weight = Plot[0][sr*NBin+sz]*Plot[1][sr*NBin+sz];
      ThickProf[sr*3+1] += Pos*Weight;
      Count[sr*3+1] += Weight;
    }
    for(int sz=0;sz<NBin;sz++){
      ThickProf[sr*3+2] += Plot[2][sr*NBin+sz];
      Count[sr*3+2] += 1.;
    }
  }
  //normalize
  double Max = 0.;
  for(int sr=0;sr<NBin;sr++){
    ThickProf[sr*3  ] /= Count[sr*3  ] > 0. ? Count[sr*3  ] : 1.;
    ThickProf[sr*3+1] /= Count[sr*3+1] > 0. ? Count[sr*3+1] : 1.;
    ThickProf[sr*3+2] /= Count[sr*3+2] > 0. ? Count[sr*3+2] : 1.;
    if(Max < ThickProf[sr*3+2]) Max = ThickProf[sr*3+2];
  }
  for(int sr=0;sr<NBin;sr++){
    if(ThickProf[sr*3+2] > .5*Max){
      ThickProf[sr*3  ] = 0.;
      ThickProf[sr*3+1] = 0.;
    }
  }
  char FName[256]; 
  sprintf(FName,"ThickDensProf.dat",cFile[NFile[0]]);
  FILE *F2Write = fopen(FName,"w");
  for(int sr=0;sr<NBin;sr++){
    double x = sr*pEdge(0)/(double)NBin;
    fprintf(F2Write,"%lf %lf %lf %lf %lf\n",x,.5*(-ThickProf[sr*3  ]+ThickProf[sr*3+1]),-ThickProf[sr*3  ],ThickProf[sr*3+1],ThickProf[sr*3+2]);
  }
  for(int t=0;t<3;t++)
    free(Plot[t]);
  free(Plot);
  free(ThickProf);
  free(Count);
}
/**
a 
*/
void ElPoly::RadDens2Thick2d(int NBin){
  double *ThickProf = (double *)calloc(NBin*NBin,sizeof(double));
  double *Count = (double *)calloc(NBin*NBin,sizeof(double));
  double Round = 0.001;
  for(int p=0;p<pNPart();p++){
    if(pType(p) != 1) continue;
    int sr = (int)((pPos(p,0)+Round)*pInvEdge(0)*NBin);
    int sz = (int)((pPos(p,1)+Round)*pInvEdge(1)*NBin);
    ThickProf[sr*NBin+sz] += pPos(p,2);
    Count[sr*NBin+sz] += 1.;
  }
  double Norm = 0.;
  for(int sr=0;sr<NBin;sr++){
    for(int sz=0;sz<NBin;sz++){
      if(Norm < ThickProf[sr*NBin+sz]) Norm = ThickProf[sr*NBin+sz];
    }
  }
  Norm = 1./Norm;
  char FName[256]; 
  sprintf(FName,"ThickDensProf2d.dat",cFile[NFile[0]]);
  FILE *F2Write = fopen(FName,"w");
  fprintf(F2Write,"%lf %lf %lf\n",0.,0.,0.);
  fprintf(F2Write,"%lf %lf %lf\n",pEdge(0),pEdge(1),1.);
  for(int sr=0;sr<NBin;sr++){
    double r = sr*pEdge(0)/(double)NBin;
    for(int sz=0;sz<NBin;sz++){
      double z = sz*pEdge(1)/(double)NBin;
      //if(Count[sr*2] <= 1. || Count[sr*2+1] <= 1.)continue;
      if(ThickProf[sr*NBin+sz]*Norm < .3) continue;
      fprintf(F2Write,"%lf %lf %lf\n",r,z,ThickProf[sr*NBin+sz]*Norm);
    }
  }
  fclose(F2Write);
  free(ThickProf);
}
/**
a 
*/
void ElPoly::ThickFromDens(int NBin){
  int NType = 3;
  double Round = 0.001;
  double InvNBin = 1./(double)NBin;
  double *Plot  = (double *)calloc(NType*NBin*NBin,sizeof(double));
  double *Line = (double *)calloc(2*NBin,sizeof(double));
  double *Count = (double *)calloc(2*NBin,sizeof(double));
  for(int p=0;p<pNPart();p++){
    int vr = (int)((pPos(p,CLat1)+Round)*pInvEdge(0)*NBin);
    if (vr<0 || vr>NBin) continue;
    int vz = (int)((pPos(p,CLat2)+Round)*pInvEdge(1)*NBin);
    if (vz<0 || vz>NBin) continue;
    //int t = Pm[p].Typ;
    for(int t=0;t<3;t++){
      Plot[(vr*NBin+vz)*NType+t] += pVel(p,t);
    }
  }
  for(int vr=0;vr<NBin;vr++){
    for(int vz=0;vz<NBin;vz++){
      double z = vz*InvNBin*pEdge(1);
      double Weight = Plot[(vr*NBin+vz)*NType+1]*Plot[(vr*NBin+vz)*NType+2];
      int t = 0;
      if(z>pCm(1)) t = 1;
      Line[vr*2+t] += z*Weight;
      Count[vr*2+t] += Weight;
    }
  }
  FILE *FLine = fopen("ThickProf.dat","w");
  for(int vr=0;vr<NBin;vr++){
    double r = vr*InvNBin*pEdge(0);
    Line[vr*2] /= Count[vr*2] > 0. ? Count[vr*2] : 1.;
    Line[vr*2+1] /= Count[vr*2+1] > 0. ? Count[vr*2+1] : 1.;
    fprintf(FLine,"%lf %lf %lf\n",r,Line[vr*2],Line[vr*2+1]);
  }  
  fclose(FLine);
  free(Plot);
  free(Count);
  free(Line);
}
/**
a 
*/
void ElPoly::SlabAngleProfs(int NBin,int NAngle,int Coord){
  NAngle = 5;
  double *ProfAngle = (double *)calloc(NAngle*NBin,sizeof(double));
  double *Count = (double *)calloc(NAngle*NBin,sizeof(double));
  double Ref[3] = {.5*pEdge(1),.5*pEdge(1),pEdge(2)};
  double BoxRad = sqrt( SQR(.5*pEdge(0)) + SQR(.5*pEdge(1)) ) ;
  double InvBoxRad = 1./BoxRad;
  int Coord1 = (Coord)%2;
  int Coord2 = (Coord+1)%2;
  FILE *Ciccia = fopen("SlabAngleSectors.dat","w");
  fprintf(Ciccia,"#l(%lf %lf %lf) v[%d] d[part]\n",pEdge(0),pEdge(1),10.,NBin);
  for(int p=0,p1=0;p<pNPart();p++){
    if(pPos(p,CNorm) < 1.) continue;
    double RelDist[2] = {pPos(p,CLat1) - Ref[CLat1],pPos(p,CLat2) - Ref[CLat2]};
    double r = hypot(RelDist[0],RelDist[1]);
    int rb = (int)(r*InvBoxRad*NBin);
    if(rb < 0 || rb >= NBin) continue;
    if(Coord == 2) RelDist[CLat2] = - RelDist[CLat2];
    if(Coord == 3) {RelDist[CLat1] = - RelDist[CLat1];RelDist[CLat2] = - RelDist[CLat2];}
    double phi = fabs(atan2(RelDist[Coord1],RelDist[Coord2]));
    int phib = (int)(phi*NAngle/(M_PI*.5));
    if(phib < 0 || phib >= NAngle) continue;
    ProfAngle[rb*NAngle+phib] += pPos(p,CNorm);
    Count[rb*NAngle+phib] += 1.;
    fprintf(Ciccia,"{t[%d 0 %d] x(%lf %lf %lf)}\n",p1++,phib,pPos(p,0),pPos(p,1),pPos(p,2));
  }
  fclose(Ciccia);
  FILE *FProf = fopen("SlabRadProfs.dat","w");
  for(int rb=0;rb<NBin;rb++){
    double r = rb*BoxRad/(double)NBin;
    fprintf(FProf,"%lf ",r);
    for(int phib=0;phib<NAngle;phib++){
      double Norm = Count[rb*NAngle+phib] > 0. ? 1./Count[rb*NAngle+phib] : 1.;
      fprintf(FProf,"%lf ",ProfAngle[rb*NAngle+phib]*Norm);
    }
    fprintf(FProf,"\n");
  }
  free(ProfAngle);
  free(Count);
}
/**
a 
*/
int ElPoly::CoreF(int NSample,int How){
  double Cm=0.;
  double Volume=0.;
  int NVolume = NSample*NSample*NSample;
  int NType = 3;//Gen->NType;
  double *Plot = (double *) calloc(NVolume*NType,sizeof(double));
  double Border[3][2];
  Border[0][0]= 0.;
  Border[0][1]= pEdge(0);
  Border[1][0]= 0.;
  Border[1][1]= pEdge(1);
  Border[2][0]= 0.;
  Border[2][1]= pEdge(2);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(How == 0){
      if(OpenRisk(cFile[f],BF_CHAIN))return 1;
    }
    else if(How == 1){
      if(OpenRisk(cFile[f],BF_NANO))return 1;
      BfDefChain();
    }
    for(int b=0;b<pNBlock();b++){
      for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
	//if(strcmp(Block[b].Name,"LIPID") ) continue;
	int vx = (int)(NSample*(pPos(p,0)-Border[0][0])/(Border[0][1]-Border[0][0]));
	if( vx < 0 || vx >= NSample) continue;
	int vy = (int)(NSample*(pPos(p,1)-Border[1][0])/(Border[1][1]-Border[1][0]));
	if( vy < 0 || vy >= NSample) continue;
	int vz = (int)(NSample*(pPos(p,2)-Border[2][0])/(Border[2][1]-Border[2][0]));
	if( vz < 0 || vz >= NSample) continue;
	int t = pType(p);
	if( CHAIN_IF_TYPE(Ch[pChain(p)].Type,CHAIN_ADDED) ) t = 2;
	if( t < 0 || t > 3) continue;
	Plot[((vx*NSample+vy)*NSample+vz)*NType+t] += 1.;
      }
      Cm += pCm(2);
    }
  }
  printf("\n");
  Cm /= NFile[1] - NFile[0];
  FILE *FileToWrite = NULL;
  FileToWrite = fopen("Core.xvl","w");
  fprintf(FileToWrite,"#l(%lf %lf %lf) v[%d] d[%s]\n",pEdge(0),pEdge(1),pEdge(2),NSample,ChooseDraw(EL_COLOR));
  double *Norm = (double *)calloc(NType,sizeof(double));
  for(int t=0;t<NType;t++){
    for(int v=0;v<NVolume;v++){
      if(Norm[t] < Plot[v*NType+t])
	Norm[t] = Plot[v*NType+t];
    }
    Norm[t] = Norm[t] <= 0. ? 1. : Norm[t];
  }
  //for(int t=0;t<NType;t++){
  for(int vx=0;vx<NSample;vx++){
    for(int vy=0;vy<NSample;vy++){
      for(int vz=NSample-1;vz>0;vz--){
	double Dens0 = Plot[((vx*NSample+vy)*NSample+vz)*NType+2]/Norm[2];
	double Dens1 = Plot[((vx*NSample+vy)*NSample+vz)*NType+0]/Norm[0];
	double Dens2 = Plot[((vx*NSample+vy)*NSample+vz)*NType+1]/Norm[1];
	if(Dens0+Dens1+Dens2 <= 0.0)continue;
	fprintf(FileToWrite,"{x(%lf %lf %lf) ",
		Border[0][0]+vx/(double)NSample*(Border[0][1]-Border[0][0]),
		Border[1][0]+vy/(double)NSample*(Border[1][1] - Border[1][0]),
		Border[2][0]+vz/(double)NSample*(Border[2][1] - Border[2][0]));
	fprintf(FileToWrite,"v(%lf %lf %lf)}\n",Dens0,Dens1,Dens2);
	Volume += 1.;
      }
    }
  }
  //}
  Volume /= NSample*NSample*NSample;
  Volume *= (Border[0][1] - Border[0][0])*(Border[1][1] - Border[1][0])*(Border[2][1] - Border[2][0]);
  fprintf(FileToWrite,"#Volume %lf\n",Volume);
  fclose(FileToWrite);
  free(Plot);
  free(Norm);
  return 0;
}
