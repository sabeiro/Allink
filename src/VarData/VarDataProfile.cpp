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
;				 
/// sum on small patches and shift the wrt the weighted average
int VarData::DensityProfile(int coord,int NBin,int NType,double *DensityAv){
  BfDefChain();
  int NPartition = 5;
  double *Density = (double *)calloc(NPartition*NPartition*NType*NBin,sizeof(double));
  double Pos[4];
  char cLine[STRSIZE];
  for(int b=0;b<Gen->NBlock;b++){
    if(!strncmp(Block[b].Name,"PEP",3) ) continue;
    if(!strcmp(Block[b].Name,"WATER") ) continue;
    for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
      int t = Pm[p].Typ;
      for(int d=0;d<3;d++)
	Pos[d] = Pm[p].Pos[d];
      int vx = (int)(NPartition*Pos[CLat1]/Gen->Edge[CLat2]);
      if(vx < 0 || vx >= NBin) continue;
      int vy = (int)(NPartition*Pos[CLat2]/Gen->Edge[CLat2]);
      if(vy < 0 || vy >= NBin) continue;
      int vz = (int)(NBin*Pos[CNorm]/Gen->Edge[CNorm]);
      if(vz < 0 || vz >= NBin) continue;
      if(CHAIN_IF_TYPE(Ch[Pm[p].CId].Type,CHAIN_ADDED))
	t = 2;
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
  return 0;
}
// Volume element in cylindrical symmetry
void VarData::VolumeCircSlab(double *VolContr,int NBin){
  if(pEdge(3) <= 0.){printf("Lateral radius not defined\n");};
  double LenMin = .5*MIN(Gen->Edge[CLat1],Gen->Edge[CLat2]);
  double LenMax = .5*MAX(Gen->Edge[CLat1],Gen->Edge[CLat2]);
  double RapCorner = LenMin/LenMax*NBin;
  double RapMin = LenMin/Gen->Edge[3]*NBin;
  double Volume2 = 0.;
  for(int v=0;v<NBin;v++){
    double Rad2 = (v+1);
    double Volume1 = PI*QUAD(Rad2);
    if(Rad2 > RapMin){
      double Pro2 = sqrt(QUAD(v+1)-QUAD(RapMin));
      double Angle2 = atan(RapMin/Pro2);
      Volume1 = 2.*QUAD(Rad2)*Angle2 + 2.*Pro2*RapMin;
      if(Rad2 > RapCorner){
	double theta = acos(RapCorner/Rad2);
	double b = Rad2*sin(theta);
	//double b = sqrt( SQR(Rad2) - SQR(RapCorner));
	Volume1 -= 4.*(Rad2*Rad2*theta - b*RapCorner);
	//printf("%lf %lf %lf %lf %lf %lf %lf\n",Rad2,RapCorner,theta,b,Volume1-Volume2,Rad2*Rad2*theta ,b*RapCorner);
      }
    }
    double Volume = (Volume1 - Volume2)*Gen->Edge[CNorm]*SQR(Gen->Edge[3])/(double)CUB(NBin);
    VolContr[v] = Volume;
    Volume2 = Volume1;
  }
}
int VarData::Worm(int Partition,int NBin,double *Border,double *dPoint){
  int cLong = CLat2;
  int cLat = CLat1;
  double *Cm;
  RETTA *r1;
  r1 = (RETTA *)calloc(Partition,sizeof(RETTA));
  Cm = (double *)calloc(3*Partition,sizeof(double));
  double* Dis0 = (double *)calloc(Gen->NPart,sizeof(double));
  double* Dis1  = (double *)calloc(Gen->NPart,sizeof(double));
  for(int pa = 0,cc=0;pa<Partition;pa++){
    for(int c=0;c<Gen->NChain;c++){
      int v = (int)(Ch[c].Pos[cLong]*Partition/Gen->Edge[cLong]);
      if(v != pa) continue;
      *(Dis0 + cc) = Ch[c].Pos[cLat];
      *(Dis1 + cc) = Ch[c].Pos[cLong];
      //printf("%lf %lf %lf\n",Ch[c].Pos[0],Ch[c].Pos[1],Ch[c].Pos[2]);
      for(int d=0;d<3;d++)
	*(Cm + 3*v + d) += Ch[c].Pos[d];
      cc++;
    }
    r1[pa] = Mat->InterRett(Dis0,Dis1,cc);
    for(int d=0;d<3;d++)
      *(Cm + 3*pa + d) /= cc;
    cc = 0;
  }
  double pos0;
  double pos1;
  double qcm;
  double qpoint = 0.;
  double dist;
  int pp0 = 0;
  int pp1 = 0;
  //FILE *Ciccia = fopen("Ciccia.dat","w");
  for(int pa=0;pa<Partition;pa++){
    //printf("m %lf Cm %lf %lf %lf \n",r1[pa].m,*(Cm+3*pa+0),*(Cm+3*pa+1),*(Cm+3*pa+2));
    for(int p=0;p<Gen->NPart;p++){
      int v = (int) (Pm[p].Pos[cLong]*Partition/Gen->Edge[cLong]);
      if(v != pa) continue;
      qcm = *(Cm + 3*pa + cLat) + *(Cm + 3*pa + cLong)/r1[pa].m;
      qpoint = Pm[p].Pos[cLat] - r1[pa].m*Pm[p].Pos[cLong];
      pos0 = (qcm-qpoint)*r1[pa].m/(QUAD(r1[pa].m) + 1.);
      pos1 = pos0*r1[pa].m + qpoint;
      dist = QUAD((pos0-Cm[3*pa+cLong])) + QUAD((pos1-Cm[3*pa+cLat]));
      if(Pm[p].Typ == 0){
	Dis0[pp0] = pow(dist,.5);
	if(Pm[p].Pos[cLat] < *(Cm+3*pa+cLat))
	  Dis0[pp0] *= -1.;
	//printf("Dis0 %lf\n",Dis0[pp0]);
	pp0++;
      }
      else if(Pm[p].Typ == 1){
	Dis1[pp1] = pow(dist,.5);
	if(Pm[p].Pos[cLat] < *(Cm+3*pa+cLat))
	  Dis1[pp1] *= -1.;
	pp1++;
      }
      if(pa == 2){
	//fprintf(Ciccia,"Ch %d %lf %lf Pos %lf %lf pos %lf %lf Dis1 %lf Dis2 %lf dist %lf\n",Pm[p].CId,Ch[Pm[p].CId].Pos[0],Ch[Pm[p].CId].Pos[2],Pm[p].Pos[0],Pm[p].Pos[2],pos0,pos1,Dis0[pp0],Dis1[pp1],pow(dist,.5));
      }
    }
    //      fclose(Ciccia);
    for(int v=0;v<NBin;v++)
      dPoint[v*2] = 0.;
    for(int i=0;i<pp0;i++){
      int v = (int) ((Dis0[i]-Border[0])/(Border[1] - Border[0])*NBin);
      if(v>= NBin) continue;
      dPoint[v*2] += 1.;
    }
    for(int v=0;v<NBin;v++){
      dPoint[v*2] *= (NBin*Partition)/(Gen->Edge[0]*Gen->Edge[1]*Gen->Edge[2]);
    }
    for(int v=0;v<NBin;v++)
      dPoint[v*2] = 0.;
    for(int i=0;i<pp1;i++){
      int v = (int) ((Dis1[i]-Border[0])/(Border[1] - Border[0])*NBin);
      if(v>= NBin) continue;
      dPoint[v*2] += 1.;
    }
    for(int v=0;v<NBin;v++){
      dPoint[v*2+1] /= (double) (NBin*Partition)/(Gen->Edge[0]*Gen->Edge[1]*Gen->Edge[2]);
      //if(pa==0)printf("%d %lf %lf\n",v,dPoint[v],dPoint1[v]);
    }
    pp0=0;
    pp1=0;
  }
  for(int v=0;v<NBin;v++){
    dPoint[v*2] /= (double) NBin;
    dPoint[v*2+1] /= (double)NBin;
  }
  return 1;
}
void VarData::Stalk(int NBin,int NLevel,double **Plot,double Threshold){
  double **Norma = (double **)calloc(NLevel,sizeof(double));
  for(int l=0;l<NLevel;l++)
    Norma[l] = (double *)calloc(NBin*NBin,sizeof(double));
  for(int b=0;b<Gen->NBlock;b++){
    for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
      for(int v=0;v<NBin;v++){
	// if(Pm[p].Typ == 1 && Pm[p-1].Typ == 0)
	  {
	  int v = (int)(Pm[p].Pos[CLat1]/Gen->Edge[CLat1]*NBin);
	  if( v < 0 || v >= NBin) continue;
	  int vv = (int)(Pm[p].Pos[CLat2]/Gen->Edge[CLat2]*NBin);
	  if( vv < 0 || vv >= NBin) continue;
	  Plot[b][v*NBin + vv] += Pm[p].Pos[CNorm];
	  Norma[b][v*NBin + vv] += 1.;
	}
      }
    }
  }
  for(int v=0;v<NBin;v++)
    for(int vv=0;vv<NBin;vv++)
      for(int l=0;l<NLevel;l++){
	if(Norma[l][v*NBin+vv] > Threshold)
	  Plot[l][v*NBin+vv] /= Norma[l][v*NBin+vv];
	else
	  Plot[l][v*NBin+vv] = 0.;
      }
  for(int l=0;l<NLevel;l++)
    free(Norma[l]);
  free(Norma);
}
