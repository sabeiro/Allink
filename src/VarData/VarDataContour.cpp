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

int VarData::DensityProfile(int coord,int Values,int NType,double *DensityAv){
  BfDefChain();
  double VolumeInv = (double) (Values)/(Gen->Edge[0]*Gen->Edge[1]*Gen->Edge[2]);
  int NPartition = 10;
  double NPartInv = 1./(double)(NPartition*NPartition);
  double *Density = (double *)calloc(NPartition*NPartition*NType*Values,sizeof(double));
  double Pos[4];
  char cLine[STRSIZE];
  int t=0;
  for(int b=0;b<Gen->NBlock;b++){
    if(!strncmp(Block[b].Name,"PEP",3) ) continue;
    if(!strcmp(Block[b].Name,"WATER") ) continue;
    for(int p=Block[b].InitIdx,link=0;p<Block[b].EndIdx;p++){
      t = Pm[p].Typ;
      for(int d=0;d<3;d++)
	Pos[d] = Pm[p].Pos[d];
      int vx = (int)(NPartition*Pos[CLat1]/Gen->Edge[CLat2]);
      if(vx < 0 || vx >= Values) continue;
      int vy = (int)(NPartition*Pos[CLat2]/Gen->Edge[CLat2]);
      if(vy < 0 || vy >= Values) continue;
      int vz = (int)(Values*Pos[CNorm]/Gen->Edge[CNorm]);
      if(vz < 0 || vz >= Values) continue;
      if(CHAIN_IF_TYPE(Ch[Pm[p].CId].Type,CHAIN_ADDED))
	t = 2;
      Density[((vz*NType+t)*NPartition+vy)*NPartition+vx] += 1.;
    }
  }
  for(int vx=0;vx<NPartition;vx++){
    for(int vy=0;vy<NPartition;vy++){
      double Media = 0.;
      double Weight = 0.;
      for(int v=0;v<Values;v++){
	Weight += Density[((v*NType+0)*NPartition+vy)*NPartition+vx];
	Media += v*Density[((v*NType+0)*NPartition+vy)*NPartition+vx];
      }
      int vMedia = Weight > 0 ? (int)(Media/Weight) : Values/2;
      int vOffset = vMedia - Values/2;
      for(int t=0;t<NType;t++){
	for(int v=0;v<Values;v++){
	  int vOther = v + vOffset;
	  if(vOther < 0)
	    vOther = Values - vOther;
	  if(vOther >= Values)
	    vOther = vOther - Values;
	  DensityAv[v*NType+t] += Density[((vOther*NType+t)*NPartition+vy)*NPartition+vx]*VolumeInv;
	}
      }
    }
  }
  free(Density);
  return 0;
}
void VarData::Stalk(int Values,int NLevel,double **Plot,double Threshold){
  double **Norma = (double **)calloc(NLevel,sizeof(double));
  for(int l=0;l<NLevel;l++)
    Norma[l] = (double *)calloc(Values*Values,sizeof(double));
  for(int b=0;b<Gen->NBlock;b++){
    for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
      for(int v=0;v<Values;v++){
	// if(Pm[p].Typ == 1 && Pm[p-1].Typ == 0)
	  {
	  int v = (int)(Pm[p].Pos[CLat1]/Gen->Edge[CLat1]*Values);
	  if( v < 0 || v >= Values) continue;
	  int vv = (int)(Pm[p].Pos[CLat2]/Gen->Edge[CLat2]*Values);
	  if( vv < 0 || vv >= Values) continue;
	  Plot[b][v*Values + vv] += Pm[p].Pos[CNorm];
	  Norma[b][v*Values + vv] += 1.;
	}
      }
    }
  }
  for(int v=0;v<Values;v++)
    for(int vv=0;vv<Values;vv++)
      for(int l=0;l<NLevel;l++){
	if(Norma[l][v*Values+vv] > Threshold)
	  Plot[l][v*Values+vv] /= Norma[l][v*Values+vv];
	else
	  Plot[l][v*Values+vv] = 0.;
      }
  for(int l=0;l<NLevel;l++)
    free(Norma[l]);
  free(Norma);
}
