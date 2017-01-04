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

int VarData::Arrange(int **Triangle,int Vertex){
  BfDefChain();
  int *nChain = (int *)calloc(Vertex,sizeof(int));
  double *dChain= (double *)calloc(Vertex,sizeof(double));
  double MaxLength = 1.5*sqrt(QUAD(Gen->Edge[0])+QUAD(Gen->Edge[1]));
  for(int c=0;c<Gen->NChain;c++){
    for(int v=0;v<Vertex;v++)
      dChain[v] = MaxLength;
    int NVertex = 0;
    if(!CHAIN_IF_TYPE(Ch[c].Type,NChType))continue;
    for(int cc=0;cc<Gen->NChain;cc++){
      if(!CHAIN_IF_TYPE(Ch[cc].Type,NChType))continue;
      if(c == cc) continue;
      double Shift[3] = {0.,0.,0.};
      double PosDiff[3] = {0.,0.,0.};
      for(int d=0;d<3;d++){
	PosDiff[d] = Ch[cc].Pos[d] - Ch[c].Pos[d];
	Shift[d] = -floor((PosDiff[d])/Gen->Edge[d]+.5)*Gen->Edge[d];
      }
      double Dist = QUAD(PosDiff[CLat1]+Shift[CLat1]) + QUAD(PosDiff[CLat2]+Shift[CLat2]);// + QUAD(PosZ);
      //if(Dist > MaxLength) continue;
      //printf("%d %d %lf %lf\n",c,cc,Dist,dChain[0]);
      for(int v=0;v<Vertex;v++){
	if(dChain[v] > Dist){
	  //printf("%d %d/%d %lf %lf\n",c,cc,v,Dist,dChain[v]);
	  for(int vv=Vertex-1;vv>v;vv--){
	    dChain[vv] = dChain[vv-1];
	    nChain[vv] = nChain[vv-1];
	  }
	  dChain[v] = Dist;
	  nChain[v] = cc;
	  NVertex++;
	  //for(int v=0;v<Vertex;v++)printf("%d %d\n",v,nChain[v]);
	  break;
	}
      }
    }
    for(int v=0;v<Vertex;v++)
      Triangle[v][c] = nChain[v];
    Triangle[Vertex][c] = NVertex<Vertex ? NVertex : Vertex;
    //Triangle[Vertex-1][c] = c;
    //for(int v=0;v<Vertex+1;v++) printf("%d ",Triangle[v][c]); printf("\n");
  }
  free(nChain);
  free(dChain);
  return 1;
}
int VarData::OrderPos(){
  // for(int c=0;c<Gen->NChain;c++){
  //   Ch[c].nPos = CalcnPos(Ch[c].Pos);
  //   //    printf("%d %d\n",c,Ch[c].nPos);
  //   for(int cc=c-1;cc>=0;cc--){
  //     int nProl = cc==0 ? 0 : Ch[cc-1].nPos;
  //     //printf("%d<%d<%d | %d<=%d<%d\n",cc,c,cc-1,Ch[cc].nPos,Ch[c].nPos,nProl);
  //     if(Ch[c].nPos < Ch[cc].nPos && Ch[c].nPos >= nProl){
  // 	  //printf("--%d/%d - %d/%d \n",c,cc+1,Ch[c].nPos,nProl);
  // 	  for(int ccc=c;ccc>cc;ccc--) SwapChain(ccc,ccc-1);
  // 	  break;
  //     }
  //   }
  //     //for(int cc=0;cc<=c;cc++)printf("%d %d\n",cc,Ch[cc].nPos);
  // }
  // return 1;
}
int VarData::CalcnPos(double *Pos){
  double BlurDist[3]={Gen->Edge[0]/(double)SubDiv[0],Gen->Edge[1]/(double)SubDiv[1],Gen->Edge[2]/(double)SubDiv[2]};
  int v[3];
  for(int d=0;d<3;d++){
    double InvSubDiv = 1./(double)SubDiv[d];
    v[d] = (int)(Pos[d]/Gen->Edge[d]*SubDiv[d]);
    if(v[d] > SubDiv[d]){printf("Chain out of range %d>%d\n",v[d],SubDiv[d]); continue;}
    if(Pos[d] - v[d]*Gen->Edge[d]*InvSubDiv < .25*BlurDist[d])
      v[d] = 2*v[d] + 1;
    else if(Pos[d] - v[d]*Gen->Edge[d]*InvSubDiv > .75*BlurDist[d])
      v[d] = 2*v[d] + 3;
    else 
      v[d] = 2*v[d];
    //printf("%d\n",v[d]);
  }
  int nPos =  (v[0]*SubDiv[1] + v[1])*SubDiv[2]+v[2];
  return nPos;
}
int VarData::PosVectInt(double *Pos){
  int vx = (int)(Pos[CLat1]/Gen->Edge[CLat1]*NEdge);
  int vy = (int)(Pos[CLat2]/Gen->Edge[CLat2]*NEdge);
  int vz = (int)(Pos[CNorm]/Gen->Edge[CNorm]*NEdge);
  return (vx*NEdge + vy)*NEdge+vz;
}
