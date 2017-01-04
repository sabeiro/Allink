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

int VarData::PairCorrelation(double *dPoint,int NSample,int How,int Type){
  BfDefChain();
  //printf("%d %d\n",NSample,How);
  double dNSample = 1./(double)NSample;
  for(int v=0;v<NSample;v++)
    dPoint[v] = 0.;
  if(How == 0){//Monomer
    for(int p=0;p<Gen->NPart;p++){
      if(Pm[p].Typ != Type)continue;
      Pm[p].Pos[3] = 0.;
      for(int pp=0;pp<Gen->NPart;pp++){
	if(Pm[pp].Typ != Type)continue;
	if( p == pp)continue;
	double Dist2 = QUAD(Pm[p].Pos[0] - Pm[pp].Pos[0]);
	Dist2 += QUAD(Pm[p].Pos[1] - Pm[pp].Pos[1]);
	Dist2 += QUAD(Pm[p].Pos[2] - Pm[pp].Pos[2]);
	Dist2 = pow(Dist2,.5);
	int v = (int)(Dist2*pInvEdge(3)*NSample);
	if( v < 0 || v >= NSample) continue;
	dPoint[v] += 1.;
      }
    }
    //      printf("Processing: %d/%d %.1f %%             \r",p,HowMany,p/(double)HowMany/100.);
    for(int v=0;v<NSample;v++){
      dPoint[v] /= Gen->NPart*(Gen->NPart-1);
    }
  }
  //printf("\n");
  if(How == 1){//Chain
    for(int c=0;c<Gen->NChain;c++){
      Ch[c].Pos[3] = 0.;
      if(Ch[c].Type != Type)continue;
      for(int cc=0;cc<Gen->NChain;cc++){
      if(Ch[cc].Type != Type)continue;
      if( c == cc)continue;
      Ch[c].Pos[3] += QUAD((Ch[c].Pos[CLat1] - Ch[cc].Pos[CLat1]));
      Ch[c].Pos[3] += QUAD((Ch[c].Pos[CLat2] - Ch[cc].Pos[CLat2]));
      }
      Ch[c].Pos[3] = pow(Ch[c].Pos[3]/(double)Gen->NChain,.5);
      int v = (int)(Ch[c].Pos[3] / Gen->Edge[3]*NSample);
      if( v < 0 || v >= NSample) continue;
      dPoint[v] += 1.;
      }
  }
  for(int v=0;v<NSample;v++){
    dPoint[v] /= DUE_PI*( QUAD((Gen->Edge[3]*(v+1)*dNSample)) - QUAD((Gen->Edge[3]*v*dNSample)) );
  }
  return 0;
}
int VarData::PairCorrelationSquare(double **dPoint,int NSample,int Type){
  BfDefChain();
  double dNSample = 1./(double)NSample;
  double InvNChain=1./(double)Gen->NChain;
  int vRef = (int)(NSample/2.);
  for(int c=0;c<Gen->NChain;c++){
    if(!CHAIN_IF_TYPE(Ch[c].Type,NChType)) continue;
    double ChX1 = (Ch[c].Pos[CLat1] - Gen->Edge[CLat1]*.5);
    double ChY1 = (Ch[c].Pos[CLat2] - Gen->Edge[CLat2]*.5);
    int vx1 = (int)(ChX1 / (Gen->Edge[CLat1])*NSample);
    int vy1 = (int)( (ChY1) / (Gen->Edge[CLat2])*NSample);
    for(int cc=0;cc<Gen->NChain;cc++){
      if(!CHAIN_IF_TYPE(Ch[cc].Type,NChType)) continue;
      if( c == cc) continue;
      double ChX2 = (Ch[cc].Pos[CLat1] - Gen->Edge[CLat1]*.5);
      double ChY2 = (Ch[cc].Pos[CLat2] - Gen->Edge[CLat2]*.5);
      int vx2 = (int)(ChX2 / (Gen->Edge[CLat1])*NSample);
      int vy2 = (int)( (ChY2) / (Gen->Edge[CLat2])*NSample);
      int vx = vx1 - vx2;
      // if(vx < -vRef) vx += vRef;
      // if(vx > vRef ) vx -= vRef;
      if( vx + vRef < 0 || vx+vRef >= NSample) continue;
      int vy = vy1 - vy2;
      // if(vy < -vRef) vy += vRef;
      // if(vy > vRef)  vy -= vRef;
      if( vy + vRef < 0 || vy+vRef >= NSample) continue;
      //      printf("%d %d %d %lf %lf\n",c,vx,vy,ChX1,ChX2);
      dPoint[vx+vRef][vy+vRef] += InvNChain;
    }
  }
  return 0;
}
int VarData::PairCorrelationPep(double **dPoint,int NSample,int Type){
  BfDefChain();
  double dNSample = 1./(double)NSample;
  double InvNChain=1./(double)Gen->NChain;
  int vRef = (int)(NSample/2.);
  for(int c=0;c<Gen->NChain;c++){
    if(!CHAIN_IF_TYPE(Ch[c].Type,NChType))continue;
    double ChX = remainder(Ch[c].Pos[CLat1] - pNanoPos(0,CLat1),Gen->Edge[CLat1]);
    double ChY = remainder(Ch[c].Pos[CLat2] - pNanoPos(0,CLat2),Gen->Edge[CLat2]);
    int vx = (int)(ChX/Gen->Edge[CLat1]*NSample);
    int vy = (int)(ChY/Gen->Edge[CLat2]*NSample);
    if( vx + vRef < 0 || vx+vRef >= NSample) continue;
    if( vy + vRef < 0 || vy+vRef >= NSample) continue;
    dPoint[vx+vRef][vy+vRef] += InvNChain;
  }
  return 0;
}
int VarData::PairCorrelationRound(double **dPoint,int NSample,int Type){
  BfDefChain();
  double dNSample = 1./(double)NSample;
  double InvNChain=1./(double)Gen->NChain;
  for(int c=0;c<Gen->NChain;c++){
    double ChRad=0.;
    double ChAngle=0.;
    if(!CHAIN_IF_TYPE(Ch[c].Type,NChType))continue;
    for(int cc=0;cc<Gen->NChain;cc++){
      if(!CHAIN_IF_TYPE(Ch[cc].Type,NChType))continue;
      if( c == cc)continue;
      double ChX = remainder(Ch[c].Pos[CLat1] - Ch[cc].Pos[CLat1],Gen->Edge[CLat1]);
      // if(ChX > .5*Gen->Edge[CLat1]) 
      // 	ChX = Ch[c].Pos[CLat1] + Ch[cc].Pos[CLat1] - Gen->Edge[CLat1];
      // else if (ChX < -.5*Gen->Edge[CLat1]) 
      // 	ChX = -Ch[c].Pos[CLat1] - Ch[cc].Pos[CLat1] + Gen->Edge[CLat1];
      double ChY = remainder(Ch[c].Pos[CLat2] - Ch[cc].Pos[CLat2],Gen->Edge[CLat2]);
      // if(ChY > .5*Gen->Edge[CLat2]) 
      // 	ChY = Ch[c].Pos[CLat2] + Ch[cc].Pos[CLat2] - Gen->Edge[CLat2];
      // else if (ChY < -.5*Gen->Edge[CLat2]) 
      // 	ChY = -Ch[c].Pos[CLat2] - Ch[cc].Pos[CLat2] + Gen->Edge[CLat2];
      ChRad = sqrt( QUAD((ChX)) + QUAD((ChY)) );
      ChAngle = acos(ChX / ChRad);
      if(ChY < 0)
	ChAngle = DUE_PI - ChAngle;
      int v = (int)(ChRad / (Gen->Edge[3])*NSample);
      if( v < 0 || v >= NSample) continue;
      int vv = (int)( (ChAngle) / (DUE_PI)*NSample);
      if( vv < 0 || vv >= NSample) continue;
      //printf("%d %d %lf %lf %lf %lf\n",v,vv,ChRad,ChX,ChY,ChAngle);
      dPoint[v][vv] += InvNChain;
    }
  }
  return 0;
}
int VarData::Scattering2d(double **Plot,int NSample,int Type){
  BfDefChain();
  double dNSample = 1./(double)NSample;
  double InvNChain=1./(double)Gen->NChain;
  int vRef = (int)(NSample/2.);
  double *CosSin = (double *)calloc(2*SQR(NSample),sizeof(double));
  for(int p=0;p<pNPart();p++){
    //if(!CHAIN_IF_TYPE(Ch[Pm[p].CId].Type,NChType))continue;
    // double ChX = (Ch[c].Pos[CLat1] - Nano->PosBf[CLat1]);
    // double ChY = (Ch[c].Pos[CLat2] - Nano->PosBf[CLat2]);
    double ChX = (Pm[p].Pos[CLat1] - pNanoPos(0,CLat1));
    double ChY = (Pm[p].Pos[CLat2] - pNanoPos(0,CLat2));
    ChX -= floor(ChX*pInvEdge(CLat1))*pEdge(CLat1);
    ChY -= floor(ChY*pInvEdge(CLat2))*pEdge(CLat2);
    ChX *= pInvEdge(CLat1);
    ChY *= pInvEdge(CLat2);
    //	printf("%lf %lf %lf %lf\n",vvx,ChX,vvx*ChX*DUE_PI,cos(vvx*ChX*DUE_PI));
    for(int vx=0;vx<NSample;vx++){
      for(int vy=0;vy<NSample;vy++){
	double qx = vx*dNSample;
	double qy = vy*dNSample;
	CosSin[(vx*NSample+vy)*2  ] += cos( (qx*ChX+qy*ChY)*DUE_PI );
	CosSin[(vx*NSample+vy)*2+1] += sin( (qx*ChX+qy*ChY)*DUE_PI );
      }
    }
  }
  for(int vx=0;vx<NSample;vx++){
    for(int vy=0;vy<NSample;vy++){
      Plot[vx][vy] = (SQR(CosSin[(vx*NSample+vy)*2  ]) + SQR(CosSin[(vx*NSample+vy)*2+1]))*SQR(InvNChain);
    }
  }
  return 0;
}
int VarData::Scattering2D(double **dPoint,int NSample,int Type){
  BfDefChain();
  double dNSample = 1./(double)NSample;
  double InvNChain=1./(double)Gen->NChain;
  for(int c=0;c<Gen->NChain;c++){
      if(!CHAIN_IF_TYPE(Ch[c].Type,NChType))continue;
    for(int cc=0;cc<Gen->NChain;cc++){
      if(!CHAIN_IF_TYPE(Ch[cc].Type,NChType))continue;
      if( c == cc)continue;
      double ChX = (Ch[c].Pos[CLat1] - Ch[cc].Pos[CLat1]);
      if(ChX > .5*Gen->Edge[CLat1]) ChX = Ch[c].Pos[CLat1] + Ch[cc].Pos[CLat1] - Gen->Edge[CLat1];
      else if (ChX < -.5*Gen->Edge[CLat1]) ChX = -Ch[c].Pos[CLat1] - Ch[cc].Pos[CLat1] + Gen->Edge[CLat1];
      double ChY = (Ch[c].Pos[CLat2] - Ch[cc].Pos[CLat2]);
      if(ChY > .5*Gen->Edge[CLat2]) ChY = Ch[c].Pos[CLat2] + Ch[cc].Pos[CLat2] - Gen->Edge[CLat2];
      else if (ChY < -.5*Gen->Edge[CLat2]) ChY = -Ch[c].Pos[CLat2] - Ch[cc].Pos[CLat2] + Gen->Edge[CLat2];
      int v = NSample/2 + (int)( 2.*(ChX / Gen->Edge[CLat1])*(NSample/2) );
      int vv = NSample/2 + (int)(2.*ChY / (Gen->Edge[CLat2])*(NSample/2));
      //if(ChX > Gen->Edge[CLat1]*.5 || ChY > Gen->Edge[CLat2]*.5 ) printf("%d %d %lf %lf\n",v,vv,ChX,ChY);
      //      printf("%d %d %lf %lf\n",v,vv,ChX,ChY);
      if( vv < 0 || vv >= NSample) continue;
      if( v < 0 || v >= NSample) continue;
      dPoint[vv][v] += InvNChain;//*sin(ChRec*ChRad)/ChRec;
     }
  }
  return 0;
}
void VarData::Spettro2d(double *Points,int NSample,int Type){
  double *Plot = (double *)calloc(SQR(NSample),sizeof(double));
  double *InPoints = (double *)calloc(NSample*NSample,sizeof(double));
  SampleSurface(Plot,NSample,Type);
  for(int v=0;v<SQR(NSample);v++){
    InPoints[v] = Plot[v];
  }
  Mat->Spettro2d(InPoints,Points,NSample);
  free(Plot);
  free(InPoints);
}
void VarData::Spettro2d(double *Plot,int NSample){
  double *Points = (double *)calloc(NSample*NSample,sizeof(double));
  for(int v=0;v<SQR(NSample);v++){
    Points[v] = Plot[v];
  }
  Mat->Spettro2d(Points,Plot,NSample);
  free(Points);
}
