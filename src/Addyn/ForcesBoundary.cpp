#include "Forces.h"

void Forces::PullBead(){
  Pm[Bead2Move].Pos[2] += IncrDist;
  if(VAR_IF_TYPE(SysShape,SYS_ROD))
  Pm[Bead2Move-1].Pos[2] += IncrDist;
}
void Forces::PushBead(){
  Pm[Bead2Move].Pos[2] -= IncrDist;
  if(VAR_IF_TYPE(SysShape,SYS_ROD))
  Pm[Bead2Move-1].Pos[2] -= IncrDist;
}
void Forces::SelectBead(int p){
  if(p<0 && p > pNPart()) return;
  Bead2Move = p;
  if(pType(p) != 2)
    Pm[p].Typ = 1;
  if(pType(Old2Move) == 1)  
    Pm[Old2Move].Typ = 0;
  Old2Move = p;
}
void Forces::Wave(){
  for(int s=0;s<NEdge;s++){
    for(int ss=0;ss<NEdge;ss++){
      Pm[s*NEdge+ss].Pos[2] = .5+.1*sin(s*Dx*60.*INV_DUE_PI+Time)*sin(ss*Dx*30.*INV_DUE_PI+Time);
    }
  }
}
void Forces::AddCircle(int n){
  for(int p=0;p<pNPart();p++){
    double CosTh = (Nano[n].Pos[0] - pPos(p,0))/Nano[n].Rad;
    double SenTh = sqrt(1 - QUAD(CosTh));
    double Sign = 1.;
    for(int d=0;d<3;d++){
      Pm[p].Vel[d] += (Fm[p].Dir[d] + Fm[p].Ext[d])*pDeltat();
    }
    if(pType(p) != 2) Pm[p].Typ = 0;
    if(p <= NEdge){
      Sign *= -1.;
      if(Pm[p].Pos[2] + .5*Pm[p].Vel[2]*pDeltat() >= Sign*(Nano[n].Rad)*SenTh + Nano[n].Pos[2]){
	Pm[p].Pos[2] = Sign*(Nano[n].Rad)*SenTh + Nano[n].Pos[2];
	Pm[p].Vel[2] = 0.;
	Sign = 1.;
	if(Pm[p].Pos[0] > Nano[n].Pos[0])
	  Sign = -1.;
	if(p > 1)
	  Fm[p-1].Dir[0] -= Kf.Cont*SenTh*Sign;
	if(p < NEdge - 1)
	  Fm[p+1].Dir[0] += Kf.Cont*SenTh*Sign;
	if(Pm[p].Typ != 2)
	  Pm[p].Typ = 1;
      }
      else 
	Fm[p].Dir[0] = 0.;
    }
    else if( p> NEdge){
      if(Pm[p].Pos[2] + .5*Pm[p].Vel[2]*pDeltat() <= Sign*(Nano[n].Rad)*SenTh + Nano[n].Pos[2]){
	Pm[p].Pos[2] = Sign*(Nano[n].Rad)*SenTh + Nano[n].Pos[2];
	Pm[p].Vel[2] = 0.;
	Sign = 1.;
	if(Pm[p].Pos[0] > Nano[n].Pos[0])
	  Sign = -1.;
	if(p > NEdge)
	  Fm[p-1].Dir[0] -= Kf.Cont*SenTh*Sign;
	if(p < pNPart() - 1)
	  Fm[p+1].Dir[0] += Kf.Cont*SenTh*Sign;
	if(Pm[p].Typ != 2)
	  Pm[p].Typ = 1;
      }
      else 
	Fm[p].Dir[0] = 0.;
    }
  }
}
void Forces::AddCylinder(int n){
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ == 2) continue;
    double Dist = sqrt(SQR(Pm[p].Pos[0]-Nano[n].Pos[0])+SQR(Pm[p].Pos[1]-Nano[n].Pos[1]));
    if(Dist > Nano[n].Rad) continue;
    double Zed = (Nano[n].Rad-Dist)*tan(Nano[n].Hamaker*DUE_PI/360.);
    Pm[p].Typ = 1;
    Pm[p].Pos[2] = .5*Nano[n].Height + Nano[n].Pos[2] + Zed;
  }
}
/**
   The boundary condition depends on the lateral distance from the border.
 */
double Forces::HeightBoundary(double *Pos,int dir){
  return Pos[2];
  int nNano = 0;
  int dir2 = (dir+1)%2;
  for(int n=0;n<pNNano();n++){
    if(fabs(Pos[dir] - Nano[n].Pos[dir]) < Nano[n].Rad){
      nNano = n;
      break;
    }
  }
  double yNano = Pos[dir2] - Nano[nNano].Pos[dir2];
  double Lat1 = Nano[nNano].Pos[dir] - sqrt(SQR(Nano[nNano].Rad)-SQR(yNano) ) - Pos[dir];
  double Lat2 = Nano[nNano].Pos[dir] + sqrt(SQR(Nano[nNano].Rad)-SQR(yNano) ) - Pos[dir];
  double Dist = MIN(fabs(Lat1),fabs(Lat2));
  //double Dist = sqrt(SQR(Pos[0]-Nano[nNano].Pos[0])+SQR(Pos[1]-Nano[nNano].Pos[1]));
  //double Zed = (Nano[nNano].Rad-Dist)*tan(Nano[nNano].Hamaker*DUE_PI/360.);
  double Zed = Dist*tan(Nano[nNano].Hamaker*DUE_PI/360.);
  Pos[2] = .5*Nano[nNano].Height + Nano[nNano].Pos[2] + Zed;
  return .5*Nano[nNano].Height + Nano[nNano].Pos[2] + Zed;
}
void Forces::AddPore(int n){
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ == 2) continue;
    double Dist = sqrt(SQR(Pm[p].Pos[0]-Nano[n].Pos[0])+SQR(Pm[p].Pos[1]-Nano[n].Pos[1]));
    if(Dist > Nano[n].Rad) continue;
    Pm[p].Typ = 1;
    Pm[p].Pos[2] = Nano->Pos[2];
  }
}
void Forces::AddRigid(){
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ == 2) continue;
    Pm[p].Typ = 0;
  }
  for(int n=0;n<pNNano();n++){
    if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_SPH)){
      AddCircle(n);
    }
    else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CYL)){
      AddCylinder(n);
    }
    else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_PORE)){
      AddPore(n);
    }
  }
}
