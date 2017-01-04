#include "../include/Matematica.h"

Piano::Piano(Vettore *P1e,Vettore *P2e,Vettore *P3e){
  P1.Copy(P1e);
  P2.Copy(P2e);
  P3.Copy(P3e);
  Dir21 = P2 - P1;
  Dir31 = P3 - P1;
  Dir23 = P2 - P3;
  for(int d=0;d<3;d++){
    Dir21[d] = P2[d] - P1[d];
    Dir31[d] = P3[d] - P1[d];
    Dir23[d] = P2[d] - P3[d];
  }
  Vettore P4e = P1 + Dir21 + Dir31;
  P4.Copy(&P4e);
  Norm = Dir21 ^ Dir31;
  Norm[0] = Dir21.x[1]*Dir31.x[2] - Dir21.x[2]*Dir31.x[1];
  Norm[1] = Dir21.x[2]*Dir31.x[0] - Dir21.x[0]*Dir31.x[2];
  Norm[2] = Dir21.x[0]*Dir31.x[1] - Dir21.x[1]*Dir31.x[0];
  Norm.Normalize();
  // printf("------\n");
  // P1.Print();
  // P2.Print();
  // Dir21.Print();
  // Norm.Print();
  for(int d=0;d<3;d++){
    if( fabs(Norm[d]) > 0.){
      InvNorm[d] = 1./Norm[d];
      IsInf[d] = 0;
    }
    else{
      InvNorm[d] = 0.;
      IsInf[d] = 1;
    }
  }
  for(int d=0;d<3;d++){
    Bound[d*2    ] = MIN(P1[d],MIN(P2[d],P3[d]));
    //Bound[d*2    ] = MIN(P1[d],MIN(P2[d],MIN(P3[d],P4[d])));
    Bound[d*2+1] = MAX(P1[d],MAX(P2[d],P3[d]));
    //Bound[d*2+1] = MAX(P1[d],MAX(P2[d],MAX(P3[d],P4[d])));
  }
  dPar = - P1[0]*Norm[0] - P1[1]*Norm[1] - P1[2]*Norm[2];
  mxy[0] = (P1[1] - P2[1])*Inv(P1[0] - P2[0]);
  qxy[0] = P1[1] - mxy[0]*P1[0];
  mxz[0] = (P1[2] - P2[2])*Inv(P1[0] - P2[0]);
  qxz[0] = P1[2] - mxz[0]*P1[0];
  mxy[1] = (P1[1] - P3[1])*Inv(P1[0] - P3[0]);
  qxy[1] = P1[1] - mxy[1]*P1[0];
  mxz[1] = (P1[2] - P3[2])*Inv(P1[0] - P3[0]);
  qxz[1] = P1[2] - mxz[1]*P1[0];
  mxy[2] = (P3[1] - P2[1])*Inv(P3[0] - P2[0]);
  qxy[2] = P1[1] - mxy[2]*P1[0];
  mxz[2] = (P3[2] - P2[2])*Inv(P3[0] - P2[0]);
  qxz[2] = P1[2] - mxz[2]*P1[0];
  Rad = 1.;
}
Piano::~Piano(){}
//--------------------------------------------------------------
double Piano::Distance(Vettore *Pos) {
  double Dist = fabs(Norm.x[0]*Pos->x[0] + Norm.x[1]*Pos->x[1] + Norm.x[2]*Pos->x[2] + dPar);
  return Dist;
}
Vettore Piano::ProjOnSurf(Vettore *Pos){
  double Dist = Distance(Pos);
  Vettore Pos1(Pos->NDim);
  double Sign = 0.;
  for(int d=0;d<3;d++){
    Sign += (P1[d] - Pos->x[d])*Norm[d];
  }
  if(Sign > 0) Sign = -1.;
  else Sign = 1.;
  for(int d=0;d<3;d++){
    Pos1[d] = Pos->x[d] - Sign*Dist*Norm.x[d];
  }
  return Pos1;
}
Vettore Piano::ProjOnNorm(Vettore *v){
  double Len = 0.;
  for(int d=0;d<v->NDim;d++){
    Len += v->x[d]*Norm[d];
  }
  Vettore Scal(Len*Norm[0],Len*Norm[1],Len*Norm[2]);
  return Scal;
}
Vettore Piano::Reflect(Vettore *v){
  // Vettore v2(v->NDim);
  // v2.Copy(v);
  Vettore n1 = ProjOnNorm(v);
  // Vettore v1 = 2.*(n1 - v2);
  Vettore v1(3);
  for(int d=0;d<3;d++){
    v1[d] = - v->x[d] + 2.*n1[d];
  }
  return v1;
}
int Piano::Impact(Vettore *Pos,Vettore *Vel) {
  double Dist = Distance(Pos);
  if(Dist > Rad) return 0;
  Vettore PosS = ProjOnSurf(Pos);
  if(!IsOnSurf(&PosS)) return 0;
  Vettore Vel1 = Reflect(Vel);
  Vel->Copy(&Vel1);
}
Vettore Piano::GetVertex(int i) {
  if(i == 1) return P1;
  if(i == 2) return P2;
  if(i == 3) return P3;
  if(i == 4) return P4;
  return P1;
}
double Piano::Inv(double x){
  return x != 0. ? 1./x : 100000000000.;
}
int Piano::IsOnSurf(Vettore *P){
  //return IsOnSurf1(P);
  return IsOnSurf2(P);
}
int Piano::IsOnSurf1(Vettore *P){
  if( SameSide(P,&P1,&P2,&P3) && SameSide(P,&P2,&P1,&P3) && SameSide(P,&P3,&P1,&P2)) return 1;
  return 0;
}
int Piano::SameSide(Vettore *P1,Vettore *A1,Vettore *B1,Vettore *C1){
  Vettore P(3); P.Copy(P1);
  Vettore A(3); A.Copy(A1);
  Vettore B(3); B.Copy(B1);
  Vettore C(3); C.Copy(C1);
  Vettore D1 = (C-B)^(P-B);
  Vettore D2 = (C-B)^(A-B);
  double Scal = 0.;
  for(int d=0;d<3;d++){
    int NUno = (d+1)%3;
    int NDue = (d+2)%3;
    D1[d] = (C[NUno]-B[NUno])*(P[NDue]-B[NDue]) - (C[NDue]-B[NDue])*(P[NUno]-B[NUno]);
    D2[d] = (C[NUno]-B[NUno])*(A[NDue]-B[NDue]) - (C[NDue]-B[NDue])*(A[NUno]-B[NUno]);
    Scal += D1[d]*D2[d];
  }
  if( Scal >= 0) return 1;
  return 0;
  if( D1%D2 >= 0) return 1;
  return 0;
}
// TOFIX!!!!
int Piano::IsOnSurf2(Vettore *PA){
  Vettore P(3);P.Copy(PA);
  Vettore DirP1 = P - P1;
  for(int d=0;d<3;d++){
    DirP1[d] = P[d] - P1[d];
  }
  // double Inv = (Dir21%Dir21)*(Dir31%Dir31) - (Dir21%Dir31)*(Dir31%Dir21);
  // double v = ( (Dir31%Dir31)*(DirP1%Dir31) - (Dir31%Dir21)*(DirP1%Dir31) )/Inv;
  // double u = ( (Dir21%Dir21)*(DirP1%Dir21) - (Dir21%Dir31)*(DirP1%Dir21) )/Inv;
  double sq2 = 0.,sq3 = 0.;
  double pd32 = 0., pdP3 = 0., pdP2 = 0.;
  for(int d=0;d<3;d++){
    sq2 += SQR(Dir21[d]);
    sq3 += SQR(Dir31[d]);
    pd32 += Dir21[d]*Dir31[d];
    pdP3 += Dir31[d]*DirP1[d];
    pdP2 += Dir21[d]*DirP1[d];
  }
  double Inv = 1./(sq2*sq3 - SQR(pd32));
  double v = (sq3*pdP2 - pd32*pdP3)*Inv;
  double u = (sq2*pdP3 - pd32*pdP2)*Inv;
  if(u >= 0. && v >= 0. && (u+v) < 1.) return 1;
  else return 0;
}
