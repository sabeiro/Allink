#include "Forces.h"
//-------------------------SPLINES--------------------------------
/** Use different types of splines to interpolate the points and calculate the second and forth derivative */
int Forces::ForceFieldLine(){
  double Sign=1.;
  SPLINE Sp1;
  SPLINE Sp2;
  SPLINE Cub1;
  SPLINE Cub2;
  SPLINE Par1;
  SPLINE Par2;
  double fSLap=0.;
  double fLap=0.;
  //Fm[Part2Move].fExt[2] = -Fm[Part2Move].fHel[2];
  //Fm[Part2Move].fExt[2] = 0.;
  for(int p=0;p<pNPart();p++){
    Fm[p].Ext[2] = -Fm[p].Dir[2]/1.1;
    if(IfInterp==FIT_SPLINE3){
      if( p == 0)
	Sp2 = Mat->Spline3Beg(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      else if( p < pNPart() - 2 )
	Sp2 = Mat->Spline3(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      else if( p == pNPart() - 2 )
	Sp2 = Mat->Spline3End(Pm[p].Pos,Pm[p+1].Pos,0,2);
      fSLap = 0.;//atan(Sp2.a1);
      fLap = 2.*Sp2.a2; 
    }
    else if(IfInterp==FIT_SPLINE4){
      if( p == 0){
	Sp2 = Mat->Spline4Beg(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,0,2);
	//Pm[p].Typ = 2;
      }
      else if( p < pNPart() - 3 )
	Sp2 = Mat->Spline4(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,0,2);
      else if( p == pNPart() - 3 ){
	Sp2 = Mat->Spline3(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      //Sp2 = Mate->Spline4PreEnd(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,1);
	//Pm[p].Typ = 2;
      }
      else if( p == pNPart() - 2 ){
	Sp2 = Mat->Spline3End(Pm[p].Pos,Pm[p+1].Pos,0,2);
      //Sp2 = Mate->Spline4End(Pm[p].Pos,Pm[p+1].Pos,0,1);
	//Pm[p].Typ = 2;
      }
      fSLap = (24.*Sp2.a4);
      fLap = (2.*Sp2.a2);
    }
//     else if(IfInterp==FIT_PARAB){
//       if(p < Gen->NPart - 2){
// 	Par1 = Mate->Parab(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
//       }
//       fSLap  = 0.;//(2.*Par2.A*Pm[p].Pos[0]+Par2.B);
//       fLap  = 2.*Par1.A/100.; 
//     }
    else if(IfInterp==FIT_PARAB2){
      if(p == 0){
	Sp2 = Mat->Parab(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      }
      if(p < pNPart() - 1){
	Sp2 = Mat->Parab2(Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,0,2);
      }
      fSLap = 0.;//(Par2.B);
      fLap = 2.*Sp2.a2;
    }
    else if(IfInterp==FIT_CUBIC){
      if(p == 0){
	Cub2 = Mat->Parab(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      }
      if(p < pNPart() - 2){
	Cub2 = Mat->Cubica(Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      }
      if(p == pNPart() - 2){
	Cub2 = Mat->Parab2(Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,0,2);
      }
      fSLap = 0.;
      fLap = 2.*Cub2.a2;
    }
    else if(IfInterp==FIT_FORTH){
      double Ref=0.;
    if(p < 1){
      //Sp1 = Mat->Forth(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,Pm[p+4].Pos,0,1);
      Sp2 = Mat->Forth(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,Pm[p+4].Pos,0,2);
      Ref = Pm[p].Pos[0] - Pm[p+2].Pos[0]; 
    }
    else if(p < 2){
      //Sp1 = Mat->Forth(Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,0,1);
      Sp2 = Mat->Forth(Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,0,2);
      Ref = Pm[p].Pos[0] - Pm[p+1].Pos[0]; 
    }
    else if(p < pNPart()-2){
      //Sp1 = Mat->Forth(Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,1);
      Sp2 = Mat->Forth(Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      Ref = 0.; 
    }
    else if(p == pNPart()-2){
      //Sp1 = Mat->Forth(Pm[p-3].Pos,Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,0,1);
      Sp2 = Mat->Forth(Pm[p-3].Pos,Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,0,2);
      Ref = Pm[p].Pos[0] - Pm[p-1].Pos[0]; 
    }
    else if(p <= pNPart()-1){
      //Sp1 = Mat->Forth(Pm[p-4].Pos,Pm[p-3].Pos,Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,0,1);
      Sp2 = Mat->Forth(Pm[p-4].Pos,Pm[p-3].Pos,Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,0,2);
      Ref = Pm[p].Pos[0] - Pm[p-2].Pos[0]; 
    }
    fSLap = -(24.*Sp2.a4)*QUAD(QUAD(Dx));
    fLap = (2.*Sp2.a2)*QUAD(Dx);
    //printf("%d) %lf %lf\n",p,Kf.kLap*fLap,Kf.kSLap*fSLap);
    }
    if(fLap > 100.)
      fLap = 100.;
    if(fSLap < -100.)
      fSLap = -100.;
    if(fSLap > 100.)
      fSLap = 100.;
    if(fSLap < -100.)
      fSLap = -100.;
    Fm[p].Dir[2] = Kf.SLap*fSLap + Kf.Lap*fLap;
    //if(Pm[p].Pos[2] <0) Sign = -1.;
    //else Sign = 1;
    //Fm[p].fHel[2] *= Sign;
    //if(p == Part2Move)
    //if(p > Part2Move -4 && p < Part2Move +4)
    {
      //printf("%d in (%lf,%lf) ha %lf|%lf Sp %lf %lf\n",p,Pm[p].Pos[2],Pm[p].Vel[2],Fm[p].fHel[2],Fm[p].fExt[2],fSLap,fLap);
    }
  }
  return 0;
}
/** Use different types of splines to interpolate the points and calculate the second and forth derivative, the two sheets have an elastic coupling */
int Forces::ForceFieldLeaves(){
  double Sign=1.;
  SPLINE Sp1;
  SPLINE Sp2;
  PARABOLA Par1;
  PARABOLA Par2;
  double fSLap=0.;
  double fLap=0.;
  for(int c=0;c<pNChain();c++){
    for(int p=c*NEdge;p<NEdge*(c+1);p++){
      Fm[p].Ext[2] = -(Fm[p].Dir[2])/1.1;
      //      Fm[p].El[2] = - Kf.El[2]*( c*Kf.Elong[2] + .5 - Pm[p].Pos[2]);
      if( p == c*NEdge){
	Sp2 = Mat->Spline4Beg(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,0,2);
	//Pm[p].Typ = 2;
      }
      else if( p < NEdge*(c+1) - 3 )
	Sp2 = Mat->Forth(Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      //Sp2 = Mat->Spline4(Pm[p].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
      //      else Pm[p].Typ = 2;
      else if( p == NEdge*(c+1) - 3 ){
	Sp2 = Mat->Spline3(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
	//Pm[p].Typ = 2;
      }
      else if( p == pNPart() - 2 ){
	Sp2 = Mat->Spline3End(Pm[p].Pos,Pm[p+1].Pos,0,2);
	//Pm[p].Typ = 2;
      }
      fSLap = -(24.*Sp2.a4)*QUAD(QUAD(Dx));
      fLap = (2.*Sp2.a2)*QUAD(Dx);
      if(fLap > 100.)fLap = 100.;
      if(fLap < -100.)fLap = -100.;
      if(fSLap > 100.)fSLap = 100.;
      if(fSLap < -100.)fSLap = -100.;
      Fm[p].Dir[2] = +Kf.SLap*fSLap + Kf.Lap*fLap;
      continue;
      {
	double Rad2=0.;
	for(int d=0;d<3;d++){
	  Rad2 += QUAD((Pm[p].Pos[d] - Nano->Pos[d]));
	}
	if( Rad2 > QUAD((3.*Nano->Rad))) continue;
	double CosTh = (Nano->Pos[0] - Pm[p].Pos[0])/Nano->Rad;
	double SenTh = sqrt(1 - QUAD(CosTh));
	double Sign = 1.;
	double Rad = sqrt(Rad2);
	double Strength = 12.*pow( (Nano->Rad/Rad) , 13) + 6.*pow( (Nano->Rad/Rad) , 7);
	if(p <= NEdge)
	  Sign *= -1.;
	Fm[p].Dir[0] = Kf.LJ * Strength * CosTh;
	Fm[p].Dir[2] = Kf.LJ * Strength * SenTh * Sign;
      }
    }
  }
}
//----------------------HARMONIC-POTENTIAL-------------------
/** Calculate the harmonic potential along the links */
int Forces::ForceFieldBulk(){
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Fm[p].Ext[d] = -(Fm[p].Dir[d])/2.;
    }
    for(int l=0;l<Ln[p].NLink;l++){
      int link = Ln[p].Link[l];
      for(int d=0;d<3;d++){
	double Delta = ASS((Pm[link].Pos[d] - Pm[p].Pos[d]));
	if(Delta < 0.001) continue;
	//	printf("%d-%d (%d) %lf ",p,link,d,Delta);
	double Sign = 1.;
	if(Pm[p].Pos[d] < Pm[link].Pos[d] )
	  Sign = -1.;
	Fm[p].Dir[d] +=  - Sign*Kf.El[d]*(Delta - Kf.Elong[d])/1.1;
      }
      //      printf("\n");
    }
    //if(p > Part2Move -4 && p < Part2Move +4)
      //printf("%d) in (%lf,%lf) ha %lf|%lf|%lf Sp %lf %lf\n",p,Pm[p].Pos[2],Pm[p].Vel[2],Fm[p].fHel[2],Fm[p].fExt[2],Fm[p].fEl[2]);
  }
}
/** Calculate the harmonic potential along the links */
void Forces::ForceFieldRod(){
  double DistRelBA[4];
  double DistRelBC[4];
  double Pre[12];
  for(int b=0;b<pNBlock();b++){
    for(int c=0;c<pNChain(b);c++){
      for(int p=c*pNPCh(b);p<(c+1)*pNPCh(b)-1;p++){
	TwoPartDist(p+1,p,DistRelBA);
	double ForceSp = pkSpr()*(1. - pSprRest()/DistRelBA[3]);
	for(int d=0;d<3;d++){
	  Fm[p].Dir[d] += ForceSp*DistRelBA[d];
	  Fm[p+1].Dir[d] -= ForceSp*DistRelBA[d];
	}
	if(p < (c+1)*pNPCh(b)-2){
	  TwoPartDist(p+2,p+1,DistRelBC);
	  double CosAngle = 0.;
	  for(int d=0;d<3;d++){
	    DistRelBA[d] /= DistRelBA[3];
	    DistRelBC[d] /= DistRelBC[3];
	    CosAngle += DistRelBA[d]*DistRelBC[d];
	  }
	  double PreFactBA = pkBen()/DistRelBA[3];
	  double PreFactBC = pkBen()/DistRelBC[3];
	  for(int d=0;d<3;d++){
	    Fm[p+0].Dir[d] += PreFactBA*(DistRelBC[d]-DistRelBA[d]*CosAngle);
	    Fm[p+1].Dir[d] -= PreFactBA*(DistRelBC[d]-DistRelBA[d]*CosAngle);
	    Fm[p+1].Dir[d] += PreFactBC*(DistRelBA[d]-DistRelBC[d]*CosAngle);
	    Fm[p+2].Dir[d] -= PreFactBC*(DistRelBA[d]-DistRelBC[d]*CosAngle);
	  }
	}	  
      }
    }
  }
}
/* obsolete? */
int Forces::SumSomeForces(int HowMany){
  double Sign=1.;
  SPLINE Par1;
  SPLINE Par2;
  double Slope1;
  double Slope2;
  if(Bead2Move - HowMany < 1 || Bead2Move + HowMany > pNPart() - 2)
    return 1;
  //Fm[Part2Move].fExt[2] = -Fm[Part2Move].fHel[2];
  //Fm[Part2Move].fExt[2] = 0.;
  int pp = Bead2Move;
  for(int h = 0;h<HowMany;h++){
    int link1 = Ln[pp].Link[0];
    int link2 = Ln[pp].Link[1];
//     printf("(%lf %lf %lf), (%lf %lf %lf), (%lf %lf %lf)\n",
// 	   Pm[link1].Pos[0],Pm[link1].Pos[1],Pm[link1].Pos[2],
// 	   Pm[pp].Pos[0],Pm[pp].Pos[1],Pm[pp].Pos[2],
// 	   Pm[link2].Pos[0],Pm[link2].Pos[1],Pm[link2].Pos[2]);
//    Par1 = Mat->Parab2(Pm[link1].Pos,Pm[pp].Pos,Pm[link2].Pos,0,2);
    Par2 = Mat->Parab2(Pm[link1].Pos,Pm[pp].Pos,Pm[link2].Pos,1,2);
    Fm[pp].Dir[2] = - Kf.Lap*QUAD(2.*Par1.a2);
    Fm[pp].Dir[2] += - Kf.Lap*QUAD(2.*Par2.a2);
    //printf("%d-%d|%d) %lf %lf %lf %lf  ",pp,link1,link2,Slope1,Par1.A,Slope2,Par2.A);
    link1 = Ln[pp].Link[2];
    link2 = Ln[pp].Link[3];
    Par1 = Mat->Parab2(Pm[link1].Pos,Pm[pp].Pos,Pm[link2].Pos,0,2);
    //Par2 = Mat->Parab2(Pm[link1].Pos,Pm[pp].Pos,Pm[link2].Pos,1,2);
    Fm[pp].Dir[2] += - Kf.Lap*QUAD(2.*Par1.a2);
    Fm[pp].Dir[2] += - Kf.Lap*QUAD(2.*Par2.a2);
    //printf("%d-%d|%d) %lf %lf %lf %lf\n",pp,link1,link2,Slope1,Par1.A,Slope2,Par2.A);
    if(Pm[pp].Pos[2] <0) Sign = -1.;
    else Sign = 1;
    Fm[pp].Dir[2] *= Sign;
    //printf("%d) Pos %lf Force %lf %lf\n",pp,Pm[pp].Pos[2],Fm[pp].fHel[2],Fm[pp].fExt[2]);
  }
  return 0;
}
//-----------------------------RIGID-------------------------
void Forces::ForceFieldRigid(){
  double Pot[2];
  for(int n=0;n<pNNano();n++){
    for(int d=0;d<3;d++){
      Nano[n].Force[d] = 0.;
      Nano[n].AMom[d] = 0.;
    }
  }
  for(int n=0;n<pNNano();n++){
    for(int nn=n+1;nn<pNNano();nn++){
      double dr[3] = {0.,0.,0.};
      double Dist = RigidDistanceAxis(n,nn,dr);
      //if( Dist > CutOff ) continue;
      double Cons = RigidLJ(n,Dist,Pot,0);
      for(int d=0;d<3;d++){
	Nano[n].Force[d] -= Cons*dr[d];
	Nano[nn].Force[d] += Cons*dr[d];
	Nano[n].AMom[d] -= Nano[n].AMomTemp[d]*Cons;
	Nano[nn].AMom[d] += Nano[n].AMomTemp[d]*Cons;
	//printf("%d %lf %lf\n",d,Cons*dr[d],Nano[nn].AMom[d]);
      }
    }
    for(int d=0;d<3;d++){
      double Ran = Nano[n].Zeta * (2.*Mat->Casuale() - 1.);
      double Dis = - Nano[n].Gamma*Nano[n].Vel[d];
      //Nano->Force[d] += Ran + Dis;
    }
    //printf("%d %lf %lf %lf\n",n,Nano[n].AMom[0],Nano[n].AMom[1],Nano[n].AMom[2]);
  }
}
double Forces::RigidDistanceAxis(int n,int nn,double *dr){
  Vettore Pos1(Nano[n].Pos,3);
  Vettore Pos2(Nano[nn].Pos,3);
  Vettore PosRel(Pos2[0]-Pos1[0],Pos2[1]-Pos1[1],Pos2[2]-Pos1[2]);
  Vettore Axis(Nano[n].Axis,3);
  Vettore Perp(3);
  Vettore AMomTemp(3);
  double Dist = Perp.PerpTo(&PosRel,&Axis);
  Perp.Export(dr);
  AMomTemp.VetV(&Perp,&PosRel);
  AMomTemp.Export(Nano[n].AMomTemp);
  //AMomTemp.Export(Nano[nn].AMomTemp);
  double HeiOnAxis = Perp.ProjOnAxis(&PosRel,&Axis);
   if( fabs(HeiOnAxis) > Nano[n].Height*.5){
    double Sign = HeiOnAxis > 0 ? 1. : -1.;
    Dist = 0.;
    for(int d=0;d<3;d++){
      dr[d] = Pos2[d] - Pos1[d];
      Dist = SQR(dr[d]);
    }
    Dist = sqrt(Dist);
  }
  return Dist;
}
void Forces::NanoInteraction(){
  double Pot[2];
  for(int n=0;n<pNNano();n++){
    for(int d=0;d<3;d++){
      Nano[n].Force[d] = 0.;
      Nano[n].AMom[d] = 0.;
    }
  }
  for(int n=0;n<pNNano();n++){
    Point2Shape(Nano[n].Shape);
    for(int p=0;p<pNPart();p++){
      double dr[3] = {0.,0.,0.};
      double Dist2 = NanoDist2(Pm[p].Pos,n);
      if( Dist2 > Kf.CutOff2 ) continue;
      double Cons = RigidLJ(n,Dist2,Pot,0);
      for(int d=0;d<3;d++){
	Nano[n].Force[d] -= Cons*dr[d];
	Fm[p].Dir[d] += Cons*dr[d];
	Nano[n].AMom[d] -= Nano[n].AMomTemp[d]*Cons;
      }
    }
    for(int d=0;d<3;d++){
      double Ran = Nano[n].Zeta * (2.*Mat->Casuale() - 1.);
      double Dis = - Nano[n].Gamma*Nano[n].Vel[d];
      //Nano->Force[d] += Ran + Dis;
    }
    //printf("%d %lf %lf %lf\n",n,Nano[n].AMom[0],Nano[n].AMom[1],Nano[n].AMom[2]);
  }
}
double Forces::NanoNrg(int p){
  double Nrg = NanoNrg(pPos(p),pType(p));
  return Nrg;
}
double Forces::NanoNrg(double *Pos,int t){
  double Pot[2];
  double Cons = 0.;
  double Nrg = 0.;
  for(int n=0;n<pNNano();n++){
    Point2Shape(Nano[n].Shape);
    double Dist2 = NanoDist2(Pos,n);
    if( Dist2 > SQR(Nano[n].CutOff) ) continue;
    double Dist = sqrt(Dist2);
    double Sign = -1.;
    if(t == 1) Sign = 1.;
    Cons += RigidHamaker(n,Dist,Pot,Sign);
    //Cons += RigidLJ(n,Dist2,Pot,0);
    Nrg += Pot[0];
  }
  Nrg *= Nano->Coating;
  return Nrg;
}
void Forces::DefForceParam(){
  double Pot[2];
  Kf.ForThr = 50.;
  Kf.DistThr = 0.;
  Potential(Kf.CutOff2,0,0,Pot);
  Kf.BaseLine = -Pot[0];
  double CutOff = sqrt(Kf.CutOff2);
  for(int i=10000;i>0;i--){
    double x = CutOff/10000.*i;
    double For = Potential(SQR(x),0,0,Pot)/x;
    if(For >= Kf.ForThr){
      Kf.DistThr = x;
      Kf.PotThr = Pot[0];
      break;
    }
  }
}
void Forces::DefNanoForceParam(){
  double Pot[2];
  for(int n=0;n<pNNano();n++){
    Nano[n].ForThr = 50.;
    Point2Shape(Nano[n].Shape);
    Nano[n].CutOff = Nano[n].Rad*3.;
    Nano[n].BaseLine = 0.;
    double For =  RigidLJ(n,Nano[n].CutOff,Pot,0);
    Nano[n].BaseLine = -Pot[0];
    Nano[n].DistThr = 0.;
    for(int i=10000;i>0;i--){
      double x = (Nano[n].CutOff)/10000.*i;
      For = RigidHamaker(n,x,Pot,0);
      //For = RigidLJ(n,x,Pot,0);
      if(For >= Nano[n].ForThr){
	Nano[n].DistThr = x;
	Nano[n].PotThr = Pot[0];
	break;
      }
    }
    //printf("DefForce %lf %lf %lf\n",Nano[n].DistThr,Nano[n].PotThr,Nano[n].BaseLine);
  }
}
double Forces::RigidDistanceRad(int n,int nn,double *dr){
  Vettore Pos1(Nano[n].Pos,3);
  Vettore Pos2(Nano[nn].Pos,3);
  double Dist = 0.;
  for(int d=0;d<3;d++){
    dr[d] = Pos2[d] - Pos1[d];
    Dist += SQR(dr[d]);
  }
  return sqrt(Dist);
}
double Forces::RigidLJ(int n,double Dist,double *Pot,double Sign){
  double Cons = 0.;
  double Strength = pow(Nano[n].Coating,9.)*Nano[n].Hamaker;
  double Hamaker  = pow(Nano[n].Coating,3.)*Nano[n].Hamaker;
  if( Dist <  Nano[n].DistThr){
    Pot[0] = Nano[n].PotThr - Nano[n].ForThr*(Dist - Nano[n].DistThr);
    Pot[1] = Nano[n].ForThr;
    Cons = Nano[n].ForThr;
    return Cons;
  }
  double idr = 1./(Dist-Nano[n].Rad+Nano[n].Coating);
  double idr2 = idr*idr;
  double idr4 = idr2*idr2;
  double idr6 = idr2*idr4;
  Pot[0] = idr*idr2*(Strength*idr6   +Sign*Hamaker)+Nano[n].BaseLine;
  Cons   = idr4    *(Strength*9.*idr6+Sign*3.*Hamaker);
  Pot[1] = idr4    *(Strength*9.*idr6+Sign*3.*Hamaker);
  return Cons;
}
/** Hamker Physica IV 10 p 1058 (1937) */
double Forces::RigidHamaker(int n,double Dr,double *Pot,double Sign){
  double SurfTens = 3.*Nano[n].Hamaker;
  double Cons = 0.;
  double Sigma  = 1.0;//0.04;
  double Slope = 1.00203;
  double Intercept = 0.31739;
  double Rnp = Nano[n].Rad/Slope-Intercept;
  if( Dr <  Nano[n].DistThr){
    Pot[0] = Nano[n].PotThr - Nano[n].ForThr*(Dr - Nano[n].DistThr);
    Pot[1] = Nano[n].PotThr;
    return Nano[n].ForThr;
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
  Pot[0] = PreRep*Rep + Sign*PreAttr*Attr + Nano[n].BaseLine;
  Pot[1] = PreRep*RepChem + Sign*PreAttr*AttrChem;
  Cons = PreRep*(RepP+RepM) + Sign*PreAttr*(AttrP+AttrM);
  return Cons;
}
double Forces::RigidCoulomb(int n,double Dist,double *Pot,double Sign){
  if(Dist < 2.*Nano[n].Rad)
    return 300.;
  Pot[0] = Sign*Nano[n].Hamaker/(Dist);
  return -Sign*Nano[n].Hamaker/SQR(Dist);
}
//----------------------------DEF-FORCES-POTENTIAL-------------
void Forces::PointShape(int iShape){
  Point2Shape(iShape);
  if(VAR_IF_TYPE(iShape,SHAPE_SPH))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_CYL))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_PILL))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_TILT))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_WALL))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_PORE))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_EXT)) 
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_JANUS))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_STALK))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_TIP))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_TORUS))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_HARM))
    CalcPot = &Forces::LJ39;
  if(VAR_IF_TYPE(iShape,SHAPE_UMBR))
    CalcPot = &Forces::LJ39;
  DefNanoForceParam();
}
double Forces::Harmonic(double Dist2,int t1,int t2,double *Pot){
  double Dist = sqrt(Dist2);
  double a = Kf.LJ;
  double b = -2.*a*Kf.LJMin;
  double c = -a*SQR(Kf.LJMin) - b*Kf.LJMin;
  Pot[2] = a*Dist2 + b*Dist + c;
  double For =  - 2.*a*Dist2 - b*Dist;
  if(t1 != t2){
    Pot[2] *= -1.;
    return -For;
  }
  //Pot[2] = 0.;
  return For;
}
double Forces::LJ39(double Dist2,int t1,int t2,double *Pot){
  double idr6 = CUBE(SQR(Kf.LJMin)/Dist2);
  double idr3 = CUBE(Kf.LJMin)/(Dist2*sqrt(Dist2));
  // if( Dist2 <  Kf.DistThr){
  //   Pot[0] = Kf.PotThr - Kf.ForThr*(Dist2 - Kf.DistThr);
  //   return Cons = Kf.ForThr;
  // }
  if(t1 == t2){
    Pot[2] = Kf.LJ*idr3*(idr6-1.)+Kf.BaseLine;
    return Kf.LJ*idr3*(9.*idr6-3.);
  }
  else{
    Pot[2] = Kf.LJ*idr3*(idr6+1.)+Kf.BaseLine;
    return Kf.LJ*idr3*(9.*idr6+3.);
  }
}
double Forces::LJPot(double Dist2,int t1,int t2,double *Pot){
  double idr6 = CUBE(SQR(Kf.LJMin)/Dist2);
  // if( Dist2 <  Kf.DistThr){
  //   Pot[0] = Kf.PotThr - Kf.ForThr*(Dist2 - Kf.DistThr);
  //   return Kf.ForThr;
  // }
  if(t1 == t2){
    Pot[2] = Kf.LJ*4.*(idr6*idr6-idr6) + Kf.BaseLine;
    return Kf.LJ*48.*(idr6*idr6-.5*idr6);
  }
  else{
    Pot[2] = Kf.LJ*4.*(idr6*idr6+idr6) + Kf.BaseLine;
    return Kf.LJ*48.*(idr6*idr6+.5*idr6);
  }
}
double Forces::StepPot(double Dist2,int t1,int t2,double *Pot){
  if(t1 != t2){
    Pot[2] = -Kf.LJ;
    return -Kf.LJ*sqrt(Dist2/Kf.CutOff2);
  }
  Pot[2] = Kf.LJ;
  return Kf.LJ*sqrt(Dist2/Kf.CutOff2);
}
double Forces::ElectroPot(double Dist2,int t1,int t2,double *Pot){
  if(t1 == t2)
    Pot[2] =  Kf.Ext/Dist2;
  else 
    Pot[2] =  -Kf.LJ/Dist2;
  return Pot[2];
}
void Forces::CalcForcesDensFunc(){
  if(!VAR_IF_TYPE(SysAlloc,ALL_FORCES)){
    printf("Forces not allocated\n");
    return;
  }
  int NTabF = 300;
  TabForceAlloc(NTabF);
  double *Count = (double *)calloc(NTab*pNType()*pNType(),sizeof(double));
  ClearDens();
  AddDens(0,pNPart());
  SumDens(0,pNPart());
  double Dist = 0.;
  double DistRelBA[4];
  double InvCutOff = 1./sqrt(Kf.CutOff2);
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Fm[p].Dir[d] = 0.;
    }
  }
  // non bonded
  for(int p1=0;p1<pNPart();p1++){
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      if(p2 <= p1) continue;
      Pc->Dist2Curr(DistRelBA);
      if(DistRelBA[3] > Kf.CutOff2) continue;
      double Force = 0.;
      double Dist = sqrt(DistRelBA[3]);
      for(int t=0;t<pNType();t++){
	Force += MInt->Coeff(Pm[p1].Typ,Pm[p2].Typ,t)*(Dens3[p1*pNType()+t]+Dens3[p2*pNType()+t]);
      }
      Force *= DerWei3(Dist,pWei3Par())*2./3.;
      Force += DerWei2(Dist,pWei2Par())*MInt->Coeff(Pm[p1].Typ,Pm[p2].Typ);
      Force /= -Dist;
      int x = (int)(Dist*InvCutOff*NTab);
      if( x < 0 || x >= NTab) continue;
      int t = MIN(Pm[p1].Typ,Pm[p2].Typ)*pNType()+MAX(Pm[p1].Typ,Pm[p2].Typ);
      FTab[x*pNType()*pNType()+t] += Force;
      Count[x*pNType()*pNType()+t] += 1.;
      for(int d=0;d<3;d++){
	Fm[p1].Dir[d] += Force*DistRelBA[d];
	Fm[p2].Dir[d] -= Force*DistRelBA[d];
      }
    }
  }
  for(int t=0;t<NTab;t++){
    for(int t1=0;t1<pNType();t1++){
      for(int t2=0;t2<pNType();t2++){
	int tType = MIN(t1,t2)*pNType()+MAX(t1,t2);
	double Norm = Count[t*pNType()*pNType()+tType] > 0. ? 1./Count[t*pNType()*pNType()+tType] : 1.;
	FTab[t*pNType()*pNType()+tType] *= Norm;
      }
    }
  }
  // bonded
  double DistRelBC[4];
  double Pre[12];
  for(int b=0;b<pNBlock();b++){
    for(int c=0;c<pNChain(b);c++){
      for(int p=c*pNPCh(b);p<(c+1)*pNPCh(b)-1;p++){
	TwoPartDist(p+1,p,DistRelBA);
	double ForceSp = pkSpr()*(1. - pSprRest()/DistRelBA[3]);
	for(int d=0;d<3;d++){
	  Fm[p].Dir[d] += ForceSp*DistRelBA[d];
	  Fm[p+1].Dir[d] -= ForceSp*DistRelBA[d];
	}
	if(p < (c+1)*pNPCh(b)-2){
	  TwoPartDist(p+2,p+1,DistRelBC);
	  double CosAngle = 0.;
	  for(int d=0;d<3;d++){
	    DistRelBA[d] /= DistRelBA[3];
	    DistRelBC[d] /= DistRelBC[3];
	    CosAngle += DistRelBA[d]*DistRelBC[d];
	  }
	  double PreFactBA = pkBen()/DistRelBA[3];
	  double PreFactBC = pkBen()/DistRelBC[3];
	  for(int d=0;d<3;d++){
	    Fm[p+0].Dir[d] -= PreFactBA*(DistRelBC[d]-DistRelBA[d]*CosAngle);
	    Fm[p+1].Dir[d] += PreFactBA*(DistRelBC[d]-DistRelBA[d]*CosAngle);
	    Fm[p+1].Dir[d] -= PreFactBC*(DistRelBA[d]-DistRelBC[d]*CosAngle);
	    Fm[p+2].Dir[d] += PreFactBC*(DistRelBA[d]-DistRelBC[d]*CosAngle);
	  }
	}	  
      }
    }
  }
}
void Forces::PrintForce(){
  double NPoint = 1000.;
  FILE *FWrite = fopen("Force.dat","w");
  double Pot[2];
  double Delta = sqrt(Kf.CutOff2)/NPoint;
  double CutOff = sqrt(Kf.CutOff2);
  for(double x=Delta;x<CutOff;x+=Delta){
    double Cons = LJPot(x,0,0,Pot);
    fprintf(FWrite,"%lf %lf %lf\n",x,Cons,*Pot);
  }
  fclose(FWrite);
}
void Forces::TabForceAlloc(int NTabExt){
  NTab = NTabExt;
  double Pot[2];
  double Delta = sqrt(Kf.CutOff2)/(double)NTab;
  double NrgLim = 20.;
  FTab = (double *)calloc(NTab*pNType()*pNType(),sizeof(double));
  MTab = (double *)calloc(NTab,sizeof(double));
  PTab = (double *)calloc(NTab*pNType()*pNType(),sizeof(double));
  VAR_ADD_TYPE(DynFlag,ALL_FORCES);
}
void Forces::TabPot(){
  NTab = 10000;
  double NrgLim = 10.;
  double Delta = 2.*NrgLim/(double)NTab;
  int i = 0;
  for(double Nrg=-NrgLim;Nrg<NrgLim;Nrg+=Delta,i++){
    PTab[i] = Nrg;
    MTab[i] = exp(-pTemp()*Nrg);
  }
  VAR_ADD_TYPE(DynFlag,ALL_METR);
  VAR_ADD_TYPE(DynFlag,ALL_POT);
}
double Forces::Wei3(const double r, const double b){
  if(r < 0.0001) return 0.;
  return (r < b) ? 15. / (2. * M_PI * CUBE(b)) * SQR(1. - r / b) : 0.;
}
double Forces::DerWei3(const double r, const double b){
  if(r < 0.0001) return 0.;
  return (r < b) ? 15. / (M_PI * SQR(SQR(b))) * (r / b - 1.) : 0.;
}
double Forces::Wei2(const double r, const double a){
  if(r < 0.0001) return 0.;
  const double b = 1.;
  double Ai = .5 * 15. / M_PI / (CUBE(a) + CUBE(a + b) + CUBE(b));
  return (r < a) ? Ai : (Ai / CUBE(a - b)) *
    ((((-2. * r) + 3. * (a + b)) * r - 6. * a * b) * r +
     3. * a * b * b - b * b * b);
}
double Forces::DerWei2(const double r, const double a){
  if(r < 0.0001) return 0.;
 const double b = 1.;
  double Ai = 15. / (M_PI * (((((4. * a) + 6.) * a) + 6.) * a + 4.));
  return (r > a && r < b) ? (6. * Ai / CUBE(a - b) * (-r * r + (a + b) *
						      r - a * b)) : 0.; 
}
