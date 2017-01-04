#include "Forces.h"
void Forces::AllocTens(){
  if(VAR_IF_TYPE(SysAlloc,ALL_TENS)) return;
  for(int d=0;d<3;d++){
    Tens.RefPos[d] = pNanoPos(0,d);//.5*pEdge(d);
  }
  Tens.Pre = (double **)calloc(Tens.NComp,sizeof(double));
  Tens.Dens = (double **)calloc(2,sizeof(double));
  //   2d
  if(VAR_IF_TYPE(Tens.CalcMode,CALC_2d)){
    Tens_Ref = &Forces::TensRefPol;
    SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
    Tens.Edge[0] = pEdge(3);
    Tens.Edge[1] = pEdge(CNorm);
    Tens.Edge[2] = 1.;
    Tens.Wrap[0] = 0;
    Tens.Wrap[1] = 1;
    for(int d=0;d<3;d++){
      Tens.EdgeInv[d] = 1./Tens.Edge[d];
    }
    for(int c=0;c<Tens.NComp;c++){
      Tens.Pre[c] = (double *)calloc(SQR(Tens.NSlab),sizeof(double));
    }
    for(int c=0;c<2;c++){
      Tens.Dens[c] = (double *)calloc(SQR(Tens.NSlab),sizeof(double));
    }
  }
  //  3d
  if(VAR_IF_TYPE(Tens.CalcMode,CALC_3d)){
    Tens_Ref = &Forces::TensRefCart;
    for(int d=0;d<3;d++){
      Tens.Edge[d] = pEdge(d);
      Tens.Wrap[d] = 1;
      Tens.EdgeInv[d] = pInvEdge(d);
    }
    for(int c=0;c<Tens.NComp;c++){
      Tens.Pre[c] = (double *)calloc(CUBE(Tens.NSlab),sizeof(double));
    }
    for(int c=0;c<2;c++){
      Tens.Dens[c] = (double *)calloc(CUBE(Tens.NSlab),sizeof(double));
    }
  }
  VAR_ADD_TYPE(SysAlloc,ALL_TENS);
}
void Forces::CalcDens(){
  // 3d
  if(VAR_IF_TYPE(Tens.CalcMode,CALC_3d)){
    double PosRel[3];
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
	PosRel[d] = pPos(p,d) - (pEdge(d)*.5 - Tens.RefPos[d]);
	PosRel[d] -= floor(PosRel[d]*pInvEdge(d))*pEdge(d);
      }
      int sx = (int)(PosRel[0]*pInvEdge(0)*Tens.NSlab);
      int sy = (int)(PosRel[1]*pInvEdge(1)*Tens.NSlab);
      int sz = (int)(PosRel[2]*pInvEdge(2)*Tens.NSlab);
      int v = (sx*Tens.NSlab+sy)*Tens.NSlab+sz;
      int t = pType(p);
      if(t > 1) continue;
      Tens.Dens[t][v] += 1.;
    }
  }
  // 2d
  else if(VAR_IF_TYPE(Tens.CalcMode,CALC_2d)){
    double PosRel[3];
    for(int p=0;p<pNPart();p++){
      PosRel[CLat1] = pPos(p,CLat1) - Tens.RefPos[CLat1];
      PosRel[CLat2] = pPos(p,CLat2) - Tens.RefPos[CLat2];
      PosRel[CNorm] = pPos(p,CNorm) + (pEdge(CNorm)*.5 - Tens.RefPos[CNorm]);
      for(int d=0;d<3;d++){
	PosRel[d] -= floor(PosRel[d]*pInvEdge(d))*pEdge(d);
      }
      int sr = (int)(hypot(PosRel[CLat1],PosRel[CLat2])*Tens.EdgeInv[0]*Tens.NSlab);
      int sz = (int)(PosRel[CNorm]*pInvEdge(CNorm)*Tens.NSlab);
      if(sr < 0 || sr >= Tens.NSlab) continue;
      if(sz < 0 || sz >= Tens.NSlab) continue;
      int v = (sr*Tens.NSlab+sz);
      int t = pType(p);
      if(t > 1) continue;
      Tens.Dens[t][v] += 1.;
      //printf("%d %d %d %lf\n",sr,sz,v,Tens.Dens[t][v]);
    }
  }
}
void Forces::CalcTens(){
  if(!VAR_IF_TYPE(SysAlloc,ALL_FORCES)){
    printf("Forces not allocated\n");
    return;
  }
  ClearDens();
  AddDens(0,pNPart());
  SumDens(0,pNPart());
  double Dist = 0.;
  double DistRelBA[4];
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Fm[p].Dir[d] = 0.;
    }
  }
  // non bonded
  double Pos[3];
  for(int p1=0;p1<pNPart();p1++){
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      if(p2 <= p1) continue;
      Pc->Dist2Curr(DistRelBA);
      if(DistRelBA[3] > Kf.CutOff2) continue;
      double Dist = sqrt(DistRelBA[3]);
      double Force = 0.;
      for(int t=0;t<pNType();t++){
	Force += MInt->Coeff(pType(p1),pType(p2),t)*(Dens3[p1*pNType()+t]+Dens3[p2*pNType()+t]);
      }
      Force *= DerWei3(Dist,pWei3Par())*2./3.;
      Force += DerWei2(Dist,pWei2Par())*MInt->Coeff(pType(p1),pType(p2));
      SumTens(p1,p2,Force,DistRelBA);
      Force /= -Dist;
      SigErr(Force > 5000.,"Forces over the limit %lf\n",Force);
      for(int d=0;d<3;d++){
	Fm[p1].Dir[d] += Force*DistRelBA[d];
	Fm[p2].Dir[d] -= Force*DistRelBA[d];
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
	SumTens(p,p+1,ForceSp,DistRelBA);
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
	    Pre[d  ] = DistRelBA[d]*pkBen()*(DistRelBC[d]-DistRelBA[d]*CosAngle);
	    Pre[d+6] = DistRelBC[d]*pkBen()*(DistRelBA[d]-DistRelBC[d]*CosAngle);
	  }
	  Pre[ 3] = DistRelBA[0]*pkBen()*(DistRelBC[1]-DistRelBA[1]*CosAngle);
	  Pre[ 4] = DistRelBA[0]*pkBen()*(DistRelBC[2]-DistRelBA[2]*CosAngle);
	  Pre[ 5] = DistRelBA[1]*pkBen()*(DistRelBC[2]-DistRelBA[2]*CosAngle);
	  Pre[ 9] = DistRelBC[0]*pkBen()*(DistRelBA[1]-DistRelBC[1]*CosAngle);
	  Pre[10] = DistRelBC[0]*pkBen()*(DistRelBA[2]-DistRelBC[2]*CosAngle);
	  Pre[11] = DistRelBC[1]*pkBen()*(DistRelBA[2]-DistRelBC[2]*CosAngle);
	  SumTens(p,p+1,Pre);
	  SumTens(p+1,p+2,Pre+6);
	}
      }
    }
  }
  return;
  //external
  double Pot[3];
  double PosBf[3];
  double dr[4];
  double NPos[3];
  for(int n=0;n<pNNano();n++){
    Point2Shape(Nano[n].Shape);
    for(int p=0;p<pNPart();p++){
      pPos(p,Pos);
      double Dr2 = NanoDist2(Pos,n);
      double InvDist = 1./Dr2;
      double Cons = Potential(Dr2,0,pType(p),Pot);
      for(int d=0;d<3;d++){
	dr[d] = Nano[n].Pos[d] - pPos(p,d);
	if(dr[d] >  .5*pInvEdge(d)) dr[d] -= pEdge(d);
	if(dr[d] < -.5*pInvEdge(d)) dr[d] += pEdge(d);
      }
      double Norm = SQR(Nano[n].Rad)/(SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]));
      Norm = sqrt(Norm);
      for(int d=0;d<3;d++){
	NPos[d] = NPos[d] + dr[d]*Norm;
	Fm[p].Dir[d] += Cons*dr[d]*InvDist;
	Pre[d  ] = Cons*dr[d]*dr[d]*InvDist;
      }
      Pre[3] = Cons*dr[0]*dr[1]*InvDist;
      Pre[4] = Cons*dr[0]*dr[2]*InvDist;
      Pre[5] = Cons*dr[1]*dr[2]*InvDist;
      SumTens(NPos,Pos,Pre);
    }
  }
}
double Forces::TensRefCart(double *Pos1,double *Pos2,double *PosP1,double *PosP2){
  for(int d=0;d<3;d++){
    PosP1[d] = Pos1[d] - (pEdge(d)*.5 - Tens.RefPos[d]);
    PosP2[d] = Pos2[d] - (pEdge(d)*.5 - Tens.RefPos[d]);
    PosP1[d] -= floor(PosP1[d]*pInvEdge(d))*pEdge(d);
    PosP2[d] -= floor(PosP2[d]*pInvEdge(d))*pEdge(d);
  }
}
double Forces::TensRefPol(double *PosO1,double *PosO2,double *PosP1,double *PosP2){
  double Pos1[2] = {PosO1[CLat1]-Tens.RefPos[CLat1],PosO1[CLat2]-Tens.RefPos[CLat2]};
  double Pos2[2] = {PosO2[CLat1]-Tens.RefPos[CLat1],PosO2[CLat2]-Tens.RefPos[CLat2]};
  for(int d=0;d<2;d++){
    Pos1[d] -= floor(Pos1[d]*pInvEdge(d) + .5)*pEdge(d);
    Pos2[d] -= floor(Pos2[d]*pInvEdge(d) + .5)*pEdge(d);
  }
  PosP1[0] = hypot(Pos1[0],Pos1[1]);
  PosP2[0] = hypot(Pos2[0],Pos2[1]);
  PosP1[1] = PosO1[CNorm] + (pEdge(CNorm)*.5 - Tens.RefPos[CNorm]);
  PosP2[1] = PosO2[CNorm] + (pEdge(CNorm)*.5 - Tens.RefPos[CNorm]);
  for(int d=1;d<2;d++){
    PosP1[d] -= floor(PosP1[d]*Tens.EdgeInv[d])*Tens.Edge[d];
    PosP2[d] -= floor(PosP2[d]*Tens.EdgeInv[d])*Tens.Edge[d];
  }
  //printf("%d %d %lf %lf %lf %lf \n",p1,p2,PosP1[0],PosP1[1],PosP2[0],PosP2[0]);
}
void Forces::SumTens(int p1,int p2,double *Pre){
  SumTens(Pm[p1].Pos,Pm[p2].Pos,Pre);
}
void Forces::SumTens(double *Pos1,double *Pos2,double *Pre){
  double PosP1[3];
  double PosP2[3];
  TensRef(Pos1,Pos2,PosP1,PosP2);
  int NPoint = 100;
  double PointInv = 1./(double)NPoint;
  double Deltav[3];
  double Pos[3];
  int vCurr[3];
  for(int d=0;d<Tens.NDim;d++){
    // starts from the smallest position
    Pos[d] = PosP1[d] < PosP2[d] ? PosP1[d] : PosP2[d];
    Deltav[d] = fabs((PosP2[d]-PosP1[d])*PointInv);
    // forward/backward?
    if( fabs(PosP2[d] - PosP1[d]) > Tens.Edge[d]*.5){
      Deltav[d] = -fabs((MAX(PosP1[d],PosP2[d])-Tens.Edge[d]-Pos[d])*PointInv);
    }
  }
  //printf("%d %d) [%lf %lf] [%lf %lf] (%lf %lf) (%lf %lf)\n",p1,p2,PosP1[0],PosP1[1],PosP2[0],PosP2[1],Deltav[0],Deltav[1],Tens.Edge[0],Tens.Edge[1]);
  for(int p=0;p<NPoint;p++){
    for(int d=0;d<Tens.NDim;d++){
      vCurr[d] = (int)(Pos[d]*Tens.NSlab*Tens.EdgeInv[d]);
      //assert(vCurr[d] >= 0 && vCurr[d] < Tens.NSlab);
      if(vCurr[d] < 0 || vCurr[d] >= Tens.NSlab) return;
      // update the current position forward/backward
      Pos[d] += Deltav[d];
      // wrap the position to the other end
      if(Pos[d] > Tens.Edge[d]){
	if(!Tens.Wrap[d]) return;
	Pos[d] -= floor(Pos[d]*Tens.EdgeInv[d])*Tens.Edge[d];
      }
    }
    int vTot = 0;
    for(int d=0;d<Tens.NDim;d++){
      vTot += vCurr[d];
      if(d < Tens.NDim-1)
	vTot *= Tens.NSlab;
    }
    //printf("%d %d %d (%lf %lf) %lf\n",vCurr[0],vCurr[1],vTot,Pos[0],Pos[1],Pre[0]);
    for(int c=0;c<Tens.NComp;c++){
      Tens.Pre[c][vTot] += Pre[c]*PointInv;
    }
  }
}
void Forces::SumTens(int p1,int p2,double Force,double *DistRel){
  if(fabs(Force) < 0.) return;
  if(fabs(Force) > 5000.) return;
  double Pre[6];
  double InvDist = Force/DistRel[3];
  Pre[0] = DistRel[0]*DistRel[0]*InvDist;
  Pre[1] = DistRel[1]*DistRel[1]*InvDist;
  Pre[2] = DistRel[2]*DistRel[2]*InvDist;
  Pre[3] = DistRel[0]*DistRel[1]*InvDist;
  Pre[4] = DistRel[0]*DistRel[2]*InvDist;
  Pre[5] = DistRel[1]*DistRel[2]*InvDist;
  SumTens(p1,p2,Pre);
}
void Forces::WriteTens2d(FILE *FWrite,int Comp,double InvNFile){
  int link[4] = {0,0,0,0};
  fprintf(FWrite,"# l(%.1f %.1f %.1f) r(%.2f %.2f %.2f) v[%d] d[color]\n",Tens.Edge[0],Tens.Edge[1],1.,Tens.RefPos[0],Tens.RefPos[1],Tens.RefPos[2],Tens.NSlab);
  double *VolContr = (double *)calloc(Tens.NSlab,sizeof(double));
  VolumeCircSlab(VolContr,Tens.NSlab);
  for(int sr=1,p=0,c=0;sr<Tens.NSlab;sr++){
    double r = sr*Tens.Edge[0]/(double)Tens.NSlab;
    for(int sz=0;sz<Tens.NSlab;sz++){
      double z = sz*Tens.Edge[1]/(double)Tens.NSlab;
      int v = sr*Tens.NSlab+sz;
      if(Tens.Dens[0][v] <= 0. && Tens.Dens[1][v] <= 0. && fabs(Tens.Pre[Comp][v]) <= 0.) continue;
      double Press = -Tens.Pre[Comp][v]*InvNFile/VolContr[sr];
      double Dens1 = Tens.Dens[0][v]*InvNFile/VolContr[sr];
      double Dens2 = Tens.Dens[1][v]*InvNFile/VolContr[sr];
      fprintf(FWrite,"{x(%.3f %.3f %.3f) v( %lf %.2f %.2f)}\n",r,z,0.,Press,Dens1,Dens2);continue;
      link[0] = (sr+0)*Tens.NSlab+(sz+0);
      link[1] = (sr+0)*Tens.NSlab+(sz+1);
      link[2] = (sr+1)*Tens.NSlab+(sz+0);
      link[3] = (sr+1)*Tens.NSlab+(sz+1);
      for(int lx=0;lx<2;lx++){
      	for(int ly=0;ly<2;ly++){
      	  int l = 2*lx+ly;
      	  int l1 = p + (p+1)%4;
      	  int l2 = p + (p+2)%4;
      	  int l3 = p + (p+3)%4;
      	  fprintf(FWrite,"{t[%d %d %d]",p,c,0);
      	  fprintf(FWrite," x(%.3f %.3f %.3f)",r,z,Press);
      	  fprintf(FWrite," v( %lf %.2f %.2f)",Press,Dens1,Dens2);
      		  // -Tens.Pre[Comp][v]*InvNFile,
      		  // Tens.Dens[0][v]*InvNFile,
      		  // Tens.Dens[1][v]*InvNFile);	
      	  fprintf(FWrite," l[%d] l[%d] l[%d]}\n",l1,l2,l3);
      	  p++;
      	}
      }
      c++;
    }
  }
  //free(VolContr);
}
void Forces::WriteTens(char *FTens,int Comp,double InvNFile){
  double InvValues = 1. / Tens.NSlab;
  double InvVolume = prho()*CUBE(Tens.NSlab)/pVol();
  if(VAR_IF_TYPE(Tens.CalcMode,CALC_2d)){
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,Tens.Pre[Comp],SQR(Tens.NSlab),MPI_DOUBLE,MPI_SUM,Proc->CommGrid);  
  for(int t=0;t<2;t++){
    MPI_Allreduce(MPI_IN_PLACE,Tens.Dens[t],SQR(Tens.NSlab),MPI_DOUBLE,MPI_SUM,Proc->CommGrid);
  }
 int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    FILE *FWrite = fopen(FTens,"w");
    WriteTens2d(FWrite,Comp,InvNFile);
    fclose(FWrite);
    for(int s=0;s<SQR(Tens.NSlab);s++){
      Tens.Pre[Comp][s] = 0.;
      Tens.Dens[0][s] = 0.;
      Tens.Dens[1][s] = 0.;      
    }
    #ifdef OMPI_MPI_H
        }
    #endif
  }
  if(VAR_IF_TYPE(Tens.CalcMode,CALC_3d)){
#ifdef OMPI_MPI_H
    MPI_Allreduce(MPI_IN_PLACE,Tens.Pre[Comp],CUB(Tens.NSlab),MPI_DOUBLE,MPI_SUM,Proc->CommGrid);  
    for(int t=0;t<2;t++){
      MPI_Allreduce(MPI_IN_PLACE,Tens.Dens[t],CUB(Tens.NSlab),MPI_DOUBLE,MPI_SUM,Proc->CommGrid);
    }
    int Rank=0;
    MPI_Comm_rank(Proc->CommGrid, &Rank);
    if(Rank==0){
#endif
    FILE *FWrite = fopen(FTens,"w");
    fprintf(FWrite,"# l(%.1f %.1f %.1f) r(%.2f %.2f %.2f) v[%d] d[color]\n",pEdge(0),pEdge(1),pEdge(2),Tens.RefPos[0],Tens.RefPos[1],Tens.RefPos[2],Tens.NSlab);
    for(int sx=0;sx<Tens.NSlab;sx++){
      double x = sx*InvValues*pEdge(CLat1);
      for(int sy=0;sy<Tens.NSlab;sy++){
	double y = sy*InvValues*pEdge(CLat2);
	for(int sz=0;sz<Tens.NSlab;sz++){
	  double z = sz*InvValues*pEdge(CNorm);
	  int v = (sx*Tens.NSlab+sy)*Tens.NSlab+sz;
	  if(Tens.Dens[0][v] <= 0 && Tens.Dens[1][v] <= 0 && ABS(Tens.Pre[Comp][v]) <= 0) continue;
	  fprintf(FWrite,"{x(%.3f %.3f %.3f)",x,y,z);
	  fprintf(FWrite," v( %lf %.2f %.2f)}\n",
		  -Tens.Pre[Comp][v]*InvNFile*InvVolume,
		    Tens.Dens[0][v]*InvNFile*InvVolume,
		  Tens.Dens[1][v]*InvNFile*InvVolume);
	}
      }
    }
    fclose(FWrite);
    for(int s=0;s<CUB(Tens.NSlab);s++){
      Tens.Pre[Comp][s] = 0.;
      Tens.Dens[0][s] = 0.;
      Tens.Dens[1][s] = 0.;      
    }
  }
#ifdef OMPI_MPI_H
  }
#endif
}

