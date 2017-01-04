#include "Forces.h"
/**
a
 */
void Forces::SolveLeaves(){
  AddRigid();
  double *Known = (double *) calloc(pNPCh(),sizeof(double));
  double *UnKnown = (double *) calloc(pNPCh(),sizeof(double));
  for(int c=0;c<pNChain();c++){
    int pInit = c*pNPCh();
    int pEnd = (c+1)*pNPCh();
    double Offset = 0.;
    for(int p=pInit;p<pEnd;p++)
      Known[p-pInit] = Pm[p].Pos[CNorm] - Offset;
    IntMatrix->Apply(Known,UnKnown);
    //IntMatrix->Solve(Known,UnKnown);
    for(int p=pInit;p<pEnd;p++){
      int p1 = p-pInit;
      if(Pm[p].Typ != 0){continue;}
      Pm[p].Pos[CNorm] = UnKnown[p1] + Offset;
    }
  }
  free(Known);
  free(UnKnown);
}
/**
a
 */
void Forces::SolveRod(){
  double *Known = (double *) calloc(pNPart(),sizeof(double));
  double *UnKnown = (double *) calloc(pNPart(),sizeof(double));
  double Offset = .5*pEdge(2);
  for(int p=0;p<pNPart();p++)
    Known[p] = Pm[p].Pos[CNorm] - Offset;
  IntMatrix->Apply(Known,UnKnown);
  //IntMatrix->Solve(Known,UnKnown);
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ != 0){continue;}
    Pm[p].Pos[CNorm] = UnKnown[p] + Offset;
  }
  free(Known);
  free(UnKnown);
}
/** 
    It uses Jacobi iterative method to converge to the solution.
 */
void Forces::SolveLinksIterative(){
  AddRigid();
  //SmoothGrid(NEdge);
  double Inter = pEdge(CLat1)/(double)NEdge;
  SPLINE Weight;
  int NIter = 2000;
  Weight.a0 = Kf.El[2];
  Weight.a1 = 0./Inter;
  Weight.a2 = Kf.Lap/SQR(Inter);
  Weight.a3 = 0./(Inter*SQR(Inter));
  Weight.a4 = Kf.SLap/(SQR(Inter)*SQR(Inter));
  int NDim = 2;
  Matrice *CoeffMatrix = new Matrice(Weight,NDim);
  double Offset = 0.;
  double *Known = (double *) calloc(pNPart(),sizeof(double));
  double *UnKnown = (double *) calloc(pNPart(),sizeof(double));
  double Sol = 0.;
  double Zed = 0.;
  int Link[2][5];
  CoeffMatrix->Print();
  for(int i=0;i<NIter;i++){
    for(int p=0;p<pNPart();p++)
      Known[p] = Pm[p].Pos[CNorm] - Offset;
    for(int p=0;p<pNPart();p++){
      if(Pm[p].Typ != 0){
	UnKnown[p] = Known[p];
	continue;
      }
      Sol = 0.;
      for(int d=0;d<2;d++){
	Link[d][1] = Ln[p].Link[d*2];
	Link[d][3] = Ln[p].Link[d*2+1];
	Link[d][0] = Ln[Link[d][1]].Link[d*2];
	Link[d][4] = Ln[Link[d][3]].Link[d*2+1];
	for(int l=0,l1=-2;l<5;l++,l1++){
	  if(l==2){
	    l1 = 0;
	    continue;
	  }
	  // if(!PeriodicImage[d]){
	  //   if(d == 0){
	  //     if(Link[d][l] != p+l1*nEdge[1]) continue;
	  //     if(Link[d][l] != p+l1) continue;
	  //   }
	  // }
	  Sol += CoeffMatrix->Val(2,l)*Pm[Link[d][l]].Pos[2];
	}
      }
      UnKnown[p] = .4*Known[p] - .6*Sol/CoeffMatrix->Val(2,2);
    }
    for(int p=0;p<pNPart();p++){
      if(Pm[p].Typ != 0){continue;}
      Pm[p].Pos[CNorm] = UnKnown[p] + Offset;
    }
  }
  free(Known);
  free(UnKnown);
  delete CoeffMatrix;
}
/**
a
 */
void Forces::SolveLinks(){		
  AddRigid();
  FillMatrix();
  double Offset = 0.;
  double *Known = (double *) calloc(pNPart(),sizeof(double));
  double *UnKnown = (double *) calloc(pNPart(),sizeof(double));
  for(int p=0;p<pNPart();p++)
    Known[p] = Pm[p].Pos[CNorm] - Offset;
  IntMatrix->Apply(Known,UnKnown);
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ != 0){continue;}
    Pm[p].Pos[CNorm] = UnKnown[p] + Offset;
  }
  free(Known);
  free(UnKnown);
}
// #include "compcol_double.h"
// #include "comprow_double.h"
// #include "ilupre_double.h"
// #include "icpre_double.h"
// #include "iotext_double.h"
// void Forces::SolveLinksSparse(){
//   int NEntries = 0;
//   int NColMax = pNPart();
//   int NRowMax = pNPart();
//   for(int p=0;p<pNPart();p++){
//     if(Pm[p].Typ != 0){
//       IntMatrix->Set(p,p,1.);Entries++;
//       continue;
//     }
//     int pym1 = Ln[p].Link[0];
//     int pyp1 = Ln[p].Link[1];
//     int pym2 = Ln[pym1].Link[0];
//     int pyp2 = Ln[pyp1].Link[1];
//     int pxm1 = Ln[p].Link[2];
//     int pxp1 = Ln[p].Link[3];
//     int pxm2 = Ln[pxm1].Link[2];
//     int pxp2 = Ln[pxp1].Link[3];
//     if(PeriodicImage[0]){ Entries += 4;
//     }
//     else{
//       if(pxm2 == p-2*nEdge[1]) Entries++;
//       if(pxm1 == p-nEdge[1]) Entries++;
//       if(pxp1 == p+nEdge[1]) Entries++;
//       if(pxp2 == p+2*nEdge[1]) Entries++;
//     }
//     Entries++;
//     if(PeriodicImage[1]){Entries += 4;
//     }
//     else{
//       if(pym2 == p-2) Entries++;
//       if(pym1 == p-1) Entries++;
//       if(pyp1 == p+1) Entries++;
//       if(pyp2 == p+2) Entries++;
//     }
//   }
//   double *MatVal = (double *)calloc(Entries,sizeof(double));
//   int *ColIdx = (int *)calloc(Entries,sizeof(int));
//   int *RowIdx = (int *)calloc(Entries,sizeof(int));
//   for(int p=0,i=0;p<pNPart();p++){
//     if(Pm[p].Typ != 0){
//       ColIdx[i] = p;
//       RowIdx[i] = p;
//       MatVal[i] = 1.;
//       i++;
//       continue;
//     }
//     int pym1 = Ln[p].Link[0];
//     int pyp1 = Ln[p].Link[1];
//     int pym2 = Ln[pym1].Link[0];
//     int pyp2 = Ln[pyp1].Link[1];
//     int pxm1 = Ln[p].Link[2];
//     int pxp1 = Ln[p].Link[3];
//     int pxm2 = Ln[pxm1].Link[2];
//     int pxp2 = Ln[pxp1].Link[3];
//     // printf("%d)\n",p);
//     //printf("%d %d %d %d\n",pym2,pym1,pyp1,pyp2);
//     // printf("%d %d %d %d\n",pxm2,pxm1,pxp1,pxp2);
//     if(PeriodicImage[0]){ 
//       RowIdx[i]=p;ColIdx[i]=pxm2;ValMat[i]=CoeffMatrix->Val(2,0);i++;
//       RowIdx[i]=p;ColIdx[i]=pxm1;ValMat[i]=CoeffMatrix->Val(2,1);i++;
//       RowIdx[i]=p;ColIdx[i]=pxp1;ValMat[i]=CoeffMatrix->Val(2,3);i++;
//       RowIdx[i]=p;ColIdx[i]=pxp2;ValMat[i]=CoeffMatrix->Val(2,4);i++;
//     }
//     else{
//       if(pxm2 == p-2*nEdge[1])
// 	RowIdx[i]=p;ColIdx[i]=pxm2;ValMat[i]=CoeffMatrix->Val(2,0);i++;
//       if(pxm1 == p-nEdge[1]) 
// 	RowIdx[i]=p;ColIdx[i]=pxm1;ValMat[i]=CoeffMatrix->Val(2,1);i++;
//       if(pxp1 == p+nEdge[1]) 
// 	RowIdx[i]=p;ColIdx[i]=pxp1;ValMat[i]=CoeffMatrix->Val(2,3);i++;
//       if(pxp2 == p+2*nEdge[1])
// 	RowIdx[i]=p;ColIdx[i]=pxp2;ValMat[i]=CoeffMatrix->Val(2,4);i++;
//     }
//     IntMatrix->Add(p,p,CoeffMatrix->Val(2,2));Entries++;
//     if(PeriodicImage[1]){
//       RowIdx[i]=p;ColIdx[i]=pym2;ValMat[i]=CoeffMatrix->Val(2,0);i++;
//       RowIdx[i]=p;ColIdx[i]=pym1;ValMat[i]=CoeffMatrix->Val(2,1);i++;
//       RowIdx[i]=p;ColIdx[i]=pyp1;ValMat[i]=CoeffMatrix->Val(2,3);i++;
//       RowIdx[i]=p;ColIdx[i]=pyp2;ValMat[i]=CoeffMatrix->Val(2,4);i++;
//     }
//     else{
//       if(pym2 == p-2) 
// 	RowIdx[i]=p;ColIdx[i]=pym2;ValMat[i]=CoeffMatrix->Val(2,0);i++;
//       if(pym1 == p-1) 
// 	RowIdx[i]=p;ColIdx[i]=pym1;ValMat[i]=CoeffMatrix->Val(2,1);i++;
//       if(pyp1 == p+1) 
// 	RowIdx[i]=p;ColIdx[i]=pyp1;ValMat[i]=CoeffMatrix->Val(2,3);i++;
//       if(pyp2 == p+2) 
// 	RowIdx[i]=p;ColIdx[i]=pyp2;ValMat[i]=CoeffMatrix->Val(2,4);i++;
//     }
//   }
//   Coord_Mat_double C(NRowMax,NColMax,Entries,MatVal,RowIdx,ColIdx);
//   free(ColIdx);
//   free(RowIdx);
//   free(MatVal);
// }
// obsolete?
// int Forces::Update(){
//   for(int i=0;i<IntMax;i++){
//     double Diff = Fm[Part2Move].Dir[2] - Fm[Part2Move].Ext[2];
//     if(QUAD(Diff) < DiffForce ){
//       printf("%d,%d) Force %lf= %lf - %lf Pos %lf \n",Part2Move,i,Diff,Fm[Part2Move].Dir[2],Fm[Part2Move].Ext[2],Pm[Part2Move].Pos[2]);
//       break;
//     }
//     double Sign=1.;
//     if( Diff < 0. ){
//       Sign = -1.;
//       Pm[Part2Move].Pos[2] -= IncrDist;      
//     }
//     else 
//       Pm[Part2Move].Pos[2] += IncrDist;
//     for(int p=0,pp=Part2Move-1,ppp=Part2Move+1;p<pNPart();p++)
//       ;
//     if( Diff < 0. ){
//       Pm[Part2Move].Pos[2] -= IncrDist;
//     }
//     else if(Diff > 0.){
//       Pm[Part2Move].Pos[2] += IncrDist;
//     }
//     ForceFieldLine();
//     Diff = Fm[Part2Move].Dir[2] - Fm[Part2Move].Ext[2];
//     if(QUAD(Diff) < DiffForce ){
//       printf("%d,%d) Force %lf= %lf - %lf Pos %lf \n",Part2Move,i,Diff,Fm[Part2Move].Dir[2],Fm[Part2Move].Ext[2],Pm[Part2Move].Pos[2]);
//       break;
//     }
//     if(Part2Move>0){
//       if( Diff < 0. ){
// 	Pm[Part2Move-1].Pos[2] += IncrDist;
//       }
//       else if(Diff > 0.){
// 	Pm[Part2Move-1].Pos[2] -= IncrDist;
//       }
//     }
//     if(Part2Move<pNPart()-1){
//       if( Diff < 0. ){
// 	Pm[Part2Move+1].Pos[2] += IncrDist;
//       }
//       else if(Diff > 0.){
// 	Pm[Part2Move+1].Pos[2] -= IncrDist;
//       }
//     }
//     ForceFieldLine();
//   }
// }
// // obsolete?
// int Forces::MinHelfrich(){
//   SPLINE Par;
//   SPLINE Sp2;
//   double a1=0.;
//   double a2=0.;
//   double a3=0.;
//   double a4=0.;
//   double Ref=0.;
//   double Dx = pEdge(0)/(double)(NEdge-1);
//   for(int p=0;p<pNPart()-1;p++){
//     if(p < 1){
//       Sp2 = Mat->Forth(Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,Pm[p+4].Pos,0,2);
//       Ref = Pm[p].Pos[0] - Pm[p+2].Pos[0]; 
//     }
//     else if(p < 2){
//       Sp2 = Mat->Forth(Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,Pm[p+3].Pos,0,2);
//       Ref = Pm[p].Pos[0] - Pm[p+1].Pos[0]; 
//     }
//     else if(p < pNPart()-2){
//       Sp2 = Mat->Forth(Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,Pm[p+2].Pos,0,2);
//       Ref = 0.; 
//     }
//     else if(p == pNPart()-2){
//       Sp2 = Mat->Forth(Pm[p-3].Pos,Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,Pm[p+1].Pos,0,2);
//       Ref = Pm[p].Pos[0] - Pm[p-1].Pos[0]; 
//     }
//     else if(p <= pNPart()-1){
//       Sp2 = Mat->Forth(Pm[p-4].Pos,Pm[p-3].Pos,Pm[p-2].Pos,Pm[p-1].Pos,Pm[p].Pos,0,2);
//       Ref = Pm[p].Pos[0] - Pm[p-2].Pos[0]; 
//     }
//     if(Pm[p].Typ == 1 || Pm[p].Typ == 2)
//       continue;
//     a2 = 2.*Sp2.a2 + 6.*Sp2.a3*Ref + 12.*Sp2.a4*Ref*Ref;
//     a4 = 24.*Sp2.a4;
//     a3 = a4*(12.*Dx*Dx - 24.*Kf.Lap/Kf.SLap) - 2.*a2;
//     //a3 = a4*(12.*Dx*Dx - 24.*.1) - 2.*a2;
//     a3 /= Dx;
//     a1 = Pm[p+1].Pos[2] - Pm[p].Pos[2] - a2*Dx*Dx - a3*Dx*Dx*Dx - a4*Dx*Dx*Dx*Dx;
//     a1 /= Dx;
//     double x = Pm[p+1].Pos[0] - Pm[p].Pos[0];
//     printf("%d) %lf %lf %lf %lf -> %lf\n",p,a1,a2,a3,a4,Pm[p].Pos[2] + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x);
//     Pm[p+1].Pos[2] = Pm[p].Pos[2] + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x;
//   }
// }
//---------------------------------RIGID--------------------------
/** Update the position of a rotating cylinder */
void Forces::VelVerletRigid(){
  for(int n=0;n<pNNano();n++){
    Vettore Axis(Nano[n].Axis,3);
    Vettore AMom(Nano[n].AMom,3);
    Vettore BackBone(3);
    Vettore VelTang(3);
    for(int d=0;d<3;d++){
      BackBone.x[d] = Nano[n].Axis[d]*Nano->Height*.5;
      double Ran = Nano[n].Zeta  * (2.*Mat->Casuale() - 1.);
      double Dis = - Nano[n].Gamma*Nano[n].AVel[d];
      Nano[n].AMom[d] += Ran + Dis;
    }
    Matrice InTensor(3,3);
    Matrice RotT(3,3);
    Matrice Resp(3,3);
    Matrice Resp2(3,3);
    Vettore AVel(3);
    Vettore Normal(0.,0.,1.);
    Vettore Ax(3);
    Axis.VetV(&Axis,&Normal);
    double Angle = Ax.Angle(&Axis,&Normal);
    Quadri q(Ax.x,Angle);
    Matrice Rot(q,3);
    Rot.CopyOn(&RotT);
    RotT.Transpose();
    InTensor.Set(0,0,1./(Nano->Mass*(.25*SQR(Nano->Rad)+1./12.*SQR(Nano->Height) ) ) );
    InTensor.Set(1,1,1./(Nano->Mass*(.25*SQR(Nano->Rad)+1./12.*SQR(Nano->Height) ) ) );
    InTensor.Set(2,2,1./(.5*Nano->Mass*SQR(Nano->Rad)) );
    // Calculating r^t I^-1 r L = w
    Resp.Mult(RotT,InTensor);
    Resp2.Mult(Resp,Rot);
    //    MatrVect(Resp2,Nano->AMom,AVel);
 
    //printf("%d) %d\n",n,Gen->Step);
    for(int d=0;d<3;d++){
      double tmp = .5*Nano[n].Force[d]*pDeltat()/Nano[n].Mass;
      Nano[n].Vel[d] += tmp;
      Nano[n].AVel[d] += Nano[n].AMom[d]/Nano[n].Mass*.01;
      Nano[n].Pos[d] += Nano[n].Vel[d]*pDeltat();
      Nano[n].Pos[d] -= floor(Nano[n].Pos[d]/pEdge(d))*pEdge(d);
    }
    if(Nano[n].Shape != 2) continue;
    Vettore Omega(Nano[n].AVel,3);
    for(int d=0;d<3;d++)
      BackBone.x[d] = Nano[n].Height*.5*Nano[n].Axis[d];
    VelTang.VetV(&Omega,&BackBone);
    for(int d=0;d<3;d++){
      BackBone.x[d] += VelTang.x[d]*pDeltat();
      Axis.x[d] = BackBone.x[d];
    }
    Axis.Normalize();
    Axis.Export(Nano[n].Axis);
    //printf("%lf %lf %lf\n",Nano[n].Vel[0],Nano[n].Vel[1],Nano[n].Vel[2]);
    //printf("%lf %lf %lf\n",Nano[n].Force[0],Nano[n].Force[1],Nano[n].Force[2]);
    //printf("%lf %lf %lf\n",Nano[n].Pos[0],Nano[n].Pos[1],Nano[n].Pos[2]);
    //printf("%lf %lf %lf\n",Nano[n].AMom[0],Nano[n].AMom[1],Nano[n].AMom[2]);
  }
}
/**
a
 */
void Forces::VelVerletRigid2(){
  for(int n=0;n<pNNano();n++) {
    for(int d=0;d<3;d++){
      Nano[n].AVel[d] += .5*Nano[n].AMom[d];
      Nano[n].Vel[d] += .5*Nano[n].Force[d]*pDeltat()/Nano[n].Mass;
    }
  }
}
//------------------------------------------------------------
//-------------------------MONTE-CARLO------------------------
//------------------------------------------------------------
/** The general version of the Metropolis algorythm */
int Forces::IfMetropolis(double Arg,double Weight){
  if(exp(pBeta()*Arg)*Weight > 1. ) return 1;
  double Ran = Mat->Casuale();
  if(exp(pBeta()*Arg)*Weight > Ran ){
    return 1;
  }
  return 0;
}
//-----------------Single-particle----------------------------
/** Insert a particle in the system and in the cell list */
int Forces::InsertBead(int p){
  Pm[p].Pos[0] = pEdge(0)*Mat->Casuale();
  Pm[p].Pos[1] = pEdge(1)*Mat->Casuale();
  Pm[p].Pos[2] = pEdge(2)*Mat->Casuale();
  Pc->AddPart(p,Pm[p].Pos);
}
/**
a
 */
int Forces::TryInsert(){
  int p = pNPart();
  ReSetNPart(pNPart()+1);
  InsertBead(p);
  double Pot[3];
  double NrgDiff = CalcNrgBead(p,Pot);
  double Arg = ChemPotId + ChemPotEx - NrgDiff;
  double Weight = pVol()/(double)(pNPart()+1);
  if( !IfMetropolis(Arg,Weight) ){
    Pc->RemPart(p,Pm[p].Pos);
    ReSetNPart(pNPart()-1);
    return 0;
  }
  //printf("Added a particle %d\n",p);
  return 1;
}
/**
a
 */
int Forces::TryRemove(){
  printf("remove\n");
  int p = (int)(Mat->Casuale()*pNPart());
  double Pot[3];
  double NrgDiff = CalcNrgBead(p,Pot); 
  double Arg = - ChemPotId - ChemPotEx - NrgDiff;
  if( !IfMetropolis(Arg,pNPart()/pVol()) ){
    return 0;
  }
  SwapPart(p,pNPart()-1);
  Pc->SwapPart(p,Pm[p].Pos,pNPart()-1,Pm[pNPart()-1].Pos);
  Pc->RemPart(pNPart()-1,Pm[pNPart()-1].Pos);
  ReSetNPart(pNPart()-1);
  return 1;
}
/** Give a new position for the particle p */
int Forces::MoveBead(int p){
  for(int d=0;d<3;d++){
    Pm[p].Pos[d] += Mat->Gaussiano(0.,GaussVar);
    Pm[p].Pos[d] -= floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);
  }
}
/**
a
 */
int Forces::TryMove(){
  int p = (int)(Mat->Casuale()*pNPart());
  if(Pm[p].Typ >= 2) return 1;
  //old situa
  double Pos[3] = {Pm[p].Pos[0],Pm[p].Pos[1],Pm[p].Pos[2]};
  double Pot[3] = {0.,0.,0.};
  double Nrg0 = CalcNrgBead(p,Pot);
  //RemDens(p,p+1);
  //double Nrg0 = OldNrgPm[p*3] + OldNrgPm[p*3+1] + OldNrgPm[p*3+2];
  //new situa
  MoveBead(p);
  Pc->MovePart(p,Pos,Pm[p].Pos);
  //AddDens(p,p+1);
  double Nrg1 = CalcNrgBead(p,Pot);
  //like the new situa?
  double NrgDiff = Nrg1 - Nrg0;
  if(!IfMetropolis(-NrgDiff,1.)){
    //RemDens(p,p+1);
    Pc->MovePart(p,Pm[p].Pos,Pos);
    for(int d=0;d<3;d++) Pm[p].Pos[d] = Pos[d];
    //AddDens(p,p+1);
    return 0;
  }
  //OldNrgSys = SumDens(0,pNPart());
  OldNrgSys -= Pot[2];//DensFuncNrgCh(c);
  return 1;
}
//-----------------------------Chains---------------------
/**
a
 */
void Forces::StudySys(){
  //getting the distribution
  double Norm = 0.;
  double *LineS = new double[NBin];
  for(int i=0;i<NBin;i++){
    FirstBeadDistr[i] = 0.;
  }
  for(int c=0;c<pNChain();c++){
    int p1 = c*pNPCh();
    int i = (int)(Pm[p1].Pos[CNorm]*NBin*pInvEdge(CNorm));
    if(i < 0 || i >= NBin) continue;
    FirstBeadDistr[i] += 1.;
    Norm += 1.;
  }
  //normalizing and smoothing
  for(int v=0;v<NBin;v++) LineS[v] = FirstBeadDistr[v];
  InterBSpline1D(LineS,FirstBeadDistr,NBin,NBin);
  for(int v=0;v<NBin;v++) LineS[v] = FirstBeadDistr[v];
  InterBSpline1D(LineS,FirstBeadDistr,NBin,NBin);
  Norm = Norm >= 0. ? 1./Norm : 1.;
  for(int i=0;i<NBin;i++){
    FirstBeadDistr[i] *= Norm;
  }
  //cumulative
  for(int i=1;i<NBin;i++){
    FirstBeadDistr[i] += FirstBeadDistr[i-1];
    if(FirstBeadDistr[i] > .25 && FirstBeadDistr[i-1] < .25)
      BorderBias[0] = i/NBin*pEdge(CNorm);
    else if(FirstBeadDistr[i] > .75 && FirstBeadDistr[i-1] < .75)
      BorderBias[1] = i/NBin*pEdge(CNorm);
  }
  delete [] LineS;
}
/**
a
 */
void Forces::ConsiderCh(int c){
  int p1 = c*pNPCh();
  for(int p=0;p<pNPCh();p++) Pc->AddPart(p+p1,Pm[p+p1].Pos);
  AddDens(p1,p1+pNPCh());
}
/**
a
 */
void Forces::IgnoreCh(int c){
  int p1 = c*pNPCh();
  RemDens(p1,p1+pNPCh());
  for(int p=0;p<pNPCh();p++) Pc->RemPart(p+p1,Pm[p+p1].Pos);
}
/**
a
 */
void Forces::RemChFromSys(int c){
  int c2 = pNChain()-1;
  int p1 = c*pNPCh();
  int p2 = c2*pNPCh();
  RemDens(p1,p1+pNPCh());
  if(c != c2){
    for(int e=0;e<3;e++)
      OldNrgCh[c*3+e] = OldNrgCh[c2*3+e];
    RemDens(p2,p2+pNPCh());
    // for(int p=0;p<pNPCh();p++)
    //   Pc->RemPart(p2+p,Pm[p2+p].Pos);
    // for(int p=0;p<pNPCh();p++)
    //   Pc->RemPart(p1+p,Pm[p1+p].Pos);
    // SwapChain(c,pNChain()-1);
    // for(int p=0;p<pNPCh();p++)
    //   Pc->AddPart(p1+p,Pm[p1+p].Pos);
    for(int p=0;p<pNPCh();p++)
      Pc->SwapPart(p+p1,Pm[p+p1].Pos,p2+p,Pm[p2+p].Pos);
    SwapChain(c,pNChain()-1);
    AddDens(p1,p1+pNPCh());
  }
  for(int p=0;p<pNPCh();p++)
    Pc->RemPart(p2+p,Pm[p2+p].Pos);
  ReSetNPart(pNPart()-pNPCh());
  ReSetNChain(pNChain()-1);
}
/**
a
 */
void Forces::SaveCh(int c){
  int p1 = c*pNPCh();
  for(int p=0;p<pNPCh();p++){
    for(int d=0;d<3;d++){
      OldPos[p][d] = Pm[p+p1].Pos[d];
    }
  }
}
/**
a
 */
void Forces::ReInsertCh(int c){
  IgnoreCh(c);
  int p1 = c*pNPCh();
  //old situa
  for(int p=0;p<pNPCh();p++){
    for(int d=0;d<3;d++){
      Pm[p+p1].Pos[d] = OldPos[p][d];
    }
  }
  ConsiderCh(c);
}
/** Insert a chain in the system and in the cell list */
double Forces::InsertCh(int c){
  int p1 = c*pNPCh();
  double Weight = 1.;
  ReSetNPart(pNPart()+pNPCh());
  ReSetNChain(pNChain()+1);
  Pm[p1].Pos[CLat1] = Mat->Casuale()*pEdge(CLat1);
  Pm[p1].Pos[CLat2] = Mat->Casuale()*pEdge(CLat2);
  Pm[p1].Pos[CNorm] = Mat->Casuale()*pEdge(CNorm);
  if(VAR_IF_TYPE(CalcMode,CALC_BIL_BIAS)){
    Pm[p1].Pos[CNorm] = Mat->Casuale()*(BorderBias[1] - BorderBias[0]) + BorderBias[0];      
      //Mat->RandDiscrProb(FirstBeadDistr,NBin)*pEdge(CNorm);
    Weight *= (BorderBias[1]-BorderBias[0])*pEdge(CLat1)*pEdge(CLat2)/pVol();
  }
  else if(VAR_IF_TYPE(CalcMode,CALC_SPH_BIAS)){
    double Incr = 0.2;
    double z = 2.0 * Mat->Casuale() - 1.0;
    double t = 2.0 * M_PI * Mat->Casuale();
    double w = sqrt( 1 - z*z );
    double x = w * cos( t );
    double y = w * sin( t );
    double Rad = Mat->Casuale()*Incr + Nano->Rad;
    Pm[p1].Pos[CLat1] = x*Rad + Nano->Pos[0];
    Pm[p1].Pos[CLat2] = y*Rad + Nano->Pos[1]; 
    Pm[p1].Pos[CNorm] = z*Rad + Nano->Pos[2];
    double Vol = 4./3.*M_PI*( CUBE(Nano->Rad + Incr) - CUBE(Nano->Rad));
    Weight *= Vol/pVol();
  }
  Weight *= InsertRest(p1,0);
  // for(int p=0;p<pNPCh();p++){
  //   fprintf(StatFile1,"%lf %lf %lf 0. 0. 0. %d\n",Pm[p1+p].Pos[0],Pm[p1+p].Pos[1],Pm[p1+p].Pos[2],pType(p+p1));
  // }
  // fflush(StatFile1);
  ConsiderCh(c);
  return Weight;
}
/** Insert the remaining particles of the new created chain */
double Forces::InsertRest(int pCurr,int StartPos){
  double Weight = 1.;
  for(int p=StartPos+1;p<pNPCh();p++){
    for(int d=0;d<3;d++){
      double Ran = Mat->Gaussiano(0.,GaussVar);
      Pm[p+pCurr].Pos[d] = Pm[p-1+pCurr].Pos[d] + Ran;
      Pm[p+pCurr].Pos[d] -= floor(Pm[p+pCurr].Pos[d]*pInvEdge(d))*pEdge(d);
    }
  }
  for(int p=StartPos-1;p>=0;p--){
    for(int d=0;d<3;d++){
      double Ran = Mat->Gaussiano(0.,GaussVar);
      Pm[p+pCurr].Pos[d] = Pm[p+1+pCurr].Pos[d] + Ran;
      Pm[p+pCurr].Pos[d] -= floor(Pm[p+pCurr].Pos[d]*pInvEdge(d))*pEdge(d);
    }
  }
  return Weight;
}
/**
a
 */
int Forces::TryMoveCh(){
  int c = (int)(Mat->Casuale()*pNChain());
  int p1 = c*pNPCh();
  double Pot[3];
  // old situa
  SaveCh(c);
  IgnoreCh(c);
  double NrgOld = CalcNrgCh(c,Pot);
  // new situa
  ReSetNPart(pNPart()-pNPCh());
  ReSetNChain(pNChain()-1);
  InsertCh(c);
  double NrgNew = CalcNrgCh(c,Pot);
  //like the new situa?
  double NrgDiff = NrgNew - NrgOld;
  if(!IfMetropolis(-NrgDiff,1.)){
    ReInsertCh(c);
    return 0;
  }
  OldNrgCh[c*3+2] = NrgNew;
  OldNrgSys += NrgDiff;
  return 1;
}
/**
a
 */
int Forces::TryRemoveCh(){
  if(pNPart() <= 0) return 0;
  int c = (int)(Mat->Casuale()*pNChain());
  int p1 = c*pNPCh();
  double Pot[3];
  //like the new situa?
  double NrgDiff = CalcNrgCh(c,Pot);
  //double NrgDiff = OldNrgCh[c*3+2];
  //fprintf(StatFile1,"%lf\n",NrgDiff);
  //fflush(StatFile1);
   double Arg = -ChemPotId -ChemPotEx +NrgDiff;
  double Weight = pNChain()/pVol();
  if(!IfMetropolis(Arg,Weight) ){
    return 0;
  }
  RemChFromSys(c);
  OldNrgSys -= NrgDiff;
  return 1;
}
/**
a
 */
int Forces::TryInsertCh(){
  //old situa
  int c = pNChain();
  int p1 = c*pNPCh();
  double Pot[3];
  //new situa
  InsertCh(c);
  //like the new situa?
  double NrgDiff = CalcNrgCh(c,Pot);
  //fprintf(StatFile2,"%lf\n",NrgDiff);
  fflush(StatFile2);
  //printf("add %lf sys %lf\n",NrgDiff,OldNrgSys);
  double Arg = ChemPotId + ChemPotEx - NrgDiff;
  double Weight = pVol()/(double)(pNChain()+1);
  //printf("Ins %lf\n",NrgDiff);
  if(!IfMetropolis(Arg,Weight)){
    RemChFromSys(c);
    return 0;
  }
  //printf("accepted\n");
  //printf("Inserted %lf\n",NrgDiff);
  OldNrgCh[c*3+2] = NrgDiff;
  OldNrgSys += NrgDiff;
  return 1;
}
//---------------------------Bias-------------------------------
/** Remove a chain in the system calculating the weight of the super detailed balance */
double Forces::RemoveChBias(int c){
  int p1 = c*pNPCh();
  RemDens(p1,p1+pNPCh());
  for(int p=0;p<pNPCh();p++) Pc->RemPart(p+p1,Pm[p+p1].Pos);
  double NrgDiff = DensFuncNrgGhost(Pm[p1].Pos,p1,Pm[p1].Typ);
  //fprintf(StatFile1,"%lf %lf\n",NrgDiff,NanoNrg(p1));
  NrgDiff += CalcBendingGhost(Pm[p1].Pos,p1);
  NrgDiff += NanoNrg(p1);
  double Weight = exp(-NrgDiff+NrgPBead);
  Pc->AddPart(p1,Pm[p1].Pos);
  AddDens(p1,p1+1);
  for(int p=p1+1;p<p1+pNPCh();p++){
    Weight *= WeightSetBond(p,Pm[p].Typ);
    Pc->AddPart(p,Pm[p].Pos);
    AddDens(p,p+1);
  }
  return Weight;
}
/** Calculate the weights of a set of bonds created for the existing monomer of the chain to remove */
double Forces::WeightSetBond(int p,int Type){
  double NrgDiff = DensFuncNrgGhost(Pm[p].Pos,p,Pm[p].Typ);
  NrgDiff += CalcBendingGhost(Pm[p].Pos,p);
  NrgDiff += NanoNrg(p);
  CumProbBias[0] = exp(-NrgDiff+NrgPBead);
  for(int t=1;t<NTrialBias;t++){
    for(int d=0;d<3;d++){
      BondPosBias[t][d] = Pm[p-1].Pos[d] + Mat->Gaussiano(0.,GaussVar);
      BondPosBias[t][d] -= floor(BondPosBias[t][d]*pInvEdge(d))*pEdge(d);
    }
    NrgDiff = DensFuncNrgGhost(BondPosBias[t],p,Type);
    //fprintf(StatFile1,"%lf %lf\n",NrgDiff,NanoNrg(BondPosBias[t],Type));
    NrgDiff += CalcBendingGhost(BondPosBias[t],p);
    NrgDiff += NanoNrg(BondPosBias[t],Type);
    CumProbBias[t] = exp(-NrgDiff+NrgPBead);
  }
  double Weight = 0.;
  for(int t=0;t<NTrialBias;t++){
    Weight += CumProbBias[t];
  }
  return Weight/(double)NTrialBias;
}
/** Create a chain choosing among a set of different bond vectors with a probability given by their Boltzmann factors */
double Forces::InsertChBias(int c){
  int p1 = c*pNPCh();
  double Weight = 1.;
  ReSetNPart(pNPart()+pNPCh());
  ReSetNChain(pNChain()+1);
  Pm[p1].Pos[CLat1] = Mat->Casuale()*pEdge(CLat1);
  Pm[p1].Pos[CLat2] = Mat->Casuale()*pEdge(CLat2);
  Pm[p1].Pos[CNorm] = Mat->Casuale()*pEdge(CNorm);
  if(VAR_IF_TYPE(CalcMode,CALC_BIL_BIAS)){
    Pm[p1].Pos[CNorm] = Mat->Casuale()*(BorderBias[1] - BorderBias[0]) + BorderBias[0];      
      //Mat->RandDiscrProb(FirstBeadDistr,NBin)*pEdge(CNorm);
    Weight *= (BorderBias[1]-BorderBias[0])*pEdge(CLat1)*pEdge(CLat2)/pVol();
  }
  else if(VAR_IF_TYPE(CalcMode,CALC_SPH_BIAS)){
    double Incr = 0.2;
    double z = 2.0 * Mat->Casuale() - 1.0;
    double t = 2.0 * M_PI * Mat->Casuale();
    double w = sqrt( 1 - z*z );
    double x = w * cos( t );
    double y = w * sin( t );
    double Rad = Mat->Casuale()*Incr + Nano->Rad;
    Pm[p1].Pos[CLat1] = x*Rad + Nano->Pos[0];
    Pm[p1].Pos[CLat2] = y*Rad + Nano->Pos[1];
    Pm[p1].Pos[CNorm] = z*Rad + Nano->Pos[2];
    double Vol = 4./3.*M_PI*( CUBE(Nano->Rad + Incr) - CUBE(Nano->Rad));
    Weight *= Vol/pVol();
    //fprintf(StatFile2,"%lf %lf %lf %d\n",Pm[p1].Pos[0],Pm[p1].Pos[1],Pm[p1].Pos[2],pType(p1));
    //fflush(StatFile2);
  }
  double NrgDiff = DensFuncNrgGhost(Pm[p1].Pos,p1,Pm[p1].Typ);
  fprintf(StatFile2,"%lf %lf\n",NrgDiff,NanoNrg(p1));
  NrgDiff += CalcBendingGhost(Pm[p1].Pos,p1);
  NrgDiff += NanoNrg(p1);
  Weight *= exp(-NrgDiff+NrgPBead);
  // printf("Weight1 %lf %lf\n",NrgDiff,Weight);
  Pc->AddPart(p1,Pm[p1].Pos);
  AddDens(p1,p1+1);
  for(int p=p1+1;p<(c+1)*pNPCh();p++){
    double Weight1 = CreateSetBond(p,Pm[p].Typ);
    Weight *= Weight1;
    Pc->AddPart(p,Pm[p].Pos);
    AddDens(p,p+1);
  }
  // for(int p=0;p<pNPCh();p++){
  //   fprintf(StatFile1,"%lf %lf %lf %d\n",Pm[p1+p].Pos[0],Pm[p1+p].Pos[1],Pm[p1+p].Pos[2],pType(p+p1));
  // }
  // fflush(StatFile1);
  return Weight;
}
/** Create a set of bond vectors and choose among them one with a probability given by its Boltzmann weight */
double Forces::CreateSetBond(int p,int Type){
  int tChoosen = 0;
  for(int t=0;t<NTrialBias;t++){
    for(int d=0;d<3;d++){
      BondPosBias[t][d] = Pm[p-1].Pos[d] + Mat->Gaussiano(0.,GaussVar);
      BondPosBias[t][d] -= floor(BondPosBias[t][d]*pInvEdge(d))*pEdge(d);
    }
    double NrgDiff = DensFuncNrgGhost(BondPosBias[t],p,Type);
    fprintf(StatFile2,"%lf %lf\n",NrgDiff,NanoNrg(BondPosBias[t],Type));
    NrgDiff += CalcBendingGhost(BondPosBias[t],p);
    NrgDiff += NanoNrg(BondPosBias[t],Type);
    CumProbBias[t] = exp(-NrgDiff+NrgPBead);
  }
  double Weight = 0.;
  for(int t=0;t<NTrialBias;t++){
    Weight += CumProbBias[t];
  }
  CumProbBias[0] /= Weight;
  double Ran = Mat->Casuale();
  for(int t=1;t<NTrialBias;t++){
    CumProbBias[t] /= Weight;
    CumProbBias[t] += CumProbBias[t-1];
  }
  for(int t=0;t<NTrialBias;t++){
    if(Ran < CumProbBias[t]){
      tChoosen = t;
      break;
    }
  }
  for(int d=0;d<3;d++){
    Pm[p].Pos[d] = BondPosBias[tChoosen][d];
    Pm[p].Pos[d] -= floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);
  }
  return Weight/(double)NTrialBias;
}
/**
a
 */
int Forces::TryRemoveChBias(){
  if(pNPart() <= 0) return 0;
  int c = (int)(Mat->Casuale()*pNChain());
  int p1 = c*pNPCh();
  double Weight1 = RemoveChBias(c);
  double Arg = - ChemPotId - ChemPotEx ;
  double Weight = pNChain()/(Weight1*pVol());
  //printf("rem %lf\n",Weight1);
  //fprintf(TempFile,"1 %lf\n",NrgDiff);
  if(!IfMetropolis(Arg,Weight)){
    return 0;
  }
  //fprintf(TempFile,"2 %lf\n",NrgDiff);
  RemChFromSys(c);
  return 1;
}
/**
a
 */
int Forces::TryInsertChBias(){
  //add one chain
  int c = pNChain();
  int p1 = c*pNPCh();
  double Weight1 = InsertChBias(c);
  double Arg = ChemPotId + ChemPotEx;
  double Weight = Weight1*pVol()/(double)(pNChain()+1);
  if(!IfMetropolis(Arg,Weight)){
    RemChFromSys(c);
    return 0;
  }
  return 1;
}
//-------------------WIDOM--------------------------
/**
a
 */
void Forces::WidomInsert(double *NrgDiff){
  int p1 = pNPart();
  InsertBead(p1);
  *NrgDiff = DensFuncNrgGhost(Pm[p1].Pos,p1,Pm[p1].Typ);
  *NrgDiff += CalcBendingGhost(Pm[p1].Pos,p1);
  *NrgDiff += NanoNrg(p1);
  Pc->RemPart(p1,Pm[p1].Pos);
}
/**
a
 */
void Forces::WidomRemove(double *NrgDiff,int p){
  double Pot[3];
  *NrgDiff = CalcNrgBead(p,Pot);
}
/**
a
 */
void Forces::WidomInsertCh(double *NrgDiff){
  int c = pNChain();
  int p1 = c*pNPCh();
  InsertCh(c);
  CalcNrgCh(c,NrgDiff);
  RemChFromSys(c);
}
/**
a
 */
void Forces::WidomRemoveCh(double *NrgDiff,int c){
  NrgDiff[0] = OldNrgCh[c*3  ];
  NrgDiff[1] = OldNrgCh[c*3+1];
  NrgDiff[2] = OldNrgCh[c*3+2];
}
/**
a
 */
void Forces::WidomBiasChIn(double *Weight){
  int c = pNChain();
  int p1 = c*pNPCh();
  *Weight = InsertChBias(c);
  RemChFromSys(c);
}
/**
a
 */
void Forces::WidomBiasChOut(double *Weight,int c){
  *Weight = RemoveChBias(c);
}
//----------------------MOLECULAR-DYNAMICS----------------------
/**
a
 */
void Forces::VelVerlet1(){
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Typ != 0) continue;
    for(int d=0;d<3;d++){
      //if(fabs(Fm[p].Dir[d]) > 500.) continue;
      Pm[p].Vel[d] += .5*(Fm[p].Dir[d] + Fm[p].Ext[d])*pDeltat();
      Pm[p].Pos[d] += Pm[p].Vel[d]*pDeltat();
      Pm[p].Pos[d] -= floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);
      //printf("%d %d) %lf %g %lf\n",p,d,Pm[p].Pos[d],Fm[p].Dir[d],Pm[p].Vel[d]);
    }
  }
}
/**
a
 */
void Forces::LangevinTherm(){
  double Gamma = Viscosity;//3.*3.14*10.*Viscosity;
  double Zeta = sqrt(12.*2.*1.*Gamma/pDeltat());
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Fm[p].Dir[d] += Zeta*(2.*Mat->Casuale() - 1.);
      //Fm[p].Dir[d] += Zeta*Mat->Gaussiano(0.,1.);
      Fm[p].Dir[d] -= Gamma*Pm[p].Vel[d];
    }
  }
}
/**
a
 */
void Forces::AndersenTherm(){
  double Sigma = sqrt(pTemp());
  for(int p=0;p<pNPart();p++){
    if(Mat->Casuale() < pDeltat()){
      for(int d=0;d<3;d++){
	Pm[p].Vel[d] = Mat->Gaussiano(0.,Sigma);
      }
    }
  }
}
/**
a
 */
void Forces::BerendsenTherm(){
  double Temp = 0.;
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Pm[p].Vel[d] += .5*Fm[p].Dir[d]*pDeltat();
      Temp += SQR(Pm[p].Vel[d]);
    }
  }
  Temp = Temp/(3*pNPart());
  double Norm = 1./sqrt(Temp);
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Pm[p].Vel[d] *= Norm;
    }
  }
}
/**
a
 */
void Forces::ChooseThermostat(int Mode){
  printf("Thermostat: ");
  if(Mode==THERM_LANG){
    CalcTherm = &Forces::LangevinTherm;
    printf("Langevin\n");
  }
  else if((Mode==THERM_AND)){
    CalcTherm = &Forces::AndersenTherm;
    printf("Andersen\n");
  }
  else if((Mode==THERM_BERE)){
    CalcTherm = &Forces::BerendsenTherm;
    printf("Berendsen\n");
  }
  else if((Mode==THERM_NO)){
    CalcTherm = &Forces::NoTherm;
    printf("no thermostat\n");
  }
  else{
   printf("mode not present\n");
   exit(1);
  }
}
/**
a
 */
void Forces::VelVerlet2(){
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Pm[p].Vel[d] += .5*(Fm[p].Dir[d] + Fm[p].Ext[d])*pDeltat();
      Fm[p].Dir[d] = 0.;
    }
  }
}
