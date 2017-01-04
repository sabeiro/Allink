#include "Forces.h"
//-----------------------ENERGY-BONDED-----------------------
double Forces::CalcSpring(int p1){
  double Nrg = 0.;
  double Dist[4];
  int NLink = 0;
  int Link[2] = {-1,-1};
  if(p1 > 0){
    if( Pm[p1].CId == Pm[p1-1].CId ){
      Link[NLink++] = p1 - 1;
    }
  }
  if(p1 <= pNPart()){
    if( Pm[p1].CId == Pm[p1+1].CId ){
      Link[NLink++] = p1 + 1;
    }
  }
  for(int l=0;l<NLink;l++){
    TwoPartDist(p1,Link[l],Dist);
    Nrg += .5*pkSpr()*SQR(Dist[3] - pSprRest());
  }
  return Nrg;
}
/** Calculates the spring and bending force between
       B
      / \
     A   C
*/
double Forces::CalcBonded(int p,double *Pot){
  double NrgSpr = 0.;
  double NrgBen = 0.;
  int pA = p-1;
  int pB = p;
  int pC = p+1;
  int IfBefore = 1;
  int IfAfter = 1;
  if(p > 0){
    if(Pm[p].CId -1 == Pm[p-1].CId){
      pA = p;
      pB = p+1;
      pC = p+2;
      IfBefore = 0;
    }
  }
  if(p <= pNPart() ){
    if(Pm[p].CId + 1 == Pm[p+1].CId || p == pNPart()-1){
      pC = p;
      pB = p-1;
      pA = p-2;
      IfAfter = 0;
    }
  }
  double DistBA[4];
  double DistCB[4];
  double CosAngle = 0.;
  for(int d=0;d<3;d++){
    DistBA[d] = Pm[pB].Pos[d] - Pm[pA].Pos[d];
    DistCB[d] = Pm[pC].Pos[d] - Pm[pB].Pos[d];
  }
  DistBA[3] = (SQR(DistBA[0])+SQR(DistBA[1])+SQR(DistBA[2]));
  DistCB[3] = (SQR(DistCB[0])+SQR(DistCB[1])+SQR(DistCB[2]));
  for(int d=0;d<3;d++)
    CosAngle += DistBA[d]*DistCB[d];
  CosAngle /= (DistBA[3]*DistCB[3]);
  if(IfAfter)
    NrgSpr += .5*pkSpr()*SQR(DistBA[3] - pSprRest());
  if(IfBefore)
    NrgSpr += .5*pkSpr()*SQR(DistCB[3] - pSprRest());
  NrgBen += pkBen()*(1.-CosAngle);
  Pot[0] = NrgSpr;
  Pot[1] = NrgBen;
  return NrgSpr + NrgBen;
}
/** Calculates the bending force between
       B
      / \
     A   C
if the ghost particle C is at least the third in the chain
*/
double Forces::CalcBendingGhost(double *Pos,int p){
  if(pkBen() <= 0.) return 0.;
  double NrgBen = 0.;
  int pA = p;
  int pB = p-1;
  int pC = p-2;
  if(pC < pChain(p)*pNPCh()){
    return 0.;
  }
  double DistBA[3];
  double DistCB[3];
  double DistBA2  = 0.;
  double DistCB2  = 0.;
  double CosAngle = 0.;
  for(int d=0;d<3;d++){
    DistBA[d] = pPos(pB,d) - pPos(pA,d);
    DistBA[d] -= floor(DistBA[d]*pInvEdge(d) + .5)*pEdge(d);
    DistCB[d] = pPos(pC,d) - pPos(pB,d);
    DistCB[d] -= floor(DistCB[d]*pInvEdge(d) + .5)*pEdge(d);
    DistBA2 += SQR(DistBA[d]);
    DistCB2 += SQR(DistCB[d]);
    CosAngle += DistBA[d]*DistCB[d];
  }
  DistCB2 = sqrt(DistCB2);
  DistBA2 = sqrt(DistBA2);
  CosAngle /= (DistBA2*DistCB2);
  NrgBen += pkBen()*(1.-CosAngle);
  return NrgBen;
}
double Forces::CalcBondedCh(int c,double *Pot){
  double DistBA[3];
  double DistCB[3];
  double NrgBend = 0.;
  double NrgSpr = 0.;
  for(int p=c*pNPCh();p<(c+1)*pNPCh()-1;p++){
    double DistBA2 = 0.;
    double DistCB2 = 0.;
    double CosAngle = 0.;
    for(int d=0;d<3;d++){
      DistCB[d] = pPos(p+1,d) + pPos(p,d);
      DistCB[d] -= floor(DistCB[d]*pInvEdge(d) + .5)*pEdge(d);
      DistCB2 += SQR(DistCB[d]);
    }
    DistCB2 = sqrt(DistCB2);
    NrgSpr += .5*pkSpr()*SQR(DistCB2 - pSprRest());
    if(p == c*pNPCh()) continue;
    for(int d=0;d<3;d++){
      DistBA[d] = pPos(p,d) - pPos(p-1,d);
      DistBA[d] -= floor(DistBA[d]*pInvEdge(d) + .5)*pEdge(d);
      DistBA2 += SQR(DistBA[d]);
      CosAngle += DistBA[d]*DistCB[d];
    }
    DistBA2 = sqrt(DistBA2);
    CosAngle /= (DistBA2*DistCB2);
    NrgBend += pkBen()*(1.-CosAngle);
  }
  Pot[0] = NrgSpr;
  Pot[1] = NrgBend;
  return NrgBend + NrgSpr;
}
double Forces::CalcBending(int p1){
  if(Ln[p1].NLink < 1) return 0.;
  if( !(p1%pNPCh()) ) return 0.;
  int l1 = p1+1;
  int l2 = p1-1;
  Vettore v1(pPos(l1,0) - pPos(p1,0),pPos(l1,1) - pPos(p1,1),pPos(l1,2) - pPos(p1,2));
  Vettore v2(pPos(p1,0) - pPos(l2,0),pPos(p1,1) - pPos(l2,1),pPos(p1,2) - pPos(l2,2));
  return pkBen()*(1. - v1.CosAngle(&v2));
}
//%%%%%%%%%%%%%%%%%%%%%%%%%MC%NON%BONDED%%%%%%%%%%%%%%%%%%%%%
/** 
Step potential
 */
double Forces::NrgStep(int p1){
  double DistRel[4];
  double Nrg = 0.;
  for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
    int p2 = Pc->p2Curr;
    Pc->Dist2Curr(DistRel);
    int t2 = pType(p2);
    if(DistRel[3] > Kf.CutOff2) continue;
    if(Pm[p1].Typ == Pm[p2].Typ)
      Nrg += Kf.LJ;
    else
      Nrg += -Kf.LJ;
  }
  return Nrg;
}
double Forces::NrgStepCh(int c,double *Pot){
  Pot[0] = 0.;Pot[1] = 0.;Pot[2] = 0.;
  int p1 = c*pNPCh();
  double DistRel[4];
  double Pos[3];
  //loop in the chain
  for(int p=p1;p<p1+pNPCh();p++){
    int t1 = pType(p);
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      int t2 = pType(p2);
      // avoid self interactions in the chains
      if(p2 > p1 && p2 < p1 + pNPCh() && p2 < p) continue;
      Pc->Dist2Curr(DistRel);
      if(DistRel[3] > Kf.CutOff2) continue;
      // if(DistRel[3] < .1) Pot[2] += 55.;
	// else Pot[2] += -.1;
      Pot[2] += 1.;
    }
  }
  return Pot[2];
}
//--------------------------Dens----------------------
int Forces::RemDens(int pInit,int pEnd){
  double DistRel[4] = {0.,0.,0.,0.};
  int NCutOff = 0;
  double Count = 0.;
  double Pos[3];
  for(int p1=pInit;p1<pEnd;p1++){
    int t1 = pType(p1);
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      int t2 = pType(p2);
      // avoid self interactions in the chains
      if(p2>=pInit && p2<pEnd && p2<=p1) continue;
      Pc->Dist2Curr(DistRel);
      if(DistRel[3] > Kf.CutOff2) continue;
      double Dist = sqrt(DistRel[3]);
      double w2 = Wei2(Dist,pWei2Par());
      double w3 = Wei3(Dist,pWei3Par());
      Dens2[p2*pNType()+t1] -= w2;
      Dens2[p1*pNType()+t2] -= w2;
      Dens3[p2*pNType()+t1] -= w3;
      Dens3[p1*pNType()+t2] -= w3;
      NCutOff++;
    }
  }
  return NCutOff;
}
//Bottleneck!!!
int Forces::AddDens(int pInit,int pEnd){
  double DistRel[4] = {0.,0.,0.,0.};
  int NCutOff = 0;
  double Pos[3];
  for(int p1=pInit;p1<pEnd;p1++){
    int t1 = pType(p1);
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      int t2 = pType(p2);
      // avoid self interactions in the chains
      if(p2>=pInit && p2<pEnd && p2<=p1)continue;
      Pc->Dist2Curr(DistRel);
      if(DistRel[3] > Kf.CutOff2) continue;
      double Dist = sqrt(DistRel[3]);
      double w2 = Wei2(Dist,pWei2Par());
      double w3 = Wei3(Dist,pWei3Par());
      Dens2[p2*pNType()+t1] += w2;
      Dens2[p1*pNType()+t2] += w2;
      Dens3[p2*pNType()+t1] += w3;
      Dens3[p1*pNType()+t2] += w3;
      NCutOff++;
    }
  }
  return NCutOff;
}
double Forces::SumDens(int pInit,int pEnd){
  double Nrg = 0.;
  double OneThird = 1./3.;
  for(int p=pInit;p<pEnd;p++){
    int t1 = pType(p);
    for(int t2=0;t2<pNType();t2++){
      Nrg += .5*Dens2[p*pNType()+t2]*MInt->Coeff(t1,t2);
      for(int t3=0;t3<pNType();t3++){
	Nrg += Dens3[p*pNType()+t2]*Dens3[p*pNType()+t3]*OneThird*MInt->Coeff(t1,t2,t3);
      }
    }
  }
  return Nrg;
}
double Forces::DensFuncNrgGhost(double *Pos,int p1,int t1){
  double DistRel[4] = {0.,0.,0.,0.};
  double Nrg = 0.;
  int NCutOff = 0;
  double OneThird = 1./3.;
  memset(LocDens2,0,pNType()*sizeof(double));
  memset(LocDens3,0,pNType()*sizeof(double));
  for(Pc->SetCurrGhost(Pos);Pc->IfCurrGhost();Pc->NextCurrGhost()){
    int p2 = Pc->p2Curr;
    int t2 = pType(p2);
    if(p1 == p2) continue;
    Pc->Dist2CurrGhost(DistRel);
    if(DistRel[3] > Kf.CutOff2) continue;
    double Dist = sqrt(DistRel[3]);
    double w2 = Wei2(Dist,pWei2Par());
    double w3 = Wei3(Dist,pWei3Par());
    LocDens2[t2] += w2*2.;
    LocDens3[t2] += w3;
    double W3Add[3] = {0.,0.,0.};
    W3Add[t1] = w3;
    for(int t3=0;t3<pNType();t3++){
      for(int t4=0;t4<pNType();t4++){
	double Fact = Dens3[p2*pNType()+t3]*W3Add[t4];
	Fact += W3Add[t3]*Dens3[p2*pNType()+t4];
	Fact += W3Add[t3]*W3Add[t4];
	//Why t2?
	Nrg += Fact*OneThird*MInt->Coeff(t2,t3,t4);
      }
    }
    NCutOff++;
  }
  for(int t2=0;t2<pNType();t2++){
    Nrg += .5*LocDens2[t2]*MInt->Coeff(t1,t2);
    for(int t3=0;t3<pNType();t3++){
      Nrg += LocDens3[t2]*LocDens3[t3]*OneThird*MInt->Coeff(t1,t2,t3);
    }
  }
  return Nrg;
}
double Forces::DensFuncNrgChInternal(int c){
  double DistRel[4] = {0.,0.,0.,0.};
  double Nrg = 0.;
  int NCutOff = 0;
  double OneThird = 1./3.;
  memset(LocDens2,0,pNPCh()*pNType()*sizeof(double));
  memset(LocDens3,0,pNPCh()*pNType()*sizeof(double));
  // for(int pt=0;pt<pNPCh()*pNType();pt++){
  //   LocDens2[pt] = 0.;
  //   LocDens3[pt] = 0.;
  // }
  int pInit = Ch[c].InitBead;
  int pEnd = Ch[c].EndBead;
  double Pos[3];
  for(int p1=pInit,pc=0;p1<pEnd;p1++,pc++){ 
    int t1 = pType(p1);
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      int t2 = pType(p2);
      if(p1 == p2) continue;
      Pc->Dist2Curr(DistRel);
      if(DistRel[3] > Kf.CutOff2) continue;
      double Dist = sqrt(DistRel[3]);
      double w2 = Wei2(Dist,pWei2Par());
      double w3 = Wei3(Dist,pWei3Par());
      LocDens2[pc*pNType()+t2] += w2;
      LocDens3[pc*pNType()+t2] += w3;
      //energy change for the neighbouring chains
      if(p2 < p1 || p2 >= pEnd){
	LocDens2[pc*pNType()+t2] += w2;
	double W3Add[3] = {0.,0.,0.};
	W3Add[t1] = w3;
	//2FIX: multiple change in the 3 order density
	for(int t3=0;t3<pNType();t3++){
	  for(int t4=0;t4<pNType();t4++){
	    double Fact = Dens3[p2*pNType()+t3]*W3Add[t4];
	    Fact += W3Add[t3]*Dens3[p2*pNType()+t4];
	    Fact += W3Add[t3]*W3Add[t4];
	    //Why t2?
	    Nrg += Fact*OneThird*MInt->Coeff(t2,t3,t4);
	  }
	  }
	NCutOff++;
      }
    }
  }
  for(int pt=0;pt<pNPCh();pt++){
    int t1 = pType(pInit+pt);
    for(int t2=0;t2<pNType();t2++){
      Nrg += .5*LocDens2[pt*pNType()+t2]*MInt->Coeff(t1,t2);
      for(int t3=0;t3<pNType();t3++){
	Nrg += LocDens3[pt*pNType()+t2]*LocDens3[pt*pNType()+t3]*OneThird*MInt->Coeff(t1,t2,t3);
      }
    }
  }
  return Nrg;
}
double Forces::DensFuncNrgGhostInternal(double *Pos,int p1,int t1){
  double DistRel[4] = {0.,0.,0.,0.};
  double Nrg = 0.;
  int NCutOff = 0;
  double OneThird = 1./3.;
  memset(LocDens2,0,pNType()*sizeof(double));
  memset(LocDens3,0,pNType()*sizeof(double));
  int c = (int)(p1/pNPCh());
  int pInit = Ch[c].InitBead;
  int pEnd = Ch[c].EndBead;
  for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
    int p2 = Pc->p2Curr;
    int t2 = pType(p2);
    if(p1 == p2) continue;
    if(p2 < pInit) continue;
    if(p2 >= pEnd) continue;
    Pc->Dist2Curr(DistRel);
    if(DistRel[3] > Kf.CutOff2) continue;
    double Dist = sqrt(DistRel[3]);
    double w2 = Wei2(Dist,pWei2Par());
    double w3 = Wei3(Dist,pWei3Par());
    LocDens2[t2] += w2*2.;
    LocDens3[t2] += w3;
    double W3Add[3] = {0.,0.,0.};
    W3Add[t1] = w3;
    for(int t3=0;t3<pNType();t3++){
      for(int t4=0;t4<pNType();t4++){
	double Fact = Dens3[p2*pNType()+t3]*W3Add[t4];
	Fact += W3Add[t3]*Dens3[p2*pNType()+t4];
	Fact += W3Add[t3]*W3Add[t4];
	//Why t2?
	Nrg += Fact*OneThird*MInt->Coeff(t2,t3,t4);
      }
    }
    NCutOff++;
  }
  for(int t2=0;t2<pNType();t2++){
    Nrg += .5*LocDens2[t2]*MInt->Coeff(t1,t2);
    for(int t3=0;t3<pNType();t3++){
      Nrg += LocDens3[t2]*LocDens3[t3]*OneThird*MInt->Coeff(t1,t2,t3);
    }
  }
  return Nrg;
}
double Forces::DensFuncNrgSys(){
  CalcDens(0,pNPart());
  return SumDens(0,pNPart());
}
void Forces::ClearDens(){
  memset(Dens2,0,pNPart()*pNType()*sizeof(double));
  memset(Dens3,0,pNPart()*pNType()*sizeof(double));
}
void Forces::CalcDens(int pInit,int pEnd){
  ClearDens();
  AddDens(pInit,pEnd);
}
//------------------Calc-Ch-Nrg----------------------------
/** Fill the array OldNrgCh with the non bonded, spring and bending energy*/
double Forces::CalcTotNrgCh(){
  Shout("Calculate chains energy");
  double Pot[3] = {0.,0.,0.};
  double Nrg = 0.;
  for(int c=0;c<pNChain();c++){
    /* sustract to the densities the values of the c chain */
    int p1 = c*pNPCh();
    CalcNrgCh(c,Pot);
    OldNrgCh[c*3  ] = Pot[0];
    OldNrgCh[c*3+1] = Pot[1];
    OldNrgCh[c*3+2] = Pot[2];
    Nrg += Pot[0] + Pot[1] + Pot[2];
    //fprintf(StatFile1,"%d %d %lf\n",c,pNPCh(c),Pot[0]);
  }
  return Nrg;
}
/** Calculate the bonded and non bonded energies for the chain c */
double Forces::DensFuncNrgCh(int c,double *Pot){
  int p1 = c*pNPCh();
  Pot[2] = DensFuncNrgChInternal(c);
  for(int p=0;p<pNPCh();p++) Pot[2] += NanoNrg(p+p1);
  return Pot[2];
}
/** Calculate the bonded and non bonded energies for the chain c */
double Forces::NrgChBondDens(int c,double *Pot){
  CalcBondedCh(c,Pot);
  int p1 = c*pNPCh();
  Pot[2] = DensFuncNrgChInternal(c);
  return Pot[0] + Pot[1] + Pot[2];
}
double Forces::DensFuncNrgChAv(int Ch){
  CalcDens(Ch*pNPCh(),(Ch+1)*pNPCh());
  return SumDens(Ch*pNPCh(),(Ch+1)*pNPCh());
}
double Forces::CalcPairwiseCh(int c,double *Pot){
  Pot[0] = 0.;Pot[1] = 0.;Pot[2] = 0.;
  int p1 = c*pNPCh();
  double DistRel[4] = {0.,0.,0.,0.};;
  double Pot1[3] = {0.,0.,0.};
  //loop in the chain
  double Pos[3];
  for(int p=p1;p<p1+pNPCh();p++){
    int t1 = pType(p);
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      Pc->Dist2Curr(DistRel);
      if(p1 >= p2) continue;
      if(p2 > p1 && p2 < p1 + pNPCh() && p2 < p) continue;
      int t2 = pType(p2);
      // avoid self interactions in the chains
      if(DistRel[3] > Kf.CutOff2) continue;
      Potential(DistRel[3],pType(p1),pType(p2),Pot1);
      Pot[2] += Pot1[0];
    }
  }
  return Pot[0] + Pot[1] + Pot[2];
}
//------------------Calc-Part-Nrg----------------------------
/** Fill the array OldNrgCh with the non bonded, spring and bending energy*/
double Forces::CalcTotNrgBead(){
  Shout("Calculate chains energy");
  double Pot[3] = {0.,0.,0.};
  double Nrg = 0.;
  for(int p=0;p<pNChain();p++){
    /* sustract to the densities the values of the c chain */
    CalcNrgBead(p,Pot);
    OldNrgBead[p*3  ] = Pot[0];
    OldNrgBead[p*3+1] = Pot[1];
    OldNrgBead[p*3+2] = Pot[2];
    Nrg += Pot[0] + Pot[1] + Pot[2];
  }
  return Nrg;
}
double Forces::DensFuncNrgBead(int p1){
  double Pos[3];
  pPos(p1,Pos);
  return DensFuncNrgGhost(Pos,p1,pType(p1));
  RemDens(p1,p1+1);
  double Nrg = OldNrgSys - SumDens(0,pNPart());
  AddDens(p1,p1+1);
  return Nrg;
}
/** Calculate the bonded and non bonded energies for the chain c */
double Forces::CalcNrgBeadDensFunc(int p,double *Pot){
  CalcBonded(p,Pot);
  Pot[2] = DensFuncNrgBead(p);
  return Pot[0] + Pot[1] + Pot[2];
}
double Forces::CalcPairwise(int p1,double *Pot){
  //return CheckDomDec(p1);
  double NrgSum = 0.;
  int NCutOff = 0;
  double Pos[3];
  double DistRel[4] = {0.,0.,0.,0.};
  Pot[0] = 0.;Pot[1] = 0.;Pot[2] = 0.;
  for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
  int p2 = Pc->p2Curr;
  // for(int p2 = 0;p2<pNPart();p2++){
    if(p2 == p1) continue;
    //if(p2 >= p1) continue;
    // TwoPartDist(p1,p2,DistRel);
    Pc->Dist2Curr(DistRel);
    if(DistRel[3] > Kf.CutOff2) continue;
    Potential(DistRel[3],pType(p1),pType(p2),Pot);
    NrgSum += Pot[0] + Pot[1] + Pot[2];
    NCutOff++;
  }
  return NrgSum;
}
//----------------------MOLECULAR-DYNAMICS----------------------
void Forces::ChooseCalcMode(int Mode){
  printf("Calculation mode: ");
  if(VAR_IF_TYPE(Mode,CALC_PAIR)){
    NrgBead = &Forces::CalcPairwise;
    NrgCh = &Forces::CalcPairwiseCh;
    printf("Pairwise interactions\n");
  }
  else if(VAR_IF_TYPE(Mode,CALC_DENS)){
    NrgBead = &Forces::CalcNrgBeadDensFunc;
    NrgCh = &Forces::DensFuncNrgCh;
    Kf.CutOff2 = 1.;
    printf("Lipid model particle energy\n");
  }
  else if(VAR_IF_TYPE(Mode,CALC_DENS_CH)){
    NrgBead = &Forces::CalcNrgBeadDensFunc;
    NrgCh = &Forces::DensFuncNrgCh;
    Kf.CutOff2 = 1.;
    printf("Lipid model chain energy\n");
  }
  else{
    printf("mode not present\n");
    exit(1);
  }
}
void Forces::ChoosePot(int Mode){
  printf("Potential: ");
  if(VAR_IF_TYPE(Mode,CALC_LJ)){
    CalcPot = &Forces::LJPot;
    printf("Lennard Jones\n");
  }
  else if(VAR_IF_TYPE(Mode,CALC_LJ39)){
    CalcPot = &Forces::LJ39;
    printf("Lennard Jones 3-9\n");
  }
  else if(VAR_IF_TYPE(Mode,CALC_HARM)){
    CalcPot = &Forces::Harmonic;
    printf("Harmonic\n");
  }
  else if(VAR_IF_TYPE(Mode,CALC_STEP)){
    CalcPot = &Forces::StepPot;
    printf("Step\n");
  }
  else if(VAR_IF_TYPE(Mode,CALC_ELECTRO)){
    CalcPot = &Forces::ElectroPot;
    printf("Step\n");
  }
  else{
   printf("not present\n");
   exit(1);
  }
  DefForceParam();
}
double Forces::SumForcesMD(){
  double Pot[3];
  double DistRel[4] = {0.,0.,0.,0.};
  double Nrg = 0.;
  for(int p1=0;p1<pNPart();p1++){
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      if(p1 >= p2) continue;
      Pc->Dist2Curr(DistRel);
      if(DistRel[3] > Kf.CutOff2) continue;
      double InvDist = 1./sqrt(DistRel[3]);
      double Cons = Potential(DistRel[3],pType(p1),pType(p2),Pot);
      for(int d=0;d<3;d++){
	Fm[p1].Dir[d] += Cons*DistRel[d]*InvDist;
	Fm[p2].Dir[d] -= Cons*DistRel[d]*InvDist;
      }
      Nrg += Pot[0];
    }
  }
  return Nrg;
}
#include <sys/time.h>
void Forces::CheckPairList(){
  int NAll = pNPart();
  int *PDom = (int *)calloc(NAll*NAll,sizeof(int));
  if(PDom == NULL){printf("Could not alloc PDom\n");return;};
  int *PLoop = (int *)calloc(NAll*NAll,sizeof(int));
  if(PLoop == NULL){printf("Could not alloc PLoop\n");return;};
  int NDom = 0;
  int NLoop = 0;
  double NPairDom  = 0.;
  double NPairLoop = 0.;
  double Dist[3];
  double Nrg = 0.;
  double Pot[2];
  double DistRel[4] = {0.,0.,0.,0.};
  double TimeDomDec;
  double TimeLoop;
  timespec TimeInit;
  timespec TimeEnd;
//  clock_gettime(CLOCK_REALTIME, &TimeInit);
  double Pos[3]; 
  for(int p1=0;p1<pNPart();p1++){
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      if(p2 <= p1) continue;
      Pc->Dist2Curr(DistRel);
      NPairDom += 1.;
      if(DistRel[3] > Kf.CutOff2) continue;
      PLoop[NDom*2+0] = p2;
      PLoop[NDom*2+1] = p2;
      NDom++;
    }
  }
//  clock_gettime(CLOCK_REALTIME, &TimeEnd);
  TimeDomDec = (double)(TimeInit.tv_nsec - TimeEnd.tv_nsec);
//  clock_gettime(CLOCK_REALTIME, &TimeInit);
  for(int p1=0;p1<pNPart();p1++){
    for(int p2=p1+1;p2<pNPart();p2++){
      double Dist2 = 0.;
      for(int d=0;d<3;d++){
	Dist[d] = pPos(p1,d) - pPos(p2,d);
	Dist[d] -= floor(Dist[d]*pInvEdge(d) + .5)*pEdge(d);
	Dist2 += SQR(Dist[d]);
      }
      NPairLoop += 1.;
      if(Dist2 > Kf.CutOff2) continue;
      PLoop[NLoop*2+0] = p1;
      PLoop[NLoop*2+1] = p2;
      NLoop++;
    }
  }
//  clock_gettime(CLOCK_REALTIME, &TimeEnd);
  TimeLoop = (double)(TimeInit.tv_nsec - TimeEnd.tv_nsec);
  //Mat->Sort(PDom,NDom);
  printf("[Pairs] %d=%d Gain (pair): %lf (time): %lf\n",NDom,NLoop,NPairLoop/NPairDom,TimeLoop/TimeDomDec);
  for(int p=0;p<MAX(NDom,NLoop);p++){
    //printf("%d %d %d\n",p,PDom[p],PLoop[p]);
    //if(p >= NDom)
    //printf(" %d %d\n",PLoop[p],Pc->pCell(PLoop[p]));
  }
  free(PDom);
  free(PLoop);
}
/** Compare the pair found in the domain decomposition method with the ones found with the simple N^2 neighbours search */
double Forces::CheckDomDec(int p1){
  const int NAll = pNPart();
  int *PDom = new int[NAll];
  int *PCell = new int[NAll];
  int NDom = 0;
  int *PLoop = new int[NAll];
  int NLoop = 0;
  double NCell = 0.;
  double Pos[3];
  double DistRel[4] = {0.,0.,0.,0.};
  double Nrg=0.;
  double Pot[2];
  for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
    int p2 = Pc->p2Curr;
    if(p2 == p1) continue;
    Pc->Dist2Curr(DistRel);
    if(DistRel[3] > Kf.CutOff2) continue;
    PDom[NDom++] = p2;
    PCell[NDom] = Pc->cCurr;
  }
  for(int p2=0;p2<pNPart();p2++){
    if(p1 == p2) continue; 
    double Dist2 = 0.;
    for(int d=0;d<3;d++){
      DistRel[d] = pPos(p1,d) - pPos(p2,d);
      DistRel[d] -= floor(DistRel[d]*pInvEdge(d) + .5)*pEdge(d);
      Dist2 += SQR(DistRel[d]);
    }
    if(Dist2 > Kf.CutOff2) continue;
    PLoop[NLoop++] = p2;
  }
  //Mat->Sort(PDom,NDom);
  printf("-------------%d=%d Ratio %lf\n",NDom,NLoop,pNPart()/NCell);
  for(int p=0;p<MAX(NDom,NLoop);p++){
    printf("%d %d %d %d\n",p,PDom[p],PLoop[p],PCell[p]);
    //if(p >= NDom)
    //printf(" %d %d\n",PLoop[p],Pc->pCell(PLoop[p]));
  }
  printf("nei \n");
  delete [] PDom;
  delete [] PCell;
  delete [] PLoop;
  return 0.;
}
