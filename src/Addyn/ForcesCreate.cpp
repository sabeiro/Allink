#include "Forces.h"
void Forces::CreateInitial(){
  if(VAR_IF_TYPE(SysShape,SYS_1D)){
    Create1d();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_2D)){
    Create2d();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_3D)){
    Create3d();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_STALK)){
    CreateStalk();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_LEAVES)){
    CreateLeaves();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_PORE)){
    CreatePore();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_ROD)){
    CreateRod();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_RIGID)){
    CreateRigid();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_MD)){
    CreateMD();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_MC)){
    CreateMC();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_ELECTRO)){
    CreateElectro();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_TRIAL)){
    int NSect = 3;
    double Pos[3];
    for(int d=0;d<3;d++){
      Pos[d] = .5*pEdge(d)/(double)NSect;
    }
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
	Pm[p].Pos[d] = Pos[d];
      }
      Pm[p].Typ = 0;
      Pos[0] += pEdge(0)/(double)NSect;
      if(Pos[0] > pEdge(0)){
	Pos[0] = .5*pEdge(0)/(double)NSect;
	Pos[1] += pEdge(1)/(double)NSect;
	if(Pos[1] > pEdge(1)){
	  Pos[1] = .5*pEdge(1)/(double)NSect;
	  Pos[2] += pEdge(2)/(double)NSect;
	}
      }
    }
    // Pm[0].Pos[0] = .1;Pm[0].Pos[1] = .5;Pm[0].Pos[2]= .5;
    // Pm[1].Pos[0] = .2;Pm[1].Pos[1] = .6;Pm[1].Pos[2]= .7;
    // Pm[2].Pos[0] = .3;Pm[2].Pos[1] = .9;Pm[2].Pos[2]= .2;
    // Pm[3].Pos[0] = .4;Pm[3].Pos[1] = .3;Pm[3].Pos[2]= .5;
    // Pm[4].Pos[0] = .5;Pm[4].Pos[1] = .4;Pm[4].Pos[2]= .6;
  }
  else{
    printf("System shape not recognized %d\n",SysShape);
    return ;
  }
  Old2Move = Bead2Move;
  VAR_ADD_TYPE(SysType,VAR_SYS_TXVL);
  VAR_ADD_TYPE(SysType,VAR_EDGE);
}
void Forces::Create2d(){
  double Dx = pEdge(0)/(double)(nEdge[0]);
  double Dy = pEdge(1)/(double)(nEdge[1]);
  for(int px=0;px<nEdge[0];px++){
    for(int py=0;py<nEdge[1];py++){
      int p = px*nEdge[1]+py;
      Pm[p].Idx = p;
      Pm[p].Pos[0] = Dx*(double)px + .5*Dx;
      Pm[p].Pos[1] = Dy*(double)py + .5*Dy;
      Pm[p].Pos[2] = 0.;
      Pm[p].CId = py;
      Pm[p].Typ = 0;
      if(px == 0 && BoundCond[0]){
      	Pm[p].Typ = 2;
      }
      if(px == nEdge[0]-1 && BoundCond[1]){
      	Pm[p].Typ = 2;
      }
      if(py == 0 && BoundCond[2]){
      	Pm[p].Typ = 2;
      }
      if(py == nEdge[1]-1 && BoundCond[3]){
      	Pm[p].Typ = 2;
      }
      Ln[p].NLink = 4;
      int pym1 = p-1;
      if( py-1 < 0 ) pym1 += nEdge[1];
      if(pym1 >= pNPart() ) pym1 -= pNPart();
      Ln[p].Link[0] = pym1;
      int pyp1 = p+1;
      if( py+1 >= nEdge[1]) pyp1 -= nEdge[1];
      if(pyp1 < 0) pyp1 += pNPart();
      Ln[p].Link[1] = pyp1;
      int pxm1 = p-nEdge[1];
      if(pxm1 < 0) pxm1 += pNPart();
      Ln[p].Link[2] = pxm1;
      int pxp1 = p+nEdge[1];
      if(pxp1 > pNPart()-1) pxp1 -= pNPart();
      Ln[p].Link[3] = pxp1;
      //printf("%d) %d %d %d %d \n",p,Ln[p].Link[0],Ln[p].Link[1],Ln[p].Link[2],Ln[p].Link[3]);
    }
  }
  //Ln[0].Link[0] = nEdge[0]-1;
  AddRigid();
}
void Forces::Create3d(){
  double Dx = pEdge(0)/(double)(NEdge-1);
  double Dy = pEdge(1)/(double)(NEdge-1);
  double Dz = pEdge(2)/(double)(NEdge-1);
  Kf.Elong[1] = Dx;
  Kf.Elong[2] = Dy;
  Kf.Elong[0] = Dz;
  Kf.El[0] = 11.;//Elastic coupling
  Kf.El[1] = 11.;//Elastic coupling
  Kf.El[2] = 11.;//Elastic coupling
  for(int px=0,ppp=0;px<NEdge;px++){
    for(int py=0;py<NEdge;py++){
      for(int pz=0;pz<NEdge;pz++){
	Pm[ppp].Idx = ppp;
	Pm[ppp].Pos[0] = Dx*(double)px;
	Pm[ppp].Pos[1] = Dy*(double)py;//Mate->Casuale();//Dy*(double)p;
	Pm[ppp].Pos[2] = Dz*(double)pz;//Mate->Casuale();
	//     if( p == 20)
	//       Pm[p].Pos[2] = .01;
	int link = 0;
	if(pz != NEdge-1){	   
	  Ln[ppp].Link[link] = ppp+1;
	  link++;
	}
	if(pz != 0){
	  Ln[ppp].Link[link] = ppp-1;
	  link++;
	}
	if(py != 0){
	  Ln[ppp].Link[link] = ppp-NEdge;
	  link++;
	}
	if(py != NEdge-1){
	  Ln[ppp].Link[link] = ppp+NEdge;
	  link++;
	}
	if(px != 0){
	  Ln[ppp].Link[link] = ppp-NEdge*NEdge;
	  link++;
	}
	if(px != NEdge-1){
	  Ln[ppp].Link[link] = ppp+NEdge*NEdge;
	  link++;
	}
	Ln[ppp].NLink = link;
	if( (px == 0 || px == NEdge -1) &&
	    (py == 0 || py == NEdge -1) &&
	    (pz == 0 || pz == NEdge -1))
	  Pm[ppp].Typ = 2;
	//for(int l=0;l<Ln[ppp].NLink;l++)
	  //printf("%d %d %lf %lf %lf\n",ppp,Pm[ppp].Link[l],Pm[ppp].Pos[0] - Pm[Pm[ppp].Link[l]].Pos[0],Pm[ppp].Pos[1] - Pm[Pm[ppp].Link[l]].Pos[1],Pm[ppp].Pos[2] - Pm[Pm[ppp].Link[l]].Pos[2]);
	ppp++;
      }
    }
  }
}
void Forces::CreateLeaves(){
  double Dx = pEdge(0)/(double)(NEdge-1);
  double Dy = pEdge(1)/(double)(NEdge-1);
  double Dz = pEdge(2)/(double)(NEdge-1);
  Bead2Move = 0;
  for(int p=0;p<NEdge;p++){
    Kf.Elong[0] = Dx;
    Pm[p].Idx = p;
    Pm[p].Pos[0] = Dx*p;
    Pm[p].Pos[1] = .5*pEdge(1);
    Pm[p].Pos[2] = Kf.Elong[2];
    Pm[p].CId = 0;
    Pm[p].Typ = 0;    
    if(p == Bead2Move) Pm[p].Typ = 1;
    Ln[p].NLink = 2;
    if(p == 0){
      if(BoundCond[0])
	Ln[p].Link[0] = NEdge-1;
      else
	Ln[p].NLink = 1;
      Pm[p].Typ = 2;
    }
    else
      Ln[p].Link[0] = p-1;
    if(p == NEdge-1){
      if(BoundCond[1])
	Ln[p].Link[1] = 0;
      else
	Ln[p].NLink = 1;
      Pm[p].Typ = 2;
    }
    else
      Ln[p].Link[1] = p+1;
  }
  AddRigid();
}
void Forces::CreateRod(){
  double Dx = pEdge(0)/(double)(NEdge-1);
  double Dy = pEdge(1)/(double)(NEdge-1);
  double Dz = pEdge(2)/(double)(NEdge-1);
  for(int p=0;p<NEdge;p++){
    Pm[p].Idx = p;
    Pm[p].Pos[0] = p*Dx*.5;
    Pm[p].Pos[1] = .5*pEdge(1);
    Pm[p].Pos[2] = .5*pEdge(2);
    Pm[p].CId = 0;
    Pm[p].Typ = 0;
    if(p < NEdge-1){
      Ln[p].NLink = 1;
      Ln[p].Link[0] = p+1;
    }
  }
  for(int p=0;p<2;p++){
    Pm[p].Typ = 2;
  }
  //SetkSpr(0.);
  SetkBen(Kf.SLap);
  Bead2Move = NEdge - 1;
  Pm[Bead2Move].Typ = 1;
  Pm[Bead2Move-1].Typ = 1;
}
void Forces::CreatePore(){
  double Dx = pEdge(0)/(double)(NEdge-1);
  double Dy = pEdge(1)/(double)(NEdge-1);
  double Dz = pEdge(2)/(double)(NEdge-1);
  Kf.Elong[0] = Dx;
  Bead2Move = 0;
  int NSegment = (int)(NEdge/6.);
  double AngleS = .25*DUE_PI/(double)NSegment;
  for(int c=0;c<pNChain();c++){
    for(int p=c*NEdge;p<NEdge*(c+1);p++){
      Pm[p].Pos[0] = Dx*(double)(p-c*NEdge);
      Pm[p].Pos[1] = .5*pEdge(1);
      Pm[p].Pos[2] = Kf.Elong[2]*(c-.5)+.5*pEdge(2);
      if(p >= c*NEdge + 2*NSegment && p <= c*NEdge + 4*NSegment){
	//double pa = (double)(p-c*NEdge + 4*NSegment);
	double pa = (double)(p-2*NSegment);
	double x = Kf.Elong[2]*.5*sin(AngleS*pa);
	double z = Kf.Elong[2]*.5*cos(AngleS*pa);
	Pm[p].Pos[0] = x + Dx*2*NSegment;
	Pm[p].Pos[2] = -z + Kf.Elong[2]*(c)+.5*pEdge(2);
      }
      else if(p > c*NEdge + 4*NSegment){
	Pm[p].Pos[0] = Dx*(6*NSegment - p);
	Pm[p].Pos[2] = Kf.Elong[2]*(c+.5)+.5*pEdge(2);
      }
      Pm[p].Idx = p;
      Pm[p].CId = c;
      Pm[p].Typ = 0;
      if(p == Bead2Move) Pm[p].Typ = 1;
      Ln[p].NLink = 3;
      if(p == c*NEdge ){
	Ln[p].Link[0] = p + 1;
      }
      else
	Ln[p].Link[0] = p - 1;
      if(p == NEdge*(c+1) - 1){
	Ln[p].Link[1] = p - 1;
      }
      else 
	Ln[p].Link[1] = p + 1;
      Ln[p].Link[2] = NEdge - p - 1;
    }
  }
  int pHalf = (int)(NEdge/2.);
  Ln[pHalf].NLink = 2;
  Pm[0].Typ = 2;
  Pm[1].Typ = 2;
  Pm[NEdge-2].Typ = 2;
  Pm[NEdge-1].Typ = 2;
  Pm[pHalf-1].Typ = 2;
  Pm[pHalf].Typ = 2;
  Pm[pHalf+1].Typ = 2;
  AddRigid();
  // for(int p=0;p<pNPart();p++)
  //   printf("%d) %lf %lf %lf %d) %d %d %d\n",p,Pm[p].Pos[0],Pm[p].Pos[1],Pm[p].Pos[2],Ln[p].NLink,Ln[p].Link[0],Ln[p].Link[1],Ln[p].Link[2]);
}
void Forces::CreateStalk(){
  double Dx = pEdge(0)/(double)(NEdge-1);
  double Dy = pEdge(1)/(double)(NEdge-1);
  double Dz = pEdge(2)/(double)(NEdge-1);
  Bead2Move = 0;
  int c=0;
  int NSegment = (int)(NEdge/3.);
  double AngleS = .5*DUE_PI/(double)NSegment;
  for(int c=0;c<pNChain();c++){
    for(int p=c*NEdge;p<NEdge*(c+1);p++){
      Kf.Elong[0] = Dx;
      Pm[p].Idx = p;
      Pm[p].Pos[0] = Dx*(double)(p-c*NEdge);
      Pm[p].Pos[1] = .5*pEdge(1);
      Pm[p].Pos[2] = Kf.Elong[2]*(double)c+.5*pEdge(2)-1.5*Kf.Elong[2];
      Pm[p].CId = c;
      Pm[p].Typ = 0;
      if(c==1){
	if(p > c*NEdge + NSegment && p <= c*NEdge + 2*NSegment ){
	  double x = Kf.Elong[2]*.5*sin(AngleS*(p-c*NEdge + NSegment));
	  double z = Kf.Elong[2]*.5*cos(AngleS*(p-c*NEdge + NSegment));
	  Pm[p].Pos[0] = x + Dx*NSegment;
	  Pm[p].Pos[2] = z + Kf.Elong[2]*(double)c+.5-1.*Kf.Elong[2];//+ Kf.Elong[2]/(double)NSegment;
	}
	else if(p > c*NEdge + 2*NSegment ){
	  Pm[p].Pos[0] = Dx*(NEdge*(c+1)-p);
	  Pm[p].Pos[2] = Kf.Elong[2]*(double)(c+1)+.5-1.5*Kf.Elong[2];
	}
      }
      if(c==2){
	if(p < c*NEdge + NSegment ){
	  Pm[p].Pos[0] =  1. - Dx*(double)(p-c*NEdge);
	  Pm[p].Pos[2] = Kf.Elong[2]*(double)(c-1)+.5-1.5*Kf.Elong[2];
	}
	else if(p >= c*NEdge + NSegment && p < c*NEdge + 2*NSegment ){
	  double x = Kf.Elong[2]*.5*sin(AngleS*(p-c*NEdge + NSegment));
	  double z = Kf.Elong[2]*.5*cos(AngleS*(p-c*NEdge + NSegment));
	  Pm[p].Pos[0] = - x + 2*Dx*NSegment;
	  Pm[p].Pos[2] = z + Kf.Elong[2]*(double)(c-1)+.5-1.*Kf.Elong[2];	  }
      }
      if(p == Bead2Move) Pm[p].Typ = 1;
      Ln[p].NLink = 3;
      if(p == c*NEdge ){
	Ln[p].Link[0] = p + 1;
	Pm[p].Typ = 2;
      }
      else
	Ln[p].Link[0] = p - 1;
      if(p == NEdge*(c+1) - 1){
	Ln[p].Link[1] = p - 1;
	Pm[p].Typ = 2;
      }
      else 
	Ln[p].Link[1] = p + 1;
      if(c == 0)
	Ln[p].Link[2] = p + NEdge;
      if(c == 1)
	Ln[p].Link[2] = p - NEdge;
      if(p == c*NEdge + 1 || p == (c+1)*NEdge - 2)
	Pm[p].Typ =2;
      //	if(p == c*NEdge + 2 || p == (c+1)*NEdge - 3)
      //Pm[p].Typ =2;
    }
  }
}
void Forces::Create1d(){
  double Dx = pEdge(0)/(double)(NEdge-1);
  double Dy = pEdge(1)/(double)(NEdge-1);
  double Dz = pEdge(2)/(double)(NEdge-1);
  Bead2Move = NEdge/2;
  for(int p=0;p<pNPart();p++){
    Pm[p].Idx = p;
    Pm[p].Pos[0] = Dx*(double)p;
    Pm[p].Pos[1] = .5;//Mate->Casuale();//Dy*(double)p;
    Pm[p].Pos[2] = .45;//Mate->Casuale();
    Ln[p].NLink = 2;
    Pm[p].Typ = 0;
    if(p == Bead2Move) Pm[p].Typ = 1;
    int link=0;
    if(p == 0 || p==1){
      Ln[p].NLink--;
      Pm[p].Typ = 2;
    }
    else{
      Ln[p].Link[link] = p - 1;
      link++;
    }
    if(p == pNPart()-1 || p == pNPart()-2){
      Ln[p].NLink--;
      Pm[p].Typ = 2;
    }
    else{
      Ln[p].Link[link] = p + 1;
      link++;
    }
  }
}
void Forces::CreateRigid(){
  SetNNano(2);
  for(int n=0;n<pNNano();n++){
    for(int d=0;d<3;d++){
      Nano[n].Pos[d] = .5*pEdge(d);
      Nano[n].Vel[d] = .0;
      Nano[n].AVel[d] = 0.;
      Nano[n].Axis[d] = 0.;
    }
    Nano[n].Axis[0] = 1.;
    Nano[n].Shape = 2;
    Nano[n].Rad = .03;
    Nano[n].Height = .3;
    Nano[n].Mass = 1.;
    Nano[n].Gamma = 10.;
    Nano[n].Zeta = 30.;
  }
  Nano[0].Mass = 1.;
  Nano[1].Shape = 1;
  for(int d=0;d<3;d++){
    Nano[1].Pos[d] = .0*pEdge(d);
    Nano[1].Vel[d] = .01;
  }
}
void Forces::CreateMC(){
  int NSect[3] = {3,3,3};
  StatFile1 = fopen("StatisticsMC1.dat","w");
  StatFile2 = fopen("StatisticsMC2.dat","w");
  NInsertion = 0;
  NRemoval = 0;
  DefNanoForceParam();
  PrintForce();
  SetNChain(NEdge);
  SetNPCh(1);
  OldNrgBead = new double[pNPart()];
  OldNrgCh = new double[pNChain()];
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  for(int d=0;d<3;d++) NSect[d] = (int)(Edge[d]/(double)(2.*sqrt(Kf.CutOff2)));
  double Dens = pNPart()/pVol();
  double Lambda3 = 1.;//CUBE(sqrt(DUE_PI/SQR(hPlanck)));
  GaussVar = sqrt(SQR(pReOverCutOff())/(double)(Block[0].NPCh-1.)/3.)/2.;
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Pm[p].Pos[d] = Mat->Casuale()*pEdge(d);
    }
  }
  ChooseCalcMode(CalcMode);
  ChoosePot(CalcMode);
}
void Forces::CreateMD(){
  int NSect[3] = {3,3,3};
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  double Dens = pNPart()/pVol();
  double Lambda3 = 1.;//CUBE(sqrt(DUE_PI/SQR(hPlanck)));
  for(int d=0;d<3;d++)NSect[d] = (int)(Edge[d]/(double)(2.*sqrt(Kf.CutOff2)));
  GaussVar = sqrt(SQR(pReOverCutOff())/(double)(Block[0].NPCh-1.)/3.)/2.;
  DefNanoForceParam();
  int n3 = 2;
  while ((n3*n3*n3)<pNPart()) n3++;
  int iix=0;
  int iiy=0;
  int iiz=0;
  for(int p=0;p<pNPart();p++){
    Pm[p].Pos[0] = ((double)iix+0.5)*pEdge(0)/n3;
    Pm[p].Pos[1] = ((double)iiy+0.5)*pEdge(1)/n3;
    Pm[p].Pos[2] = ((double)iiz+0.5)*pEdge(2)/n3;
    iix++;
    if (iix==n3) {
      iix=0;
      iiy++;
      if (iiy==n3) {
	iiy=0;
	iiz++;
      }
    }
    for(int d=0;d<3;d++){
      Pm[p].Vel[d] = Mat->Gaussiano(0.,1.0);
    }
  }
  ChooseCalcMode(CalcMode);
  ChoosePot(CalcMode);
}
void Forces::CreateElectro(){
  int NSect[3] = {3,3,3};
  StatFile1 = fopen("StatisticsMC1.dat","w");
  StatFile2 = fopen("StatisticsMC2.dat","w");
  NInsertion = 0;
  NRemoval = 0;
  DefNanoForceParam();
  PrintForce();
  SetNChain(NEdge);
  SetNPCh(1);
  OldNrgBead = new double[pNPart()];
  OldNrgCh = new double[pNChain()];
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  for(int d=0;d<3;d++) NSect[d] = (int)(Edge[d]/(double)(2.*sqrt(Kf.CutOff2)));
  double Dens = pNPart()/pVol();
  double Lambda3 = 1.;//CUBE(sqrt(DUE_PI/SQR(hPlanck)));
  GaussVar = sqrt(SQR(pReOverCutOff())/(double)(Block[0].NPCh-1.)/3.)/2.;
  int NCenter = 6;
  double *PCenter = (double *)calloc(NCenter*2,sizeof(double));
  GaussVar = Kf.Elong[0];
  for(int c=0;c<NCenter;c++){
    for(int d=0;d<2;d++){
      PCenter[c*2+d] = Mat->Casuale()*pEdge(d)*.7 + pEdge(d)*.2;
    }
  }
  for(int p=0;p<NEdge;p++){
    int c = (p%NCenter);
    for(int d=0;d<2;d++){
      Pm[p].Pos[d] = Mat->Gaussiano(PCenter[c*2+d],GaussVar);
      Pm[p].Pos[d] -= floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);
      // if(Pm[p].Pos[d] > pEdge(d)) Pm[p].Pos[d] += pEdge(d) - Pm[p].Pos[d];
      // if(Pm[p].Pos[d] < 0) Pm[p].Pos[d] *= -1.;
    }
    Pm[p].Pos[2] = .2*pEdge(2)*sin(6*Pm[p].Pos[0]*pInvEdge(0))*sin(4*Pm[p].Pos[1]*pInvEdge(1)) + .5*pEdge(2);
  }
  GaussVar = Kf.Elong[1];
  for(int p=NEdge;p<NEdge+NSpline;p++){
    for(int d=0;d<2;d++){
      Pm[p].Pos[d] = (p-NEdge)/(double)NSpline*pEdge(d);
      // double Move = Mat->Gaussiano(Pm[p-1].Pos[d],GaussVar);
      // if(p == NEdge){ Pm[p].Pos[d] = 0.; Move = 0.;}
      // if(Move > pEdge(d) ) Pm[p].Pos[d] = Move - pEdge(d);
      // else if(Move < 0. ) Pm[p].Pos[d] = -Move;
      // else Pm[p].Pos[d] = Move;
      // Pm[p].Pos[d] -= floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);      
    }
    Pm[p].Pos[2] = .2*pEdge(2)*sin(6*Pm[p].Pos[0]*pInvEdge(0))*sin(4*Pm[p].Pos[1]*pInvEdge(1)) + .5*pEdge(2);
  }
  ChooseCalcMode(CalcMode);
  ChoosePot(CalcMode);
}
