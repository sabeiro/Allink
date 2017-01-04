/***********************************************************************
VarDataCreate: Functions that create an initial configuration system as read from Polymer.conf. 
Copyright (C) 2008-2010 by Giovanni Marelli <sabeiro@virgilio.it>


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#include "../include/VarData.h"

static int Tries = 0;

int VarData::DefSoft(char *nome2,char *ConfF){
  /*---------Variables--------------------*/
  char cSystem[512];
  double ReOSigma = 5.;
  int NChain = 1250;
  //  Dens = 5.;
  double ChainPArea = 1./(0.01974*QUAD(ReOSigma));
  double Thickness = 0.849*ReOSigma;
  double Volume = 1.;
  double Area =0.;
  int *arch;
  if(ReadConf(ConfF)) return 1;
  /*----------Setting----------*/
  Gen->NChain = 0;
  Gen->NPart  = 0;
  Gen->Step   = 0;
  Gen->NBlock = 0;
  //if(Gen->kappaBend  <= 0.) Gen->kappaSpring = 3.*(Gen->NPCh-1)/(2.*QUAD(Gen->ReOverCutOff));
  VAR_ADD_TYPE(SysType,VAR_EDGE);
  double sigma = sqrt(QUAD(Gen->ReOverCutOff)/(Gen->NPCh-1)/3.);
  sigma = 1./sqrt(pkSpr());
  Mat->InizializzaGaussiano(sigma,100);
  int AddNSoft = 0;
  for(int s=0;s<NSoft;s++){
    //------------------------16---------------------------
    if (Gen->NPCh == 16){
      ChainPArea = 2.6*Soft[s].Size[2];
      Thickness = 0.93*ReOSigma;
    }
    //------------------------11------------------------
    else if (Gen->NPCh == 11){
      ChainPArea = 3.19*Soft[s].Size[2];
      Thickness = 2.55;
    }
    //------------------------12------------------------
    else if (Gen->NPCh == 12){
      ChainPArea = 3.19*Soft[s].Size[2];
      Thickness = 2.55;
    }
    //------------------------13------------------------
    else if (Gen->NPCh == 13){
      ChainPArea = 3.19*Soft[s].Size[2];
      Thickness = 2.55;
    }
    //------------------------10------------------------
    else if (Gen->NPCh == 10){
      ChainPArea = 3.98*Soft[s].Size[2];
      Thickness = 0.8*ReOSigma;
    }
    //------------------------14------------------------
    else if (Gen->NPCh == 14){
      ChainPArea = 3.65*Soft[s].Size[2];
      Thickness = 0.5*ReOSigma;
    }
    //------------------------20------------------------
    else if (Gen->NPCh==21||Gen->NPCh==20){
      ChainPArea = 3.65*Soft[s].Size[2];
      Thickness = 0.6*ReOSigma;
    }
    //------------------------32---------------------
    else if (Gen->NPCh == 32){
      ChainPArea = 2.6*Soft[s].Size[2];
      Thickness = 3.5;//1.14*ReOSigma;
    }
    if(VAR_IF_TYPE(SysType,VAR_TWOTAILS)){
      if (Gen->NPCh == 21 ||Gen->NPCh == 20 ){
	ChainPArea = 3.65*Soft[s].Size[2];
	Thickness = 0.6*ReOSigma;
      }
      if (Gen->NPCh == 11||Gen->NPCh == 12||Gen->NPCh == 13){
	ChainPArea = 3.19*Soft[s].Size[2];
	Thickness = 2.55;
      }
    }
    //-------------------arch--------------------
    arch = (int *)calloc(Gen->NPCh,sizeof(int));
    for(int b=0;b<Block[0].Asym;b++){
      arch[b] = 0;
    }
    for(int b=Block[0].Asym;b<Gen->NPCh;b++){
      arch[b] = 1;
    }
    //--------------------Area---------------------
    Volume = Gen->Edge[CLat1]*Gen->Edge[CLat2]*Thickness;
    if(VAR_IF_TYPE(Soft[s].Topology,VAR_PLANAR)){
      Area = Gen->Edge[CLat1]*Gen->Edge[CLat2];
    }
    else if(VAR_IF_TYPE(Soft[s].Topology,VAR_COATING)){
      Thickness = 4.0;//1.15;
      Area = .7*DUE_PI*(Nano->Rad+Thickness)*Nano->Height;
      if(VAR_IF_TYPE(Nano->Shape,SHAPE_SPH)){
	Area = .3*2.*DUE_PI*SQR(Nano->Rad+Thickness);
	if (Gen->NPCh == 32){
	  Thickness = 6.0;//1.15;
	  Area = .1*2.*DUE_PI*SQR(Nano->Rad+Thickness);
	}
      }
    }
    else if(VAR_IF_TYPE(Soft[s].Topology,VAR_TUBE)){
      Area = DUE_PI*(Soft[s].Size[0]+Thickness)*Soft[s].Size[1];
    }
    else if(VAR_IF_TYPE(Soft[s].Topology,VAR_VESICLE)){
      ChainPArea = 3.6/Soft[s].Size[2];
      //Thickness = 3.5;
      Area  = DUE_PI*SQR(Soft[s].Size[0]+.5*Thickness);
      Area += DUE_PI*SQR(Soft[s].Size[0]-.5*Thickness);
      AddNSoft += 1;
    }
    for(int n=0;n<Gen->NNano;n++){
      if(Nano[n].Shape == SHAPE_NONE) continue;
      if(Nano[n].Shape == SHAPE_WALL) continue;
      Area -= PI*SQR(Nano[n].Rad);
    }
     //------------------NChain-------------------
    if(Nano->Shape == SHAPE_SPH){
      Volume -= 4.*PI*CUB((Nano->Rad))/3.;
      Nano->Height = 0.;
    }
    if(Nano->Shape == SHAPE_CYL || Nano->Shape == SHAPE_TILT)
      Volume -= PI*QUAD(Nano->Rad)*Nano->Height;
    if(Nano->Shape == SHAPE_CLUSTER){
      Volume -= PI*QUAD(Nano->Rad)*Nano->Height;
    }
    if(VAR_IF_TYPE(Soft[s].Topology,VAR_DISTRIBUTED)){
      Volume = Gen->Edge[CLat1]*Gen->Edge[CLat2]*Gen->Edge[CNorm];
      NChain = (int)( Volume*Soft[s].Size[2]/Gen->NPCh);
    }
    else if(VAR_IF_TYPE(Soft[s].Topology,VAR_PLANAR)){
      NChain = (int)(Area*ChainPArea);
    }
    else if(VAR_IF_TYPE(Soft[s].Topology,VAR_OBSTACLE)){
      Soft[s].NPCh = 2*Gen->NPCh;
      NChain = (int)(Area*ChainPArea*Soft[s].Size[2]);
      NChain = 200;
    }
    else if(VAR_IF_TYPE(Soft[s].Topology,VAR_COATING)){
      NChain = (int)( .53*Area*ChainPArea);    
    }
    else if(VAR_IF_TYPE(Soft[s].Topology,VAR_TUBE)){
      NChain = (int)( .69*Area*ChainPArea);    
    }
    else if(VAR_IF_TYPE(Soft[s].Topology,VAR_VESICLE)){
      ChainPArea = 3.0/Soft[s].Size[2];
      //Thickness = 3.5;
      NChain = (int)(Area*ChainPArea);
    }
    //------------DefBlock---------------
    Soft[s].NChain = NChain;
    Soft[s].NPCh = Gen->NPCh;
    if(VAR_IF_TYPE(Soft[s].Topology,VAR_PLANAR_PE)){
      Soft[s].NPCh = Gen->NPCh-1;
    }
    if(VAR_IF_TYPE(Soft[s].Topology,VAR_OBSTACLE))
      Soft[s].NPCh = 2*Gen->NPCh;
    Soft[s].NPart = Soft[s].NPCh*NChain;
    //Soft[s].NPart = NChain*Gen->NPCh;
    Gen->NChain += Soft[s].NChain;
    Gen->NPart += Soft[s].NPart;
  }
  Gen->NBlock = NSoft+AddNSoft;
  //----------Nano---------------
  for(int n=0;n<Gen->NNano;n++){
    if(Nano[n].Shape == SHAPE_CLUSTER){
      Gen->NBlock++;
      Gen->NChain += 1;
      Gen->NPart += Nano[n].NCircle*Nano[n].NHeight;
      if(Gen->NLink < 6) Gen->NLink = 6;
    }
  }
  if(NAddChain > 0) Gen->NBlock++;
  for(int s=0;s<NSoft;s++)
    if( VAR_IF_TYPE(Soft[s].Topology,VAR_PLANAR) )
      if(NAddChol  > 0) Gen->NBlock++;
  if(NStuffing > 0) Gen->NBlock++;
  if(NSolvent  > 0) Gen->NBlock++;
  //--------------Block---------------
  int DiblockLim = Block[0].Asym;
  Block = (BLOCK *)calloc(Gen->NBlock,sizeof(BLOCK));
  Soft[0].InitIdx = 0;
  for(int s=1;s<NSoft;s++){
    Soft[s].InitIdx += Soft[s-1].NPart + Soft[s-1].InitIdx;
  }
  for(int s=0,s1=0;s<NSoft;s++,s1++){
    Soft[s].EndIdx   = Soft[s].InitIdx + Soft[s].NPart;
    Block[s1].Asym = DiblockLim;
    Block[s1].InitIdx = Soft[s].InitIdx;
    Block[s1].NChain  = Soft[s].NChain;
    Block[s1].EndIdx  = Soft[s].EndIdx;
    Block[s1].NPCh = Soft[s].NPCh;
    sprintf(Block[s1].Name,"LIPID%d",s);
    if(VAR_IF_TYPE(SysType,VAR_TWOTAILS))
      sprintf(Block[s1].Name,"TT%d",s);
    if(VAR_IF_TYPE(Soft[s].Topology,VAR_VESICLE)){
      int NLayerIn = (int)(Soft[s].NChain*SQR(Soft[s].Size[0])/(SQR(Soft[s].Size[0])+SQR(Soft[s].Size[0]+Thickness)));
      int NLayerOut = Soft[s].NChain - NLayerIn;
      Block[s1].NChain = NLayerIn;
      Block[s1].EndIdx = Block[s1].InitIdx+Soft[s].NPCh*NLayerIn;
      Block[s1+1].InitIdx = Block[s1].EndIdx;
      Block[s1+1].EndIdx  = Soft[s].EndIdx;
      Block[s1+1].NChain = NLayerOut;
      Block[s1+1].NPCh = Soft[s].NPCh;
      sprintf(Block[s1].Name,"INNER%d",s);
      sprintf(Block[s1+1].Name,"OUTER%d",s);
      s1++;
    }
  }
  for(int n=0,b=NSoft+AddNSoft;n<Gen->NNano;n++){
    if(Nano[n].Shape == SHAPE_CLUSTER){
      Block[b].InitIdx = b != 0 ? Block[b-1].EndIdx : 0 ;
      Block[b].NChain  = 1;
      Block[b].EndIdx  = Block[b].InitIdx + Nano[n].NCircle*Nano[n].NHeight;
      Block[b].NPCh = Nano[n].NCircle*Nano[n].NHeight;
      Block[b].NPart = Block[b].NPCh;
      sprintf(Block[b].Name,"PEP%d",n);
      Nano[n].nBlock = b;
      b++;
    }
  }
  //-------------Alloc---------------
  Pm = (PART *)calloc(Gen->NPart,sizeof(PART));
  Ln = (LINKS *)calloc(Gen->NPart,sizeof(LINKS));
  if(Pm == NULL){printf("Non s'alloca\n"); return 1;}
  for(int p=0;p<Gen->NPart;p++){
    Ln[p].Link = (int *)calloc(Gen->NLink,sizeof(int));
    if(Ln[p].Link == NULL){printf("Non s'alloca\n"); return 1;}
  }
  SysInfo(cSystem);
  printf("%s\n",cSystem);
  SysDef(cSystem);
  printf("%s",cSystem);
  if(!HeaderSoft(cSystem))
    printf("%s",cSystem);
  /*---------Creating--------*/
  char ArchFile[60];
  int NMonCluster = 0;
  for(int s=0,b=0;s<NSoft;s++,b++){
    CreateSoft(arch,Thickness,s);
    if(VAR_IF_TYPE(Soft[s].Topology,VAR_VESICLE))
      b++;
    NMonCluster = Block[b].EndIdx;
  }
  for(int n=0;n<Gen->NNano;n++){
    if(Nano[n].Shape == SHAPE_CLINKS){
      FindNeighbours("CrossLinks.dat");
    }
    if(Nano[n].Shape != SHAPE_CLUSTER) continue;
    sprintf(ArchFile,"Architecture%d.dat",n);
    CreateProtein(n,NMonCluster);
    NMonCluster += Nano[n].NCircle*Nano[n].NHeight;
  }
  SetCoeff();
  double Norm2 = CUBE(pReOverCutOff()) / SQR(pNPCh());
  double Norm3 = CUBE(SQR(pReOverCutOff()) / pNPCh());
  MInt->Rescale(Norm2,2);
  MInt->Rescale(Norm3,3);
  Write(nome2);
  AddStuffing(nome2,NStuffing,0);
  AddChains(nome2,Thickness);
  AddSolvent(nome2,NSolvent);
  for(int s=0;s<NSoft;s++){
    if( VAR_IF_TYPE(Soft[s].Topology,VAR_PLANAR) ){
      AddCholesterol(nome2,Thickness,s);
    }
  }
  return 0;
}
bool VarData::CreateSoft(int *arch,double Thickness,int s){
  if( VAR_IF_TYPE(Soft[s].Topology,VAR_COATING) )
    CreateCoating(arch,Thickness,s);
  else if( VAR_IF_TYPE(Soft[s].Topology,VAR_VESICLE) )
    CreateVesicle(arch,Thickness,s);
  else if( VAR_IF_TYPE(Soft[s].Topology,VAR_TUBE) )
    CreateTube(arch,Thickness,s);
  else if( VAR_IF_TYPE(Soft[s].Topology,VAR_PLANAR) )
    CreatePlanar(arch,Thickness,s);
  else if( VAR_IF_TYPE(Soft[s].Topology,VAR_PLANAR_PE) )
    CreatePlanar(arch,Thickness,s);
  else if( VAR_IF_TYPE(Soft[s].Topology,VAR_OBSTACLE) )
    CreateObstacle(arch,Thickness,s);
  else if( VAR_IF_TYPE(Soft[s].Topology,VAR_DISTRIBUTED) ){
    double sigma = 1./sqrt(pkSpr());
    //sqrt(QUAD(Gen->ReOverCutOff)/(Gen->NPCh-1)/3.);
    for(int c=0,p=Soft[s].InitIdx;c<Soft[s].NChain;){
      for(int d=0;d<3;d++)
	Pm[p].Pos[d] = Mat->Casuale()*Gen->Edge[d];
      int IfContinue = 1;
      for(int i=1;i<Gen->NPCh;i++){
	for(int d=0;d<3;d++)
	  Pm[p+i].Pos[d] = Pm[p+i-1].Pos[d]+Mat->Gaussiano(0.,sigma);
	if(CheckNano(Pm[p+i].Pos,s)){IfContinue=0;break;}
      }
      if(IfContinue){
	c++;
	p += Gen->NPCh;
      }
    }
  }
  else {
    printf("System topology not recognized\n");
    return 1;
  }
  DefRest(arch,s);
  printf("Efficency: %d Tries for %d Particle\n",Tries,Soft[s].NPart);
  return 0;
}
int VarData::TrialSys(){
  int Values=32;
  double dValues = 1./(double)(Values);
  double **Plot = (double **)calloc(Values,sizeof(double));
  Gen->NPart = Values*Values;
  Gen->NChain = Gen->NPart;
  Gen->NPCh = 1;
  for(int d=0;d<3;d++){
    Gen->Edge[d] = (double)Values;
  }
  Pm = (PART *)calloc(Gen->NPart,sizeof(PART));
  Ch = (CHAIN *)calloc(Gen->NChain,sizeof(CHAIN));
  for(int i=0;i<Values;i++){
    Plot[i] = (double *)malloc(Values*sizeof(double));
    for(int j=0;j<Values;j++){
      Pm[i*Values+j].Pos[CLat1] = Ch[i*Values+j].Pos[CLat1] = (double) j;
      Pm[i*Values+j].Pos[CLat2] = Ch[i*Values+j].Pos[CLat2] = (double) i;
      Pm[i*Values+j].Pos[CNorm] = Ch[i*Values+j].Pos[CNorm] =
	Gen->Edge[CNorm]*.5 + cos(DUE_PI*(i)*dValues) + cos(DUE_PI*(i)*dValues*2.) + cos(DUE_PI*(i)*dValues*4.) +
	cos(DUE_PI*(j)*dValues) + cos(DUE_PI*(j)*dValues*2.) + cos(DUE_PI*(j)*dValues*4.);
    }
  }
  return 0;
}
void VarData::DefRest(int *arch,int s){
  double sigma = 1./sqrt(pkSpr());
  for(int c=0,p=Soft[s].InitIdx;c<Soft[s].NChain;c++){
    for(int i=0;i<Soft[s].NPCh;i++){
      for(int d=0;d <3;d++){
	Pm[p].Vel[d] = Mat->Gaussiano(0.,sigma) + Soft[s].Vel[d];
      }
      if(i == Soft[s].NPCh-1){
	Ln[p].Link[0] = p-1;
	Ln[p].NLink = 1;
      }
      else if(i == 0){
	Ln[p].Link[0] = p+1;
	Ln[p].NLink = 1;
      }
      else{
	Ln[p].NLink = 2;
	Ln[p].Link[0] = p-1;
	Ln[p].Link[1] = p+1;
      }
      if( VAR_IF_TYPE(Soft[s].Topology,VAR_ADDED) )
	Pm[p].Typ = 0;
      else 
	Pm[p].Typ = arch[i];
      Pm[p].CId = c;
      Pm[p].Idx = p;
      p++;
    }
  }
}
int VarData::PutPart(int j,int p,int HalfLim,double sigma){
  int i=0;
  if( j <= Block[0].Asym){
    i = Block[0].Asym - j;
    for(int d=0;d <3;d++){
      Pm[p+i].Pos[d] = Pm[p+1+i].Pos[d]+Mat->Gaussiano(0.,sigma);
    }
    if( VAR_IF_TYPE(SysType,VAR_TWOTAILS) ){
      if(i==HalfLim-1){
	for(int d=0;d <3;d++)
	  Pm[p+i].Pos[d] = Pm[p+Block[0].Asym].Pos[d]+Mat->Gaussiano(0.,sigma);
      }
    }
  }
  else if(j>Block[0].Asym){
    i = j;
    for(int d=0;d <3;d++){
      Pm[p+i].Pos[d] = Pm[p-1+i].Pos[d]+Mat->Gaussiano(0.,sigma);
    }
  }
  return i;
}
void VarData::CreateObstacle(int *arch,double Thickness,int s){
  double sigma = sqrt(QUAD(Gen->ReOverCutOff)/(Gen->NPCh-1)/3.);
  int HalfLim = (int)(Block[0].Asym*.5);
  double Dz = Thickness/16.;
  double Leaflet = -.5*Thickness - 2.*Dz;
  for(int c=0,p=Soft[s].InitIdx;c<Soft[s].NChain;c++,p+=Soft[s].NPCh){
    double Cas1 = Mat->Casuale(),Cas2 = Mat->Casuale();
    Pm[p].Pos[CLat1] = Cas1*Gen->Edge[CLat1]*.5;
    Pm[p].Pos[CLat2] = Cas2*Gen->Edge[CLat2]*.5;
    Pm[p].Pos[CNorm] = Soft[s].Pos[CNorm]+Leaflet;
    Pm[p].Typ = 1;
    for(int pn=1;pn<Soft[s].NPCh;pn++){
      Pm[p+pn].Pos[CLat1] = Pm[p].Pos[CLat1];
      Pm[p+pn].Pos[CLat2] = Pm[p].Pos[CLat2];
      Pm[p+pn].Pos[CNorm] = Pm[p+pn-1].Pos[CNorm] + Dz;
      Pm[p+pn].Typ = 0;
      if( (pn < 2) || (pn >= Soft[s].NPCh - 2)){
	Pm[p+pn].Typ = 1;
      }
    }
  }
}
void VarData::CreatePlanar(int *arch,double Thickness,int s){
  double sigma = sqrt(QUAD(Gen->ReOverCutOff)/(Gen->NPCh-1)/3.);
  int HalfLim = (int)(Block[0].Asym*.5);
  for(int c=0,p=Soft[s].InitIdx;c<Soft[s].NChain;){
    printf("%d %d %d %lf \r",p,c,Tries,p/(double)Soft[s].NPart);
    double Leaflet = -.5*Thickness;
    if(c >= Soft[s].NChain/2) Leaflet = .5*Thickness;
    //first part
    double Cas1 = Mat->Casuale(),Cas2 = Mat->Casuale();
    Pm[p+Block[0].Asym].Pos[CLat1] = Cas1*Gen->Edge[CLat1];
    Pm[p+Block[0].Asym].Pos[CLat2] = Cas2*Gen->Edge[CLat2];
    Pm[p+Block[0].Asym].Pos[CNorm] = Soft[s].Pos[CNorm]+Leaflet;
    int IfContinue = 1;
    //others
    for(int j=1;j<Soft[s].NPCh;j++){
      int i = PutPart(j,p,HalfLim,sigma);
      //absorbing boundary condition
      if(arch[i] == 0 ){
	if(Pm[p+i].Pos[CNorm] > Soft[s].Pos[CNorm]+.5*Thickness ||
	   Pm[p+i].Pos[CNorm] < Soft[s].Pos[CNorm]-.5*Thickness ){
	  Tries++;
	  IfContinue = 0;
	  break;
	}
      }
      else if (arch[i] == 1 ){
	if(Pm[p+i].Pos[CNorm] < Soft[s].Pos[CNorm]+.5*Thickness &&
	   Pm[p+i].Pos[CNorm] > Soft[s].Pos[CNorm]-.5*Thickness ){
	  Tries++;
	  IfContinue = 0;
	  break;
	}
      }
      if(CheckNano(Pm[p+i].Pos,s)){IfContinue=0;break;}
    }
    if(IfContinue){
      c++;
      p += Soft[s].NPCh;
    }
  }
}
void VarData::CreateVesicle(int *arch,double Thickness,int s){
  double sigma = sqrt(QUAD(Gen->ReOverCutOff)/(Gen->NPCh-1)/3.);
  int NLayerIn = (int)(Soft[s].NChain*SQR(Soft[s].Size[0])/(SQR(Soft[s].Size[0])+SQR(Soft[s].Size[0]+Thickness)) );
  int NLayerOut = Soft[s].NChain - NLayerIn;
  int HalfLim = (int)(Block[0].Asym*.5);
  double inc = M_PI * (3. - sqrt(5.));
  double NInv = 1. / (double)NLayerIn;
  double Leaflet = -.5*Thickness;
  for(int c=0,cc=0,p=Soft[s].InitIdx;c<Soft[s].NChain;){
    printf("%d %d %d %lf \r",p,c,Tries,p/(double)Soft[s].NPart);
    if(c == NLayerIn){
      cc=0;
      NInv = 1./(double)NLayerOut;
      Leaflet = .5*Thickness;
    }
    double x = cc*2.*NInv - 1. + (NInv);
    double r = sqrt(1.-x*x);
    double phi = cc*inc;
    double y = cos(phi)*r;
    double z = sin(phi)*r;
    //first part
    Pm[p+Block[0].Asym].Pos[CLat1] =
      (Soft[s].Size[0]+Leaflet)*x + Soft[s].Pos[CLat1];
    Pm[p+Block[0].Asym].Pos[CLat2] = 
      (Soft[s].Size[0]+Leaflet)*y + Soft[s].Pos[CLat2];
    Pm[p+Block[0].Asym].Pos[CNorm] = 
      (Soft[s].Size[0]+Leaflet)*z + Soft[s].Pos[CNorm];
    //others
    int IfContinue = 1;
    for(int j=1;j<Gen->NPCh;j++){
      int i = PutPart(j,p,HalfLim,sigma);
      // absorbing boundary condition
      double Dist = SQR(Pm[p+i].Pos[CLat1]-Soft[s].Pos[CLat1]) + SQR(Pm[p+i].Pos[CLat2]-Soft[s].Pos[CLat2]) + SQR(Pm[p+i].Pos[CNorm]-Soft[s].Pos[CNorm]);
      if(arch[i] == 0){
	if( Dist < SQR(Soft[s].Size[0]-.5*Thickness) || Dist > SQR(Soft[s].Size[0]+.5*Thickness) ) {
	  Tries++;
	  IfContinue = 0;
	  break;
	}
      }
      else if(arch[i] == 1){
	if(Dist < SQR(Soft[s].Size[0]+.5*Thickness) && Dist > SQR(Soft[s].Size[0] - .5*Thickness) ){
	  Tries++;
	  IfContinue = 0;
	  break;
	}
      }
      if(CheckNano(Pm[p+i].Pos,s)){IfContinue=0;break;}
    }
    if(IfContinue){
      c++;
      cc++;
      p += Gen->NPCh;
    }
  }
  // for(int p=0;p<pNPart();p++){
  //   pPos(p);
  // }
}
void VarData::CreateCoating(int *arch,double Thickness,int s){
  double sigma = sqrt(QUAD(Gen->ReOverCutOff)/(Gen->NPCh-1)/3.);
  int HalfLim = (int)(Block[0].Asym*.5);
  double NInv = 1./(double)Soft[s].NChain;
  double inc = 3.141592654 * (3. - sqrt(5.));
  double Dist = 0.;
  for(int c=0,p=Soft[s].InitIdx;c<Soft[s].NChain;){
    printf("%d %d %d %lf \r",p,c,Tries,p/(double)Soft[s].NPart);
    //first part
    double Cas1 = Mat->Casuale()*DUE_PI;
    if(VAR_IF_TYPE(Nano->Shape,SHAPE_HEI)){
      Pm[p+Block[0].Asym].Pos[CLat1] = 
	(Nano->Rad+.5*Thickness)*cos(Cas1) + Nano->Pos[CLat1];
      Pm[p+Block[0].Asym].Pos[CLat2] = 
	(Nano->Rad+.5*Thickness)*sin(Cas1) + Nano->Pos[CLat2];
      Pm[p+Block[0].Asym].Pos[CNorm] = 
	(.5 - Mat->Casuale())*Nano->Height + Nano->Pos[CNorm];
    }
    else if(VAR_IF_TYPE(Nano->Shape,SHAPE_SPH)){
      double x = c*2.*NInv - 1. + (NInv);
      double r = sqrt(1.-x*x);
      double y = cos(c*inc)*r;
      double z = sin(c*inc)*r;
      Pm[p+Block[0].Asym].Pos[CLat1] =
	(Nano->Rad+.5*Thickness)*x + Nano->Pos[CLat1];
      Pm[p+Block[0].Asym].Pos[CLat2] = 
	(Nano->Rad+.5*Thickness)*y + Nano->Pos[CLat2];
      Pm[p+Block[0].Asym].Pos[CNorm] = 
	(Nano->Rad+.5*Thickness)*z + Nano->Pos[CNorm];
    }
    int IfContinue = 1;
    //others
    for(int j=1;j<Gen->NPCh;j++){
      int i = PutPart(j,p,HalfLim,sigma);
      // absorbing boundary condition
      if(VAR_IF_TYPE(Nano->Shape,SHAPE_HEI))
	Dist = SQR(Pm[p+i].Pos[CLat1]-Nano->Pos[CLat1]) + SQR(Pm[p+i].Pos[CLat2]-Nano->Pos[CLat2]);
      else if(VAR_IF_TYPE(Nano->Shape,SHAPE_SPH))
	Dist = SQR(Pm[p+i].Pos[CLat1]-Nano->Pos[CLat1]) + SQR(Pm[p+i].Pos[CLat2]-Nano->Pos[CLat2]) + SQR(Pm[p+i].Pos[CNorm]-Nano->Pos[CNorm]);
      if(arch[i] == 0 && Dist > SQR(Nano->Rad+.5*Thickness) ){
	Tries++;
	IfContinue = 0;
	break;
      }
      else if(arch[i] == 1 && Dist < SQR(Nano->Rad+.5*Thickness) ){
	Tries++;
	IfContinue = 0;
	break;
      }
      if(CheckNano(Pm[p+i].Pos,s)){IfContinue=0;break;}
    }
    if(IfContinue){
      c++;
      p += Gen->NPCh;
    }
  }
}
void VarData::CreateTube(int *arch,double Thickness,int s){
  double sigma = sqrt(QUAD(Gen->ReOverCutOff)/(Gen->NPCh-1)/3.);
  int HalfLim = (int)(Block[0].Asym*.5);
  int NLayerIn = (int)(Soft[s].NChain/2.*SQR(Soft[s].Size[0])/(SQR(Soft[s].Size[0])+SQR(Soft[s].Size[0]+Thickness)) );
  for(int c=0,p=Soft[s].InitIdx;c<Soft[s].NChain;){
    printf("%d %d %d %lf \r",p,c,Tries,p/(double)Soft[s].NPart);
    double Leaflet = -.5*Thickness;
    if(c >= NLayerIn) Leaflet = .5*Thickness;
    //first part
    double Cas1 = DUE_PI*Mat->Casuale();
    //FIXME: this distribution is somehow not uniform
    Pm[p+Block[0].Asym].Pos[CLat1] = 
      (Soft[s].Size[0]+Leaflet)*cos(Cas1) + Soft[s].Pos[CLat1];
    Pm[p+Block[0].Asym].Pos[CLat2] = 
      (Soft[s].Size[0]+Leaflet)*sin(Cas1) + Soft[s].Pos[CLat2];
    Pm[p+Block[0].Asym].Pos[CNorm] = (Mat->Casuale() - .5)*Soft[s].Size[1] + Soft[s].Pos[CNorm];

    //others
    int IfContinue = 1;
    for(int j=1;j<Gen->NPCh;j++){
      int i = PutPart(j,p,HalfLim,sigma);
      // absorbing boundary condition
      double Dist = SQR(Pm[p+i].Pos[CLat1]-Soft[s].Pos[CLat1]) + SQR(Pm[p+i].Pos[CLat2]-Soft[s].Pos[CLat2]);
      if(arch[i] == 0)
	if( Dist < SQR(Soft[s].Size[0]-.5*Thickness) || Dist > SQR(Soft[s].Size[0]+.5*Thickness) ) {
	  Tries++;
	  IfContinue = 0;
	  break;
	}
	else if(arch[i] == 1)
	  if(Dist > SQR(Soft[s].Size[0]-.5*Thickness) || Dist < SQR(Soft[s].Size[0] + .5*Thickness) ){
	    Tries++;
	    IfContinue = 0;
	    break;
	  }
      if(CheckNano(Pm[p+i].Pos,s)){IfContinue=0;break;}
    }
    if(IfContinue){
      c++;
      p += Gen->NPCh;
    }
  }
}
int VarData::CheckNano(double *Pos,int s){
  for(int n=0;n<pNNano();n++){
    Point2Shape(Nano[n].Shape);
    double Add = .0;
    double Radius2 = 0.;
    if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_BOUND)) Add = .3;
    if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_NONE)) continue;
    else{
      Radius2 = NanoDist2(Pos,n);
    }
    // else{
    //   Vettore PosRel(3);
    //   Vettore NanoAxis(3);
    //   Vettore Distance(3);
    //   for(int d=0;d<3;d++){
    // 	PosRel.Set(Nano[n].Pos[d] - Pos[d],d);
    // 	NanoAxis.Set(Nano[n].Axis[d],d);
    //   }
    //   Distance.VetV(&NanoAxis,&PosRel);
    //   Radius2 = SQR(Distance.Norm());
    //   double HeiOnAxis = PosRel.ProjOnAxis(&NanoAxis);
    //   if( fabs(HeiOnAxis) > Nano[n].Height*.5)
    // 	Radius2 = QUAD(Nano[n].Rad+Add);
    // }
    if(Radius2 < QUAD(Nano[n].Rad+Add)){
      Tries++;
      return 1;
    }
  }
  return 0;
}
//Obsolete
void VarData::AddProtein(int NCircle,int NHeight,int nNano,char *filename){
  int NPart = NCircle*NHeight;
  PART *Pn = (PART *)calloc(NPart,sizeof(PART));
  double CirInv = 1./(double)NCircle;
  double HeiInv = Nano[nNano].Height/(double)NHeight;
  for(int c=0;c<NCircle;c++){
    double Sin = sin(c*CirInv*DUE_PI);
    double Cos = cos(c*CirInv*DUE_PI);
    double x = Nano[nNano].Rad * Cos + Nano[nNano].Pos[0];
    double y = Nano[nNano].Rad * Sin + Nano[nNano].Pos[1];
    for(int h=0;h<NHeight;h++){
      int p = c*NHeight + h;
      double z = h*HeiInv + Nano[nNano].Pos[2] - .5*Nano[nNano].Height;
      Pn[p].Pos[CLat1] = x;
      Pn[p].Pos[CLat2] = y;
      Pn[p].Pos[CNorm] = z;
    }
  }
  FILE *ReSave = fopen(filename,"a");
  fprintf(ReSave,"# n=1 N=%d name=PEP1\n",NPart);
  for(int p=0;p<NPart;p++)
    fprintf(ReSave,"%lf %lf %lf %lf %lf %lf %d\n",
	    Pn[p].Pos[0],Pn[p].Pos[1],Pn[p].Pos[2],
	    0.0,0.0,0.0,0);
  Nano[nNano].OffSet = (double)NCircle;
}
void VarData::CreateProtein(int nNano,int np){
  int NCircle = Nano[nNano].NCircle;
  int NHeight = Nano[nNano].NHeight;
  int NCyl = NCircle*NHeight;
  int NSph = (NCircle)*(NCircle/2-1) + 2;
  int NTot = NCyl;// + NSph;
  int IfDoubleSided = 0;
  int IfKink = 0;
  double AsymPhil = -.1*Nano[nNano].Rad;
  double CirInv = 1./(double)NCircle;
  double HeiInv = Nano[nNano].Height/(double)(NHeight-1);
  double HeiSph = Nano[nNano].Height;
  if(Nano[nNano].NHeight == 0 ) HeiSph + HeiInv;
  int pWrote = 0;
  double Shift = .5*sin(1*CirInv*DUE_PI);
  double SegSide = Nano[nNano].Height/(double)Nano[nNano].NHeight;
  double SegCirc = DUE_PI*Nano[nNano].Rad/(double)Nano[nNano].NCircle;
  // Part positions
  // Cylinder
  for(int c=0;c<NCircle;c++){
    for(int h=0;h<NHeight;h++){
      int NLink = 0;
      int p = c*NHeight + h;
      //pos current particle
      double Sin = sin(c*CirInv*DUE_PI);
      double Cos = cos(c*CirInv*DUE_PI);
      double Rad = Nano[nNano].Rad;
      double Weight = 1.;
      double Axes[3] = {pNanoPos(nNano,0),pNanoPos(nNano,1),pNanoPos(nNano,2)};
      if(IfKink){
	Axes[0] += .25*Nano->Height*(fabs((double)(h-NHeight/2))/(double)(NHeight));
	// Weight -= .1*Cos;
	// Weight -= .8*(1.-fabs((double)(h-NHeight/2))/(double)(NHeight));
      }
      if( ((h)%2)==0 ){
	Sin = sin((c+.5)*CirInv*DUE_PI);
	Cos = cos((c+.5)*CirInv*DUE_PI);
	Rad = Nano[nNano].Rad - .1;
      }
      Ln[p+np].NLink = 0;
      Pm[p+np].Pos[0] = Rad * Cos + Axes[0];
      Pm[p+np].Pos[1] = Rad * Sin + Axes[1];
      Pm[p+np].Pos[2] = h*HeiInv*Weight - .5*Nano[nNano].Height*Weight + Axes[2];
      Pm[p+np].Typ = 0;
      if(IfDoubleSided){
	if(Pm[p+np].Pos[0] < (Nano[nNano].Pos[0] + AsymPhil))
	  Pm[p+np].Typ = 1;
      }
      if(h < 2 || h > NHeight - 3) Pm[p+np].Typ = 1;
      pWrote++;
      //-------------Connections------------------
      // Right
      int pp = (c+1)*NHeight + h;
      if( pp >= NCyl ) pp -= NCyl;
      Ln[p+np].Link[NLink++] = pp + np;
      pp = p + 1;
      if( (pp%NHeight)!=0 ) Ln[p+np].Link[NLink++] = pp + np;
      //Left up
      pp = (c-1)*NHeight + h + 1;
      if(c==0) pp = (NCircle-1)*NHeight + h + 1;
      if( ((h)%2)==0 ){
	pp = (c+1)*NHeight + h + 1;
	if(c==NCircle-1) pp = 0 + h + 1;
      }
      if( (pp%Nano[nNano].NHeight) ) Ln[p+np].Link[NLink++] = pp + np;    
      // Side
      if( h == 0 ) {
	pp = (c)*NHeight + NHeight - 2;
	//if(pp < NCyl) Ln[p+np].Link[NLink++] = pp + np;
	// Diagonal
	pp = (c+NCircle/2)*NHeight + NHeight-2;//p + NCyl/2 + NCircle - 1;
	if(pp > NCyl) pp = (c-NCircle/2)*NHeight + NHeight-2;
	//if(pp < NCyl) Ln[p+np].Link[NLink++] = pp + np;
      }
      // Disks
      pp = p + NCircle/2*NHeight;
      //if(pp < NCyl-1) Ln[p+np].Link[NLink++] = pp + np;
      Ln[p+np].NLink = NLink;
    }
  }
  Vettore Axis(0.,0.,1.);
  Vettore Axis1(Nano[nNano].Axis[0],Nano[nNano].Axis[1],Nano[nNano].Axis[2]);
  Vettore Axis2(3);
  Axis2 = Axis1 + Axis;
  for(int d=0;d<3;d++){
    Axis2.Set(Axis1.Val(d)+Axis.Val(d),d);
  }
  Vettore Origin(Nano[nNano].Pos[0],Nano[nNano].Pos[1],Nano[nNano].Pos[2]);
  int b = nNano + NSoft;
  RotateBlock(&Axis2,&Origin,b);
  char Filename[60];
  sprintf(Filename,"Architecture%d.dat",nNano);
  FILE *CSave = fopen(Filename,"w");
  fprintf(CSave,"# Cylinder\n");
  //write the initial mutual distances
  for(int p=np;p<NTot+np;p++){
    for(int l=0;l<Ln[p].NLink;l++){
      int l2 = Ln[p].Link[l];
      double Dist = sqrt( SQR(Pm[p].Pos[0]-Pm[l2].Pos[0]) + SQR(Pm[p].Pos[1]-Pm[l2].Pos[1]) + SQR(Pm[p].Pos[2]-Pm[l2].Pos[2]) );
      double kSpr = 10000.;
      if(Dist > 2.) kSpr = 10000.;
      if(Dist > 4.) kSpr = 10000.;
      fprintf(CSave,"%d %d %lf %.0f\n",p-np,l2-np,Dist,kSpr);
    }
  }
  fclose(CSave);return;
  // FIXME: the links are not closing at the boundaries
  // Cupola
  fprintf(CSave,"# Cupola\n");
  int PType = 2;
  int NCircHalf = NCircle/2;
  // Vertical
  for(int cc=1;cc<NCircHalf-1;cc++){
    double Sin2 = sin(cc*CirInv*DUE_PI);
    double Cos2 = cos(cc*CirInv*DUE_PI);
    // Horizontal
    for(int c=0;c<NCircle;c++){
      double Sin = sin(c*CirInv*DUE_PI);
      double Cos = cos(c*CirInv*DUE_PI);
      if( (cc%2)==0 ){
	Sin = sin((c+.5)*CirInv*DUE_PI);
	Cos = cos((c+.5)*CirInv*DUE_PI);
      }
      double x = Nano[nNano].Rad * Cos * Sin2 + Nano[nNano].Pos[0];
      double y = Nano[nNano].Rad * Sin * Sin2 + Nano[nNano].Pos[1];
      double Quota = Nano[nNano].Height*.5;
      if(cc > NCircle/4) Quota = - Nano[nNano].Height*.5;
      double z = Nano[nNano].Rad * Cos2 + Nano[nNano].Pos[2] + Quota;
      pWrote++;
      //fprintf(PSave,"%lf %lf %lf %lf %lf %lf %d\n",x,y,z,0.,0.,0.,PType);
      // Connections
      double Dx = x - Nano[nNano].Rad*cos((c+1)*CirInv*DUE_PI)*Sin2 - Nano[nNano].Pos[0];
      double Dy = y - Nano[nNano].Rad*sin((c+1)*CirInv*DUE_PI)*Sin2 - Nano[nNano].Pos[1];
      double Elong = sqrt(SQR(Dx)+SQR(Dy));
      int p = NCyl + c + cc*NCircle;
      int pp = NCyl + (c+1) + cc*NCircle;
      if( c+1==NCircle ) pp = NCyl + (0) + cc*NCircle;
      if(c != 0 && c != NCircle-1) fprintf(CSave,"%d %d %lf\n",p,pp,Elong);
      if( cc != NCircle/2-1 && cc != NCircle/4){
	pp = NCyl + (c+1) + (cc+1)*NCircle;
	Elong = sqrt( SQR(Nano[nNano].Rad*1.*CirInv*DUE_PI) + SQR(.5*Elong));
	fprintf(CSave,"%d %d %lf\n",p,pp,Elong);
	pp = NCyl + (c) + (cc-1)*NCircle;
	if( (cc%2)==0 )
	  pp = NCyl + (c+1) + (cc)*NCircle;
	fprintf(CSave,"%d %d %lf\n",p,pp,Elong);
      }
    }
  }
  // Closing the extrema
  fprintf(CSave,"# Extrema\n");
  {
    double x = Nano[nNano].Pos[0];
    double y = Nano[nNano].Pos[1];
    double z = Nano[nNano].Pos[2] - HeiSph*.5 - Nano[nNano].Rad;
    pWrote++;
    //fprintf(PSave,"%lf %lf %lf %lf %lf %lf %d\n",x,y,z,0.,0.,0.,PType);
    z = Nano[nNano].Pos[2] + HeiSph*.5 + Nano[nNano].Rad;
    pWrote++;
    //fprintf(PSave,"%lf %lf %lf %lf %lf %lf %d\n",x,y,z,0.,0.,0.,PType);
    {
      int p = pWrote - 1;
      int pp = pWrote - 2;
      double Elong = Nano[nNano].Height + 2.*Nano[nNano].Rad;
      fprintf(CSave,"%d %d %lf\n",p,pp,Elong);
    }
    for(int c=0;c<NCircle;c++){
      int p = pWrote-2;
      int pp = c + NCyl + NCircle*(NCircle/2-2);
      double Elong = Nano[nNano].Rad * 1.*CirInv*DUE_PI;
      fprintf(CSave,"%d %d %lf\n",p,pp,Elong);
      p = pWrote-1;
      pp = c + NCyl;
      fprintf(CSave,"%d %d %lf\n",p,pp,Elong);
      // Conntecting to the cylinder
      p = c*NHeight;
      pp = c + NCyl + NCircle*(NCircle/4);
      Elong = HeiInv;//Actually it is a bit more
      fprintf(CSave,"%d %d %lf\n",p,pp,Elong);
      p = c*(NHeight) + NHeight - 1;
      pp = c + NCyl + NCircle*(NCircle/4-1);
      fprintf(CSave,"%d %d %lf\n",p,pp,Elong);
    }
  }
  fclose(CSave);
  //fclose(PSave);
}
void VarData::AddStuffing(char *filename,int NWater,int nNano){
  if(NWater == 0) return;
  FILE *PSave = fopen(filename,"a");
  fprintf(PSave,"# n=%d N=10 name=STUFFING\n",NWater/10);
  for(int p=0;p<NWater;p++){
    double Angle = Mat->Casuale()*DUE_PI;
    double Rad = 2.*(Mat->Casuale()-.5)*Nano[nNano].Rad;
    double x = cos(Angle)*Rad + Nano[nNano].Pos[0];
    double y = sin(Angle)*Rad + Nano[nNano].Pos[1];
    double z = (Mat->Casuale()-.5)*Nano[nNano].Height + Nano[nNano].Pos[2];
    fprintf(PSave,"%lf %lf %lf %lf %lf %lf %d\n",x,y,z,0.,0.,0.,1);
  }
  fclose(PSave);
}
void VarData::AddSolvent(char *filename,int NWater){
  if(NSolvent <= 0) return;
  FILE *PSave = fopen(filename,"a");
  double sigma = 1./sqrt(pkSpr());
  fprintf(PSave,"# n=%d N=1 name=SOLVENT\n",NWater);
  for(int p=0;p<NWater;p++){
    double x = Mat->Casuale()*Gen->Edge[CLat1];
    double y = Mat->Casuale()*Gen->Edge[CLat2];
    double z = Mat->Casuale()*Gen->Edge[CNorm]*.3;
    fprintf(PSave,"%lf %lf %lf %lf %lf %lf %d\n",x,y,z,Mat->Gaussiano(0.,sigma),Mat->Gaussiano(0.,sigma),Mat->Gaussiano(0.,sigma),2);
  }
  fclose(PSave);
}
void VarData::AddChains(char *filename,double Thickness){
  if(NAddChain <= 0) return;
  double sigma = 1./sqrt(pkSpr());
  PART *Pn = (PART *)calloc(Gen->NPCh,sizeof(PART));
  FILE *PSave = fopen(filename,"a");
  int NPCh = (int)(.5*pNPCh());
  fprintf(PSave,"# n=%d N=%d name=ADDED\n",NAddChain,NPCh);
  int s = 0;
  for(int c=0;c<NAddChain;){
    Pn[0].Pos[CLat1] = Mat->Casuale()*Gen->Edge[CLat1];
    Pn[0].Pos[CLat2] = Mat->Casuale()*Gen->Edge[CLat2];
    Pn[0].Pos[CNorm] = .2*(Mat->Casuale()-.5)*Thickness + Soft[s].Pos[CNorm];
    int IfContinue = !CheckNano(Pn[0].Pos,s);
    for(int i=1;i<NPCh;i++){
      for(int d=0;d<3;d++)
	Pn[i].Pos[d] = Pn[i-1].Pos[d]+Mat->Gaussiano(0.,sigma);
      if(Pn[i].Pos[CNorm] > Soft[s].Pos[CNorm]+.5*Thickness ||
	 Pn[i].Pos[CNorm] < Soft[s].Pos[CNorm]-.5*Thickness ){
	Tries++;
	IfContinue = 0;
	break;
      }
      //if(CheckNano(Pn[i].Pos,0)){IfContinue=0;break;}
    }
    if(IfContinue){
      for(int i=0;i<NPCh;i++){
	fprintf(PSave,"%lf %lf %lf %lf %lf %lf %d\n",Pn[i].Pos[0],Pn[i].Pos[1],Pn[i].Pos[2],Mat->Gaussiano(0.,sigma),Mat->Gaussiano(0.,sigma),Mat->Gaussiano(0.,sigma),0);
      }
      c++;
      s++;if(s==NSoft)s=0;
    }
  }
  fclose(PSave);
  free(Pn);
}
void VarData::AddCholesterol(char *filename,double Thickness,int s){
  if(NAddChol <= 0) return;
  double sigma = 1./sqrt(pkSpr());
  int HalfLim = (int)(Block[0].Asym*.5);
  int NPCh = 5;
  int DLim = NPCh-1;
  PART *Pn = (PART *)calloc(Gen->NPCh,sizeof(PART));
  FILE *PSave = fopen(filename,"a");
  fprintf(PSave,"# n=%d N=%d name=CHOL%d\n",NAddChol,NPCh,s);
  for(int c=0;c<NAddChol;){
    printf("%d %d %lf \r",c,Tries,c/(double)(NAddChol));
    double Leaflet = -.5*Thickness;
    if(c >= NAddChol/2) Leaflet = +.5*Thickness;
    //first part
    double Cas1 = Mat->Casuale(),Cas2 = Mat->Casuale();
    Pn[DLim].Pos[CLat1] = Cas1*Gen->Edge[CLat1];
    Pn[DLim].Pos[CLat2] = Cas2*Gen->Edge[CLat2];
    Pn[DLim].Pos[CNorm] = Soft[s].Pos[CNorm]+Leaflet;
    Pn[DLim].Typ = 1;
    int IfContinue = !CheckNano(Pn[DLim].Pos,s);
    //others
    for(int j=DLim-1;j>=0;j--){
      Pn[j].Typ = 0;
      for(int d=0;d<3;d++)
	Pn[j].Pos[d] = Pn[j+1].Pos[d]+Mat->Gaussiano(0.,sigma);
      if(Pn[j].Pos[CNorm] > Soft[s].Pos[CNorm]+.5*Thickness ||
	 Pn[j].Pos[CNorm] < Soft[s].Pos[CNorm]-.5*Thickness ){
	Tries++;
	IfContinue = 0;
	break;
      }
      if(CheckNano(Pn[j].Pos,s)){IfContinue=0;break;}
    }
    if(IfContinue){
      for(int i=0;i<NPCh;i++){
	fprintf(PSave,"%lf %lf %lf %lf %lf %lf %d\n",Pn[i].Pos[0],Pn[i].Pos[1],Pn[i].Pos[2],Mat->Gaussiano(0.,sigma),Mat->Gaussiano(0.,sigma),Mat->Gaussiano(0.,sigma),Pn[i].Typ);
      }
      c++;
    }
  }
  fclose(PSave);
}
#include <Cubo.h>
typedef DdLinkedList Cubo;
void VarData::FindNeighbours(char *FileName){
  double CutOff = 1.0;
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  int NeiList[27];
  int NPair = (int)(pNChain()/10.);
  int NChOffSet = 0;
  int bInner = 0;
  for(int b=0,NCh=0;b<pNBlock();NCh+=Block[b++].NChain){
    if(strcmp(Block[b].Name,"INNER0")) continue;
    NChOffSet = NCh;
    bInner = b;
  }
  int NChain = 0;
  for(int b=0;b<pNBlock();b++){
    NChain += Block[b].NChain;
  }
  SetNChain(NChain);
  for(int c=0;c<pNChain();c++){
    for(int d=0;d<3;d++){
      Ch[c].Pos[d] = Pm[c*pNPCh()+pNPCh()-1].Pos[d];
    }
    // for(int p=c*pNPCh();p<(c+1)*pNPCh();p++){
    //   for(int d=0;d<3;d++){
    // 	Ch[c].Pos[d] += Pm[p].Pos[d];
    //   }
    // }
    // for(int d=0;d<3;d++){
    //   Ch[c].Pos[d] /= (double)pNPCh();
    // }
  }
  double *cPair = (double *)calloc(3*pNChain(),sizeof(double));
  for(int c=0;c<pNChain();c++){
    cPair[c*3+0] = c;
    cPair[c*3+1] = c;
    cPair[c*3+2] = 1000.;
  }
  Cubo *Pc = new Cubo(Edge,pNChain(),CutOff);
  for(int c=0;c<pNChain();c++){
    int p = c*pNPCh()+pNPCh()-1;
    //Pc->AddPart(c,Pm[p].Pos);
    Pc->AddPart(c,Ch[c].Pos);
  }
  FILE *FWrite = fopen(FileName,"w");
  for(int c=NChOffSet;c<NChOffSet+Block[bInner].NChain-1;c++){
    cPair[c*3+0] = (double)c;
    int p1 = c*pNPCh(bInner)+pNPCh(bInner)-1;
    double MinDist = 1000.;
    int NNei = Pc->GetNei(Pm[p1].Pos,NeiList);
    for(int i=0;i<NNei;i++){
      int c1 = NeiList[i];
      for(Pc->SetCounters(c1);Pc->IfItCell(c1);Pc->IncrCurr(c1)){
	int c2 = Pc->ItCell(c1);
	if(c2 <= c) continue;
	int p2 = c2*pNPCh(bInner)+pNPCh(bInner)-1;
	double Dist2 = 0.;
	for(int d=0;d<3;d++){
	  Dist2 += SQR(Ch[c].Pos[d] - Ch[c2].Pos[d]);
	  //Dist2 += SQR(Pm[p1].Pos[d] - Pm[p2].Pos[2]);
	}
	if(MinDist > Dist2){
	  MinDist = Dist2;
	  cPair[c*3+1] = (double)c2;
	  cPair[c*3+2] = MinDist;
	}
      }
    }
  }
  double Temp[3];
  for(int c=1;c<pNChain();c++){
    for(int c1=c;c1>=0;c1--){
      if(cPair[c1*3+2] >= cPair[(c1-1)*3+2]) break;
      for(int d=0;d<3;d++){
      	Temp[d] = cPair[c1*3+d];
      	cPair[c1*3+d] = cPair[(c1-1)*3+d];
      	cPair[(c1-1)*3+d] = Temp[d];
      }
      // Mat->Swap(cPair,c1*3+2,cPair,(c1-1)*3+2,3);
      //printf("%d) %d %lf %lf\n",c,c1,cPair[(c1)*3+2],cPair[(c1-1)*3+2]);
    }
  }
  int pRef = 0;//NChOffSet*pNPCh(bInner);
  double KEl = 1.;
  for(int c=0;c<NPair;c++){
    int p1 = (int)cPair[c*3+0]*pNPCh(bInner)+pNPCh(bInner)-1;
    int p2 = (int)cPair[c*3+1]*pNPCh(bInner)+pNPCh(bInner)-1;
    double Dist = 0.;
    for(int d=0;d<3;d++){
      Dist += SQR(Ch[(int)cPair[c*3+0]].Pos[d] - Ch[(int)cPair[c*3+1]].Pos[d]);
      //Dist += SQR(Pm[p1].Pos[d] - Pm[p2].Pos[d]);
    }
    Dist = sqrt(Dist);
    //Dist = 0.;
    fprintf(FWrite,"%d %d %lf %lf\n",p1-pRef,p2-pRef,Dist,KEl);
  }
  fclose(FWrite);
}
