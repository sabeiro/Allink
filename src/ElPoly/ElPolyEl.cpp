/***********************************************************************
ElPoly:This progam provide a graphical visualisation of the data 
opend by VarData using openGL glut. The most important option are 
the possibility of changing the backfold of the polymers with 'c', 
see the subsequent file in the list with '>', see the bond with 'b'. 
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
#include "ElPoly.h"
;

ElPoly::ElPoly(int argc,char **argv,int *FilePos){
  cFile = (char **)calloc(argc,4);//sizeof(**cFile));
  if(cFile == NULL){printf("cFile not allocated\n");return;}
  for(int f=0;f<argc;f++){
    cFile[f] = (char *)calloc(STRSIZE,sizeof(*cFile));
    if(cFile[f] == NULL){printf("cFile not allocated\n");return;}
    cFile[f] = argv[FilePos[f]];
  }
  NFileTot = argc;
  NFile[0] = 0; NFile[1] = NFileTot;
  NPro =0;//?
  quando=0;
  if(!pNLink()) IfLine = 1;
  else IfLine = 0;
  // IfIntorno=0;//?
  // IfColour=0;//?
  // IfChains=0;//?
  // IfChType=0;//?
  // Diameter=0.1;//?
  // StepDiameter=0.;//?
  // ExtraDiam=0.;//?
  // Vicinanze=1.;//?
  InvScaleUn=1.;
  ScaleFact = 1.;
  LineSize = 1;
  NVisSkip = 1;
  NChType = CHAIN_EVERY;
  NBackFold = BF_PART;
  What2Draw = EL_PART;
  for(int d=0;d<3;d++) ShiftPos[d] = 0.;
  for(int d=0;d<3;d++) ScaleF[d] = 1.;
#ifdef OMPI_MPI_H 
  MPI_Init(&argc,&argv);
  int Rank=0,Size=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &Size);
  int Partition = (int)(NFileTot/(double)Size);
  NFile[0] = Partition*Rank;
  NFile[1] = Partition*(Rank+1);
  if(Rank==Size-1) NFile[1] += NFileTot%Size;
  Proc = new SingProc(Size,Rank);
#endif
  //TrialSys();
}
ElPoly::~ElPoly(){
  // for(int f=0;f<NFileTot;f++)
  //   free(cFile[f]);
  // free(cFile);
#ifdef OMPI_MPI_H 
  MPI_Finalize();
#endif
}
int ElPoly::Angle(int Values){
  int NRadici = 4;
  int NZeri = NRadici;
  double *Radici;
  Radici = (double *)malloc(NRadici*sizeof(double));
  double *Average;
  Average = (double *)malloc(NRadici*sizeof(double));
  for(int i=0;i<NRadici;i++){
    Radici[i] = 0.;
    Average[i] = 0.;
  }
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART)==0)return 1;
    Mat->Ypsilon = pEdge(2)-pCm(2)-1.;
    Mat->PreFact = pow((2.*pNPart()/10.)/DUE_PI,.33);
    Mat->Func = &Matematica::ContactAngle;
    printf("Volume  %lf Cm %lf # to Invert %lf\n",(double)pNPart()/10.,pCm(2),pow(DUE_PI/(2.*pNPart()/10.),.25));
    NZeri = Mat->Zeri(0.,DUE_PI/2.,Radici,NRadici);
    for(int i=0;i<NZeri;i++){
      Radici[i] *= 360./DUE_PI;
      Average[i] += Radici[i];
    }
    if(NZeri == 0){
      printf("The function has no real roots\n");
    }
  }
  for(int i=0;i<NZeri;i++){
    Average[i] /= NFile[1]-NFile[0];
    printf("Angle %lf\n",Average[i]);
  }
  return 0;
}
int ElPoly::CenterOfMass(int Coord){
  double SysCm;
  FILE *FileToWrite = NULL;
  FileToWrite = fopen("CenterOfMass.dat","w");
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],1)==0)return 1;
    SysCm=0.;
    if(Coord == 3){
      for(int p=0;p<pNPart();p++){
	SysCm += pVel(p,2);
      }
    }
    else{
      for(int p=0;p<pNPart();p++){
	SysCm += pPos(p,Coord);
      }
    }
      fprintf(FileToWrite,"%lf \n",SysCm/pNPart()); 
  }
  printf("\n"); 
  fclose(FileToWrite);
  return 0;
}
int ElPoly::ChangeFile(){
  printf("Introduce the first of %d file from the list: ",NFileTot);
  scanf("%d",&NFile[0]);
  if(NFile[0] < 0 || NFile[0] > NFile[1])
    NFile[0] = 0;
  quando = NFile[0];
  printf("Introduce the second of %d file from the list: ",NFileTot);
  scanf("%d",&NFile[1]);
  if(NFile[1] < NFile[0] || NFile[1] > NFileTot)
    NFile[1] = NFileTot;
  return 0;
}
int ElPoly::OpenFile(int f){
  if(f >= NFileTot || f < 0)
    return 0;
  if(Open(cFile[f],NBackFold)) return 1;
  NFile[0] = f;
  quando = f;
  return 0;
}
int ElPoly::OpenFile(char *FileName){
  if(Open(FileName,NBackFold))return 1;
  return 0;
}
double ElPoly::ContactAngle(double x){
  double Ris;
  Ris = QUAD( 2*x - QUAD(sin(2*x)) ) / (4./3. * CUB(sin(x))-2.*x*cos(x)+cos(x)*sin(2*x) ) - 5.;
  return Ris;
}
// OpenGl----------------------------------
int ElPoly::PropertiesF(){
  Properties Pr;
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART))return 0;
    Pr = Pr + SysProperties();
  }
  Pr = Pr * (1./(NFile[1] - NFile[0]));
  printf("\n");
  char *NomeFile = (char *) malloc(60*sizeof(char));
  sprintf(NomeFile,"PropertiesChi%.0fRho%.0fKappa%.0f.txt",pchiN(),prho(),pkappaN());
  FILE *FileToWrite = NULL;
  FileToWrite = fopen(NomeFile,"w");
  fprintf(FileToWrite,"RePhob %lf RePhil %lf\n",Pr.RePhob,Pr.RePhil);
  fprintf(FileToWrite,"RadGyrPhob %lf RadGyrPhil %lf\n",Pr.GyrPhob,Pr.GyrPhil);
  fprintf(FileToWrite,"StructFactPhob %lf StructFactPhil %lf\n",Pr.FactPhob,Pr.FactPhil);
  fprintf(FileToWrite,"VolPhob %lf VolPhil %lf\n",Pr.VolPhob,Pr.VolPhil);
  fprintf(FileToWrite,"ChainDiff %lf\n",Pr.ChDiff);
  Pr.Print();
  fclose(FileToWrite);
  return 0;
}
void ElPoly::DivideLayers(int How){
  int NLayer = 2;
  if(How == VAR_OPPOSED) NLayer = 4;
  int *NChStep = (int *)calloc(NLayer,sizeof(int));
  int *ChainLevel = (int *)calloc(pNChain(),sizeof(int));
  if(How == VAR_OPPOSED){
    for(int c=0,b=0;c<pNChain();c++){
      if(c*Block[b].NPCh > Block[b].EndIdx) b++;
      int Level = 0;
      if(CHAIN_IF_TYPE(Ch[c].Type,CHAIN_UP))
	Level = 1;
      if(!strcmp(Block[b].Name,"LIPID1")) Level += 2;
      NChStep[Level]++;
      ChainLevel[c] = Level;
    }
  }
  else if(How == VAR_TUBE || How == VAR_VESICLE){
    for(int c=0,b=0;c<pNChain();c++){
      if(c*Block[b].NPCh > Block[b].EndIdx) b++;
      int Level = 0;
      if(CHAIN_IF_TYPE(Ch[c].Type,CHAIN_OUTER))
	Level = 1;
      NChStep[Level]++;
      ChainLevel[c] = Level;
    }
  }
  //for(int c=0;c<pNChain();c++)printf("%lf %d %d\n",Ch[c].Pos[CNorm],c,ChainLevel[c]);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NO))return ;
    PART *Pn = (PART *)calloc(pNPart(),sizeof(PART));
    int pp=0;
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
    	Pn[p].Pos[d] = pPos(p,d);
    	Pn[p].Vel[d] = pVel(p,d);
      }
      Pn[p].Typ = pType(p);
    }
    for(int l=0;l<NLayer;l++){
      for(int p=0;p<pNPart();p++){
    	if(ChainLevel[pChain(p)] == l){
    	  for(int d=0;d<3;d++){
	    SetPos(pp,Pn[p].Pos);
	    SetVel(pp,Pn[p].Vel);
    	  }
	  SetType(pp,Pn[p].Typ);
    	  pp++;
    	}
      }
    }
    DefBlock(NChStep,How);
    Write(cFile[f]);
    free(Pn);
  }
  printf("\n");
  free(ChainLevel);
  free(NChStep);
}
void ElPoly::RemoveChains(){
  int nChain = pNChain();
  for(int b=0;b<pNBlock();b++){
    for(int c=0;c<Block[b].NChain;c++){
      if(Ch[c].Pos[CLat1] > .4*pEdge(CLat1) && Ch[c].Pos[CLat1] < .6*pEdge(CLat1) ){
	SwapChain(c,nChain-1,b);
	nChain--;
	SetNPart(pNPart()-pNPCh());
      }
    }
    Block[b].NChain = nChain;
  }
  SetNChain(nChain);
  Block[0].EndIdx = Block[0].NChain*Block[0].NPCh;
  for(int b=1;b<pNBlock();b++){
    Block[b].InitIdx = Block[0].NChain*Block[b].NPCh + Block[b-1].EndIdx;
  }
  Write("cleft.dat");
}
// int ElPoly::DensPatch(int NPatch,int Values){
//   int NType = 2;
//   double *DensProf = (double *)calloc(NPatch*Values,sizeof(double));
//   for(int f=NFile[0];f<NFile[1];f++){
//    Processing(f);
//     fprintf(stderr,"Processing: %d/%d %.1f %%\r",f,(NFile[1]-NFile[0]),(f-NFile[0])/(double)(NFile[1]-NFile[0])*100.);
//     if(OpenRisk(cFile[f],BF_PARTICLE))return 0;
    
//   }
//   char *FileName = (char *)calloc(60,sizeof(char));
//   sprintf(FileName,"ChainPArea%.0fRho%.0fChiN%.0fKappa.dat",Gen->rho,Gen->chiN,Gen->kappaN);
//   FileToWrite = fopen(FileName,"w");
//   fprintf(FileToWrite,"#%s\n",SysInfo());
// }
char *ElPoly::ChooseDraw(int ExtWhat2Draw){
  char *String = (char *)calloc(16,sizeof(char));
  if(ExtWhat2Draw == EL_PART)
    sprintf(String,"part");
  else if(ExtWhat2Draw == EL_PART)
    sprintf(String,"part");
  else if(ExtWhat2Draw == EL_CHAIN)
    sprintf(String,"chain");
  else if(ExtWhat2Draw == EL_DENS)
    sprintf(String,"dens");
  else if(ExtWhat2Draw == EL_QUAD)
    sprintf(String,"quad");
  else if(ExtWhat2Draw == EL_QUAD1)
    sprintf(String,"quad1");
  else if(ExtWhat2Draw == EL_POTENTIAL)
    sprintf(String,"pot");
  else if(ExtWhat2Draw == EL_SURF)
    sprintf(String,"surf");
  else if(ExtWhat2Draw == EL_CROSS)
    sprintf(String,"cross");
  else if(ExtWhat2Draw == EL_ISOIPSE)
    sprintf(String,"isoipse");
  else if(ExtWhat2Draw == EL_SAMPLE)
    sprintf(String,"sample");
  else if(ExtWhat2Draw == EL_COLOR)
    sprintf(String,"color");
  else if(ExtWhat2Draw == EL_SKIN)
    sprintf(String,"skin");
  else if(ExtWhat2Draw == EL_POLYGON)
    sprintf(String,"polygon");
  else if(ExtWhat2Draw == EL_ISOLEVEL)
    sprintf(String,"isolevel");
  else if(ExtWhat2Draw == EL_SQUAREMESH)
    sprintf(String,"squaremesh");
  return String;
}
void ElPoly::ChooseDraw(char *String){
  if(!strcmp(String,"part"))
    What2Draw = EL_PART;
  else if(!strcmp(String,"quad"))
    What2Draw = EL_QUAD;
  else if(!strcmp(String,"chain"))
    What2Draw = EL_CHAIN;
  else if(!strcmp(String,"dens"))
    What2Draw = EL_DENS;
  else if(!strcmp(String,"quad1"))
    What2Draw = EL_QUAD1;
  else if(!strcmp(String,"pot"))
    What2Draw = EL_POTENTIAL;
  else if(!strcmp(String,"surf"))
    What2Draw = EL_SURF;
  else if(!strcmp(String,"cross"))
    What2Draw = EL_CROSS;
  else if(!strcmp(String,"isoipse"))
    What2Draw = EL_ISOIPSE;
  else if(!strcmp(String,"sample"))
    What2Draw = EL_SAMPLE;
  else if(!strcmp(String,"color"))
    What2Draw = EL_COLOR;
  else if(!strcmp(String,"stalk"))
    What2Draw = EL_STALK;
  else if(!strcmp(String,"skin"))
    What2Draw = EL_SKIN;
  else if(!strcmp(String,"polygon"))
    What2Draw = EL_POLYGON;
  else if(!strcmp(String,"isolevel"))
    What2Draw = EL_ISOLEVEL;
  else if(!strcmp(String,"squaremesh"))
    What2Draw = EL_SQUAREMESH;
  else
    What2Draw = 0;
}
void ElPoly::SetBoundFile(int InitFile,int EndFile){
  if(InitFile < 0)
    NFile[0] = 0;
  else if(InitFile < EndFile){
    NFile[0] = InitFile;
    quando = InitFile;
  }
  else {
    printf("Wrong init file %d !< %d\n",InitFile,EndFile);
  }
  if(EndFile >= NFileTot)
    NFile[1] = NFileTot-1;
  else if(EndFile > InitFile)
    NFile[1] = EndFile;
  else {
    printf("Wrong end file %d !< %d\n",InitFile,EndFile);
  }
}
void ElPoly::Processing(int f){
  fprintf(stderr,"Processing: %s %d/%d %.1f %%\r",cFile[f],f,(NFile[1]-NFile[0]),(f-NFile[0])/(double)(NFile[1]-NFile[0])*100.);
}
void ElPoly::Shift2Center(){
  double sigma = 
    sqrt(QUAD(pReOverCutOff())/(pNPCh()-1)/3.);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    Open(cFile[f],BF_NO);
    // for(int d=0;d<3;d++){
    //   ShiftPos[d] = .5 - pCm(d)*pInvEdge(d);
    // }
    // ShiftRef(BF_CHAIN);
    double Temp = pEdge(CLat1);
    SetEdge(pEdge(CNorm),CLat1);
    SetEdge(Temp,CNorm);
    
    // FILE *FWrite = fopen("Response.dat","w");
    // fprintf(FWrite,"# L=%lf %lf %lf t=%lf blocks=%d\n",pEdge(0),pEdge(1),pEdge(2),pTime(),pNBlock());
    // HeaderInteraction(FWrite);
    // HeaderNano(FWrite);
    // fprintf(FWrite,"# n=%d N=%d name=%s\n",Block[0].NChain,Block[0].NPCh,Block[0].Name);
    // int NChain = 0;
    // for(int c=0;c<pNChain();c++){
    //   double Dist = 0.;
    //   for(int d=0;d<3;d++){
    // 	Dist += SQR(Nano->Pos[d] - Ch[c].Pos[d]);
    //   }
    //   if(Dist > SQR(2.5*Nano->Rad)) continue;
    //   for(int p=c*pNPCh();p<(c+1)*pNPCh();p++){
    // 	fprintf(FWrite,"%lf %lf %lf %lf %lf %lf %d\n",
    // 		pPos(p,0),pPos(p,1),pPos(p,2),
    // 		pVel(p,0),pVel(p,1),pVel(p,2),pType(p));
    //   }
    //   NChain++;
    // }
    // printf("Number of chains %d\n",NChain);
    // return;

    for(int p=0;p<pNPart();p++){
      double Pos[3] = {pPos(p,CNorm),pPos(p,CLat1),pPos(p,CLat2)};
      SetPos(p,Pos);
      // for(int d=0;d<3;d++){
      // 	Pm[p].Vel[d] = Mat->Gaussiano(0.,sigma);
      // }
    }
    Write(cFile[f]);
  }
}
