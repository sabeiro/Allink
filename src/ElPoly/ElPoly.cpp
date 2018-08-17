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
//#include "../include/Draw.h"

ElPoly *Pol;
void Legenda();

#ifdef __glut_h__
void MenuChoise(int option){
  //Pol->What2Draw = option;
  Pol->ElMenuChoise(option);
}
void MenuVisual(int option){
  //Pol->What2Draw = option;
  Pol->ElMenuVisual(option);
}
#endif
int main(int argc,char **argv){
  if(argc<2){
    printf("Mi serve come argomento un file!\n");
    return 0;
  }
  int *FilePos = (int *)calloc(argc,sizeof(int));
  char Comando[60];
  char ConfF[256];
  char cWhat2Draw[12];
  sprintf(cWhat2Draw,"part");
  sprintf(ConfF,"Polymers.conf");
  sprintf(Comando,"ciccia");
  int NBin = 100;
  int NSample = 32;
  int NFile = 0;
  int IfUser = 1;
  int IfCreate = 0;
  int nBlock = 0;
  int quando = 0;
  int Coord = 2;
  int CoordProfile = 0;
  int nNano = 0;
  int BackFold = BF_PART;
  int NParam = 0;
  double Param[10];
  //Number = (char *)malloc(100*sizeof(char));
  //  frame = (char *)malloc(100*sizeof(char));
  //SpecFuntore <Matematica> SpecMat(&Mat,Matematica::Eval);
  //SpecFuntore <ElPoly> SpecPol(&Pol,ElPoly::ContactAngle);
  for(int i=1;i<argc;i++){
    if(argv[i][0] != '-'){
      FilePos[NFile] = i;
      NFile++;
    }
    if(!strcmp(*(argv+i),"--help")){
      Legenda();
      return 0;
    }
    if(!strcmp(*(argv+i),"-b")){
      char sBf[12];
      if(argc < i){printf("Which backfold?\n");return 1;}
      i += 1;
    }
    else if(!strcmp(*(argv+i),"-bl")){
      if(argc < i){printf("Which block?\n");return 1;}
      sscanf(argv[i+1],"%d",&nBlock);
      i+=1;
    }
    else if(!strcmp(*(argv+i),"-c")){
      if(argc < i){printf("Which normal coordinate?\n");return 1;}
      i+=1;
    }
    else if(!strcmp(*(argv+i),"-cc")){
      if(argc < i){printf("Which coordinate for the profile?\n");return 1;}
      sscanf(argv[i+1],"%d",&CoordProfile);
      i+=1;
    }
    else if(!strcmp(*(argv+i),"-d")){
      sprintf(Comando,"d");
      if(argc < i){printf("What to visualize?\n");return 1;}
      sscanf(argv[i+1],"%s",&cWhat2Draw);
      i+=1;
      IfUser = 0;
    }
    else if(!strcmp(*(argv+i),"-f")){
      if(argc < i+1){printf("Which files?\n");return 1;}
      i+=2;
    }
    else if(!strcmp(*(argv+i),"-n")){
      if(argc < i){printf("Which nanoparticle?\n");return 1;}
      sscanf(argv[i+1],"%d",&nNano);
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--normalize")){
    }
    else if(!strcmp(*(argv+i),"-p")){
      if(argc < i){printf("Which value for the paramter?\n");return 1;}
      i+=1;
    }
    else if(!strcmp(*(argv+i),"-s")){
      if(argc < i){printf("How many samples?\n");return 1;}
      sscanf(argv[i+1],"%d",&NSample);
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--scale")){
      if(argc < i){printf("What scale?\n");return 1;}
      i+=3;
    }
    else if(!strcmp(*(argv+i),"--scalez")){
      if(argc < i){printf("What scale?\n");return 1;}
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--skip")){
      if(argc < i){printf("How many chains to skip?\n");return 1;}
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--Shift")){
      i+=3;
    }
    else if(!strcmp(*(argv+i),"-v")){
      if(argc < i){printf("How many values?\n");return 1;}
      sscanf(argv[i+1],"%d",&NBin);
      i+=1;
    }
    else if(!strcmp(*(argv+i),"-x")){
      IfCreate = 1;
    }
    else if(!strncmp(*(argv+i),"--",2)){
      sprintf(Comando,"%s",argv[i]+2);
      IfUser = 0;
    }
  }
  if(IfCreate){
    VarData *Dat = new VarData;
    Dat->DefSoft(argv[FilePos[0]],ConfF);
    return 0;
  }
  Pol = new ElPoly(NFile,argv,FilePos);
  for(int i=1;i<argc;i++){
    if(!strcmp(*(argv+i),"-b")){
      char sBf[12];
      if(argc < i){printf("Which backfold?\n");return 1;}
      sscanf(argv[i+1],"%s",sBf);
      if(!strcmp(sBf,"part"))
	BackFold = BF_PART;
      else if(!strcmp(sBf,"no"))
	BackFold = BF_NO;
      else if(!strcmp(sBf,"chain"))
	BackFold = BF_CHAIN;
      else if(!strcmp(sBf,"skip"))
	BackFold = BF_SKIP;
      else if(!strcmp(sBf,"nano"))
	BackFold = BF_NANO;
      else if(!strcmp(sBf,"tilt"))
	BackFold = BF_TILT;
      else{
	printf("Backfold not recognized set to chain\n");
	BackFold = BF_CHAIN;
      }
      Pol->SetBackFold(BackFold);
      i+=1;
    }
    else if(!strcmp(*(argv+i),"-c")){
      if(argc < i){printf("Which normal coordinate?\n");return 1;}
      sscanf(argv[i+1],"%d",&Coord);
      Pol->SetCNorm(Coord);
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--normalize")){
      Pol->SetIfNormalize(1);
    }
    else if(!strcmp(*(argv+i),"-f")){
      int InitFile = 0;
      int EndFile = 0;
      if(argc < i+1){printf("Which files?\n");return 1;}
      sscanf(argv[i+1],"%d",&InitFile);
      sscanf(argv[i+2],"%d",&EndFile);
      Pol->SetBoundFile(InitFile,EndFile);
      i+=2;
    }
    else if(!strcmp(*(argv+i),"-p")){
      if(argc < i){printf("Which value for the paramter?\n");return 1;}
      sscanf(argv[i+1],"%lf",Param);
      Pol->ExtParam = Param[0];
      NParam = 1;
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--Ref")){
      double RefPos[3]   = {.0,.0,.0};
      double ShiftPos[3] = {.0,.0,.0};
      if(i<argc-2){
	for(int d=0;d<3;d++){
	  sscanf(*(argv+d+i+1),"%lf",RefPos+d);
	  ShiftPos[d] = RefPos[d] - .5;
	}
	i+=3;
      }
      Pol->SetShiftPos(ShiftPos);
    }  
    else if(!strcmp(*(argv+i),"-s")){
      Pol->NEdge = NSample;
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--scale")){
      double Scale[3];
      if(argc < i){printf("What scale?\n");return 1;}
      for(int d=0;d<3;d++){
	sscanf(*(argv+d+i+1),"%lf",Scale+d);
      }
      Pol->SetScaleF(Scale);
      i+=3;
    }
    else if(!strcmp(*(argv+i),"--scalez")){
      double Scale;
      if(argc < i){printf("What scale?\n");return 1;}
      sscanf(argv[i+1],"%lf",&Scale);
      Pol->ScaleFact = Scale;
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--skip")){
      if(argc < i){printf("How many lipid to skip in the visualisation?\n");return 1;}
      int NSkip = 0;
      sscanf(argv[i+1],"%d",&NSkip);
      Pol->SetNVisSkip(NSkip);
      i+=1;
    }
    else if(!strcmp(*(argv+i),"--Shift")){
      double ShiftPos[3] = {.0,.0,.0};
      if(i<argc-2){
	for(int d=0;d<3;d++){
	  sscanf(*(argv+d+i+1),"%lf",ShiftPos+d);
	}
	i+=3;
      }
      Pol->SetShiftPos(ShiftPos);
    }
    else if(!strcmp(*(argv+i),"-S")){
      int Shift = SHIFT_NO;
      char sBf[12];
      if(argc < i){printf("Which shift?\n");return 1;}
      sscanf(argv[i+1],"%s",sBf);
      if(!strcmp(sBf,"cm"))
	Shift = BF_PART;
      else if(!strcmp(sBf,"nano"))
	Shift = BF_NO;
      else if(!strcmp(sBf,"cm_nano"))
	Shift = BF_CHAIN;
      else{
	printf("shift not recognized set to no\n");
	Shift = SHIFT_NO;
      }
      Pol->ShiftSys(Shift);
      i+=1;
    }
  }
  Pol->Open(argv[FilePos[0]],BackFold);
  if(!IfUser)sprintf(Pol->cWhat2Draw,"%s",cWhat2Draw);
  char cSystem1[STRSIZE];
  char cSystem2[STRSIZE];
  Pol->SysInfo(cSystem1);
  Pol->SysDef(cSystem2);
  printf("------------------------------------------------------\n");
  printf("%s NFile %d\n%s\n",cSystem1,NFile,cSystem2);
  printf("------------------------------------------------------\n");
  while(strcmp(Comando,"q")){
    if(IfUser){
      printf("ElPoly> ");
      scanf("%s",Comando);
    }
    if(!strcmp(Comando,"d")){
#ifdef __glut_h__
      Pol->Graphics(argc,argv);
#else 
      printf("Graphics libraries (GL/glut.h) not supplied\n");
#endif //__glut_h__
    }
    else if(!strcmp(Comando,"?") ){
      Legenda();
    }
    else if(!strcmp(Comando,"!") ){
      printf("Insert a shell command: ");
      scanf("%s",Comando);
      system(Comando);
    }//--------------------Dens-Func----------------------
    else if(!strcmp(Comando,"dens")){
      Pol->DensProf(NBin,NSample,Coord);
      if(!IfUser) return 0;
    }//--------------------3To2-------------------------
    else if(!strcmp(Comando,"3to2")){
      Pol->From3To2d(NSample,Param[0]);
      if(!IfUser) return 0;
    }//--------------------3To1-------------------------
    else if(!strcmp(Comando,"3to1")){
      Pol->From3To1d(Coord);
      if(!IfUser) return 0;
    }//--------------------2To1-------------------------
    else if(!strcmp(Comando,"2to1")){
      Pol->From2To1d(Coord);
      if(!IfUser) return 0;
    }//--------------------Projection----------------------
    else if(!strcmp(Comando,"pro")){
      Pol->ProjectionF(NBin,Coord);
      if(!IfUser) return 0;
    }//--------------------Core-Part----------------------
    else if(!strcmp(Comando,"core")){
      Pol->CoreF(NSample,0);
      if(!IfUser) return 0;
    }//--------------------Core-Nano----------------------
    else if(!strcmp(Comando,"coreNano")){
      Pol->CoreF(NSample,1);
      if(!IfUser) return 0;
    }//--------------------Radial-Shell--------------------
    else if(!strcmp(Comando,"RadShell")){
      Pol->RadialShell(NBin);
      if(!IfUser) return 0;
    }//--------------------Contact-Angle-------------------
    else if(!strcmp(Comando,"angle")){
      Pol->Angle(NBin);
      if(!IfUser) return 0;
    }//--------------------Slab-Prof-----------------------
    else if(!strcmp(Comando,"SlabProf")){
      Pol->SlabProf(NBin,nNano,CoordProfile);
      if(!IfUser) return 0;
    }//--------------------Slab-Prof-----------------------
    else if(!strcmp(Comando,"SlabAngleProfs")){
      Pol->SlabAngleProfs(NBin,10,CoordProfile);
      if(!IfUser) return 0;
    }//--------------------Cart-Dens-----------------------
    else if(!strcmp(Comando,"CartDens")){
      Pol->CartDens(NBin,nNano);
      if(!IfUser) return 0;
    }//--------------------Rad-Nano-----------------------
    else if(!strcmp(Comando,"radNano")){
      Pol->RadDistrF(NBin,1,nNano);
      if(!IfUser) return 0;
    }//--------------------Rad-Cm----------------------------
    else if(!strcmp(Comando,"radCm")){
      Pol->RadDistrF(NBin,0,nNano);
      if(!IfUser) return 0;
    }//--------------------Rad-CmN----------------------------
    else if(!strcmp(Comando,"radCmN")){
      Pol->RadDistrF(NBin,2,nNano);
      if(!IfUser) return 0;
    }//--------------------Rad-CmN----------------------------
    else if(!strcmp(Comando,"radStalk")){
      Pol->RadDistrF(NBin,4,nNano);
      if(!IfUser) return 0;
    }//--------------------Rad-CmN----------------------------
    else if(!strcmp(Comando,"radPore")){
      Pol->RadDistrF(NBin,5,nNano);
      if(!IfUser) return 0;
    }//----------Cartesian-around-the-pep---------------------
    else if(!strcmp(Comando,"cartNano")){
      Pol->RadDistrF(NBin,6,nNano);
      if(!IfUser) return 0;
    }//----------Distribution-of-the-bond-length-----------------
    else if(!strcmp(Comando,"BondDistr")){
      Pol->BondDistr(NSample);
      if(!IfUser) return 0;
    }//----------Distribution-of-the-splay-angle-----------------
    else if(!strcmp(Comando,"SplayDistr")){
      Pol->SplayDistr(NSample);
      if(!IfUser) return 0;
    }//--------------------Center-Of-Mass------------------
    else if(!strcmp(Comando,"cm")){
      Pol->CenterOfMass(Coord);
      if(!IfUser) return 0;
    }//-------------------Thick-prof----------------------
    else if(!strcmp(Comando,"Dens2Thick")){
      Pol->RadDens2Thick(NBin);
      if(!IfUser) return 0;
    }//-------------------Thick-prof----------------------
    else if(!strcmp(Comando,"Dens2Thick2d")){
      Pol->RadDens2Thick2d(NBin);
      if(!IfUser) return 0;
    }//-------------------Thick-prof----------------------
    else if(!strcmp(Comando,"ThickFromDens")){
      Pol->ThickFromDens(NBin);
      if(!IfUser) return 0;
    }//--------------------Temperature----------------------
    else if(!strcmp(Comando,"temp")){
      Pol->Temperature(NBin,Coord);
      if(!IfUser) return 0;
    }//--------------------Nano-Particle---------------------
    else if(!strcmp(Comando,"nano")){
      Pol->NanoParticle(NBin);
      if(!IfUser) return 0;
    }//--------------------Surface------------------------
    else if(!strcmp(Comando,"surf")){
      Pol->Surface(NBin,1);
      if(!IfUser) return 0;
    }//--------------------Diffusivity------------------------
    else if(!strcmp(Comando,"diff")){
      Pol->Diffusivity();
      if(!IfUser) return 0;
    }//--------------------DiffusivitySlab------------------------
    else if(!strcmp(Comando,"DiffSlab")){
      Pol->DiffSlab(NSample);
      if(!IfUser) return 0;
    }//--------------------Pair-Corr--------------------
    else if(!strcmp(Comando,"PairCorr")){
      Pol->PairCorr(NBin,1);
      if(!IfUser) return 0;
    }//--------------------Pair-Correlation--------------------
    else if(!strcmp(Comando,"pairMon")){
      Pol->PairCorrelationF(NBin,0);
      if(!IfUser) return 0;
    }//--------------------Pair-Correlation--------------------
    else if(!strcmp(Comando,"pairChain")){
      Pol->PairCorrelationF(NBin,1);
      if(!IfUser) return 0;
    }//--------------------Pair-Correlation-----------------
    else if(!strcmp(Comando,"pairRound")){
      Pol->PairCorrelationF(NBin,2);
      if(!IfUser) return 0;
    }//--------------------Pair-Correlation------------------
    else if(!strcmp(Comando,"pairSquare")){
      Pol->PairCorrelationF(NBin,3);
      if(!IfUser) return 0;
    }//--------------------Pair-Correlation---------------
    else if(!strcmp(Comando,"pairPep")){
      Pol->PairCorrelationF(NBin,4);
      if(!IfUser) return 0;
    }//--------------------Average-snapshots---------------
    else if(!strcmp(Comando,"AvSnap")){
      Pol->AvSnap();
      if(!IfUser) return 0;
    }//--------------------Scattering--------------------
    else if(!strcmp(Comando,"scatt")){
      Pol->ScatteringF(NBin,0);
      if(!IfUser) return 0;
    }//--------------------Scattering--------------------
    else if(!strcmp(Comando,"scatt2")){
      Pol->ScatteringF(NBin,1);
      if(!IfUser) return 0;
    }//--------------------Worm--------------------
    else if(!strcmp(Comando,"worm")){
      Pol->WormF(20,NBin);
      if(!IfUser) return 0;
    }//--------------------Header--------------------
    else if(!strcmp(Comando,"header")){
      Pol->HeaderAverage(nNano);
      if(!IfUser) return 0;
    }//--------------------Spectrum--------------------
    else if(!strcmp(Comando,"spe")){
      Pol->SpectrumF(NSample);
      return 0;
      if(!IfUser) return 0;
    }//--------------------RadNormPos--------------------
    else if(!strcmp(Comando,"RadNormPos")){
      Pol->RadNormPos(NBin,NSample);
      if(!IfUser) return 0;
    }//--------------------PoreDistr--------------------
    else if(!strcmp(Comando,"AreaDistr")){
      Pol->AreaDistrF(NBin);
      if(!IfUser) return 0;
    }//--------------------Midplane----------------
    else if(!strcmp(Comando,"mid")){
      Pol->Midplane(NSample);
      if(!IfUser) return 0;
    }//--------------------Midplane----------------
    else if(!strcmp(Comando,"midSpe")){
      Pol->SpectrumMidplane(NSample);
      if(!IfUser) return 0;
    }//-------------------chain-properties--------------------
    else if(!strcmp(Comando,"prop")){
      Pol->PropertiesF();
      if(!IfUser) return 0;
    }//-------------------Divide-Layers--------------------
    else if(!strcmp(Comando,"divOpp")){
      Pol->DivideLayers(VAR_OPPOSED);
      if(!IfUser) return 0;
    }//-------------------Divide-Layers--------------------
    else if(!strcmp(Comando,"divTube")){
      Pol->DivideLayers(VAR_TUBE);
      if(!IfUser) return 0;
    }//-------------------Divide-Layers--------------------
    else if(!strcmp(Comando,"divVes")){
      Pol->DivideLayers(VAR_VESICLE);
      if(!IfUser) return 0;
    }//-------------------ChainPArea--------------------
    else if(!strcmp(Comando,"AreaCompr")){
      Pol->AreaCompr(NSample);
      if(!IfUser) return 0;
    }//-------------------Pre-Trace--------------------
    else if(!strcmp(Comando,"PTrace")){
      Pol->PressTrace();
      if(!IfUser) return 0;
    }//-------------------Press-Trace--------------------
    else if(!strcmp(Comando,"PRadial")){
      Pol->PressRadial();
      if(!IfUser) return 0;
    }//-------------------TensProf--------------------
    else if(!strcmp(Comando,"Tens")){
      Pol->SurfTens(NBin);
      if(!IfUser) return 0;
    }//-------------------TensProf--------------------
    else if(!strcmp(Comando,"TensCartRad")){
      Pol->Tens2dCartRad();
      if(!IfUser) return 0;
    }//-------------------TensSum--------------------
    else if(!strcmp(Comando,"SumTens")){
      Pol->SumTens();
      //Pol->Prova();
      if(!IfUser) return 0;
    }//-------------------Prova--------------------
    else if(!strcmp(Comando,"Prova")){
      Pol->Prova();
      if(!IfUser) return 0;
    }//-------------------Stalk--------------------
    else if(!strcmp(Comando,"stalk")){
      Pol->StalkF(32);
      if(!IfUser) return 0;
    }//-------------------Stalk--------------------
    else if(!strcmp(Comando,"StalkLine")){
      Pol->StalkLineProfF(NBin);
      if(!IfUser) return 0;
    }//-------------------BackBone--------------------
    else if(!strcmp(Comando,"BackBone")){
      double *Line = new double[NBin];
      Pol->BackBone(Line,NBin);
      if(!IfUser) return 0;
    }//-------------------Remove--------------------
    else if(!strcmp(Comando,"rem")){
      Pol->RemoveChains();
      if(!IfUser) return 0;
    }//-----------------Chain-position--------------------
    else if(!strcmp(Comando,"ChPos")){
      FILE *FWrite = fopen("ChainPos.dat","w");
      for(int c=0;c<Pol->pNChain();c++){
	fprintf(FWrite,"%lf %lf %lf\n",Pol->Ch[c].Pos[0],Pol->Ch[c].Pos[1],Pol->Ch[c].Pos[2]);
      }
      fclose(FWrite);
      if(!IfUser) return 0;
    }//-------------------Transform--------------------
    else if(!strcmp(Comando,"tra")){
      Pol->Transform(nBlock);
      Pol->Write("transformed.dat");
      if(!IfUser) return 0;      
    }//-------------------WriteXyz--------------------
    else if(!strcmp(Comando,"xyz")){
      Pol->WriteXyz("resume.dat");
      if(!IfUser) return 0;
    }//-------------------WriteXyz--------------------
    else if(!strcmp(Comando,"ConvLattice")){
      Pol->ConvLattice(NSample,"LatticePoints.dat");
      if(!IfUser) return 0;
    }//-------------------Write--------------------
    else if(!strcmp(Comando,"write")){
      Pol->Write("resume.dat");
      if(!IfUser) return 0;
    }//-------------------WidomOut--------------------
    else if(!strcmp(Comando,"WidomOut")){
      //Pol->WidomOut("NrgLipid.dat",NBin);
      Pol->WidomOut();
      if(!IfUser) return 0;
    }//-------------------WidomIn--------------------
    else if(!strcmp(Comando,"WidomIn")){
      //Pol->WidomIn("NrgLipid.dat",NBin);
      Pol->WidomIn();
      if(!IfUser) return 0;
    }//-------------------End2EndDistr--------------------
    else if(!strcmp(Comando,"End2End")){
      Pol->End2EndDistr("End2EndDistr");
      if(!IfUser) return 0;
    }//-----------------ElasticCoupling--------------------
    else if(!strcmp(Comando,"ElCoup")){
      Pol->ElasticCoupling(NSample);
      if(!IfUser) return 0;
    }//-----------------ElasticCoupling--------------------
    else if(!strcmp(Comando,"ElCoupNVT")){
      Pol->ElasticCouplingNVT();
      if(!IfUser) return 0;
    }//-----------------E2EDecoupling--------------------
    else if(!strcmp(Comando,"E2EDec")){
      Pol->Decoupling(1);
      if(!IfUser) return 0;
    }//-----------------SmoothGrid------------------------
    else if(!strcmp(Comando,"SmoothGrid")){
      Pol->SmoothGrid(NBin,"GridSmoothed.dat");
      if(!IfUser) return 0;
    }//-----------------IsoSurface------------------------
    else if(!strcmp(Comando,"IsoSurf")){
      Pol->IsoSurf(NSample,Param,NParam);
      if(!IfUser) return 0;
    }//-----------------IsoLine------------------------
    else if(!strcmp(Comando,"IsoLine")){
      NParam = 10;
      double OldParam = Param[0];
      for(int p=0;p<NParam;p++){
	Param[p] = OldParam*p/(double)(NParam);
      }
      Pol->IsoLine(NSample,Param,NParam,0);
      if(!IfUser) return 0;
    }//-----------------IsoLineDens------------------------
    else if(!strcmp(Comando,"IsoLineDens")){
      Pol->IsoLine(NSample,Param,NParam,1);
      if(!IfUser) return 0;
    }//-----------------Follow-Pore------------------------
    else if(!strcmp(Comando,"FPore")){
      Pol->FetchPore();
      if(!IfUser) return 0;
    }//-----------------Follow-Stalk------------------------
    else if(!strcmp(Comando,"FStalk")){
      Pol->FetchStalk();
      if(!IfUser) return 0;
    }//-----------------Area-Stalk------------------------
    else if(!strcmp(Comando,"AStalk")){
      Pol->StalkArea();
      if(!IfUser) return 0;
    }//-----------------Shift-to-center------------------------
    else if(!strcmp(Comando,"Shift2Center")){
      Pol->Shift2Center();
      if(!IfUser) return 0;
    }//-----------------AvSnap------------------------
    else if(!strcmp(Comando,"AvSnap")){
      Pol->AvSnap();
      if(!IfUser) return 0;
    }//-----------------DirDecoupling--------------------
    else if(!strcmp(Comando,"DirDec")){
      Pol->Decoupling(0);
      if(!IfUser) return 0;
    }//-----------------ElasticCoupling--------------------
    else if(!strcmp(Comando,"BilDist")){
      Pol->BilayerDistance("BilayerDistance.dat",NSample);
      if(!IfUser) return 0;
    }//-------------------BondDistr--------------------
    else if(!strcmp(Comando,"Bond")){
      Pol->BondDistr("BondDistr.dat",NBin);
      if(!IfUser) return 0;
    }//-------------------Difference,-two-files---------------
    else if(!strcmp(Comando,"Diff2Dens")){
      Pol->Diff2Files(NSample,1);
      if(!IfUser) return 0;
    }//-------------------Difference,-two-files---------------
    else if(!strcmp(Comando,"Diff2Pre")){
      Pol->Diff2Files(NSample,0);
      if(!IfUser) return 0;
    }//-------------------Difference,-two-files---------------
    else if(!strcmp(Comando,"RestPre")){
      Pol->RestPress(NBin);
      if(!IfUser) return 0;
    }//-------------------Conv-to-tecplot-----------------
    else if(!strcmp(Comando,"tecplotDens")){
      Pol->Conv2Tecplot(NBin,0);
      if(!IfUser) return 0;
    }//-------------------Conv-to-tecplot-----------------
    else if(!strcmp(Comando,"tecplotPre")){
      Pol->Conv2Tecplot(NBin,1);
      if(!IfUser) return 0;
    }//-------------------Conv-to-tecplot-----------------
    else if(!strcmp(Comando,"tecplotHei")){
      Pol->Conv2Tecplot(NBin,2);
      if(!IfUser) return 0;
    }//-------------------Conv-to-rzd-----------------
    else if(!strcmp(Comando,"rzd")){
      Pol->Conv2rzd(NBin);
      if(!IfUser) return 0;
    }//-------------------Conv-to-rzd-----------------
    else if(!strcmp(Comando,"xyzd")){
      Pol->Conv2xyzd(NBin);
      if(!IfUser) return 0;
    }//-------------------Conv-to-vmd--------------
    else if(!strcmp(Comando,"vmd")){
      Pol->Conv2Vmd();
      if(!IfUser) return 0;
    }//-------------------povray--------------
    else if(!strcmp(Comando,"povray")){
      Pol->Conv2Povray();
      if(!IfUser) return 0;
    }//-----------------------Coord------------------
    else if(!strcmp(Comando,"coord")){
      printf("Enter coordinate number [0 2] or radius [3]: ");
      scanf("%d",&Coord);
      if(Coord < 0 || Coord > 3){
	printf("Value not valid, coord set to 0\n");
	Coord = 0;
      }
      if(!IfUser) return 0;
    }//-----------------------Normal------------------
    else if(!strcmp(Comando,"norm")){
      printf("Enter the direction of the normal: ");
      int Temp = 0;
      scanf("%d",&Temp);
      if(Temp >= 0 && Temp < 3){
	Pol->CNorm = Temp;
	Pol->CLat1 = (Pol->CNorm+1)%3;
	Pol->CLat2 = (Pol->CNorm+2)%3;
	printf("n %d c1 %d c2 %d\n",Pol->CNorm,Pol->CLat1,Pol->CLat2);
      }
      if(!IfUser) return 0;
    }//--------------------Chain-Type------------------
    else if(!strcmp(Comando,"type")){
      int NChType = CHAIN_EVERY;
      printf("up %d down %d flabby %d stretch %d added %d every %d: ",CHAIN_UP,CHAIN_DOWN,CHAIN_FLABBY,CHAIN_STRETCH,CHAIN_ADDED,CHAIN_EVERY);
      scanf("%d",&NChType);
      if(NChType != CHAIN_UP && NChType != CHAIN_DOWN && NChType != CHAIN_FLABBY && NChType !=  CHAIN_STRETCH && NChType != CHAIN_ADDED && NChType != CHAIN_EVERY){
	printf("Value not valid, type set to %d\n",CHAIN_EVERY);
	NChType = CHAIN_EVERY;
      }
      Pol->NChType = NChType;
      if(!IfUser) return 0;
    }//------------------Change-Files------------------------
    else if(!strcmp(Comando,"file")){
      Pol->ChangeFile();
      if(!IfUser) return 0;
    }//------------------Change-NBin---------------------
    else if(!strcmp(Comando,"val")){
      printf("Number bins (NBin) :");
      scanf("%d",&NBin);
      if(!IfUser) return 0;
    }//------------------Find-neighbours---------------------
    else if(!strcmp(Comando,"Nei")){
      Pol->FindNeighbours("CrossLinks.dat");
      if(!IfUser) return 0;
    }//------------------Change-NSample---------------------
    else if(!strcmp(Comando,"sample")){
      Pol->Sample(NSample);
      if(!IfUser) return 0;
    }//-------------------Open-Another-File-----------------------------
    else if(!strcmp(Comando,"open")){
      char InFile[256];
      printf("New file name: ");
      scanf("%s",Comando);
      sprintf(InFile,Comando);
      Pol->OpenFile(InFile);
      if(!IfUser) return 0;
    }//---------------Next-File-in-the-List----------------------
    else if(!strcmp(Comando,">")){
      quando++;
      if(quando >= NFile){
	quando = 0;
      }
      if(quando >=0 && quando < NFile){
	Pol->OpenFile(quando);
	printf("Opening: %s\n",argv[FilePos[quando]]);
      }
      if(!IfUser) return 0;
    }//---------------Previous--------------------------
    else if(!strcmp(Comando,"<")){
      quando--;
      if(quando < 0){
	quando = NFile-1;
      }
      if(quando >=0 && quando < NFile){
	Pol->OpenFile(quando);
	printf("Opening: %s\n",argv[FilePos[quando]]);
      }
      if(!IfUser) return 0;
    }//----------------Open-f-File-----------------
    else if(!strcmp(Comando,"f")){
      printf("Open file num 0-%d: ",NFile);
      scanf("%d",&quando);
      if( quando >= 0 && quando < NFile){
	Pol->OpenFile(quando);
      }
      if(!IfUser) return 0;
    }//-----------------System-Info---------------------------
    else if(!strcmp(Comando,"info")){
      printf("------------------------------------------------------\n");
      char cSystem[STRSIZE];
      Pol->SysInfo(cSystem);
      printf("%s\n",cSystem);
      Pol->SysDef(cSystem);
      printf("%s\n",cSystem);
    }
    else if(!strcmp(Comando,"q") ){
      if(!IfUser) return 0;
    }
    else {
      printf("Comando non valido, scrivere ?\
 per la lista dei comandi\n");
      IfUser = 1;
    }
  }
  printf("Te se qe te ve be te ne?\n");
  if(Pol) delete Pol;
  return 0;
}
#ifdef __glut_h__
void Slide(){
  Pol->ESlide();
}
void ParticleList(){
  Pol->RenderPart();
}
void ParticleRealTime(){
  return;
  Pol->DrRunTime();
}
void keyboard(unsigned char key,int x, int y){
  Pol->keyboard(key,x,y);
}
void mouse(int button, int state,int x,int y){
  Pol->ElDrawMouse(button,state,x,y);
}
void Menu(){
  Pol->Menu();
}
#endif // __glut_h__
void Legenda(){
  printf("\n\
*****************************************************************\n\
(                                                               )\n\
)    Program that reads and elaborate a spefic file format      (\n\
(    the option are the following, Needs a input file           )\n\
)                                                               (\n\
(-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-)\n\
(   d     draw the 3D position of the particle                  )\n\
(   ! execute a shell command                                   )\n\
(   q    quit                                                   )\n\
(   dens density of the [coord] coordinate for the [file]       )\n\
(   pro calculate the projection long [coord]                   )\n\
(   core spatial average sample on NBin bin                     )\n\
(   surf ratio between the actual surcafe and the circular      )\n\
(   core three dimentional sampled average of the system        )\n\
(   diff diffusivity of the extra particle                      )\n\
(   > <  open the file form the list position                   )\n\
(   cm   movement  of the center of mass of the system          )\n\
(   file number of file to use from the list                    )\n\
(   info ...                                                    )\n\
(   val  number of values for the density profile               )\n\
(   radShell average the outer shell of a drop                  )\n\
(   radCm normal and radial density wrt the center of mass      )\n\
(   radCmN normal wrt to the center of mass, radial to the nano )\n\
(   radNano normal and radial density wrt the Nano              )\n\
(   temp Temperature profile along a direction                  )\n\
(   nano density and diffisivity of the nano in the membrane    )\n\
(   surf sampling and averaging a surface                       )\n\
(   pairMon pair correlation function of the monomers           )\n\
(   pairChain pair correlation function of the chains           )\n\
(   pairRound radial pair correlation function                  )\n\
(   scatt scattering og the membrane                            )\n\
(   scatt2 alternative formulation                              )\n\
(   worm projecting the density of the system on the local norma)\n\
(   header averaging the values written on the header           )\n\
(   spe calculating the spectrum                                )\n\
(   mid sampling the membrane in midplanes                      )\n\
(   prop calculating some characteristic properties of the chain)\n\
(   Area calculating the area compressibility                   )\n\
(                                                               )\n\
(                                                               )\n\
(                                                               )\n\
*****************************************************************\n");
}
//  else if(!strcmp(Comando,"coreNano")){
//  else if(!strcmp(Comando,"radShell")){
//  else if(!strcmp(Comando,"angle")){
//  else if(!strcmp(Comando,"radNano")){
// else if(!strcmp(Comando,"radCm")){
//  else if(!strcmp(Comando,"radCmN")){
//  else if(!strcmp(Comando,"cm")){
//  else if(!strcmp(Comando,"temp")){
//  else if(!strcmp(Comando,"nano")){
//  else if(!strcmp(Comando,"surf")){
//  else if(!strcmp(Comando,"diff")){
//  else if(!strcmp(Comando,"pairMon")){
//  else if(!strcmp(Comando,"pairChain")){
//  else if(!strcmp(Comando,"pairRound")){
//  else if(!strcmp(Comando,"pairSquare")){
//  else if(!strcmp(Comando,"scatt")){
//  else if(!strcmp(Comando,"scatt2")){
//  else if(!strcmp(Comando,"worm")){
//  else if(!strcmp(Comando,"header")){
//  else if(!strcmp(Comando,"spe")){
//  else if(!strcmp(Comando,"mid")){
//  else if(!strcmp(Comando,"prop")){
//   else if(!strcmp(Comando,"Area")){
//  else if(!strcmp(Comando,"coord")){
//  else if(!strcmp(Comando,"norm")){
//  else if(!strcmp(Comando,"type")){
//  else if(!strcmp(Comando,"file")){
//  else if(!strcmp(Comando,"val")){
//  else if(!strcmp(Comando,"sample")){
//  else if(!strcmp(Comando,"open")){
//  else if(!strcmp(Comando,">")){
//  else if(!strcmp(Comando,"<")){
// else if(!strcmp(Comando,"f")){
