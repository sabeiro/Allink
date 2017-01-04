/***********************************************************************
VarData: This  Program reads and writes a specific file format 
storing all the information relative to a set of equal structure
polymers to the CHAIN, PART and GENERAL structures. It provides 
two different ways to backfold the coordinates, a fanction that 
creates an initial system with different option and some function
for the data analisys. The first calculate the distribution of the
monomer in the box, the second the distribution of the bonds.
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
#include "../include/VarData.h"
;
bool VarData::ReadConf(char *InFile){
  FILE *FileToRead;
  VarMessage("ReadConf");
  if((FileToRead = fopen(InFile,"r"))==0){
    printf("The file %s is missing\n",InFile);
    return 1;
  }
  //double *buff = (double *)malloc(sizeof(double));
  double buff[1];
  char cLine[STRSIZE];
  char Topology[20];
  int NCircle = 0;
  int NHeight = 0;
  int DiblockLim = 0;
  SysFormat = 0;
  //  fgets(cLine,256,FileToRead);
  for(int k=0;!(fgets(cLine,STRSIZE,FileToRead)==NULL);k++){
    if(cLine[0] == '#') continue;
    if(ReadString("NSoft",cLine,buff)==1){
      NSoft = (int)*buff;
      if(NSoft != 0) 
	Soft = (SOFT *)realloc(Soft,NSoft*sizeof(SOFT));
    }
    //sytem type
    if(ReadString("IfSystem",cLine,buff)==1){
      if( (int)*buff == 1 ){
	SysFormat = VAR_SYS_TXVL;
      }
      else{
	SysFormat = VAR_SYS_XVT;
      }
    }
    //
    if(ReadString("CNorm",cLine,buff)==1){
      CNorm = (int)*buff;
      CLat1 = (CNorm+1)%3;
      CLat2 = (CNorm+2)%3;
    }
    if(ReadString("NAddChain",cLine,buff)==1)
      NAddChain = (int)*buff;
    if(ReadString("NAddChol",cLine,buff)==1)
      NAddChol = (int)*buff;
    if(ReadString("NSolvent",cLine,buff)==1)
      NSolvent = (int)*buff;
    if(ReadString("NStuffing",cLine,buff)==1)
      NStuffing = (int)*buff;
    if(ReadString("DiblockLim",cLine,buff)==1)
      DiblockLim = (int)*buff;
    if(ReadString("IfTwoTails",cLine,buff)==1)
      if( (int)*buff == 1 ){
	VAR_ADD_TYPE(SysType,VAR_TWOTAILS);
      }      
    if(ReadString("NPartPChain",cLine,buff)==1)
      SetNPCh((int)*buff);
    if(ReadString("NNano",cLine,buff)==1){
      Gen->NNano = (int)*buff;
      if(Gen->NNano != 0) 
	Nano = (NANO *)realloc(Nano,Gen->NNano*sizeof(NANO));
    }
    if(ReadString("NCircle",cLine,buff)==1){
      NCircle = (int)*buff;
    }
    if(ReadString("NHeight",cLine,buff)==1){
      NHeight = (int)*buff;
    }
    if(ReadString("rho",cLine,buff)==1)
      Gen->rho = *buff;
    if(ReadString("chiN",cLine,buff)==1)
      Gen->chiN = *buff;
    if(ReadString("kappaN",cLine,buff)==1)
      Gen->kappaN = *buff;
    if(ReadString("kappaBend",cLine,buff)==1)
      Gen->kappaBend = *buff;
    if(ReadString("kappaSpring",cLine,buff)==1)
      Gen->kappaSpring = *buff;
    if(ReadString("ReOverCutOff",cLine,buff)==1)
      Gen->ReOverCutOff = *buff;
    if(ReadString("WFuncStraight2",cLine,buff)==1)
      Gen->WFuncStraight2 = *buff;
    if(ReadString("WFuncStraight3",cLine,buff)==1)
      Gen->WFuncStraight3 = *buff;
    if(ReadString("vBB",cLine,buff)==1)
      Gen->vBB = *buff;
    if(ReadString("Lx",cLine,buff)==1)
      Gen->Edge[0] = *buff;
    if(ReadString("Ly",cLine,buff)==1)
      Gen->Edge[1] = *buff;
    if(ReadString("Lz",cLine,buff)==1)
      Gen->Edge[2] = *buff;
  }
  ReadSoft(FileToRead);
  ReadNano(FileToRead,NCircle,NHeight);
  SetNBlock(1);
  Block[0].Asym = DiblockLim;
  fclose(FileToRead);
  // char cSystem[512];
  // SysDef(cSystem);
  // printf("%s\n",cSystem);
  // HeaderNano(cSystem);
  // printf("%s",cSystem);
  // HeaderSoft(cSystem);
  // printf("%s",cSystem);
  return 0;
}
int VarData::NanoString(char *cLine,int n){
  double Pos[3];
  double Vel[3];
  double Axis[3];
  double Char[5];
  char Shape[20];
  char Name[20];
  if( !Fetch(cLine,"x","%lf %lf %lf",Pos,Pos+1,Pos+2)){
    for(int d=0;d<3;d++)
      Nano[n].Pos[d] = pEdge(d)*.5;
  }
  else{
    for(int d=0;d<3;d++){
      if(isnan(Pos[d])) Pos[d] = .5*pEdge(d);
      Nano[n].Pos[d] = Pos[d]*ScaleF[d];
      Nano[n].Bkf[d] = - floor(Nano[n].Pos[d]/pEdge(d))*pEdge(d);
    }
  }    
  if( !Fetch(cLine,"a","%lf %lf %lf",Axis,Axis+1,Axis+2)){
    Nano[n].Axis[CLat1] = 0.;
    Nano[n].Axis[CLat2] = 0.;
    Nano[n].Axis[CNorm] = 1.;
  }
  else{
    if(isnan(Axis[0])) Axis[0] = 0.;
    if(isnan(Axis[1])) Axis[1] = 0.;
    if(isnan(Axis[2])) Axis[2] = 1.;
    Nano[n].Axis[0] = Axis[0];
    Nano[n].Axis[1] = Axis[1];
    Nano[n].Axis[2] = Axis[2];
  }
  if( Fetch(cLine,"c","%lf %lf %lf %lf",Char,Char+1,Char+2,Char+3)){
    double r = sqrt(.5*SQR(ScaleF[0])+.5*SQR(ScaleF[1]));
    Nano[n].Rad = Char[0]*r;
    Nano[n].Viscosity = 0.1;
    Nano[n].Hamaker = Char[1];
    Nano[n].Height = Char[2]*ScaleF[2];
    Nano[n].OffSet = Char[3];
    Nano[n].Coating = Char[3];
    Nano[n].Gamma = 3.*3.14*Nano[n].Rad*Nano[n].Viscosity;
    //Nano[n].Zeta = sqrt(12. * 2. * 1. * Nano[n].Gamma/Nano[n].Mass/dt);
  }
  // if( ret == 8){
  if( Fetch(cLine,"s","%s",Shape)){
    Nano[n].Shape = ShapeId(Shape);
  }
  Nano[n].Mass = 100.;
  return 0;
}
int VarData::ShapeId(char *Shape){
  int iShape=0;
  if(!strcmp(Shape,"no"))
    VAR_ADD_TYPE(iShape,SHAPE_NONE);
  else if(!strcmp(Shape,"sph")) 
    VAR_ADD_TYPE(iShape,SHAPE_SPH);
  else if(!strcmp(Shape,"tip")) 
    VAR_ADD_TYPE(iShape,SHAPE_TIP);
  else if(!strcmp(Shape,"dip")) 
    VAR_ADD_TYPE(iShape,SHAPE_SPH);
  else if(!strcmp(Shape,"cyl")){
    VAR_ADD_TYPE(iShape,SHAPE_CYL);
    VAR_ADD_TYPE(iShape,SHAPE_HEI);
  }
  else if(!strcmp(Shape,"tilt")){ 
    VAR_ADD_TYPE(iShape,SHAPE_TILT);
    VAR_ADD_TYPE(iShape,SHAPE_HEI);
  }
  else if(!strcmp(Shape,"pill")){
    VAR_ADD_TYPE(iShape,SHAPE_PILL);
    VAR_ADD_TYPE(iShape,SHAPE_HEI);
  }
  else if(!strcmp(Shape,"wall")) 
    VAR_ADD_TYPE(iShape,SHAPE_WALL);
  else if(!strcmp(Shape,"cluster")) 
    VAR_ADD_TYPE(iShape,SHAPE_CLUSTER);
  else if(!strcmp(Shape,"harm")) 
    VAR_ADD_TYPE(iShape,SHAPE_HARM);
  else if(!strcmp(Shape,"clinks")) 
    VAR_ADD_TYPE(iShape,SHAPE_CLINKS);
  else if(!strcmp(Shape,"pore")) 
    VAR_ADD_TYPE(iShape,SHAPE_PORE);
  else if(!strcmp(Shape,"ext")) 
    VAR_ADD_TYPE(iShape,SHAPE_EXT);
  else if(!strcmp(Shape,"janus")) 
    VAR_ADD_TYPE(iShape,SHAPE_JANUS);
  else if(!strcmp(Shape,"stalk")) 
    VAR_ADD_TYPE(iShape,SHAPE_STALK);
  else if(!strcmp(Shape,"torus"))
    VAR_ADD_TYPE(iShape,SHAPE_TORUS);
  else if(!strcmp(Shape,"umbr"))
    VAR_ADD_TYPE(iShape,SHAPE_UMBR);
  else if(!strcmp(Shape,"bound"))
    VAR_ADD_TYPE(iShape,SHAPE_BOUND);
  else 
    printf("Nano type %s not recognized\n",Shape);
  return iShape;
}
void VarData::ReadNano(FILE *ConfFile,int NCircle,int NHeight){
  char cLine[STRSIZE];
  rewind(ConfFile);
  for(int k=0,n=0;!(fgets(cLine,STRSIZE,ConfFile)==NULL);k++){
    if(n == Gen->NNano) break;
    if(cLine[0] == '#') continue;
    if(strstr(cLine, "Rigid") == cLine) {
      NanoString(cLine,n);
      for(int d=0;d<3;d++){
	Nano[n].Pos[d] *= pEdge(d);
	Nano[n].Vel[d] = 0.;
	Nano[n].AVel[d] = 0.;
      }
      Nano[n].NCircle = (int)(NCircle*Nano[n].Rad);
      if( (Nano[n].NCircle%2) != 0)
	Nano[n].NCircle++;
      Nano[n].NHeight = (int)(NHeight*Nano[n].Height);
      if( (Nano->NHeight%2) != 0)
	Nano->NHeight++;
      n++;
    }
  }
}
void VarData::SubNanoHeader(char *cFile){
  FILE *FRead = fopen(cFile,"r+");
  char cLine[STRSIZE];
  fpos_t fPosNano;
  do {
    fgetpos(FRead,&fPosNano);
    fgets(cLine, sizeof(cLine),FRead);
    if(strstr(cLine, "# Rigid") == cLine){
      break;
    }
    else if(strstr(cLine, "#") != cLine){
      printf("No # Rigid present!\n");
      break;
    }
  } while (1==1);
  fgets(cLine, sizeof(cLine),FRead);
  fsetpos(FRead,&fPosNano);
  HeaderNano(FRead);
  fprintf(FRead,"%s",cLine);
  fclose(FRead);
}
int VarData::ReadSoft(FILE *ConfFile){
  double Pos[3];
  double Vel[3];
  double Char[5];
  int iShape=0;
  char Shape[20];
  char Name[20];
  char cLine[STRSIZE];
  rewind(ConfFile);
  for(int k=0,n=0;!(fgets(cLine,STRSIZE,ConfFile)==NULL);k++){
    if(cLine[0] == '#') continue;
    if(n == NSoft) break;
    if (strstr(cLine, "Soft") == cLine) {
      if( !Fetch(cLine,"x","%lf %lf %lf",Pos,Pos+1,Pos+2)){
	Soft[n].Pos[0] = 0.;
	Soft[n].Pos[1] = 0.;
	Soft[n].Pos[2] = 0.;
      }
      else{
	Soft[n].Pos[0] = Pos[0]*pEdge(0);
	Soft[n].Pos[1] = Pos[1]*pEdge(1);
	Soft[n].Pos[2] = Pos[2]*pEdge(2);
      }    
      if( !Fetch(cLine,"v","%lf %lf %lf",Vel,Vel+1,Vel+2)){
	Soft[n].Vel[0] = 0.;
	Soft[n].Vel[1] = 0.;
	Soft[n].Vel[2] = 0.;
      }
      else{
	Soft[n].Vel[0] = Vel[0];
	Soft[n].Vel[1] = Vel[1];
	Soft[n].Vel[2] = Vel[2];
      }    
      if( !Fetch(cLine,"c","%lf %lf %lf",Char,Char+1,Char+2)){
	printf("Rigid characteristic are not specified\n");
	return 1;
      }
      Soft[n].Size[0] = Char[0];
      Soft[n].Size[1] = Char[1];
      Soft[n].Size[2] = Char[2];
      if( !Fetch(cLine,"s","%s",Shape)){
	printf("Soft shape is not specified\n");
	return 1;
      }
      if(!strcmp(Shape,"planar")){
	VAR_ADD_TYPE(Soft[n].Topology,VAR_PLANAR);
      }
      else if(!strcmp(Shape,"planarPE")){
	VAR_ADD_TYPE(Soft[n].Topology,VAR_PLANAR_PE);
      }
      else if(!strcmp(Shape,"tube"))
	VAR_ADD_TYPE(Soft[n].Topology,VAR_TUBE);
      else if(!strcmp(Shape,"obstacle"))
	VAR_ADD_TYPE(Soft[n].Topology,VAR_OBSTACLE);
      else if(!strcmp(Shape,"coating"))
	VAR_ADD_TYPE(Soft[n].Topology,VAR_COATING);
      else if(!strcmp(Shape,"distributed"))
	VAR_ADD_TYPE(Soft[n].Topology,VAR_DISTRIBUTED);
      else if(!strcmp(Shape,"vesicle"))
	VAR_ADD_TYPE(Soft[n].Topology,VAR_VESICLE);
      else{
	printf("Soft type not recognized %s\n",Shape);
	return 1;
      }
      sprintf(Soft[n].Name,"LIPID");
      // if( !Fetch(cLine,"n","%s",Name)){
      // 	printf("Soft name not specified\n");
      // 	return 1;
      // }
      //      sprintf(Soft[n].Name,"%s",Name);
      n++;
    }
  }
  return 0;
}
//#################SYS#INFO##################################
void VarData::ReadHeader(FILE *FileToRead){
  VarMessage("ReadSysInfo");
  SysFormat = 0;
  VAR_REM_TYPE(SysType,VAR_CHAIN_DEF);
  char cLine[STRSIZE];
  double Val[6];
  for(int d=0;d<3;d++){
    Gen->Cm[d] = 0.;
  }
  fgets(cLine,STRSIZE,FileToRead);
  if(Fetch(cLine,"l",3,Val)){
    SysFormat = VAR_SYS_TXVL;
    VAR_ADD_TYPE(SysType,VAR_EDGE);
    SetEdge(Val[0]*ScaleF[0],0);
    SetEdge(Val[1]*ScaleF[1],1);
    SetEdge(Val[2]*ScaleF[2],2);
  }
  //else if(!strcspn(cLine,"# L=\n")){
  else if(!strncmp(cLine,"# L=",4)){
    SysFormat = VAR_SYS_XVT;
  }
  else{
    do{
      fgets(cLine,STRSIZE,FileToRead);
      char *pLine = strchr(cLine,'#');
      if(pLine == NULL){
	int NPar = sscanf(cLine,"%lf %lf %lf %lf\n",Val,Val+1,Val+2,Val+3);
	if(NPar == 3){
	  SysFormat = VAR_SYS_XYZ;
	  SetNChain(1);
	  SetNBlock(1);
	}
	else if(NPar == 4){
	  SysFormat = VAR_SYS_XYZT;
	  SetNChain(1);
	  SetNBlock(1);
	}
	break;
      }
    } while(1==1);
  }
  rewind(FileToRead);
  if( SysFormat == VAR_SYS_TXVL ){
    ReadHeaderTxvl(FileToRead);
  }
  else if ( SysFormat == VAR_SYS_XVT ){
    ReadHeaderXvt(FileToRead);
  }
  rewind(FileToRead);
}
void VarData::ReadHeaderTxvl(FILE *FileToRead){
  double Val[6];
  char cLine[STRSIZE];
  fgets(cLine,STRSIZE,FileToRead);
  Gen->Time += 1.;
  Gen->NNano = 0;
  //Energy
  //#Chain
  if(Fetch(cLine,"c",1,Val)){
    SetNChain((int)Val[0]);
  }
  //What to draw
  if(Fetch(cLine,"d","%s",cWhat2Draw)); 
  //Delta t
  if(Fetch(cLine,"D",1,Val)) SetDeltat(Val[0]);
  if(Fetch(cLine,"e",3,Val)){
    Gen->Energy[0] = Val[0];
    Gen->Energy[1] = Val[1];
    Gen->Energy[2] = Val[2];
  }
  //inputs
  if(Fetch(cLine,"i",3,Val)){
    Gen->rho    = Val[0];
    Gen->chiN   = Val[1];
    Gen->kappaN = Val[2];
  }
  //#link
  if(Fetch(cLine,"L",1,Val)){
    SetNLink((int)Val[0]);
  }
  //#Part
  if(Fetch(cLine,"n",1,Val)){
    //SetNLink(2);
    SetNPart((int)Val[0]);
    VAR_ADD_TYPE(SysType,VAR_OPEN_TRUST);
  }
  //Nano
  if(Fetch(cLine,"N",4,Val)){
    Gen->NNano = 1;
    //    Nano = (NANO *) realloc(Nano,Gen->NNano*sizeof(NANO));
    Nano->Rad = Val[0];
    Nano->Hamaker = Val[1];
    Nano->Viscosity = Val[2];
    Nano->Height = Val[3];
    if(Nano->Hamaker > 0.1){
      Nano->Shape = SHAPE_SPH;
      if(Nano->Height > 0.1){
	Nano->Shape = SHAPE_CYL;
      }
    }
  }
  if(Fetch(cLine,"P",1,Val)){
    if( *Val < 0. || *Val > 3.);
    else{
      CNorm = (int)*Val;
      CLat1 = (CNorm+1)%3;
      CLat2 = (CNorm+2)%3;
    }
  }
  if(Fetch(cLine,"r",4,Val)){
    Gen->NNano = 1;
    //    Nano = (NANO *) realloc(Nano,Gen->NNano*sizeof(NANO));
    Nano->Pos[0] = Val[0];
    Nano->Pos[1] = Val[1];
    Nano->Pos[2] = Val[2];
    Nano->Axis[0] = 0.;
    Nano->Axis[1] = 0.;
    Nano->Axis[2] = 1.;
    Nano->Shape = SHAPE_NONE;
  }
  //step
  if(Fetch(cLine,"s",1,Val)) SetStep((int) Val[0]);
  //Temperature
  if(Fetch(cLine,"T",1,Val)) SetTemp(Val[0]);
  //# Values
  if(Fetch(cLine,"v",1,Val)){
    NEdge = (int)Val[0];
    if(!strcmp(cWhat2Draw,"color")){
      SetNLink(0);
      SetNPart(NEdge*NEdge*NEdge);
      SetNChain(1);
      SetNPCh(pNPart());
      VAR_ADD_TYPE(SysType,VAR_OPEN_TRUST);
    }
  }
  SetNBlock(1);
  //----------inclusion------------------------
  // how many Nano
  //Gen->NNano = 0;
  fpos_t PosTemp;
  fgetpos(FileToRead,&PosTemp);
  do {
    fgets(cLine, sizeof(cLine),FileToRead);
    if(strstr(cLine, "# Rigid") == cLine)
      Gen->NNano++;
    else if(strstr(cLine, "# Pep") == cLine)
      Gen->NNano++;
    else 
      break;
  } while (1==1);
  fsetpos(FileToRead,&PosTemp);
  if(Gen->NNano != 0)
    Nano = (NANO *) realloc(Nano,Gen->NNano*sizeof(NANO));
  for(int n=0;n<Gen->NNano;n++){
    fgetpos(FileToRead,&PosTemp);
    fgets(cLine,STRSIZE,FileToRead);
    if( strncmp(cLine,"# Rigid",7)){// + strncmp(cLine,"# Pep",5)){
      fsetpos(FileToRead,&PosTemp);
      break;
    }
    NanoString(cLine,n);
  }
}
void VarData::ReadHeaderXvt(FILE *FileToRead){
  double Val[6];
  char cLine[STRSIZE];
  fgets(cLine,STRSIZE,FileToRead);
  fpos_t PosTemp;
  //-----------edge-nblock------------
  fgetpos(FileToRead,&PosTemp);
  double Time = 0.;
  int NBlock = 0;
  if( sscanf(cLine,"# L=%lf %lf %lf t=%lf blocks=%d",&Val[0],&Val[1],&Val[2],&Time,&NBlock) == 5 ){
    for(int d=0;d<3;d++) SetEdge(Val[d]*ScaleF[d],d);
    SetTime(Time);
    SetStep((int)(Time/0.05));
    SetNBlock(NBlock);
    VAR_ADD_TYPE(SysType,VAR_EDGE);
    for(int d=0;d<3;d++)
      Nano->Pos[d] = (.5 - ShiftPos[d])*pEdge(d);
    SetNanoBkf(0);
    Nano->Axis[0] = 0.;Nano->Axis[1] = 0.;Nano->Axis[2] = 1.;
  }
  //---------------virial-coefficients---------------
  //fgets(cLine,STRSIZE,FileToRead);
  double v[6];
  double w[10];
  int IfTwoType = 1;
  int IfThreeType = 1;
  fscanf(FileToRead,"# v=");
  for(int i=0;i<6;i++){
    if(fscanf(FileToRead,"%lf ",v+i) != 1){
      //fsetpos(FileToRead,&PosTemp);
      if(i == 2) IfThreeType = 0;
      break;
    }
    if(i >  2) IfTwoType = 0;
  }
  int NThird = 10;
  Gen->vBB = v[3];
  if(IfTwoType) NThird = 4;
  fscanf(FileToRead,"w=");
  for(int i=0;i<NThird;i++){
    if(fscanf(FileToRead,"%lf ",w+i) != 1){
      break;
    }
  }
  Gen->rho = 3.*(- .5*v[0] + sqrt(QUAD(v[0])*.25-8.*w[0]/3.))/(4.*w[0]);
  Gen->kappaN = - v[0] * Gen->rho*.5-3.;
  if(IfTwoType)
    Gen->chiN = (v[1] - .5*(v[0]+v[2]))*Gen->rho;
  else if(IfThreeType)
    Gen->chiN = (v[1] - .5*(v[0]+v[3]))*Gen->rho;
  //---------------intra-forces-------------------
  fgets(cLine,STRSIZE,FileToRead);
  double Nigot1,Nigot2,Nigot3;
  fgetpos(FileToRead,&PosTemp);
  int ret = sscanf(cLine,"# a2=%lf a3=%lf Re=%lf N=%lf ks=%lf kb=%lf l0=%lf",&Gen->WFuncStraight2,&Gen->WFuncStraight3,&Gen->ReOverCutOff,&Nigot3,&Gen->kappaSpring,&Gen->kappaBend,&Gen->SpringRest);
  if(ret == 7);
  else fsetpos(FileToRead,&PosTemp);
  SetCoeff(v,w);
  //----------inclusion------------------------
  double Pos[3];
  double Char[5];
  int NSide[2];
  double Axis[3];
  char Shape[10];
  char FileName[60];
  int NNanoTemp = 0;
  // how many Nano
  fgetpos(FileToRead,&PosTemp);
  Gen->NNano = 0;
  for(int t=0;t<50;t++){
    if(NULL == fgets(cLine,sizeof(cLine),FileToRead)) break;
    if(strstr(cLine, "# Rigid") == cLine)
      Gen->NNano++;
    else if(strstr(cLine, "# Pep") == cLine)
      Gen->NNano++;
    else 
      break;
  }
  fsetpos(FileToRead,&PosTemp);
  if(Gen->NNano != 0)
    Nano = (NANO *) realloc(Nano,Gen->NNano*sizeof(NANO));
  for(int n=0;n<Gen->NNano;n++){
    fgetpos(FileToRead,&PosTemp);
    fgets(cLine,STRSIZE,FileToRead);
    if( strncmp(cLine,"# Rigid",7)){
      fsetpos(FileToRead,&PosTemp);
      continue;
    }
    NanoString(cLine,n);
    NNanoTemp++;
  }
  //printf("%lf %lf %lf %s\n",Pos[0],Axis[0],Char[0],Shape);
  // Soft
  for(int n=NNanoTemp;n<Gen->NNano;n++){
    fgetpos(FileToRead,&PosTemp);
    fgets(cLine,STRSIZE,FileToRead);
    if( strncmp(cLine,"# Pep",4) ){
      fsetpos(FileToRead,&PosTemp);
      break;
    }
    if( Fetch(cLine,"g","%lf %lf %lf",Char,Char+1,Char+2)){
      Nano[n].Rad = Char[0];
      Nano[n].Height = Char[1];
      Nano[n].Hamaker = Char[2];
    }
    if( Fetch(cLine,"d","%d %d",NSide,NSide+1)){
      Nano[n].NCircle = NSide[0];
      Nano[n].NHeight = NSide[1];
    }
    if( Fetch(cLine,"fn","%s",FileName)){
      sprintf(Nano[n].ArchFile,"%s",FileName);
    }
    Nano[n].Shape = SHAPE_CLUSTER;
  }
}
//###############################PASS#THRU########################
//If the information for the allocation are missing
int VarData::ReadPassThru(FILE *FileToRead){
  VarMessage("PassThru");
  int Val[3];
  char cLine[STRSIZE];
  char cVar[STRSIZE];
  char cVal[STRSIZE];
  int Paren[2];
  int NChain = 0;
  int NPart = 0;
  int NPCh = 0;
  int NType = 0;
  int NLink = 0;
  if( SysFormat == VAR_SYS_TXVL ){
    for(int k=0;!(fgets(cLine,STRSIZE,FileToRead)==NULL);k++){
      if(cLine[0] == '#') continue;
      int iLen = (int) (strlen(cLine));
      if( Fetch(cLine,"t","%d %d %d",Val,Val+1,Val+2)){
	int p = Val[0];
	int c = Val[1];
	int t = Val[2];
	if(c >= NChain) NChain++;
	if(t >= NType) NType++;
      }
      if(cLine[0] == '{') NPart++;
      if(NChain==0) NPCh++;
      for(int i=0,link=0,part=0;i<iLen;i++){
	if(cLine[i] == 'l' && cLine[i+1] == '['){
	  link++;
	  if(link >= NLink) NLink++;
	}
      }
    }
    if(NChain==0) NChain = 1;
    Block[0].InitIdx = 0;
    Block[0].NChain = NChain;
    Block[0].NPCh = NPCh;
    Block[0].NPart = NPart;
    Block[0].EndIdx = NPart;
  }
  else if ( SysFormat == VAR_SYS_XVT ){
    for(int k=0,b=0;!(fgets(cLine,STRSIZE,FileToRead)==NULL);k++){
      if(3 == sscanf(cLine,"# n=%d N=%d name=%s",&Val[0],&Val[1],&Block[b].Name)){
	Block[b].InitIdx = NPart;
	Block[b].NChain = Val[0];
	Block[b].NPCh = Val[1];
	Block[b].NPart = Block[b].NChain*Block[b].NPCh;
	NChain += Block[b].NChain;
	NPart  += Block[b].NPart;
	Block[b].EndIdx = NPart;
	Block[b].Arch = ARCH_LINES;
	NLink = 1;
	if(strcasestr(Block[b].Name, "TT") == Block[b].Name){
	  Block[b].Arch = ARCH_TWOTAILS;
	  NLink = 2;
	}
	if(strcasestr(Block[b].Name, "PEP") == Block[b].Name){
	  Block[b].Arch = ARCH_CLUSTER;
	}
	//printf("Found block # %d name %s #chain %d #part %d #partPchain%d from %d to %d\n",b,Block[b].NChain,Block[b].Name,Block[b].NPart,Block[b].NPCh,Block[b].InitIdx,Block[b].EndIdx);
	b++;
	if(Gen->NBlock == b) break;
      }
    }
    for(int b=0,nNano=0;b<pNBlock();b++){
      if(VAR_IF_TYPE(Block[b].Arch,ARCH_CLUSTER)){
	for(int n=0;n<pNNano();n++,nNano=n){
	  if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CLUSTER)){
	    Nano[n].nBlock = b;
	  }
	}
      }
    }
    NPCh = Block[0].NPCh;
  }
  else{
    for(int k=0,b=0;!(fgets(cLine,STRSIZE,FileToRead)==NULL);k++){
      if(cLine[0] == '#' || cLine[0] == '$'){continue;}
      NPart++;
    }
    if(SysFormat == VAR_SYS_XYZ) NPCh = NPart;
  }
  SetNLink(NLink);
  SetNPart(NPart);
  SetNChain(NChain);
  SetNPCh(NPCh);
  SetNType(NType);
  rewind(FileToRead);
 return 0;
}
//####################READ#PART#################################3
int VarData::ReadPart(FILE *FileToRead){
  for(int d=0;d<3;d++) Gen->Cm[d] = 0.;  
  if( SysFormat == VAR_SYS_TXVL ){
    if(ReadPartTxvl(FileToRead)) return 1;
    Block[0].NPart = Gen->NPart;
    Block[0].NChain = Gen->NChain;
    Block[0].NPCh = Gen->NPCh;
    Block[0].InitIdx = 0;
    Block[0].EndIdx = Gen->NPart;
  }
  else if( SysFormat == VAR_SYS_XVT ){
    if(ReadPartXvt(FileToRead)) return 1;
  }
  else if( SysFormat == VAR_SYS_XYZT ){
    if(ReadPartXyzt(FileToRead)) return 1;
  }
  else {
    if(ReadPartXyz(FileToRead)) return 1;
  }
  for(int d=0;d<3;d++) Gen->Cm[d] /= (double)pNPart();
  if(pNNano() == 0){
    for(int d=0;d<3;d++){
      Nano->Pos[d] = .5*pEdge(d);
    }
    SetNanoBkf(0);
    Nano->Axis[0] = 0.;Nano->Axis[1] = 0.;Nano->Axis[2] = 1.;
  }
  return 0;
}
int VarData::ReadPartTxvl(FILE *FileToRead){
  VarMessage("ReadPartSys");
  int NChain = 0;
  int NPart = 0;
  Gen->NPCh = 0;
  int NType = 0;
  int NLink = 0;
  double Val[6];
  int Char[6];
  char cLine[STRSIZE];
  for(int d=0;d<4;d++) Gen->Vel[d] = 0.000001;
  for(int c=0;c<Gen->NChain;c++){
    for(int d=0;d<3;d++)
      Ch[c].Pos[d] = 0.;
  }
  Block[0].InitIdx = 0;
  Block[0].EndIdx = Gen->NPart;
  int NPCh = 0;
  for(int p=0,c=0;!(fgets(cLine,STRSIZE,FileToRead)==NULL);p++){
    //if(p >= Gen->NPart) break;
    int IfContinue = 1;
    int iLen = strlen(cLine);
    for(int k=0;k<iLen;k++){
      if(cLine[k] == '#'){
	p--;
	IfContinue = 0;
	break;
      }
    }
    if(!IfContinue) continue;
    if(p >= Gen->NPart){
      //printf("More particles than expected %d > %d\n",p,Gen->NPart);
      //return 1;
    }
    Pm[p].Idx=p;Pm[p].Typ=0;Pm[p].CId=0;
    int sPos = Fetch(cLine,"t","%d %d %d",Char,Char+1,Char+2);
    if(sPos){
      Pm[p].Idx = Char[0];
      Pm[p].CId = Char[1];
      Pm[p].Typ = Char[2];
    }
    if(p>0 && Pm[p].CId != Pm[p-1].CId){
      Ch[c].NPCh = NPCh;
      Ch[c].InitBead = p-NPCh;
      Ch[c].EndBead = p;
      NPCh = 0;
    }
    sPos += Fetch(cLine+sPos,"x","%lf %lf %lf",Pm[p].Pos,Pm[p].Pos+1,Pm[p].Pos+2);
    sPos += Fetch(cLine+sPos,"v","%lf %lf %lf",Pm[p].Vel,Pm[p].Vel+1,Pm[p].Vel+2);
    for(int d=0;d<3;d++){
      Pm[p].Pos[d] *= ScaleF[d];
      Ch[c].Pos[d] += Pm[p].Pos[d];
      Gen->Cm[d] += Pm[p].Pos[d];
      //Ch[c].Vel[d] += Pm[p].Vel[d];
      if(Gen->Vel[d] < Pm[p].Vel[d])
	Gen->Vel[d] = Pm[p].Vel[d];
    }
    if(Pm[p].CId > c){
      c++;
      if(c >= Gen->NChain){
	printf("More chains than expected %d > %d\n",c,Gen->NChain);
      }
    }
    if(Pm[p].CId == 0){Gen->NPCh++;}
    if(Gen->NType < Pm[p].Typ ){
      NType++;
      if(Pm[p].CId == 0)
	Block[0].Asym = p;
    }
    //Ch[c].NPart++;
    Pm[p].Vel[3] = sqrt( QUAD((Pm[p].Vel[0])) + QUAD((Pm[p].Vel[1])) +  QUAD((Pm[p].Vel[2])) );
    Gen->Vel[3] += sqrt( QUAD((Pm[p].Vel[0])) + QUAD((Pm[p].Vel[1])) +  QUAD((Pm[p].Vel[2])) );
    for(int l=0;l<Gen->NLink;l++){
      Ln[p].Link[l] = 0;
      int sPosOld = Fetch(cLine+sPos,"l",1,Val);
      //int sPosOld = Fetch(cLine+sPos,"l","%d",Val);
      if(sPosOld){
	Ln[p].Link[l] = (int)Val[0];
	Ln[p].NLink = l+1;
	sPos += sPosOld;
      }
    }
    //printf("%s{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf) l[%d] l[%d]}\n\n",p,cLine,Pm[p].Idx,Pm[p].CId,Pm[p].Typ,Pm[p].Pos[0],Pm[p].Pos[1],Pm[p].Pos[2],Pm[p].Vel[0],Pm[p].Vel[1],Pm[p].Vel[2],Pm[p].Link[0],Pm[p].Link[1]);
  }
  for(int c=0;c<Gen->NChain;c++)
    for(int d=0;d<3;d++)
      Ch[c].Pos[d] /= Gen->NPCh;
  Gen->Vel[3] /= (double)NPart;
  NType++;
  //  Gen->NChain++;
  if(Gen->NPCh == 0) Gen->NPCh = Gen->NPart;
  return 0;
}
int VarData::ReadPartXvt(FILE *FileToRead){
  VarMessage("ReadPartNoSys");
  double *buff = (double *) malloc(sizeof(double));
  char cLine[STRSIZE];
  sprintf(cLine,"               ");
  int NChain = 0;
  int NPart = 0;
  int NPCh = 0;
  int NType = 0;
  int Asym = 0;
  int Val[3];
  int NLink=0;
  fpos_t PosTemp;
  for(int b=0,NCh=0;b<pNBlock();NCh+=Block[b++].NChain){
    for(int t=0;t<50;t++){
      fgets(cLine,STRSIZE,FileToRead);
      if(3 == sscanf(cLine,"# n=%d N=%d name=%s",&Val[0],&Val[1],&Block[b].Name)){
	Block[b].InitIdx = NPart;
	Block[b].NChain = Val[0];
	Block[b].NPCh = Val[1];
	Block[b].NPart = Block[b].NChain*Block[b].NPCh;
	NChain += Block[b].NChain;
	NPart  += Block[b].NPart;
	Block[b].EndIdx = NPart;
	Block[b].Arch = ARCH_LINES;
	NLink = 1;
	if(strcasestr(Block[b].Name, "TT") == Block[b].Name){
	  Block[b].Arch = ARCH_TWOTAILS;
	  NLink = 2;
	}
	if(strcasestr(Block[b].Name, "PEP") == Block[b].Name){
	  Block[b].Arch = ARCH_CLUSTER;
	}
	break;
      }
    }
    // // while(strncmp(cLine,"# n=",3)){fgets(cLine,STRSIZE,FileToRead);printf("%s\n",cLine);}
    // printf("%s\n",cLine);
    for(int c=NCh;c<NCh+Block[b].NChain;c++){
      int pCurr = Block[b].InitIdx + (c-NCh)*Block[b].NPCh;
      Ch[c].NPCh = Block[b].NPCh;
      Ch[c].InitBead = pCurr;
      Ch[c].EndBead = pCurr + Block[b].NPCh;
      for(int ppc=0;ppc<Block[b].NPCh;ppc++){
	if(NULL == fgets(cLine,STRSIZE,FileToRead)) break;
	if(cLine[0] == '#'){ppc--;continue;}
	int p =  pCurr + ppc;
	if(p >= pNPart()){
	  printf("Reading more particles than the ones allocated\n");
	  exit(1);
	}
	Pm[p].Idx = p;
	Pm[p].CId = c;
	//sscanf(cLine,"%lf %lf %lf %lf %lf %lf %d",&Pm[p].Pos[0],&Pm[p].Pos[1],&Pm[p].Pos[2],&Pm[p].Vel[0],&Pm[p].Vel[1],&Pm[p].Vel[2],&Pm[p].Typ);
	int Incr = ReadVal(cLine  ,buff);Pm[p].Pos[0] = *buff*ScaleF[0];
	Incr += ReadVal(cLine+Incr,buff);Pm[p].Pos[1] = *buff*ScaleF[1];
	Incr += ReadVal(cLine+Incr,buff);Pm[p].Pos[2] = *buff*ScaleF[2];
	Incr += ReadVal(cLine+Incr,buff);Pm[p].Vel[0] = *buff;
	Incr += ReadVal(cLine+Incr,buff);Pm[p].Vel[1] = *buff;
	Incr += ReadVal(cLine+Incr,buff);Pm[p].Vel[2] = *buff;
	Incr += ReadVal(cLine+Incr,buff);Pm[p].Typ    = (int)*buff;
	for(int d=0;d<3;d++) Gen->Cm[d] += Pm[p].Pos[d];
	if(NType < Pm[p].Typ ) NType++;
      }
    }
  }
  for(int b=0;b<pNBlock();b++){
    int p1 = Block[b].InitIdx;
    for(int ppc=1;ppc<Block[b].NPCh;ppc++){
      //if(Pm[p1+ppc].Typ == 1 && Pm[p1+ppc-1].Typ == 0)
      if(Pm[p1+ppc].Typ != Pm[p1+ppc-1].Typ)
	Block[b].Asym = ppc;
    }
    if(Block[b].Asym == 0) Block[b].Asym = pNPCh();
  }
  // links
  if(pNLink() > 0){
    for(int b=0;b<pNBlock();b++){
      int NLink = 1;
      if(Block[b].Arch == ARCH_CLUSTER) NLink = 0;
      for(int c=0;c<Block[b].NChain;c++){ 
	int pCurr = Block[b].InitIdx + c*Block[b].NPCh;
  	for(int ppc = 0;ppc<Block[b].NPCh-1;ppc++){
	  Ln[pCurr+ppc].NLink = NLink;
	  Ln[pCurr+ppc].Link[0] = pCurr+ppc+1;
	  if(Block[b].Arch == ARCH_TWOTAILS){
	    if(ppc == Block[b].Asym/2 - 1){
	      Ln[pCurr+ppc].NLink = 1;
	      Ln[pCurr+ppc].Link[0] = pCurr + Block[b].Asym - 1;
	    }
	    // if(ppc == Block[b].Asym){
	    //   Ln[pCurr+ppc].NLink = 1;
	    //   Ln[pCurr+ppc].Link[1] = pCurr + Block[b].Asym + 2;
	    // }
	  }
	}
      }
    }
  }
  //  for(int c=0;c<Gen->NChain;c++)printf("%d %d\n",c,Ch[c].NPart);
  double Norm2 = CUBE(pReOverCutOff()) / SQR(pNPCh());
  double Norm3 = CUBE(SQR(pReOverCutOff()) / pNPCh());
  MInt->Rescale(Norm2,2);
  MInt->Rescale(Norm3,3);
  free(buff);
  NType++;
  if(NPCh == 0) NPCh = Gen->NPart;
  SetNType(NType);
  if(0 == pNPart()){
    Block[0].InitIdx = 0;
    Block[0].NChain = 0;
    Block[0].NPCh = 0;
    Block[0].NPart = 0;
    Block[0].EndIdx = 0;
  }
  return 0;
}
int VarData::ReadLineXvt(char *cLine,double *Pos,int *Type){
  double Temp = 0.;
  int Incr = ReadVal(cLine  ,&Temp);Pos[0] = Temp*ScaleF[0];
  Incr += ReadVal(cLine+Incr,&Temp);Pos[1] = Temp*ScaleF[1];
  Incr += ReadVal(cLine+Incr,&Temp);Pos[2] = Temp*ScaleF[2];
  Incr += ReadVal(cLine+Incr,&Temp);
  Incr += ReadVal(cLine+Incr,&Temp);
  Incr += ReadVal(cLine+Incr,&Temp);
  Incr += ReadVal(cLine+Incr,&Temp);*Type   = (int)Temp;
  return 0;
}
int VarData::ReadPartXyz(FILE *FileToRead){
  double *buff = (double *) malloc(sizeof(double));
  char cLine[STRSIZE];
  double Pos[3];
  for(int p=0,c=-1;!(fgets(cLine,STRSIZE,FileToRead)==NULL);p++){
    if(p >= pNPart()-1) continue;
    Pm[p].Idx=p;Pm[p].Typ=0;Pm[p].CId=0;
    if(cLine[0] == '#' || cLine[0] == '$'){p--;continue;}
    sscanf(cLine,"%lf %lf %lf",Pos,Pos+1,Pos+2);
    for(int d=0;d<3;d++)
      Pm[p].Pos[d] = Pos[d];
    // int Incr = ReadVal(cLine,buff)  ;Pm[p].Pos[0] = *buff;
    // Incr += ReadVal(cLine+Incr,buff);Pm[p].Pos[1] = *buff;
    // Incr += ReadVal(cLine+Incr,buff);Pm[p].Pos[2] = *buff;
  }
  BfEdge();
  if(IfNormalize){
    for(int d=0;d<3;d++){
      ScaleF[d] = pInvEdge(d);
      SetEdge(1.,d);
    }
  }
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Pm[p].Pos[d] *= ScaleF[d];
    }
  }
  SetNNano(0);
  Block[0].InitIdx = 0;
  Block[0].EndIdx = pNPart();
  Block[0].NChain = 1;
  Block[0].NPCh = pNPart();
  free(buff);
  return 0;
}
int VarData::ReadPartXyzt(FILE *FileToRead){
  double *buff = (double *) malloc(sizeof(double));
  char cLine[STRSIZE];
  double Pos[4];
  int NChain = 0;
  int NPCh = 0;
  int OldTyp = 0;
  for(int p=0,c=-1;!(fgets(cLine,STRSIZE,FileToRead)==NULL);p++){
    if(p >= pNPart()-1) continue;
    Pm[p].Idx=p;Pm[p].Typ=0;Pm[p].CId=0;
    if(cLine[0] == '#' || cLine[0] == '$'){p--;continue;}
    sscanf(cLine,"%lf %lf %lf %lf",Pos,Pos+1,Pos+2,Pos+3);
    for(int d=0;d<3;d++)
      Pm[p].Pos[d] = Pos[d];
    Pm[p].Typ = (int)Pos[3];
    if(OldTyp > Pm[p].Typ){
      NChain++;
    }
    Pm[p].CId = NChain;
    if(NChain == 0) NPCh++;
    OldTyp = Pm[p].Typ;
  }
  BfEdge();
  if(IfNormalize){
    for(int d=0;d<3;d++){
      ScaleF[d] = pInvEdge(d);
      SetEdge(1.,d);
    }
  }
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Pm[p].Pos[d] *= ScaleF[d];
    }
  }
  if(NChain == 0) NChain = 1;
  SetNChain(NChain);
  SetNPCh(NPCh);
  Block[0].InitIdx= 0;
  Block[0].EndIdx = pNPart();
  Block[0].NChain = pNChain();
  Block[0].NPCh   = pNPCh();
  free(buff);
  return 0;
}
