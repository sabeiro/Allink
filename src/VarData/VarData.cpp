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
#include <stdarg.h>

void VarData::VarMessage(const char * s, ...)
{
#ifdef VAR_DEBUG
  va_list args;
  va_start(args, s);
  fprintf(stderr,"VarData] ");
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
#else
  return;
#endif
}
VarData::VarData(){
  VarMessage("Constructing");
  Mat = new Matematica();
  MInt = new MatInt(3,3);
  Gen = (GENERAL *) calloc(1,sizeof(GENERAL));
  Nano = (NANO *) calloc(1,sizeof(NANO));
  Soft = (SOFT *) calloc(1,sizeof(SOFT));
  Gen->NBlock = 1;
  Block = (BLOCK *)calloc(Gen->NBlock,sizeof(BLOCK));
  if( !Gen | !Mat ){printf("Il costruttore non alloca\n");return;}
  sprintf(cWhat2Draw,"part");
  CNorm = 2;
  CLat1 = (CNorm+1)%3;
  CLat2 = (CNorm+2)%3;
  SetNPCh(32); //# beads
  for(int d=0;d<3;d++) ShiftPos[d] = 0.;
  for(int d=0;d<3;d++) ScaleF[d] = 1.;
  Gen->NLink = 0;//max # of links
  NAddChain = 0;//added homopolymer chains
  NAddChol = 0;
  NSolvent = 0;
  NStuffing = 0;
  Gen->Deltat = 0.01;
  Gen->chiN = 30.;
  Gen->kappaN = 100.;
  Gen->rho = 50.;
  Gen->vBB = -0.1;
  Gen->NChain = 1;
  Gen->Temp = 1.;
  Gen->Edge[0] = 16.;
  Gen->Edge[1] = 16.;
  Gen->Edge[2] = 32.;
  Gen->NBlock = 1;
  Gen->WFuncStraight2 = .9;
  Gen->WFuncStraight3 = 1.0;
  SysType = 0;
  SysFormat = 0;
  SysCreate = 0;
  VAR_ADD_TYPE(SysCreate,VAR_ABSORBING);
  NPartNearSphere=0;
  NChType = CHAIN_EVERY;
  NPType = BEAD_EVERY;
  NEdge = 200;
  IfNormalize = 0;
  IfPlotMem = 0;
  Gen->NNano = 1;
  NSoft = 0;
  Pm = NULL;Ch = NULL;
  Point2Shape(SHAPE_NONE);
}
VarData::~VarData(){
  VarMessage("Destroying");
  delete Mat;
  free(Gen);
  free(Nano);
  free(Block);
  free(Soft);
  free(Ch);
  for(int p=0;p<Gen->NPart;p++)
    free(Ln[p].Link);
  free(Pm);
  free(Ln);
}
bool VarData::Open(char *InFile,int Bf){
  VarMessage("Open");
  Gen->NType = 0;
  for(int d=0;d<3;d++)
    Gen->Edge[d]=0.;
  VAR_REM_TYPE(SysType,VAR_EDGE);
  // VAR_REM_TYPE(SysType,VAR_CHAIN_DEF);
  FILE *FileToRead;
  if((FileToRead = fopen(InFile,"r"))==0){
    printf("The file is missing %s \n",InFile);
    return 1;
  }
  ReadHeader(FileToRead);
  if( !VAR_IF_TYPE(SysType,VAR_OPEN_TRUST) ) 
    if(ReadPassThru(FileToRead)) return 1;
  ReadPart(FileToRead);
  //if(AllocPart()) return 1;
  ShiftRef(Bf);
  fclose(FileToRead);
  return 0;
}
bool VarData::OpenRisk(char *InFile,int Bf){
  FILE *FileToRead;  
  if((FileToRead = fopen(InFile,"r"))==0){
    printf("The file is missing %s \n",InFile);
    return 1;
  }
  ReadHeader(FileToRead);
  ReadPart(FileToRead);
  ShiftRef(Bf);
  fclose(FileToRead);
  return 0;
}
Properties Properties::operator+ (const Properties& Pr) const{
  Properties Result;
  Result.RePhob = (this->RePhob + Pr.RePhob);
  Result.RePhil = (this->RePhil + Pr.RePhil);
  Result.VolPhob = (this->VolPhob + Pr.VolPhob);
  Result.VolPhil = (this->VolPhil + Pr.VolPhil);
  Result.FactPhob = (this->FactPhob + Pr.FactPhob);
  Result.FactPhil = (this->FactPhil + Pr.FactPhil);
  Result.GyrPhob = (this->GyrPhob + Pr.GyrPhob);
  Result.GyrPhil = (this->GyrPhil + Pr.GyrPhil);
  Result.ChDiff = (this->ChDiff + Pr.ChDiff);
  return Result;
}
Properties Properties::operator* (const double& c) const{
  Properties Result;
  Result.RePhob = (this->RePhob * c);
  Result.RePhil = (this->RePhil * c);
  Result.VolPhob = (this->VolPhob * c);
  Result.VolPhil = (this->VolPhil * c);
  Result.FactPhob = (this->FactPhob * c);
  Result.FactPhil = (this->FactPhil * c);
  Result.GyrPhob = (this->GyrPhob * c);
  Result.GyrPhil = (this->GyrPhil * c);
  Result.ChDiff = (this->ChDiff * c);
  return Result;
}
void VarData::SysInfo(char *cSystem){
  sprintf(cSystem,"# NPart %d NType %d NChain %d NPCh %d Asymmetry %d Edge: [%.3g %.3g %.3g] rho %.0f chiN %.0f kappaN %.0f Nano %d (# %d) (%lf %lf %lf) NBlock %d",Gen->NPart,Gen->NType,Gen->NChain,Gen->NPCh,Block[0].Asym,Gen->Edge[0],Gen->Edge[1],Gen->Edge[2],Gen->rho,Gen->chiN,Gen->kappaN,Nano->Shape,Gen->NNano,Nano->Rad,Nano->Hamaker,Nano->Height,Gen->NBlock);
}
void VarData::SysDef(char *cSystem){
  sprintf(cSystem,"# SysType) Edge %d Sys (Txvl %d Xvt %d) \n",
	  VAR_IF_TYPE(SysType,VAR_EDGE),SysFormat==VAR_SYS_TXVL,SysFormat==VAR_SYS_XVT);
}
void VarData::SetCoeff(){
  double v2[6];
  double v3[10];
  v2[0] = -2. * (Gen->kappaN+3.)/Gen->rho;//00
  v2[3] = Gen->vBB;//11
  v2[1] = Gen->chiN/Gen->rho + .5*(v2[0]+v2[2]);//01
  v2[2] = v2[1];//02
  v2[4] = v2[3];//12
  if( 1 == 1){//Hydrophobic
    v2[2] = 1.*v2[0];//02
    v2[4] = v2[1];//12
  }
  v2[5] = 0.;//22
  v3[0] = 1.5 * (Gen->kappaN+2.)/QUAD(Gen->rho);//000
  v3[1] = v3[0];//001
  v3[2] = v3[0];//002
  v3[3] = v3[0];//011
  v3[4] = v3[0];//012
  v3[5] = v3[0];//022
  v3[6] = 0.;//111
  v3[7] = v3[0];//112
  v3[8] = v3[0];//122
  v3[9] = 0.;//222
  SetCoeff(v2,v3);
}
void VarData::SetCoeff(double *v2,double *v3){
  MInt->SetCoeff(v2[0],0,0);
  MInt->SetCoeff(v2[1],0,1);
  MInt->SetCoeff(v2[2],0,2);
  MInt->SetCoeff(v2[3],1,1);
  MInt->SetCoeff(v2[4],1,2);
  MInt->SetCoeff(v2[5],2,2);
  MInt->SetCoeff(v3[0],0,0,0);
  MInt->SetCoeff(v3[1],0,0,1);
  MInt->SetCoeff(v3[2],0,0,2);
  MInt->SetCoeff(v3[3],0,1,1);
  MInt->SetCoeff(v3[4],0,1,2);
  MInt->SetCoeff(v3[5],0,2,2);
  MInt->SetCoeff(v3[6],1,1,1);
  MInt->SetCoeff(v3[7],1,1,2);
  MInt->SetCoeff(v3[8],1,2,2);
  MInt->SetCoeff(v3[9],2,2,2);
}
MatInt::MatInt(int NTypeExt,int NOrdExt){
  NOrd = NOrdExt;
  NType = NTypeExt;
  if(NOrd > 1){
    int NMax = IntType(NType-1,NType-1)+1;
    IntMatr2 = new double[NMax];
  }
  if(NOrd > 2){
    int NMax = IntType(NType-1,NType-1,NType-1)+1;
    IntMatr3 = new double[NMax];
  }
}
void MatInt::FillEntries(double *Matr,int Ord){
  if(Ord == 2){
    for(int t=0;t<NType;t++){
      for(int tt=t;tt<NType;tt++){
	int ttt = IntType(t,tt);
	IntMatr2[ttt] = Matr[ttt];
      }
    }
  }
  if(Ord == 3){
    for(int t=0;t<NType;t++){
      for(int tt=t;tt<NType;tt++){
	for(int ttt=tt;ttt<NType;ttt++){
	  int tf = IntType(t,tt,ttt);
	  IntMatr3[tf] = Matr[tf];
	}
      }
    }
  }
}
double MatInt::Coeff(int t1,int t2){
  return IntMatr2[IntType(t1,t2)];
}
double MatInt::Coeff(int t1,int t2,int t3){
  return IntMatr3[IntType(t1,t2,t3)];
}
void MatInt::SetCoeff(double Co,int t1,int t2){
  IntMatr2[NType*t1+t2] = Co;
  IntMatr2[NType*t2+t1] = Co;
  //IntMatr2[IntType(t1,t2)] = Co;
}
void MatInt::SetCoeff(double Co,int t1,int t2,int t3){
  IntMatr3[(t1*NType+t2)*NType+t3] = Co;
  IntMatr3[(t1*NType+t3)*NType+t2] = Co;
  IntMatr3[(t2*NType+t1)*NType+t3] = Co;
  IntMatr3[(t2*NType+t3)*NType+t1] = Co;
  IntMatr3[(t3*NType+t2)*NType+t1] = Co;
  IntMatr3[(t3*NType+t1)*NType+t2] = Co;
  //  IntMatr3[IntType(t1,t2,t3)] = Co;
}
int MatInt::IntType(int t1,int t2){
  return NType*t1+t2;
  return NType*MIN(t1,t2)+MAX(t1,t2);
  int t = 0;
  for(int ta=0;ta<MIN(t1,t2);ta++){
    for(int tb=ta;tb<MAX(t1,t2);tb++){
      t++;
    }
  }
  return t;
}
int MatInt::IntType(int t1,int t2,int t3){
  return (t1*NType+t2)*NType+t3;
  int tTemp = 0;
  int tSort[3] = {t1,t2,t3};
  for(int d=1;d<3;d++){
    if(tSort[d] < tSort[d-1]){
      tTemp = tSort[d];
      tSort[d] = tSort[d-1];
      tSort[d-1] = tTemp;
    }
  }
  return (NType*tSort[0]+tSort[1])*NType+tSort[2];
}
void MatInt::Print(){
  printf(" |_");
  for(int tc=0;tc<NType;tc++)
    printf("_____%d_____",tc);
  printf("|\n");
  for(int tr=0;tr<NType;tr++){
    printf("%d| ",tr);
    for(int tc=0;tc<NType;tc++){
      printf("%lf ",Coeff(tr,tc));
    }
    printf("|\n");
  }
  printf("\n");
  for(int t=0;t<NType;t++){
    for(int tt=t;tt<NType;tt++){
      for(int ttt=tt;ttt<NType;ttt++){
	int tf = IntType(t,tt,ttt);
	printf("%d%d%d) %lf\n",t,tt,ttt,Coeff(t,tt,ttt));
      }
    }
  }
}
void MatInt::Rescale(double SFactor,int Order){
  if(Order == 2){
    for(int t1=0;t1<NType;t1++)
      for(int t2=0;t2<NType;t2++)
	IntMatr2[(t1*NType+t2)] *= SFactor;
  }
  if(Order == 3){
    for(int t1=0;t1<NType;t1++)
      for(int t2=0;t2<NType;t2++)
	for(int t3=0;t3<NType;t3++)
	  IntMatr3[(t1*NType+t2)*NType+t3] *= SFactor;
  }
}
//memchr in wich position a caracter is
//memcmp how distant are to pointers
//memcpy copy n bytes from the second pointer to the first 
//memmove the same but uses a buffer and the memory can overlap
//memset set the value of the n character in a string
//strcat append a value to a string
//strchr says to a pointer where is a character
//strcmp compare two strings
//strcoll compare two strings avoiding null-character
//strcopy copy two strings
//strcspn find the position in wich one char of the second string appear in the first
//strncat append n char to the string
//strncmp compare n char from the begining
