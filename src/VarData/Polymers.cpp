/***********************************************************************
Polymers: program that uses VarData's function for creating an initial 
system of bonded monomers or for modifing a given system. 
Every parameter of the VarData class are here set.
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
#define COMM_STRETCH 1
#define COMM_REMOVE 2

int main(int argc,char **argv){
  char nome1[120];
  char nome2[120];
  char ConfF[60];
  sprintf(ConfF,"Polymers.conf");
  VarData *Dat = new VarData();
  int *FilePos = (int *)calloc(argc,sizeof(int));
  int NFile = 0;
  double StretchFact = 1.;
  int Comando = 0;
  double RefPos[3] = {.0,.0,.0};
  sprintf(nome1,"ciccia.dat");
  sprintf(nome2,"resume.dat");
  for(int i=1;i<argc;i++){
    if(argv[i][0] != '-'){
      FilePos[NFile] = i;
      NFile++;
    }
    else if(!strcmp(*(argv+i),"--Stretch")){
      if(i<argc){
	sscanf(*(argv+i+1),"%lf",&StretchFact);
	i+=1;
	Comando = COMM_STRETCH;
      }
    }
    else if(!strcmp(*(argv+i),"--Rem")){
      if(i<argc){
	Comando = COMM_REMOVE;
      }
    }
    else if(!strcmp(*(argv+i),"--Shift")){
      if(i<argc-2){
	for(int d=0;d<3;d++){
	  sscanf(*(argv+d+i+1),"%lf",RefPos+d);
	}
	i+=3;
      }
    }
    else if(!strcmp(*(argv+i),"-c")){
      if(i<argc){
	sscanf(*(argv+i+1),"%s",ConfF);
	i+=1;
      }
    }    
  }
  if(NFile == 1)
    sprintf(nome2,argv[FilePos[0]]);
  if(NFile == 2){
    sprintf(nome1,argv[FilePos[0]]);    
    sprintf(nome2,argv[FilePos[1]]);
  }
  printf("Apro %s e %s\n",nome1,nome2);
  if(NFile == 2){
    Dat->Open(nome1,BF_CHAIN);
    char cSystem[STRSIZE];
    Dat->SysInfo(cSystem);
    if(Comando == COMM_STRETCH){
      for(int d=0;d<2;d++)
	Dat->Gen->Edge[d] = StretchFact*Dat->Gen->Edge[d];
      for(int p=0;p<Dat->Gen->NPart;p++){
	for(int d=0;d<3;d++)
	  Dat->Pm[p].Pos[d] = StretchFact*Dat->Pm[p].Pos[d];
	Dat->Pm[p].Pos[2] = Dat->Pm[p].Pos[2] - 1.5;
      }
      for(int n=0;n<Dat->NNano;n++){
	for(int d=0;d<2;d++)
	  Dat->Nano[n].Pos[d] *= StretchFact;
	Dat->Nano[n].Rad *= StretchFact;
	Dat->Nano[n].Height *= StretchFact;
      }
    }
    if(Comando == COMM_REMOVE){
      int NChain = Dat->Block[0].NChain;
    }
    Dat->Write(nome2);  
  }
  else {
    Dat->DefSoft(nome2,ConfF);
  }
  printf("Te se qe te ve be ne?\n");
  return 0;
}

