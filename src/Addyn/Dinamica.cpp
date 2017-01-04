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
#include "Forces.h"
#include <signal.h>
#include <fenv.h>
void Legenda();
void Particle();
void DrawParticle();
void keyboard(unsigned char key,int x, int y);
void InitConstant();
Forces *Fm;
static volatile int IfExit = 0;
static int IfUser = 1;

static void StopLoop(int sig)
{
  IfExit = 1;
  printf("Loop stopped\n");
}
int main(int argc,char **argv){
  char Comando[50];
  char ConfFile[256];
  signal(SIGTSTP,StopLoop);
  signal(SIGSTOP,StopLoop);
  signal(SIGSTOP,StopLoop);
  //  frame = (char *)malloc(100*sizeof(char));
  //SpecFuntore <Matematica> SpecMat(&Mat,Matematica::Eval);
  //SpecFuntore <ElPoly> SpecPol(&Pol,ElPoly::ContactAngle);
  if(argc <2){ printf("ConfFile not specified\n"); return 1;}
  int NEdge = 20+1;
  //Fm->Solve();
  int *FilePos = (int *)calloc(argc,sizeof(int));
  int NFile = 0;
  char Snapshot[60];
  int InitFile = 0;
  //  ScaleUn = Fm->Gen->Edge[0];
  //Principal(argc,argv);
  for(int i=1;i<argc;i++){
    if(argv[i][0] != '-'){
      FilePos[NFile] = i;
      NFile++;
    }
    else if(!strcmp(*(argv+i),"-d")){
      sprintf(Comando,"d");
      IfUser = 0;
    }
    else if(!strcmp(*(argv+i),"-r")){
      sprintf(Comando,"r");
      IfUser = 0;
    }
    else if(!strcmp(*(argv+i),"-s")){
      sprintf(Comando,"s");
      IfUser = 0;
    }
    else if(!strcmp(*(argv+i),"-c")){
      if(argc > i){
	sprintf(ConfFile,"%s",argv[i+1]);
	i++;
      }
      else {
	printf("Missing config file\n");
      }
    }
    else if(!strcmp(*(argv+i),"-i")){
      if(argc > i){
	sprintf(Snapshot,"%s",*(argv+i+1));
	InitFile = i + 1;
      }
      else{
	printf("Snapshot file missing\n");
      }
    }
    else if(!strncmp(*(argv+i),"--",2)){
      sprintf(Comando,"%s",argv[i]+2);
      IfUser = 0;
    }
  }
  if(InitFile)
    Fm = new Forces(argc,argv,ConfFile,argv[InitFile]);
  else
    Fm = new Forces(argc,argv,NEdge,ConfFile);
  fedisableexcept(FE_ALL_EXCEPT);
  char cSystem1[STRSIZE];
  char cSystem2[STRSIZE];
  Fm->SysInfo(cSystem1);
  Fm->SysDef(cSystem2);
  printf("------------------------------------------------------\n");
  printf("%s NFile %d\n%s\n",cSystem1,NFile,cSystem2);
  printf("------------------------------------------------------\n");
  while(strcmp(Comando,"q")){
    if(IfUser){
      printf("Dinamica> ");
      scanf("%s",Comando);
    }
    if(!strcmp(Comando,"d")){
#ifdef __glut_h__
      Fm->Graphics(argc,argv);
#else 
      printf("Graphics libraries (GL/glut.h) not supplied\n");
#endif //__glut_h__
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"s")){
      Fm->RunDynamics();
      if(!IfUser) return 0;      
    }
    else if(!strcmp(Comando,"r")){
      Fm->RunDynamics();
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"nrg")){
      for(int f=0;f<NFile;f++){
	fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
	Fm->CalcTotNrg(argv[FilePos[f]],f);
      }
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"PepNrg")){
      for(int f=0;f<NFile;f++){
	fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
	Fm->CalcNrgPep(argv[FilePos[f]],f);
      }
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"ExpPep")){
      Fm->ExplorePepSize();
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"ExpPep2d")){
      Fm->ExplorePepSize2d();
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"ExpDist")){
      Fm->ExploreDoubleMin();
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"Widom")){
      for(int f=0;f<NFile;f++){
	fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
	Fm->RunWidom(argv[FilePos[f]],f);
      }
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"WidomChIn")){
      for(int f=0;f<NFile;f++){
	fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
	Fm->RunWidomChIn(argv[FilePos[f]],f);
      }
      printf("\n");
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"WidomChOut")){
      for(int f=0;f<NFile;f++){
	fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
	Fm->RunWidomChOut(argv[FilePos[f]],f);
      }
      printf("\n");
      if(!IfUser) return 0;
    }
    else if(!strcmp(Comando,"Rosenbluth")){
      char File2Write[60];
      for(int f=0;f<NFile;f++){
	sprintf(File2Write,"Rosenbluth%05d.dat",f);
	FILE *WidomIn = fopen(File2Write,"w");
	fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
	Fm->ReOpen(argv[FilePos[f]],BF_PART);
	Fm->RosenIn(WidomIn);
	fclose(WidomIn);
      }
      printf("\n");
      if(!IfUser) return 0;
   }
    else if(!strcmp(Comando,"RosenOut")){
      char File2Write[60];
      for(int f=0;f<NFile;f++){
	sprintf(File2Write,"RosenOut%05d.dat",f);
	FILE *WidomIn = fopen(File2Write,"w");
	fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
	Fm->ReOpen(argv[FilePos[f]],BF_PART);
	Fm->RosenOut(WidomIn);
	fclose(WidomIn);
      }
      printf("\n");
      if(!IfUser) return 0;
   }
   else if(!strcmp(Comando,"Tens")){
     Fm->CalcTens(argv,FilePos,NFile);
     if(!IfUser) return 0;
   }
   else if(!strcmp(Comando,"AvForces")){
     Fm->AvForces(argv,FilePos,NFile);
     if(!IfUser) return 0;
   }
   else if(!strcmp(Comando,"Trial")){
     Fm->Trial();
     if(!IfUser) return 0;
   }
    else if(!strcmp(Comando,"?") ){
      Legenda();
    }
    else if(!strcmp(Comando,"!") ){
      printf("Insert a shell command: ");
      scanf("%s",Comando);
      system(Comando);
    }//-----------------System-Info---------------------------
    else if(!strcmp(Comando,"info")){
      printf("------------------------------------------------------\n");
      //printf("%s\n",Fm->SysInfo());
    }
    else if(!strcmp(Comando,"write") ){
      Fm->Write("ciccia.dat");
    }
    else if(!strcmp(Comando,"q") );
    else {
      IfUser = 1;
      printf("Comando non valido, scrivere ? per la lista dei comandi\n");
    }
  }
  printf("Te se qe te ve be te ne?\n");
  return 0;
}
#ifdef __glut_h__
void ParticleList(){
  Fm->DrawScene();
}
// void ParticleRealTime(){
//   //Fm->DrawNano();
// }
// void InitConstant(){
//   Fm->InitConstant();
// }
void keyboard(unsigned char key,int x, int y){
  Fm->keyboard(key,x,y);
}
void Menu(){
  Fm->Menu();
}
void DynamicsMotion(){
  Fm->DynamicsView();
}
#endif
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
(   threshold set how many particle per cell will be considered )\n\
(   surf ratio between the actual surcafe and the circular      )\n\
(   core three dimentional sampled average of the system        )\n\
(   diff diffusivity of the extra particle                      )\n\
(   > <  open the file form the list position                   )\n\
(   cm   movement  of the center of mass of the system          )\n\
(   file number of file to use from the list                    )\n\
(   info ...                                                    )\n\
(   chbf change backfold                                        )\n\
(   val  number of values for the density profile               )\n\
*****************************************************************\n");
}
