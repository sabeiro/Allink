#include "Segnali.h"
  int NSegn=1000;
  int Periodi=10000;
  double Passo = .5;
  double MediaSu=1.;
  double ScartoSu=1.;
  double MediaGiu=30.;
  double ScartoGiu=1.;
  double Incr=.01;
int main(int argc,char **argv){
  char *NomeMedia;
  char *Comando;
  NomeMedia = new char[30];
  Comando = new char[50];
  Legenda();
  Variabili *v1 = new Variabili(NSegn,Periodi,MediaSu,ScartoSu,MediaGiu,ScartoGiu,Passo,Incr);
  for(int i=0;i<argc;i++){
    if(!strcmp(*(argv+i),"-help")){
      return 0;
    }
    if(!strcmp(*(argv+i),"-sim")){
      (void)time(&ti);
      v1->Simula();
      (void)time(&tf);
      printf("Il tutto dopo %d [s] ovvero %.1f [min]\n",
	     (int)(tf-ti),(double)(tf-ti)/60.);
      //      delete v1;
      printf("Te ve be te e`?\n");
      return 0;
    }    
  }
  while(strcmp(Comando,"esci") ){
    printf("Segnali> ");
    scanf("%s",Comando);
    if(!strcmp(Comando,"num")){
      printf("Introdurre il numero di segnali: ");
      scanf("%d",&NSegn);
      v1->Cambia(Comando,NSegn);
    }
    else if(!strcmp(Comando,"?") ){
      Legenda();
    }
    else if(!strcmp(Comando,"!") ){
      printf("Introdurre il comando per la shell");
      scanf("%s",Comando);
      system(Comando);
    }
    else if(!strcmp(Comando,"per")){
      printf("Introdurre i periodi: ");
      scanf("%d",&Periodi);
      v1->Cambia(Comando,Periodi);
    }
    else if(!strcmp(Comando,"pas")){
      printf("Introdurre il passo di campionamento: ");
      scanf("%lf",&Passo);
      v1->Cambia(Comando,Passo);
    }
    else if(!strcmp(Comando,"mSu") ){
      printf("Introdurre la media: ");
      scanf("%lf",&MediaSu);
      v1->Cambia(Comando,MediaSu);
    }
    else if(!strcmp(Comando,"sSu") ){
      printf("Introdurre lo scarto: ");
      scanf("%lf",&ScartoSu);
      v1->Cambia(Comando,ScartoSu);
    }
    else if(!strcmp(Comando,"mGiu") ){
      printf("Introdurre la media: ");
      scanf("%lf",&MediaGiu);
      v1->Cambia(Comando,MediaGiu);
    }
    else if(!strcmp(Comando,"sGiu") ){
      printf("Introdurre lo scarto: ");
      scanf("%lf",&ScartoGiu);
      v1->Cambia(Comando,ScartoGiu);
    }
    else if(!strcmp(Comando,"incr") ){
      printf("Introdurre l'incremento: ");
      scanf("%lf",&Incr);
      v1->Cambia(Comando,Incr);
    }
    else if(!strcmp(Comando,"sim") ){
      (void)time(&ti);
      v1->Simula();
      (void)time(&tf);
      printf("Il tutto dopo %d [s] ovvero %.1f [min]\n",
	     (int)(tf-ti),(double)(tf-ti)/60.);
    }
    else if(!strcmp(Comando,"esci") );
    else {
      printf("Comando non valido, scrivere ?\
 per la lista dei comandi\n");
    }
  }
  //  delete v1;
  printf("Te ve be te e`?\n");
  return 0;
}
void Legenda(){
  printf("\n\
*****************************************************************\n\
(                                                               )\n\
)    Qui inizia il programma che pasticcia con i segnali        (\n\
(    Mettetevi comodi che saran ore di conti                    )\n\
)                                                               (\n\
(-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-)\n\
(   num %%d imposta in numero di segnali  (pred=%d)              )\n\
)   per %%d impsta per il numero di periodi    (pred=%d)         (\n\
(   pas %%lf imposta il numero di passi (%.1f)                   )\n\
(   mSu %%lf imposta la media dei segnali  (%.1f)                )\n\
)   sSu %%lf imposta lo scarto su cui son distribuiti (%.1f)     (\n\
(   mGiu%%lf imposta la media dei segnali  (%.1f)                )\n\
)   sGiu%%lf imposta lo scarto su cui son distribuiti (%.1f)     (\n\
)   incr%%lf imposta lo l'incremento (%.1f)                      (\n\
(   sim simula                                                  )\n\
(   ! esegue un comando (come cambiare cartella)                )\n\
(   esci ...                                                    )\n\
*****************************************************************\n",
	 NSegn,Periodi,Passo,MediaSu,ScartoSu,MediaGiu,ScartoGiu,Incr);
}
    
