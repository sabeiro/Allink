#include "Rinato.h"
// int *leArray(int & numEls) {
//   int tam=0; numEls=0;
//   int valorEnt;
//   int *arrayN=NULL;
//   cout << "Escreva valores inteiros: ";
//   while (cin >> valorEnt) {
//       if (numEls==tam) {
//              int *original = arrayN;
//              arrayN = new int[tam*2+1];
//              for (int i=0; i<tam; i++)
//                     arrayN[i] = original[i];
//              delete [] original;
//              tam = tam*2+1;
//       }
//       arrayN[numEls++] = valorEnt;   }
//   return arrayN;
// }
int main(int argc,char **argv){
  int Lato[3]={128,128,128};
  int Partizioni[3]={16,16,16};
  int Tipo=0;
  int ATipo=0;
  int Quale=0;
  int modo=0;
  unsigned long int tempo=1;
  FILE *Posizioni;
  char *file = new char[40];
  Particella *n[2];
  time_t ti=0,tf=0;
  (void)time(&ti);
  n[0] = new Particella("Interst.dat","ValInterst.dat",10000);
  n[1] = new Particella("Vacanze.dat","ValVacanze.dat",10000);
  Cubo *c1 = new Cubo(Lato,Partizioni);
  for(int i=0;i<argc;i++){
    if(!strcmp(*(argv+i),"-comm")){
      if(argc>=i+1){
	if(!strcmp(*(argv+i+1),"1"))
	  comm=1;
	else if(!strcmp(*(argv+i+1),"2"))
	  comm=2;
	else if(!strcmp(*(argv+i+1),"3"))
	  comm=3;
	else comm=0;
      }
    }
  }
  printf("Va la be?\n");
  c1->Posiziona(*n[0]);
  c1->Posiziona(*n[1]);
  Tipo=c1->ScegliTipo(*n,2);//Funzione da sistemare
  system("rm -f ./Posizioni/*.dat");
  for(unsigned long int k=0;!(n[0]->MNEffet()==0||n[1]->MNEffet()==0);k++){
    if(comm==1)printf("%d)  ",k);
    if(k==-614){
      n[0]->Stato();
      n[1]->Stato();
      c1->Stato();
    }
    Tipo=c1->ScegliT(*n[0],*n[1]);
    ATipo=(Tipo+1)%2;
    if(( Quale=n[Tipo]->ScegliP() )>0){
      n[ATipo]->AggiornaTempo();
      modo=n[Tipo]->ScegliModo();
      if(c1->Sposta(*n[Tipo],modo,Quale)==1){
	c1->Controlla(*n[Tipo],Quale,*n[ATipo]);
      }
    }
    if(t>tempo){
      sprintf(file,"Posizioni/Pos_t%d.dat",tempo);
      if((Posizioni = fopen(file,"w"))!=NULL){
	fprintf(Posizioni,"%d\t%d\n",n[0]->MNEffet()+n[1]->MNEffet(),tempo);
	for(int it=0;it<2;it++){
	  for(int i=0;i<n[it]->MNEffet();i++){
	    fprintf(Posizioni,"%d\t%d\t%d\t%d\t%d\n",n[it]->MTocc(i),n[it]->MTipo(),n[it]->MPos(i,0),n[it]->MPos(i,1),n[it]->MPos(i,2));
	  }
	}
	fclose(Posizioni);
      }
      else printf("Non s'apre %s\n",file);
      printf("\nAl tempo %.0f ps sono rimaste %d n{%d} e %d n{%d} particelle\n",t,n[0]->nEffet,n[0]->Tipo,n[1]->nEffet,n[1]->Tipo);
        tempo*=10;
    }
  }
  (void)time(&tf);
  printf("-------------------------------------------------------------\n\
Alla fine della simulazione dopo %d s per %.0f ps sono rimaste %d particelle di tipo %d e %d di tipo %d\n\
----------------------------------------------------------------\n",(int)(tf-ti),t,n[0]->nEffet,n[0]->Tipo,n[1]->nEffet,n[1]->Tipo);
  return 0;
}
