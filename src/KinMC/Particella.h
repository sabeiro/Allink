#ifndef PARTICELLA_H
#define PARTICELLA_H
#include "Variabili.h"
#define quad(a) ((a)*(a))//Quadrato
#define DisRet 2//DistanzaReticolo
#define Min(a,b) a<b?a:b
using namespace std;
extern float t;
extern int comm;
class Particella{
private:
public:
  Particella(char *cDis,char *cVal,int NPart);
  int ScegliP();
  void AggiornaTempo();
  void CambiaTemperatura(float DeltaT);
  void Elimina(int Tocc);
  int  ScegliModo();
  int ScegliSup();
  void Stato();
  void Stato(int Tocc);
  void fStato(int temp);
  //  ~Particella();
  FILE   *FDis;
  FILE   *FVal;
  FILE *Posizioni;
  float Norm;
  vector <float> PrDis;
  float   PrMoto[7];
  float   PrSup[4];
  float   Pr0,Pr,E0;
  int nDis;
  int nEffet;
  int nPart;
  int Tipo;
  PART *Part;
  float Dis;
  float PrCum;
  float T;
  float Deltat;
  int   MTocc(int i);
  int MPos(int i,int l);
  int   MTipo(void);
  int   MNEffet(void);
};
#endif //VARIABILI_H
