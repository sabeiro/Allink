#ifndef CUBO_H
#define CUBO_H
#include "Particella.h"
#include "Variabili.h"
#define quad(a) ((a)*(a))//Quadrato
#define DisRet 2//DistanzaReticolo
#define Min(a,b) a<b?a:b
using namespace std;
extern int comm;
extern float t;
class Particella;
class Cubo{
private:
  int yMod10,zMod10,nMass;
  int nCelle[3];
  int nPartizioni[3];
  int Lato[3];
public:
  //Costruttore, definisce le dimensioni delle celle
  Cubo(int LatoCubo[3],int Partizioni[3]);
  //Posiziona le particelle
  bool Posiziona(Particella &p);
  //Fornisce la posizione relativa nella cella
  int PosRet(const int Pos[3]);
  //Da le nformazioni contenute in ogni cella
  void Stato();
  //Informazione della i-esima cella
  void Stato(int nCella);
  //sceglie il tipo di particella da muovere
  int ScegliTipo(const Particella p[],const int tipi);
  //uguale ma funziona
  int ScegliT(const Particella &p1,const Particella &p2);
  //sposta la particella scelta, rit. 0 se fallisce
  bool Sposta(Particella &p1,const int modo,const int Quale);
  //controlla che non sia vicina ad altre celle, rit.0 se fallisce
  bool Controlla(Particella &p1,const int Quale,Particella &p2);
  //Elimina le particelle che si ricombinano, rit.1 se ricombinano
  bool nmTocc(Particella &p1,Particella &p2,const int dove1,const int dove2,const int Quale);
  CELLA *Cella;
  //fornisce il valore (x,y,z)->n
  int nCella(const int pos[3]);
  //fornisce la posizione nella cella della particella Tocc
  int Questa(const int Tipo,const int Tocc,const int dove);
  //ritorna iRet->(x,y,z)
  COORD Coord(const int Ret);
  //ritorna nCella->(x,y,z)
  COORD NCella(const int mCella);
  //File sulle particelle uscite
  FILE *Uscite,*Ricombinate;
};
#endif //CUBO_H
