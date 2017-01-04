#ifndef VARIABILI_H
#define VARIABILI_H
#include "../include/Matematica.h"

typedef struct {bool **Segnale;int NMass;int NSegn;}SEGNALE;
class Variabili:Matematica{
 private:
  FILE *MEDIA;
  FILE *SEGNALE;
 public:
  Variabili(int NSegn,int Potenza2,double MediaSu,double ScartoSu,double MediaGiu,double ScartoGiu,double Passo,double Incr);
  Variabili(char *NomeMedia,int Valori);
  ~Variabili();
  void Verificare(int cosa,int Valore);
  bool Simula();
  bool Sifula();
  void ScriviFile();
  void Cambia(char*,int);
  void Cambia(char*,double);
  bool *st;
  double *sMedia;
  int NMass;
  int Periodi;
  int NSegn;
  double MediaSu;
  double ScartoSu;
  double MediaGiu;
  double ScartoGiu;
  double Passo;
  double Incr;
  Matematica *mSu;
  Matematica *mGiu;
};
#endif //VARIABILI_H
