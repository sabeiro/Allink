#ifndef VARIABILI_H
#define VARIABILI_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector> 
#include <list>
#include <iostream>
#include <time.h>
#define quad(a) ((a)*(a))//Quadrato
#define DisRet 2//DistanzaReticolo
#define Min(a,b) a<b?a:b
using namespace std;
double ran1(void);
//Coordinate Posizione e nome particella
typedef struct {int Pos[3]; int Tocc;}PART;
//nome, tipo e disposizione nella cella
typedef struct {int Tocc;int Tipo;int Ret;}OCCUP;
//struttura contenente le infomazioni della particella
typedef struct {vector <OCCUP> Part;}CELLA;
//Coordinate
typedef struct{int x;int y;int z;}COORD;
#endif //VARIABILI_H
