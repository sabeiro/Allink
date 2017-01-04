#ifndef DISEGNA_H
#define DISEGNA_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>

#define QUAD(A) ((A*A))

#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4
#endif


typedef struct {float *pos0;float *pos1;float *pos2;int *tocc;
  int *tipo;int Num;int lato[3];float time;int Type;} PART1;
PART1 *P;

static GLfloat spin =  0.,angolo = 0., dspin = .1 ,
  xa=-15. , ya=0. , za= .0,
  xf=90. , yf=0. , zf= 0.,
  xp=-.5 , yp=-.5 , zp= -.1;
static int pr = 1,//prospettiva
  sc = 0,la=1,gr=0;//sfera-cubo
GLuint SferaB,SferaV,Griglia,SferaVi;
int arc = 0;
int quando=0;
int Passo=4;
int menu,submenu;
int Diap=0,tDiap,tDiapBase=0;
bool ChangeMouse = 0;
float Diameter = .005,StepDiameter = .005,ExtraDiam = .0;
int Detail = 20;
int xRem=0, yRem=0;
int CountT0=0, CountT1=0,CountT2 =0;
int MainWindow,SubWindow1,SubWindow2;
int width=600,height=600;
char *info;

static int lu = 1;
GLfloat ambientMaterial[4]     = { 0.0, 1.0, 1.0,1.};
GLfloat diffuseMaterial[4]     = { .8, .8, .0,1.};
GLfloat specularMaterial[4]    = { .0, .0, 0.0,1.0};
GLfloat shininessMaterial      = 10.0;
GLfloat emissionMaterial[4]    = { 0.0, 0.0, 0.0,1.};
static int il = 1;
GLfloat luceAmbiente[4]        = { .1, .1, .1,1.};
GLfloat luceDiffusa[4]         = { 1., 1., 1.,1.};
GLfloat luceSpeculare[]        = { 1.0, 1.0, 1.0,1.};
GLfloat direzione[4]           = { 1.0, 1.0, 0.0,0.};
GLfloat posizioneLight0[4]     = { .0, 1.0, 1.0,1.};
static int sp = 1;
GLfloat direzioneSpot[3]       = { .0, 1.0, 1.0};
GLfloat riflettivitaSpeculare[]= { 0.0, .0, .0,1.};
GLfloat coloreNebbia[]         = { 0.2, 0.2, 0.2,1.};
GLfloat angoloSpot             = 120.;
static int ne = 0;
GLfloat zIniNebbia             = 10.;
GLfloat zFinNebbia             = -10.;
GLfloat esponenteSpot          = 1.5;
GLfloat DensitaNebbia          = 1.;

GLint   MaterialShininess      = 10;

GLfloat nColore[] = {0.,0.,1.,1.};
GLfloat mColore[] = {0.,1.,0.,1.};


#endif//DISEGNA_H
