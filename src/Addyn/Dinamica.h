#include "../include/VarData.h"
#include <GL/glut.h>

extern float ExtraDiam;
extern float ScaleUn;
extern int Values;
extern GLuint Particles,SferaB,SferaV,Griglia,Quad,Point,SferaR,Chains;
extern GLfloat spin,angolo, dspin,
  xa, ya, za,
  xf, yf, zf,
  xp, yp, zp,
  xi, yi, zi,
  scale, dscale,tscale;
extern char *info;
extern GENERAL *Gen1;
extern float NanoRad,ExtraDiam;
extern float Diameter,StepDiameter;
extern int Passo;//Griglia
extern int IfVideo;
extern int MainWindow,SubWindow1,SubWindow2;
extern int la;
extern int gr;
extern int IfShift;
extern int IfPoint;//point or Sphere
extern int pr;
void Legenda();
//ElPoly Draw
extern void Principal(int argc,char** argv);
extern void Lista();
extern void Particle(void);
extern void processEvent(int value);
extern void Menu();
extern void keyboard(unsigned char key,int x, int y);
extern void InitConstant();
extern void keyboardDraw(unsigned char key);
extern void Timer(int v);
extern int Picture();
PART *Pn;
static int NVisChainMin=0;
static int NVisChainMax=0;
static int IfIntorno =0;//Intorno
static int IfColour=0;
static int IfWhite=0;
static int IfLine=0;//disegna i legami
int menu,submenu;

Forces *Fm;
int Part2Move;
