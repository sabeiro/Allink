#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>

static GLfloat spin =  0.,angolo = 0., dspin = 1 ,
  xa=-15. , ya=0. , za= .0,
  xp=-.5 , yp=-.5 , zp= -.1;
static int pr = 1,//prospettiva
  sc = 0,la=0,gr=0;//sfera-cubo
GLuint SferaB,SferaV,Griglia,SferaVi;
int Passo=4;

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

void Sfera(void){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if(lu != 0){
        glMaterialfv(GL_FRONT,GL_AMBIENT,ambientMaterial);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuseMaterial);
    glMaterialfv(GL_FRONT,GL_SPECULAR,specularMaterial);
    glMaterialf(GL_FRONT,GL_SHININESS,shininessMaterial);
    glMaterialfv(GL_FRONT,GL_EMISSION,emissionMaterial);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,riflettivitaSpeculare);
  }
  glPushMatrix();
  glTranslatef(xp,yp,zp);


  glRotatef(xa,1.,0.,0.);
  glRotatef(ya,0.,1.,0.);
  glRotatef(za,0.,0.,1.);
  glRotatef(spin,.1,.2,.3);
  //1
  glTranslatef(0.,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);//
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);//
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  //2
  glTranslatef(-1.,(float)1/Passo,0.);
  glCallList(SferaVi);//
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);//
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  //3
  glTranslatef(-1.,(float)1/Passo,0.);
  glCallList(SferaV);//
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  //4
  glTranslatef(-1.,(float)1/Passo,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  //5
  glTranslatef(-1.,(float)1/Passo,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
  glTranslatef((float)1/Passo,0.,0.);
  glCallList(SferaVi);
 
   glTranslatef(-.12,-.12,0.);
   glCallList(SferaB);
//    glTranslatef(-.25,0,0.);
//    glCallList(SferaB);
//    glTranslatef(-.5,0.,0.);
//    glCallList(SferaB);
//    glTranslatef(.5,-.25,0.);
//    glCallList(SferaB);
//    glTranslatef(.25,0,0.);
//    glCallList(SferaB);

  


  glPopMatrix();

  glPushMatrix();
  glTranslatef(xp,yp,zp);
  glRotatef(xa,1.,0.,0.);
  glRotatef(ya,0.,1.,0.);
  glRotatef(za,0.,0.,1.);
  glCallList(Griglia);
  glTranslatef(0.,0.,-.25);
  //  glCallList(Griglia);

  glPopMatrix();
  glutSwapBuffers();
}
void lista(void){
  SferaB = glGenLists(1);
  glNewList(SferaB,GL_COMPILE);
    //1
    glPushMatrix();
    glColor4f(.0,.0,1.,1.);
    glutSolidSphere(.045,20,20);
    glPopMatrix();
    glEndList();
  SferaV = glGenLists(1);
  glNewList(SferaV,GL_COMPILE);
    //2
    glPushMatrix();
    glColor4f(.0,1.0,.0,1.);
    glutSolidSphere(.045,20,20);
    glPopMatrix();
    glEndList();
    SferaVi = glGenLists(1);
  glNewList(SferaVi,GL_COMPILE);
    //3
    glPushMatrix();
    glColor4f( .4,.4,1.,1.);
    glutSolidSphere(.045,20,20);
    glPopMatrix();
    glEndList();
  Griglia = glGenLists(1);
  glNewList(Griglia,GL_COMPILE);
    glBegin(GL_LINE_STRIP);
      glColor3f(.0,.0,.0);
      glVertex3f(0.,0.,0.);
      glVertex3f(1.,0.,0.);
      glVertex3f(1.,1.,0.);
      glVertex3f(0.,1.,0.);
      glVertex3f(0.,0.,0.);
    glEnd();
    glBegin(GL_LINES);
    glColor3f(.0,.0,.0);
    for(int i=1;i<Passo;i++){
      glVertex3f((GLfloat)i/Passo,1,0);
      glVertex3f((GLfloat)i/Passo,0,0);
    }
    for(int i=1;i<5;i++){
      glVertex3f(1,(GLfloat)i/Passo,0);
      glVertex3f(0,(GLfloat)i/Passo,0);
    }
    glEnd();
  glEndList();
}

  
void spinDisplay(void){
  spin += dspin;
  if(spin > 360.)
    spin -= 360.;
  glutPostRedisplay();
}
void init (void){
  glClearColor(1.0,1.0,1.0,.0);//colore sfondo

}
void Illuminazione(void){
  glEnable(GL_LIGHTING);// Illuminazione
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);  
  glEnable(GL_DEPTH_TEST);//Controllo di profondit`a

  glEnable(GL_LIGHT0);//Fonte di raggi paralleli
  glLightfv(GL_LIGHT0,GL_AMBIENT,luceAmbiente);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,luceDiffusa);
  glLightfv(GL_LIGHT0,GL_SPECULAR,luceSpeculare);
  glLightfv(GL_LIGHT0,GL_POSITION,posizioneLight0);
}
void Spot(void){
  glLightf(GL_LIGHT0,GL_SPOT_CUTOFF,angoloSpot);
    glLightfv(GL_LIGHT0,GL_SPOT_DIRECTION,direzioneSpot);
  glLightf(GL_LIGHT0,GL_SPOT_EXPONENT,esponenteSpot);
}
void Nebbia(void){
  glEnable(GL_FOG);
  glFogf(GL_FOG_START,zIniNebbia);
  glFogf(GL_FOG_END,zFinNebbia);
  glFogi(GL_FOG_MODE,GL_EXP);
  glFogf(GL_FOG_DENSITY,DensitaNebbia);
}
void reshape(int w,int h){
  glViewport(0,0,(GLsizei)  w, (GLsizei) h);
  glLoadIdentity();
  if(pr == 1){
    glMatrixMode(GL_PROJECTION);
    gluPerspective(60,(GLfloat) w/(GLfloat) h,1.,2.);//ang,rapp,zmin,zmax
  glTranslatef(0.,0.,-1.);
  }
  if(pr == 0){
    glOrtho(0.,1.,0.,1.,0.,1.);
  }
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  //gluLookAt( 0., 0., 1.,0., 0., 0.,0., 1., 0.);//osservatore,centro,orientazion
  //glFrustum(-3,3.,-3.,3.,-3.,3.);//Dimensione massima dell'immagine
}
static void MouseMove(int x,int y){
  xa = y/2;
  ya = x/2;
  glutPostRedisplay();
}
void mouse(int button, int state,int x,int y){
  switch (button){
  case GLUT_LEFT_BUTTON:
    if(state==GLUT_DOWN){
      glutIdleFunc(spinDisplay);
    }
    break;
  case GLUT_RIGHT_BUTTON:
    if(state == GLUT_DOWN)
      glutIdleFunc(NULL);
    break;
  default:
    break;
  }
}
static void keyboard(unsigned char key,int x, int y){
  switch (key){
  case 's':
    dspin += .1;
    glutPostRedisplay();
    break;
  case 'S':
    dspin -= .1;
    glutPostRedisplay();
    break;
  case 'C':
    printf("Coordinate (%f,%f,%f) Angolo (%f,%f,%f)\n"
	   ,xp,yp,zp,xa,ya,za);
    break;
  case 'p':
    pr += 1;
    if(pr == 2) pr =0;
    glutPostRedisplay();
    break;
  case 'w':
    il += 1    ;
    if(il == 2) il =0;
    glutPostRedisplay();
    break;
  case 'l':
    la += 1    ;
    if(la == 2) la =0;
    glutPostRedisplay();
    break;
  case 'r':    
    glutPostRedisplay();
    break;
  case 'R':
    glutIdleFunc(Sfera);
    glutPostRedisplay();
    break;
  case 't':
    sp += 1    ;
    if(sp == 2) sp =0;
    glutPostRedisplay();
    break;
  case 'g':
    gr += 1    ;
    if(gr == 2) gr =0;
    glutPostRedisplay();
    break;
  case 'u':
    sc =!sc;
    break;
  case 'U':
    break;
  case 'z':
    zp += .1;
    glutPostRedisplay();
    break;
  case 'Z':
    zp -= .1;
    glutPostRedisplay();
    break;
  case 'x':
    xp += .1;
    glutPostRedisplay();
    break;
  case 'X':
    xp -= .1;
    glutPostRedisplay();
    break;
  case 'y':
    yp += .1;
    glutPostRedisplay();
    break;
  case 'Y':
    yp -= .1;
    glutPostRedisplay();
    break;
  case 'q':
    exit(0);
    break;
  case 27:
    exit(0);
    break;
  default:
    break;
  }
}
static void
special(int k, int x, int y)
{
  switch (k) {
  case GLUT_KEY_UP:
    xa += 5.0;
    break;
  case GLUT_KEY_DOWN:
    xa -= 5.0;
    break;
  case GLUT_KEY_LEFT:
    ya += 5.0;
    break;
  case GLUT_KEY_RIGHT:
    ya -= 5.0;
    break;
  default:
    return;
  }
  glutPostRedisplay();
}
void Timer(int v){
  glutPostRedisplay();
  glutTimerFunc(1,Timer, 1);
}
void ChangeSize(GLsizei w,GLsizei h){
  GLfloat aspectRatio;
  if(h==0)h=1;
  aspectRatio=(GLfloat)w/(GLfloat)h;
  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if(w<=h)
    glOrtho(-100.,100.,-100./aspectRatio,100.
	    /aspectRatio,1.,-1.);
  else
    glOrtho(-100*aspectRatio,-100*aspectRatio,
	    -100,100,1.,-1.);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  }
int main(int argc,char** argv)
{
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
  glutInitWindowSize(400,400);//Dimensione finestra
  glutInitWindowPosition(000,300);//Posizione finestra
  glutCreateWindow("Diffusione");//Titolo

  lista();
  init();
  if (il != 0)Illuminazione();
  if (lu != 0);
  if (ne != 0)Nebbia();
  if (sp != 0)Spot();

  //  glutTimerFunc(1,Timer,1);// tempo tra due diapositive
  //  glutReshapeFunc(ChangeSize);//Abilita la ridimensione per mouse
  glutDisplayFunc(Sfera);//Disegna la figura
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special);
  //   glutPassiveMotionFunc( MouseMove );
   glutMainLoop();//Mantiene la finestra
  return 0;
}
