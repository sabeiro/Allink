#include "Disegna.h"

void lista(void){
  SferaB = glGenLists(1);
  glNewList(SferaB,GL_COMPILE);
    //1
    glPushMatrix();
    glColor4f(.0,.0,1.,1.);
    glutSolidSphere(Diameter,Detail,Detail);
    glPopMatrix();
    glEndList();
  SferaV = glGenLists(1);
  glNewList(SferaV,GL_COMPILE);
    //2
    glPushMatrix();
    glColor4f(.0,1.0,.0,1.);
    glutSolidSphere(Diameter,Detail,Detail);
    glPopMatrix();
    glEndList();
    SferaVi = glGenLists(1);
  glNewList(SferaVi,GL_COMPILE);
    //3
    glPushMatrix();
    glColor4f( .4,.4,1.,1.);
    glutSolidSphere(Diameter+ExtraDiam,Detail,Detail);
    glPopMatrix();
    glEndList();
  Griglia = glGenLists(1);
  glNewList(Griglia,GL_COMPILE);
    glBegin(GL_LINE_STRIP);
      glColor3f(1.0,1.0,1.0);
      glVertex3f(0.,0.,0.);
      glVertex3f(1.,0.,0.);
      glVertex3f(1.,1.,0.);
      glVertex3f(0.,1.,0.);
      glVertex3f(0.,0.,0.);
    glEnd();
    glBegin(GL_LINES);
    glColor3f(1.0,1.0,1.0);
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
static void Particle(void){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  //  glutTimerFunc(tDiap, Particle, 0);
  Diap++;
  tDiap=glutGet(GLUT_ELAPSED_TIME);
  if (tDiap - tDiapBase > 1000) {
    //    sprintf(s,"FPS:%4.2f",frame*1000.0/(time-timebase));
    tDiapBase = Diap;		
    Diap = 0;
  }
  if(lu != 0){
    glMaterialfv(GL_FRONT,GL_AMBIENT,ambientMaterial);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuseMaterial);
    glMaterialfv(GL_FRONT,GL_SPECULAR,specularMaterial);
    glMaterialf(GL_FRONT,GL_SHININESS,shininessMaterial);
    glMaterialfv(GL_FRONT,GL_EMISSION,emissionMaterial);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,riflettivitaSpeculare);
  }
  if(P[quando].Num <= 100){
    Detail = 20;
    lista();
  }
  else if(P[quando].Num > 100 && P[quando].Num <= 1000){
    Detail =10;
    lista();
  }
  else if(P[quando].Num > 1000){
    Detail =4;
    lista();
  }
  glPushMatrix();
  glTranslatef(xp,yp,zp);
  glTranslatef(.5,.5,.0);
  glRotatef(xa,1.,0.,0.);
  glRotatef(ya,0.,1.,0.);
  glRotatef(za,0.,0.,1.);
  glRotatef(spin,.1,.2,.3);
  glColor4f(1.,1.,1.,1.);
  if(la!=0)glutWireCube((GLfloat)1);
  for(int j=0;j<P[quando].Num;j++){
    glPushMatrix();
    glTranslatef((GLfloat)(P[quando].pos0[j]/P[quando].lato[0]-.5),
		 (GLfloat)(P[quando].pos1[j]/P[quando].lato[1]-.5),
		 -(GLfloat)(P[quando].pos2[j]/P[quando].lato[2]-.5));
    //printf("%d) [%d] (%.3f,%.3f,%.3f)\n",j,P[quando].tocc[j],P[quando].pos0[j]/P[quando].lato[0],P[quando].pos1[j]/P[quando].lato[1],P[quando].pos2[j]/P[quando].lato[2]);
    if(P[quando].tipo[j] == 0){
      glCallList(SferaV);
    }
    else if(P[quando].tipo[j] == 1){
      glCallList(SferaB);
    }
    else if(P[quando].tipo[j] == 2){
      glCallList(SferaVi);
    }
    glPopMatrix();
  }
  glPopMatrix();
  glPushMatrix();
  glTranslatef(-1.1,-1.1,.0);
  glColor3f(1.0,1.0,1.0);
  glRasterPos3f(0., 0., 0.);
  sprintf(info,"There are %d particle with %d different types",P[quando].Num,P[quando].Type+1);
  int len = (int)strlen(info);
  for (int i = 0; i < len; i++) {
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, info[i]);
  }
  glPopMatrix();
  
  glutSwapBuffers();
}
void figura(void){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   //  GLUquadricObj *sfera;
   //  sfera = gluNewQuadric();
  glPushMatrix();
  
  glRotatef(xf,1.,0.,0.);
  glRotatef(yf,0.,1.,0.);
  glTranslatef(0.,0.,zf);
  //  gluScaleImage
  glColor4f(1.,1.,1.,1.);
  if(la!=0)glutWireCube((GLfloat)1);
  glPopMatrix();
 
  glutSwapBuffers();
}
void Mostra(void){
  for(int j =0;j<P[quando].Num;j++){
    printf("%d) [%d] (%.3f,%.3f,%.3f)\n",j,P[quando].tocc[j],P[quando].pos0[j]/P[quando].lato[0],P[quando].pos1[j]/P[quando].lato[1],P[quando].pos2[j]/P[quando].lato[2]);
  }
}
void spinDisplay(void){
  spin += dspin;
  if(spin > 360.)
    spin -= 360.;
  glutPostRedisplay();
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
void init (void){
  glClearColor(.0,.0,.0,.0);//colore sfondo

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
    glLoadIdentity();
    gluPerspective(60,(GLfloat) w/(GLfloat) h,1.,200.);//ang,rapp,zmin,zmax
  glTranslatef(0.,0.,-2.);
  }
  if(pr == 0){
    glOrtho(0.,1.,0.,1.,0.,1.);
  }
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  //gluLookAt( 0., 0., 1.,0., 0., 0.,0., 1., 0.);//osservatore,centro,orientazion
  //glFrustum(-3,3.,-3.,3.,-3.,3.);//Dimensione massima dell'immagine
}
void processEvent(int value){
  switch(value){
  case '1':
    glutIdleFunc(spinDisplay);
    glutPostRedisplay();
  case '2':
    glutIdleFunc(NULL);
    glutPostRedisplay();
  }
}
void Menu(){
  submenu = glutCreateMenu(processEvent);
  glutAddMenuEntry("SottoPrima",3);

  menu = glutCreateMenu(processEvent);
  glutAddMenuEntry("Move",1);
  glutAddMenuEntry("Stop",2);
  glutAddSubMenu("Altro",submenu);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}
static void MouseMove(int x,int y){
  if(y > yRem)
    xa += 2;
  else if(y < yRem)
    xa -= 2;
  if(x > xRem)
    ya += 2;
  else if(x < xRem)
    ya -= 2;
  xRem = x;
  yRem = y;
  //  printf("(%f,%f)\n",xa,ya);
  glutPostRedisplay();
}
void mouse(int button, int state,int x,int y){
  switch (button){
  case GLUT_LEFT_BUTTON:
    if(state==GLUT_DOWN){
      if(ChangeMouse == 0){
	//	glutPassiveMotionFunc( MouseMove );
	//glutIdleFunc(spinDisplay);
      }
      else {
	//	  glutPassiveMotionFunc(NULL);
      }
      ChangeMouse = !ChangeMouse;
    }
    break;
  case GLUT_WHEEL_DOWN:
    zp -= .1;
    glutPostRedisplay();
    break;
  case GLUT_WHEEL_UP:
    zp += .1;
    glutPostRedisplay();
    break;    
  default:
    break;
  }
}
static void keyboard(unsigned char key,int x, int y){
  switch (key){
  case '>':
    if(quando<arc-1)
      quando++;
    glutPostRedisplay();
    break;
  case '<':
    if(quando>0)
      quando--;
    glutPostRedisplay();
    break;
  case 'C':
    printf("Coordinate (%f,%f,%f) Angolo (%f,%f,%f)\n"
	   ,xp,yp,zp,xa,ya,za);
    break;
  case 'd':
    Diameter += StepDiameter;
    lista();
    glutPostRedisplay();
    break;
  case 'D':
    Diameter -= StepDiameter;
    lista();
    glutPostRedisplay();
    break;
  case 'e':
    ExtraDiam += StepDiameter;
    lista();
    glutPostRedisplay();
    break;
  case 'E':
    ExtraDiam -= StepDiameter;
    lista();
    glutPostRedisplay();
    break;
  case 'f':
    glutFullScreen();
    glutPostRedisplay();
    break;
  case 'F':

    glutPostRedisplay();
    break;
  case 'g':
    gr += 1    ;
    if(gr == 2) gr =0;
    glutPostRedisplay();
    break;
  case 'l':
    la += 1    ;
    if(la == 2) la =0;
    glutPostRedisplay();
    break;
  case 'm':
    Menu();
    glutPostRedisplay();
    break;
  case 'M':
    glutDestroyMenu(menu);
    glutPostRedisplay();
    break;
  case 'p':
    pr += 1;
    if(pr == 2) pr =0;
    glutPostRedisplay();
    break;
  case 'q':
    exit(0);
    break;
  case 'r':
    spin =  0.;angolo = 0.; dspin = .1 ;
    xa=-15. ; ya=0. ; za= .0;
    xf=90. ; yf=0. ; zf= 0.;
    xp=-.5 ; yp=-.5 ; zp= -.1;
    glutPostRedisplay();
    break;
  case 't':
    sp += 1    ;
    if(sp == 2) sp =0;
    glutPostRedisplay();
    break;
  case 's':
    dspin += .01;
    glutPostRedisplay();
    break;
  case 'S':
    dspin -= .01;
    glutPostRedisplay();
    break;
  case 'u':
    sc =!sc;
    break;
  case 'U':
    break;
  case 'w':
    il += 1    ;
    if(il == 2) il =0;
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
  case 'z':
    zp += .1;
    glutPostRedisplay();
    break;
  case 'Z':
    zp -= .1;
    glutPostRedisplay();
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
  width = (int) w;
  height = (int) h;
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

int main(int argc,char **argv){
  arc=argc-1;
  printf("\n------------------------Disegna------------------------------\n\
*        Questo programma richiede in entrata un file      *\n\
*                formattato nel seguente modo:             *\n\
* %%d\t%%d\t%%f\t%%f\t%%f                                 *\n\
* num\ttype\tx\ty\tz                                  *\n\
* Richiede come argomento almeno un file                   *\n\
------------------------------------------------------------\n");
  time_t ti,tf;
  (void) time(&ti);
  FILE *Pos;
  char *nome;
  int *dec;
  int Tipo;
  double *pos0,*pos1,*pos2;
  dec  = (int *)malloc(sizeof(int));
  pos0 = (double *)malloc(sizeof(double));
  pos1 = (double *)malloc(sizeof(double));
  pos2 = (double *)malloc(sizeof(double));  
  nome = (char *)malloc(30*sizeof(char));
  P  = (PART1 *)malloc(arc*sizeof(PART1));
  info = (char *)malloc(100*sizeof(char));
  if(argc<2){
    printf("Mi serve come argomento un file!\n");
    return 0;
  }
  for(int i=0;i<arc;i++){
    sprintf(nome,"%s",argv[i+1]);
    printf("Apro %s su %d\n",nome,arc);
    if((Pos=fopen(nome,"r"))==0)printf("Non s'apre %s\n",nome);
    fscanf(Pos,"%d",dec);
    P[i].Num = *dec;
    fscanf(Pos,"%d",dec);
    P[i].time = *dec;
    P[i].tocc = (int *)malloc(P[i].Num*sizeof(int));
    P[i].tipo = (int *)malloc(P[i].Num*sizeof(int));
    P[i].pos0 = (float *)malloc(P[i].Num*sizeof(float));
    P[i].pos1 = (float *)malloc(P[i].Num*sizeof(float));
    P[i].pos2 = (float *)malloc(P[i].Num*sizeof(float));
    for(int l=0;l<3;l++)
      P[i].lato[l] = 0;
    /* reads and store number type and position */
    CountT0 = CountT1 = 0;
    for(int j=0;!(fscanf(Pos,"%d",dec)==0||j>=P[i].Num);j++){
      P[i].tocc[j]=*dec;
      fscanf(Pos,"%d",dec);
      if(*dec==0){
	P[i].tipo[j] = 0;
	CountT0++;
      }
      else if(*dec==1){
	P[i].tipo[j] = 1;      
	CountT1++;
      }
      else if(*dec==2){
	P[i].tipo[j] = 2;      
	CountT2++;
      }
      fscanf(Pos,"%lf",pos0);
      P[i].pos0[j]=*pos0;
      if(P[i].pos0[j]>P[i].lato[0])
	P[i].lato[0] = P[i].pos0[j];
      fscanf(Pos,"%lf",pos1);
      P[i].pos1[j]=*pos1;
      if(P[i].pos1[j]>P[i].lato[1])
	P[i].lato[1] = P[i].pos1[j];
      fscanf(Pos,"%lf",pos2);
      P[i].pos2[j]=*pos2;
      if(P[i].pos2[j]>P[i].lato[2])
	P[i].lato[2] = P[i].pos2[j];
    }
    P[i].Num = CountT0 + CountT1 + CountT2;
    if(CountT0 != 0)
      P[i].Type = 0;
    if(CountT1 != 0)
      P[i].Type = 1;
    if(CountT2 != 0)
      P[i].Type = 2;    
    printf("At the time %d we have %d green, %d blue and %d violet particles\n",i,CountT0,CountT1,CountT2);
  }
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
  glutInitWindowSize(400,400);//Dimensione finestra
  glutInitWindowPosition(000,300);//Posizione finestra
  MainWindow = glutCreateWindow("Disegna");//Titolo

  lista();
  init();
  Menu();
//   SubWindow1 = glutCreateSubWindow(MainWindow, 0, 100,
// 				  (width), (height-100));
  if (il != 0)Illuminazione();
  if (lu != 0);
  if (ne != 0)Nebbia();
  if (sp != 0)Spot();
  //glutTimerFunc(tDiap,Particle,1);// tempo tra due diapositive
  //  glutReshapeFunc(ChangeSize);//Abilita la ridimensione per mouse
  glutDisplayFunc(Particle);//Disegna la figura
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutKeyboardFunc(keyboard);
  glutMotionFunc( MouseMove );
  glutSpecialFunc(special);
  glutMainLoop();//Mantiene la finestra

//   SubWindow2 = glutCreateSubWindow(MainWindow, 0, 0,
// 				  (width), 100);
//   glClearColor(.0,.0,.0,.0);//colore sfondo
//   glutDisplayFunc(figura);


  printf("Te ve be te e?\n");
 
  return 0;
}
