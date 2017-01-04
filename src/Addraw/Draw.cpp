/***********************************************************************
Draw: Creator and initialisation of the graphical framework
Copyright (C) 2008-2010 by Giovanni Marelli <sabeiro@virgilio.it>


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#include "../../include/Draw.h"

#ifdef __glut_h__

Draw::Draw(){
  info = (char *)malloc(256*sizeof(char));
  Number = (char *)malloc(100*sizeof(char));
  InvScaleUn = 1.;
  for(int d=0;d<3;d++){
    Edge[d] = 1.;
  }
  ImHeight = 800;
  ImWidth = 800;
  pixel = NULL;
  pixel = (GLubyte*)malloc(sizeof(GLubyte));
  NLevel = 4;
  InitConstant();
  xRem = 0;yRem=0;
  IncrVisDxSx = 2.;
  IncrVisSuGiu = 2.;
}
void Draw::Shout(const char * s, ...)
{
#ifdef DRAW_DEBUG
  va_list args;
  va_start(args, s);
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
#else
  return;
#endif
}
// #include <SDL.h>
// #include "SDL/SDL_opengl.h" 
// int Draw::WindowSDL(){
//   int ret = SDL_Init(SDL_INIT_EVERYTHING);
//   if(  ret < 0 ) return false;
//   //ret = SDL_SetVideoMode(gScreenWidth, gScreenHeight, 0, SDL_OPENGL | SDL_RESIZABLE | SDL_ASYNCBLIT | SDL_HWSURFACE | SDL_ANYFORMAT);
//   SDL_Surface *screen;
//   screen = SDL_SetVideoMode(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, SDL_OPENGL);
//   if( screen == NULL ) return false;
//   SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
//   SDL_MapRGB(main_surface->format, 255, 255, 255);
//   InitConstant();
//   Illuminazione();
//   if(ne) Nebbia();
//   if(sp) Spot();
//   ParticleList();
//   ReadScript();
//   ChooseBlend(1);
//   SDL_WM_SetCaption("Disegna",NULL);

//   SDL_Event event;
  
//   SDL_WaitEvent(&event);
  
//   switch (event.type) {
//   case SDL_KEYDOWN:
//     printf("Il tasto %s e' stato premuto!\n",
// 	   SDL_GetKeyName(event.key.keysym.sym));
//     break;
//   case SDL_QUIT:
//     exit(0);
//   }
//   {
//     SDL_Event event;
    
//     while ( SDL_PollEvent(&event) ) {
//       switch (event.type) {
//       case SDL_MOUSEMOTION:
// 	printf("Il Mouse e' stato mosso da %d,%d "	\
// 	       "a (%d,%d)\n", 
// 	       event.motion.xrel, event.motion.yrel,
// 	       event.motion.x, event.motion.y);
// 	break;
//       case SDL_MOUSEBUTTONDOWN:
// 	printf("Il pulsante %d del Mouse %d e' stato "	\
// 	       "premuto a (%d,%d)\n",
// 	       event.button.button, event.button.x, 
// 	       event.button.y);
//                 break;
//       case SDL_QUIT:
// 	exit(0);
//       }
//     }
//   }

//   return true;
// }
void Draw::Window(int argc,char **argv){
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
  glutInitWindowSize(ImWidth,ImHeight);//Dimensione finestra
  //glutInitWindowPosition(100,300);//Posizione finestra
  MainWindow = glutCreateWindow("Disegna");//Titolo
  //  sprintf(info,"There are %d particle with %d different types",Gen1->NPart,Gen1->NType);
  InitConstant();
  ReadConf();
  Lista(Values);
  //   SubWindow1 = glutCreateSubWindow(MainWindow, 0, 100,
  // 				  (width), (height-100));
  ParticleList();
  ReadScript();
  // cout << "GL_VENDOR      : " << glGetString(GL_VENDOR) << endl;
  // cout << "GL_RENDERER    : " << glGetString(GL_RENDERER) << endl;
  // cout << "GL_VERSION     : " << glGetString(GL_VERSION) << endl;
  // //cout << "GL_EXTENSION   : " << glGetString(GL_EXTENSIONS) << endl;
  glutTimerFunc(0,Timer,0);// tempo tra due diapositive
  glutDisplayFunc(Figure);//Disegna la figura
  glutReshapeFunc(reshape);
  //glutReshapeFunc(ChangeSize);//Abilita la ridimensione per mouse
  glutMouseFunc(mouse);
  glutKeyboardFunc(keyboard);
  glutMotionFunc( MouseMove );
  glutSpecialFunc(special);
  ChooseBlend(1);
  Illuminazione();
  if(ne) Nebbia();
  if(sp) Spot();
  //glutMainLoop();//Mantiene la finestra
  //   SubWindow2 = glutCreateSubWindow(MainWindow, 0, 0,
  // 				  (width), 100);
  //   gl
  glClearColor(.0,.0,.0,.0);//colore sfondo
}
void Draw::DTimer(int v){
  glutPostRedisplay();
  glutTimerFunc(30,Timer, 0);
}
// void Draw::Illuminazione(){

// }
void Draw::Illuminazione(void){
  //------------------Depth------------------
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  //glDepthFunc(GL_LEQUAL);
  glClearDepth( 1.0 ); 
  //glClearDepth(1.0);
  //glDepthMask(GL_TRUE);  // make Z-buffer writable
  //  glDepthMask(GL_FALSE);
  //  glDepthRange(0,-10.);
  //----------------Blend---------------------------
  glEnable(GL_BLEND); 
  glShadeModel(GL_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_ALPHA_TEST);
  glAlphaFunc(GL_GREATER,0.1);
  //------------Miscellaneus--------------------------
  //glEnable(GL_STENCIL_TEST); /// Enable stencil test
  //glStencilOp(GL_KEEP, GL_KEEP, GL_INCR); ///
  //glStencilFunc(GL_ALWAYS, uniqueStencilValue, ~0);
  //glEnable(GL_TEXTURE_2D);
  glReadBuffer(GL_FRONT);/* so glReadPixel() always get the right image */
  //glEnable(GL_NORMALIZE);
  //---------------Material----------------------
  // glEnable(GL_COLOR_MATERIAL);// have materials set by curr color
  //glEnable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//GL_LINE
  glCullFace(GL_FRONT_AND_BACK);
  glFrontFace(GL_CCW);//Senso antiorario
  //glEnableClientState(GL_VERTEX_ARRAY);
  //glEnableClientState(GL_NORMAL_ARRAY);
  //glEnable(GL_VERTEX_ARRAY);
  GLfloat ambientMaterial[4]     = { .5, .5, .5,1.};
  GLfloat diffuseMaterial[4]     = { .0, .0, .0,1.};
  GLfloat specularMaterial[4]    = { .0, .0, .0,1.0};
  GLfloat shininessMaterial      = 25.0;//50.0;
  GLfloat emissionMaterial[4]    = { .0, .0, .0,1.};
  GLfloat riflettivitaSpeculare[]= { 0.0, .0, .0,1.};
  //glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  //glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
  //glColorMaterial(GL_FRONT_AND_BACK,GL_EMISSION);
  //glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  //glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
  //glDisable(GL_POLYGON_SMOOTH);// make sure not to antialias polygons
  glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,DrAmbientWhite);
  glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,DrDiffuseWhite);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,DrSpecularWhite);
  glMaterialf(GL_FRONT,GL_SHININESS,shininessMaterial);
  glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,emissionMaterial);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,riflettivitaSpeculare);
  glMaterialfv(GL_BACK,  GL_AMBIENT,   DrAmbientGreen); 
  glMaterialfv(GL_BACK,  GL_DIFFUSE,   DrDiffuseGreen); 
  glMaterialfv(GL_FRONT, GL_AMBIENT,   DrAmbientBlue); 
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   DrDiffuseBlue); 
  glMaterialfv(GL_FRONT, GL_SPECULAR,  DrSpecularWhite); 
  glMaterialf( GL_FRONT, GL_SHININESS, 25.0); 
  //--------------------Light--------------------------
  GLfloat DrPropertiesAmbient [] = {0.50, 0.50, 0.50, 1.00}; 
  GLfloat DrPropertiesDiffuse [] = {0.75, 0.75, 0.75, 1.00}; 
  GLfloat DrPropertiesSpecular[] = {1.00, 1.00, 1.00, 1.00}; 

  glEnable(GL_LIGHTING);

  //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1.0);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0); 

  glEnable( GL_LIGHT0 );
  glLightfv( GL_LIGHT0, GL_AMBIENT,  DrPropertiesAmbient); 
  glLightfv( GL_LIGHT0, GL_DIFFUSE,  DrPropertiesDiffuse); 
  glLightfv( GL_LIGHT0, GL_SPECULAR, DrPropertiesSpecular); 
  GLfloat LightPos0[]={xl0,yl0,zl0, 0.0f };
  GLfloat LightPos1[]={xl1,yl1,zl1, 0.0f };
  glLightfv(GL_LIGHT0, GL_POSITION,LightPos0);
  glEnable(GL_LIGHT1);
  glLightfv(GL_LIGHT1, GL_AMBIENT, DrPropertiesAmbient);		
  glLightfv(GL_LIGHT1, GL_DIFFUSE, DrPropertiesDiffuse);		
  glLightfv(GL_LIGHT1, GL_POSITION,LightPos1);
  glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 1.5);
  glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.5);
  glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.2);
  glLightf(GL_LIGHT1, GL_SPOT_CUTOFF, 45.0);
  glLightf(GL_LIGHT1, GL_SPOT_EXPONENT, 2.0);
}
void Draw::Dreshape(int w,int h){
  glViewport(0,0,(GLsizei)  w, (GLsizei) h);
  WinWidth = (int) w;
  WinHeight = (int) h;
  XCenter = WinWidth / 2;
  YCenter = WinHeight / 2;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if(pr == 1){
    //glLoadIdentity();
    gluPerspective(30,(GLfloat) w/(GLfloat) h,.1,200.);//ang,rapp,zmin,zmax
    glFrustum(-1.,1.,-Edge[1]*InvScaleUn,Edge[1]*InvScaleUn,-Edge[2]*InvScaleUn,Edge[2]*InvScaleUn);//Sx,Dx,giu,su,vicino,lontano//Dimensione massima dell'immagine   
    //glTranslatef(0.,0.,-2.);
  }
  if(pr == 0){
    //glOrtho(-1.,1.,(float)(w-1),(float)(h-1),.1,-200.);//Sx,Dx,giu,su,vicino,lontano
    //    glTranslatef(0.,0.,1.0);
    glOrtho(-1.,1.,-WinHeight/(float)WinWidth,WinHeight/(float)WinWidth,.1,20);//Sx,Dx,giu,su,vicino,lontano  
    glEnable(GL_DEPTH_TEST);
  }
  glMatrixMode(GL_MODELVIEW);
  //gluLookAt( -1., -1., -1.,.0, .0, .0,1., 1., 1.);//osservatore,centro,orientazion
  glLoadIdentity();
  gluLookAt( 0., 0.,1.5/* 1.5*Edge[2]*InvScaleUn*/,
	     0., 0., 0.,
	     0., 1., 0.);//osservatore,centro,orientazion
  //  glLoadIdentity();
  //glFrustum(-3,3.,-3.,3.,-3.,3.);//Sx,Dx,giu,su,vicino,lontano//Dimensione massima dell'immagine
  //glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);// Really Nice Perspective Calculations
}
void Draw::Spot(void){
  GLfloat direzioneSpot[3]       = { .0, 1.0, 1.0};
  GLfloat coloreNebbia[]         = { 0.2, 0.2, 0.2,1.};
  GLfloat esponenteSpot          = 1.5;
  GLfloat angoloSpot             = 120.;
  glLightf(GL_LIGHT0,GL_SPOT_CUTOFF,angoloSpot);
  glLightfv(GL_LIGHT0,GL_SPOT_DIRECTION,direzioneSpot);
  glLightf(GL_LIGHT0,GL_SPOT_EXPONENT,esponenteSpot);
}
void Draw::Nebbia(void){
  GLfloat zIniNebbia             = 10.;
  GLfloat zFinNebbia             = -10.;
  GLfloat DensitaNebbia          = 1.;
  glEnable(GL_FOG);
  glFogf(GL_FOG_START,zIniNebbia);
  glFogf(GL_FOG_END,zFinNebbia);
  glFogi(GL_FOG_MODE,GL_EXP);
  glFogf(GL_FOG_DENSITY,DensitaNebbia);
}
void Draw::ChangeSize(GLsizei w,GLsizei h){
  GLfloat aspectRatio;
  WinWidth = (int) w;
  WinHeight = (int) h;
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
void Draw::Numera(double *Pos,int n){
  char *Number = (char *)calloc(50,sizeof(char));
  sprintf(Number,"%d",n);
  glRasterPos3f((Pos[0]*InvScaleUn),
		(Pos[1]*InvScaleUn),
		(Pos[2]*InvScaleUn));
  for (int l = 0; l < strlen(Number); l++) {
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, Number[l]);
  }
  glCallList(Point);
  free(Number);
}

#endif// __glut_h__
