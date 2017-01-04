/***********************************************************************
DrawDefinition: Create some useful lists
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
void Draw::InitConstant(){
  spin =  0.;angolo = 0.; dspin = .1 ;
  xp=0.  ; yp=0.  ; zp= 0.; zw = -1.;
  xa=0.  ; ya=0.  ; za= 0.;
  xf=55. ; yf=5.  ; zf= 0.;
  xi=-.4 ; yi=-.4 ; zi= -0.;
  xLeg=-1.4; yLeg=-.3 ; zLeg= 0.;
  dxLeg=.4;dyLeg=.8;
  xl0=1. ; yl0=1. ; zl0=1.;
  xl1=1. ; yl1=0. ; zl1=1.;
  Rback = 1.0; Gback = 1.0;Bback = 1.0;Aback = 0.0;
  GridStep = 2.;
  for(int d=0;d<3;d++)
    GridEdge[d] = (int)(GridStep*Edge[d]*InvScaleUn);
  Step = 0;
  pr = 1;//prospettiva
  IfPoint = 1;//point or Sphere
  IfInfo = 1;
  IfImage = 0;
  IfScript = 1;
  IfBlend = 1;
  la = 1;
  gr = 0;
  lu = 1;
  sp = 1;
  ne = 0;
  AngleDxSx = 0.;
  AngleSuGiu = 0.;
  IfMaterial = 1;
  ChooseDepthMap(0);
  sprintf(info,"");
}
int Draw::DefQuad(int NSquare){
  glDeleteLists(Quad,1);
  Quad = glGenLists(1);
  double xEdge = 1./(double)NSquare*Edge[0]*InvScaleUn;
  double yEdge = 1./(double)NSquare*Edge[1]*InvScaleUn;
  glNewList(Quad,GL_COMPILE);
  glBegin(GL_POLYGON);
  glNormal3f(0.,0.,1.);
  glVertex3f(0,0,0);
  glVertex3f(xEdge,0,0);
  glVertex3f(xEdge,yEdge,0);
  glVertex3f(0,yEdge,0);
  glEnd();
  glEndList();
  return Quad;
}
int Draw::DefHexagon(){
  glDeleteLists(Hexagon,1);
  Hexagon = glGenLists(1);
  glNewList(Hexagon,GL_COMPILE);
  double Lato=0.01;
  glBegin(GL_POLYGON);
  //  glBegin(GL_POINTS);
  //  glVertex3f(0.,0.,0.);
  glNormal3f(0.,0.,1.);
  glVertex3f(Lato*.5,-Lato*.87,0.);
  glVertex3f(Lato,0.,0.);
  glVertex3f(Lato*.5,Lato*.87,0.);
  glVertex3f(-Lato*.5,Lato*.87,0.);
  glVertex3f(-Lato,0.,0.);
  glVertex3f(-Lato*.5,-Lato*.87,0.);
  glEnd();
  glEndList();
  return Hexagon;
}
int Draw::DefCube(int NSquare){
  glDeleteLists(Cube,1);
  Cube = glGenLists(1);
  double xEdge = 1./(double)NSquare*Edge[0]*InvScaleUn;
  double yEdge = 1./(double)NSquare*Edge[1]*InvScaleUn;
  double zEdge = 1./(double)NSquare*Edge[2]*InvScaleUn;
  glNewList(Cube,GL_COMPILE);
  glBegin(GL_QUADS);
  for(int f=0;f<2;f++){
    glNormal3d(0.,0.,1.);
    glVertex3d(0,0,f*zEdge);
    glVertex3d(xEdge,0,f*zEdge);
    glVertex3d(xEdge,yEdge,f*zEdge);
    glVertex3d(0,yEdge,f*zEdge);
    glNormal3d(1.,0.,0.);
    glVertex3d(f*xEdge,0,0);
    glVertex3d(f*xEdge,0,zEdge);
    glVertex3d(f*xEdge,yEdge,zEdge);
    glVertex3d(f*xEdge,yEdge,0);
    glNormal3d(0.,1.,0.);
    glVertex3d(0,f*yEdge,0);
    glVertex3d(0,f*yEdge,zEdge);
    glVertex3d(xEdge,f*yEdge,zEdge);
    glVertex3d(xEdge,f*yEdge,0);
  }
  glEnd();
  glEndList();
  return Cube;
}
int Draw::DefWall(){
  glDeleteLists(GlWall,1);
  GlWall = glGenLists(1);
  glEnable(GL_TEXTURE_2D);
  glNewList(GlWall,GL_COMPILE);
  DefTexture();
  //  GL_DECAL,GL_MODULATE,GL_BLEND,GL_REPLACE
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,
		  GL_LINEAR_MIPMAP_LINEAR);
  glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
  // glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,
  //                                      GL_NEAREST_MIPMAP_NEAREST);
  // glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glBegin(GL_QUADS);
  glNormal3d(0.,0.,1.);
  glColor4f(1.,.0,.0,.5);
  glVertex3d(0.,0.,0.);
  glColor4f(1.,.4,.0,.5);
  glVertex3d(Edge[0]*InvScaleUn,0.,0.);
  glColor4f(1.,.4,.4,.5);
  glVertex3d(Edge[0]*InvScaleUn,Edge[1]*InvScaleUn,0.);
  glColor4f(1.,.0,.4,.5);
  glVertex3d(0.,Edge[1]*InvScaleUn,0.);
  glColor4f(1.,.2,.2,.5);
  glTexCoord2f(0.8, 0.8);glVertex3d(0.,0.,0.);
  glColor4f(1.,.2,.0,.5);
  glTexCoord2f(0.2, 0.8);glVertex3d(Edge[0]*InvScaleUn,0.,0.);
  glColor4f(1.,.0,.0,.5);
  glTexCoord2f(0.2, 0.2);glVertex3d(Edge[0]*InvScaleUn,Edge[1]*InvScaleUn,0.);
  glColor4f(1.,.0,.2,.5);
  glTexCoord2f(0.8, 0.2);glVertex3d(0.,Edge[1]*InvScaleUn,0.);
  glEnd();
  glEndList();
  glDisable(GL_TEXTURE_2D);
  return GlWall;
}
int Draw::DefArrow(){
  glDeleteLists(Arrow,1);
  Arrow = glGenLists(1);
  glNewList(Arrow,GL_COMPILE);
  glNormal3f(0.,0.,-1.);
  glBegin(GL_LINES);
  glVertex3f(0.,0.,0.);
  glVertex3f(1.,0.,0.);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f(1.,0.,0.);
  glVertex3f(.95,.02,0.);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f(1.,0.,0.);
  glVertex3f(.95,-.02,0.);
  glEnd();
  glEndList();
  return Arrow;
}
int Draw::DefArrowThin(){
  glDeleteLists(Arrow,1);
  Arrow = glGenLists(1);
  glNewList(Arrow,GL_COMPILE);
  glNormal3f(0.,0.,-1.);
  glBegin(GL_LINES);
  glVertex3f(0.,0.,0.);
  glVertex3f(1.,0.,0.);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f(1.,0.,0.);
  glVertex3f(.95,.002,0.);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f(1.,0.,0.);
  glVertex3f(.95,-.002,0.);
  glEnd();
  glEndList();
  return Arrow;
}
int Draw::DefPoint(){
  glDeleteLists(Point,1);
  Point = glGenLists(1);
  glNewList(Point,GL_COMPILE);
  glBegin(GL_POINTS);
  glNormal3f(0.,0.,1.);
  glVertex3f(0,0,0);
  glEnd();
  glEndList();
  return Point;
}
double Draw::Normal(double *v,double *u,double *w,double *n){
  double Norma = 0.;
  for(int d=0;d<3;d++){
    int Coord1 = (d+1)%3;
    int Coord2 = (d+2)%3;
    double Uno = (v[Coord1] - w[Coord1])*(u[Coord2] - w[Coord2]); 
    double Due = (v[Coord2] - w[Coord2])*(u[Coord1] - w[Coord1]); 
    n[d] = Uno - Due;
    Norma += SQR(n[d]);
  }
  Norma = sqrt(Norma);
  for(int d=0;d<3;d++)
    n[d] /= Norma;
  return Norma;
}
int Draw::DefCylinder(double CyEdge,double CyHeight){
  int ListRef = glGenLists(1);
  glEnable(GL_LIGHTING);
  glEnable( GL_LIGHT0 );
  glMaterialfv(GL_BACK,  GL_AMBIENT,   DrAmbientRed); 
  glMaterialfv(GL_BACK,  GL_DIFFUSE,   DrDiffuseRed); 
  glMaterialfv(GL_FRONT, GL_AMBIENT,   DrAmbientRed); 
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   DrDiffuseRed); 
  glMaterialfv(GL_FRONT, GL_SPECULAR,  DrSpecularWhite); 
  glMaterialf( GL_FRONT, GL_SHININESS, 25.0); 
  glNewList(ListRef,GL_COMPILE);
  int NEdges = 40;
  double Angle = DUE_PI/(double)(NEdges);
  CyEdge *= InvScaleUn;
  CyHeight *= InvScaleUn;
  //    glEnableClientState(GL_NORMAL_ARRAY);
  glBegin(GL_POLYGON);//Kork
  glNormal3f(0.,0.,-1.);
  for(int l=0;l<NEdges+1;l++){
    glVertex3d(CyEdge*cos(Angle*l),CyEdge*sin(Angle*l),-CyHeight*.5);
  }
  glEnd();//Kork
  glBegin(GL_POLYGON);//Kork
  glNormal3f(0.,0.,1.);
  for(int l=0;l<NEdges+1;l++){
    glVertex3d(CyEdge*cos(Angle*l),CyEdge*sin(Angle*l),CyHeight*.5);
  }
  glEnd();//Kork
  for(int l=0;l<NEdges;l++){
    double v1[3] = {CyEdge*cos(Angle*l),CyEdge*sin(Angle*l),-CyHeight*.5};
    double v2[3] = {CyEdge*cos(Angle*l),CyEdge*sin(Angle*l),CyHeight*.5};
    double v3[3] = {CyEdge*cos(Angle*(l+1)),CyEdge*sin(Angle*(l+1)),CyHeight*.5};
    double v4[3] = {CyEdge*cos(Angle*(l+1)),CyEdge*sin(Angle*(l+1)),-CyHeight*.5};
    double vn[3];
    glBegin(GL_POLYGON);//Side
    Normal(v2,v3,v1,vn);
    //printf("%lf %lf %lf\n",vn[0],vn[1],vn[2]);
    //printf("%lf %lf %lf\n",-sin(Angle*(l)),cos(Angle*(l)),0);
    glNormal3f(vn[0],vn[1],vn[2]);
    glNormal3f(-CyEdge*sin(Angle*(l)),CyEdge*cos(Angle*(l)),0);
    glVertex3d(v1[0],v1[1],v1[2]);
    Normal(v3,v4,v2,vn);
    glNormal3f(vn[0],vn[1],vn[2]);
    glNormal3f(-CyEdge*sin(Angle*(l)),CyEdge*cos(Angle*(l)),0);
    glVertex3d(v2[0],v2[1],v2[2]);
    Normal(v4,v1,v3,vn);
    glNormal3f(vn[0],vn[1],vn[2]);
    glNormal3f(-CyEdge*sin(Angle*(l+1)),CyEdge*cos(Angle*(l+1)),0);
    glVertex3d(v3[0],v3[1],v3[2]);
    Normal(v1,v2,v4,vn);
    glNormal3f(vn[0],vn[1],vn[2]);
    glNormal3f(-CyEdge*sin(Angle*(l+1)),CyEdge*cos(Angle*(l+1)),0);
    glVertex3d(v4[0],v4[1],v4[2]);
    glEnd();//Side
  }
  glEndList();
  glDisable(GL_LIGHTING);
  return ListRef;
  //Cupola
  for(int cc=0;cc<NEdges/2;cc+=1){
    for(int c=0;c<NEdges;c+=1){
      double Quota = CyHeight*.5;
      if(cc > NEdges/4) Quota = - CyHeight*.5;
      double v1[3] = {CyEdge*cos(Angle*c)*sin(cc*Angle),
		      CyEdge*sin(Angle*c)*sin(cc*Angle),
		      CyEdge*cos(cc*Angle) + Quota};
      double v2[3] = {CyEdge*cos(Angle*(c+1))*sin(cc*Angle),
		      CyEdge*sin(Angle*(c+1))*sin(cc*Angle),
		      CyEdge*cos(cc*Angle) + Quota};
      double v3[3] = {CyEdge*cos(Angle*(c+1))*sin((cc+1)*Angle),
		      CyEdge*sin(Angle*(c+1))*sin((cc+1)*Angle),
		      CyEdge*cos((cc+1)*Angle) + Quota};
      double v4[3] = {CyEdge*cos(Angle*c)*sin((cc+1)*Angle),
		      CyEdge*sin(Angle*c)*sin((cc+1)*Angle),
		      CyEdge*cos((cc+1)*Angle) + Quota};
      double vn[3];
      glBegin(GL_POLYGON);//Kork
      glNormal3f(vn[0],vn[1],vn[2]);
      Normal(v1,v2,v3,vn);
      glVertex3d(v1[0],v1[1],v1[2]);
      glNormal3f(vn[0],vn[1],vn[2]);
      Normal(v2,v3,v4,vn);
      glVertex3d(v2[0],v2[1],v2[2]);
      glNormal3f(vn[0],vn[1],vn[2]);
      Normal(v3,v4,v1,vn);
      glVertex3d(v3[0],v3[1],v3[2]);
      glNormal3f(vn[0],vn[1],vn[2]);
      Normal(v4,v1,v2,vn);
      glVertex3d(v4[0],v4[1],v4[2]);
      glEnd();//Kork
    }
  }
  glEndList();
  return ListRef;
      // 	static GLUquadricObj *quadObj;	  
      // 	quadObj = gluNewQuadric();
      // 	gluQuadricDrawStyle(quadObj, GLU_FILL);
      // 	gluQuadricNormals(quadObj, GLU_SMOOTH);
      // 	gluCylinder(quadObj,  Nano->Rad*InvScaleUn,  Nano->Rad*InvScaleUn,  Nano->Height*InvScaleUn, 20, 20);
}
#include <vector> 
int Draw::DefMetalCylinder(double CyEdge,double CyHeight){
  // GLfloat Light1[4] = {1.,1.,1.,1.};
  // glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,Light1);
  // GLfloat Light2[4] = {1.,1.,1.,1.};
  // glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Light2);
  // GLfloat Light3[4] = {0.,0.,0.,1.};
  // glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Light3);
  int MetalCylinder = 0;
  MetalCylinder = glGenLists(1);
  glNewList(MetalCylinder,GL_COMPILE);
  int NEdges = 20;
  double Angle = DUE_PI/(double)(NEdges);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnable(GL_VERTEX_ARRAY);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glFrontFace(GL_CCW);   //Tell OGL which orientation shall be the front face
  glShadeModel(GL_SMOOTH);
  glEnableClientState(GL_NORMAL_ARRAY);
  double Normal[3]={.1,.1,.1};
  glNormalPointer(GL_FLOAT,sizeof(double),&Normal[0]);  //Pointer to the first color*/
  GLfloat *Kork = new GLfloat [2*3*(NEdges)];
  vector <GLuint> VecIndices;
  for(int l=0;l<NEdges;l++){
    Kork[(l)*3+0] = CyEdge*cos(Angle*l);
    Kork[(l)*3+1] = CyEdge*sin(Angle*l);
    Kork[(l)*3+2] = CyHeight*.5;
    VecIndices.push_back(l);
  }
  VecIndices.push_back(0);
  int NIndex = VecIndices.size();
  //    GLuint *Indices = new GLuint [NIndex];
  //    for(int i=0;i<NIndex;i++)Indices[i] = VecIndices[i];
  GLuint *Indices = &VecIndices[0];
  glVertexPointer(3,GL_FLOAT,3*sizeof(GLfloat),Kork);
  glDrawElements(GL_POLYGON,NIndex,GL_UNSIGNED_INT,Indices);
  for(int l=0;l<NEdges;l++)
    Kork[(l)*3+2] = -CyHeight*.5;
  glDrawElements(GL_POLYGON,NIndex,GL_UNSIGNED_INT,Indices);
  delete [] Kork;
  VecIndices.clear();
  GLfloat *Side = new GLfloat [3*4];   
  for(int l=0;l<NEdges+1;l++){
    Side[(0)*3+0] = CyEdge*cos(Angle*l);
    Side[(0)*3+1] = CyEdge*sin(Angle*l);
    Side[(0)*3+2] = -CyHeight*.5;
    VecIndices.push_back(0);
    Side[(1)*3+0] = CyEdge*cos(Angle*(l));
    Side[(1)*3+1] = CyEdge*sin(Angle*(l));
    Side[(1)*3+2] = CyHeight*.5;
    VecIndices.push_back(1);
    Side[(2)*3+0] = CyEdge*cos(Angle*(l+1));
    Side[(2)*3+1] = CyEdge*sin(Angle*(l+1));
    Side[(2)*3+2] = CyHeight*.5;
    VecIndices.push_back(2);
    Side[(3)*3+0] = CyEdge*cos(Angle*(l+1));
    Side[(3)*3+1] = CyEdge*sin(Angle*(l+1));
    Side[(3)*3+2] = -CyHeight*.5;
    VecIndices.push_back(3);
    Indices = &VecIndices[0];
    glVertexPointer(3,GL_FLOAT,3*sizeof(GLfloat),Side);
    glDrawElements(GL_QUADS,4,GL_UNSIGNED_INT,Indices);
  }
  delete [] Side;
  VecIndices.clear();
  glEndList();
  return MetalCylinder;
}

// Written by Chris Halsall (chalsall@chalsall.com) for the 
// O'Reilly Network on Linux.com (oreilly.linux.com).
// May 2000.
int Draw::DefTexture(){
  GLuint Texture_ID = 0;
  GLenum gluerr;
  GLubyte tex[128][128][4];
  int x,y,t;
  int hole_size = 3300; // ~ == 57.45 ^ 2.
  // Generate a texture index, then bind it for future operations.
  glGenTextures(1,&Texture_ID);
  glBindTexture(GL_TEXTURE_2D,Texture_ID);
  // Iterate across the texture array.
  for(y=0;y<128;y++) {
    for(x=0;x<128;x++) {
      // A simple repeating squares pattern.
      // Dark blue on white.
      if ( ( (x+4)%32 < 8 ) && ( (y+4)%32 < 8)) {
	tex[x][y][0]=tex[x][y][1]=0; tex[x][y][2]=120;
      } else {
	tex[x][y][0]=tex[x][y][1]=tex[x][y][2]=240;
      }
      // Make a round dot in the texture's alpha-channel.
      // Calculate distance to center (squared).
      t = (x-64)*(x-64) + (y-64)*(y-64) ;
      if ( t < hole_size) // Don't take square root; compare squared.
	tex[x][y][3]=255; // The dot itself is opaque.
      else if (t < hole_size + 100)
	tex[x][y][3]=128; // Give our dot an anti-aliased edge.
      else
	tex[x][y][3]=0;   // Outside of the dot, it's transparent.
    }
  }
  // The GLU library helps us build MipMaps for our texture.
  if ((gluerr=gluBuild2DMipmaps(GL_TEXTURE_2D, 4, 128, 128, GL_RGBA,
				GL_UNSIGNED_BYTE, (void *)tex))) {
    fprintf(stderr,"GLULib%s\n",gluErrorString(gluerr));
    exit(-1);
  }
  // Some pretty standard settings for wrapping and filtering.
  glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
  glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
  // We start with GL_DECAL mode.
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL);
  return Texture_ID;
}
int Draw::DefGriglia(){
  glDisable(GL_LIGHTING);
  glDeleteLists(Griglia,1);
  Griglia = glGenLists(1);
  glNewList(Griglia,GL_COMPILE);
  glBegin(GL_LINES);
  glColor4f(1.-Rback,1.-Gback,1.-Bback,1.0);
  double Step[3];
  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);
  for(int d=0;d<3;d++){
    //GridEdge[d] = (int)(GridStep*Edge[d]*InvScaleUn);
    Step[d] = Edge[d]*InvScaleUn/GridEdge[d];
  }
  //--------------------------grid
  if(gr == 1){
    for(int i=0;i<=GridEdge[0];i++){
      for(int j=0;j<=GridEdge[1];j++){
	glVertex3f(i*Step[0],j*Step[1],Edge[2]*InvScaleUn);
	glVertex3f(i*Step[0],j*Step[1],0);
      }
    }
    for(int j=0;j<=GridEdge[1];j++){
      for(int k=0;k<=GridEdge[2];k++){
	glVertex3f(0,j*Step[1],k*Step[2]);
	glVertex3f(Edge[0]*InvScaleUn,j*Step[1],k*Step[2]);
      }
    }
    for(int i=0;i<=GridEdge[0];i++){
      for(int k=0;k<=GridEdge[2];k++){
	glVertex3f(i*Step[0],0.,k*Step[2]);
	glVertex3f(i*Step[0],Edge[1]*InvScaleUn,k*Step[2]);
      }
    }
  }
  else if (gr == 2){
    for(int i=0;i<=GridEdge[0];i++){
      glVertex3f(i*Step[0],0.,Edge[2]*InvScaleUn);
      glVertex3f(i*Step[0],0.,0.);
      glVertex3f(i*Step[0],Edge[1]*InvScaleUn,0.);
      glVertex3f(i*Step[0],0.,0.);
    }
    for(int j=0;j<=GridEdge[1];j++){
      glVertex3f(0,j*Step[1],0.);
      glVertex3f(Edge[0]*InvScaleUn,j*Step[1],0.);
      glVertex3f(0.,j*Step[1],0.);
      glVertex3f(0.,j*Step[1],InvScaleUn*Edge[2]);
    }
    for(int k=0;k<=GridEdge[2];k++){
      glVertex3f(0.,0.,k*Step[2]);
      glVertex3f(0.,Edge[1]*InvScaleUn,k*Step[2]);
      glVertex3f(0.,0.,k*Step[2]);
      glVertex3f(Edge[0]*InvScaleUn,0.,k*Step[2]);
    }
  }
  //glDisable(GL_LINE_SMOOTH);
  //----------------------label
  for(int i=0;i<=GridEdge[0];i++){
    //glPushMatrix();//info
    sprintf(Number,"%.2f",i/(double)GridEdge[0]);
    glRasterPos3f(i*Step[0],0.,0.);
    for (int l = 1; l < strlen(Number); l++) {
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, Number[l]);
    }
    glCallList(Point);
    //glPopMatrix();//info
  } 
  for(int j=0;j<=GridEdge[1];j++){
    //glPushMatrix();//info
    sprintf(Number,"%.2f",j/(double)GridEdge[1]);
    glRasterPos3f(0.,j*Step[1],0.);
    for (int l = 1; l < strlen(Number); l++) {
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, Number[l]);
    }
    glCallList(Point);
    //glPopMatrix();//info
  } 
  for(int k=0;k<=GridEdge[2];k++){
    //glPushMatrix();//info
    sprintf(Number,"%.2f",k/(double)GridEdge[2]);
    glRasterPos3f(0.,0.,k*Step[2]);
    for (int l = 1; l < strlen(Number); l++) {
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, Number[l]);
    }
    glCallList(Point);
    //glPopMatrix();//info
  }
  glColor4f(.2,.1,.7,1.);
  glRasterPos3f(-.1,0.,0.);
  glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,'0');
  glRasterPos3f(Edge[0]*InvScaleUn+.1,0.,0.);
  glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,'x');
  glRasterPos3f(0.,Edge[1]*InvScaleUn+.1,0.);
  glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,'y');
  glRasterPos3f(0.,0.,Edge[2]*InvScaleUn+.1);
  glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,'z');
  glEnd();
  glPopAttrib(); 
  glEndList();
  return Griglia;
}
void Draw::Lista(int NSquare){
  // glEnable(GL_LIGHTING);
  // glEnable(GL_COLOR_MATERIAL);
  // glEnable(GL_NORMALIZE);
  // //glEnable(GL_CULL_FACE);
  // glCullFace(GL_FRONT_AND_BACK);//Disegna anche il lato dietro
  // //glFrontFace(GL_CCW);//Senso antiorario
  // //glColorMaterial(GL_FRONT_AND_BACK,GL_EMISSION);
  // GLfloat Light1[4] = {0.,0.,0.,1.};
  // glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,Light1);
  // GLfloat Light2[4] = {1.,1.,1.,1.};
  // glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Light2);
  // GLfloat Light3[4] = {0.,0.,0.,1.};
  // glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Light3);
  // GLfloat Light4[4] = {1.,1.,1.,1.};
  // glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,Light3);
  // GLfloat Light5[4] = {1.,1.,1.,1.};
  // glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,Light3);
  // glEnable(GL_NORMALIZE);
  DefPoint();
  DefGriglia();
  DefLegend();
  DefArrowThin();
}
static double WeightFuncSpline(const double x){
  if(x < 0.)
    return 1.;
  if(fabs(x) < 1.)
    return 2.*fabs(x)*x*x - 3.*x*x + 1.;
  else
    return 0.;
}
static double Parab(const double x){
  double Resp = -4.*x*x+4.;
  if(Resp < 0) Resp = 0.;
  return Resp;
}
void Draw::DepthMap1(double Val,GLfloat *Color){
  Color[0] = (GLfloat)WeightFuncSpline((1.-Val)*.7);
  Color[1] = (GLfloat)Parab((Val-.5)*3.);
  Color[2] = (GLfloat)WeightFuncSpline((Val)*.7);
  Color[3] = 1.;
}
void Draw::ChooseDepthMap(int i){
  if(i == 0)
    Depth_Map = &Draw::DepthMap1;
}
int Draw::DefLegend(){
  glDisable(GL_LIGHTING);
  glDeleteLists(DrLegend,1);
  DrLegend = glGenLists(1);
  glNewList(DrLegend,GL_COMPILE);
  int NSlab = 100;
  GLfloat InvNSlab = 1./(GLfloat)NSlab;
  double dx = dxLeg*InvNSlab;
  double dy = dyLeg*InvNSlab;
  GLfloat Color[4];
  char Number[10];
  for(int i=0;i<NSlab;i++){
    GLfloat x = xLeg;// + Size[0]*i*dx;
    GLfloat y = yLeg + dyLeg*i*dy;
    glBegin(GL_POLYGON);
    glNormal3d(0.,0.,1.);
    DepthMap(i*InvNSlab,Color);
    for(int d=0;d<3;d++)Color[d] -= .2;
    glColor4fv(Color);
    glVertex3d(x,y,zLeg);
    glVertex3d(x+dx,y,zLeg);
    glVertex3d(x+dx,y+dy,zLeg);
    glVertex3d(x,y+dy,zLeg);
    glEnd();
    if(!(i%10)){
      glPushMatrix();//Info
      sprintf(Number,"%.2f",i*InvNSlab);
      glRasterPos3f(xLeg+4.*dx,y,zLeg);
      for(int l=1;l<strlen(Number);l++){
	//glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,Number[l]);
      	glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,Number[l]);
      }
      glPopMatrix();
    }
  }
  glPushMatrix();//Info
  sprintf(Number,"1.00");
  glRasterPos3f(xLeg+4.*dx,yLeg+dyLeg,zLeg);
  for(int l=1;l<strlen(Number);l++){
    //glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,Number[l]);
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,Number[l]);
  }
  glPopMatrix();
  glEndList();
  return DrLegend;
}
#endif // __glut_h__
