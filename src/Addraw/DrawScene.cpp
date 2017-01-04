#include <Draw.h>

#ifdef __glut_h__
#include <vector>
int Detail = 20;
GLfloat luceAmbiente[4]        = { .1, .1, .1,1.}; 
GLfloat luceDiffusa[4]         = { 1., 1., 1.,1.};
GLfloat luceSpeculare[]        = { 1.0, 1.0, 1.0,1.};
GLfloat direzione[4]           = { 1.0, 1.0, 0.0,1.};
GLfloat posizioneLight0[4]     = { 1., 1., 1.0,0.};
GLfloat posizioneLight1[4]     = { 1., 0., 1.0,0.};
//int width=600,height=600;
GLint   MaterialShininess      = 50;

void Draw::DMinimal(void){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor(Rback,Gback,Bback,Aback);//colore sfondo 
  glPushMatrix();//Whell
  glTranslatef(0,0,zw);
  glPushMatrix();//Rotation
  glRotatef(xa,1.,0.,0.);
  glRotatef(ya,0.,1.,0.);
  glRotatef(za,0.,0.,1.);
  glPushMatrix();//Translation
  //glTranslatef(xp,yp,zp);
  //glTranslatef(-.5*Edge[0]*InvScaleUn,-.5*Edge[1]*InvScaleUn,-.5*Edge[2]*InvScaleUn);
  glCallList(Particles);
  ParticleRealTime();
  glPopMatrix();//Translation
  glPopMatrix();//Rotation
  glPopMatrix();//Wheel
  //glFlush();			//Finish rendering
  glutSwapBuffers();
}
void Draw::DFigure(void){
  glMatrixMode( GL_MODELVIEW );
  //  pixel = (GLuint *)realloc(pixel,4*width*height*sizeof(*pixel));
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor(Rback,Gback,Bback,Aback);//colore sfondo
  glPushMatrix();//Whell
  glTranslatef(0,0,zw);
  glPushMatrix();//Rotation
  glRotatef(xa,1.,0.,0.);
  glRotatef(ya,0.,1.,0.);
  glRotatef(za,0.,0.,1.);
  glPushMatrix();//Translation
  glTranslatef(xp,yp,zp);
  glTranslatef(-.5*Edge[0]*InvScaleUn,-.5*Edge[1]*InvScaleUn,-.5*Edge[2]*InvScaleUn);
  if(IfBlend) glEnable(GL_DEPTH_TEST); 
  else glDisable(GL_DEPTH_TEST); 
  glBlendFunc(BlendSource,BlendDest);
  if(lu != 0){
    glPushMatrix();//Light
    glRotatef(xf,1.,0.,0.);
    glRotatef(yf,0.,1.,0.);
    glRotatef(zf,0.,0.,1.);
    glTranslatef(1.5*Edge[0]*InvScaleUn,1.5*Edge[1]*InvScaleUn,1.5*Edge[2]*InvScaleUn);
    GLfloat PosLight[4] = {xl0,yl0,zl0,1.};
    glLightfv(GL_LIGHT0,GL_POSITION,PosLight);
    //glLightfv(GL_LIGHT1,GL_POSITION,posizioneLight1);
    //glutSolidSphere(0.05,20,20);
    glPopMatrix();//Light
  }
  int PreChain = 0;
  glCallList(Particles);
  if(IfScript)
    glCallList(ScriptList);
  ParticleRealTime();
  if(IfImage){
    glPushMatrix();
    // glTranslatef(.5,.7,1.6);
    glDrawPixels(ImWidth,ImHeight,GL_RGBA,GL_UNSIGNED_BYTE,pixel);
    glPopMatrix();
  }
  if(la!=0){
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glColor4f(1.-Rback,1.-Gback,1.-Bback,1.0);
    glPushMatrix();//Cube
    glTranslated(.5*Edge[0]*InvScaleUn,.5*Edge[1]*InvScaleUn,.5*Edge[2]*InvScaleUn);
    glScalef((GLfloat)Edge[0]*InvScaleUn,
    	     (GLfloat)Edge[1]*InvScaleUn,
    	     (GLfloat)(Edge[2])*InvScaleUn);
    glutWireCube((GLfloat)1.);
    glPopMatrix();//Cube
    glPopAttrib(); 
  }
  if(gr != 0){
    glColor4f(1.-Rback,1.-Gback,1.-Bback,1.0);
    glPushMatrix();//Griglia
    //glTranslated(0.,0.,Edge[2]*InvScaleUn);
    glPopMatrix();//Griglia
    glCallList(Griglia);
  }
  glPopMatrix();//Translation
  glPopMatrix();//Rotation
  glPopMatrix();//Wheel
  //Text
  // glDisable(GL_TEXTURE_2D);
  // glDisable(GL_LIGHTING);
  // glDisable(GL_DEPTH_TEST); 
  // glLoadIdentity();
  // glMatrixMode(GL_PROJECTION);
  // glPushMatrix();
  // glLoadIdentity();
  // glOrtho(0,WinWidth,0,WinHeight,-1.0,1.0);
  if(IfInfo){
    glPushMatrix();//Info
    glColor4f(1.-Rback,1.-Gback,1.-Bback,1.0);
    glRasterPos3d(xi, yi, zi);
    int len = (int)strlen(info);
    for (int i = 0; i < len; i++) {
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
      //glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, info[i]);
    }
    glPopMatrix();//Info
  }
  glCallList(DrLegend);
  glFlush();			//Finish rendering
  glutSwapBuffers();
}
void Draw::Draw1(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();//Translation
  glTranslatef(xp,yp,zp);
  glPushMatrix();//Rotation
  glRotatef(xa,1.,0.,0.);
  glRotatef(ya,0.,1.,0.);
  glRotatef(za,0.,0.,1.);
  glTranslatef(-.5*Edge[0]*InvScaleUn,-.5*Edge[1]*InvScaleUn,-.5*Edge[2]*InvScaleUn);
  if(lu != 0){
    glPushMatrix();//Ligth
    glRotatef(xf,1.,0.,0.);
    glRotatef(yf,0.,1.,0.);
    glRotatef(zf,0.,0.,1.);
    glTranslatef(1.,1.,Edge[2]*InvScaleUn);
    glLightfv(GL_LIGHT0,GL_AMBIENT,luceAmbiente);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,luceDiffusa);
    glLightfv(GL_LIGHT0,GL_SPECULAR,luceSpeculare);
    glLightfv(GL_LIGHT0,GL_POSITION,posizioneLight0);
    glutSolidSphere(0.05,20,20);
    glPopMatrix();//Light
  }
  glPushMatrix();//Cube
  glTranslatef(.5*Edge[0]*InvScaleUn,.5*Edge[1]*InvScaleUn,.5*Edge[2]*InvScaleUn);
  glScalef((GLfloat)Edge[0]*InvScaleUn,
	   (GLfloat)Edge[1]*InvScaleUn,
	   (GLfloat)(Edge[2])*InvScaleUn);
  glColor4f(1.,1.,1.,1.);
  glutWireCube((GLfloat)1);
  glPopMatrix();//Cube
//   glPushMatrix();
//   glTranslatef(.5,.7,8.);
//   glDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,pixel);
//   glPopMatrix();
  if(1==0)
  {
    GLfloat Normal[3]={.1,.1,.1};
    glNormalPointer(GL_FLOAT,sizeof(GLfloat),&Normal[0]);  //Pointer to the first color*/
    int NPoint = 20;
    GLfloat Deltax = 1./(GLfloat) NPoint;
    double *Surface = new double [3*NPoint*NPoint];
    glColor4f(1.0,.0,1.,1.);
    for(int i=0;i<NPoint;i++){
      for(int j=0;j<NPoint;j++){
	Surface[(i*NPoint+j)*3+0] = Deltax*i+drand48()*.01;
	Surface[(i*NPoint+j)*3+1] = Deltax*j+drand48()*.01;
	Surface[(i*NPoint+j)*3+2] = drand48()*.05;
      }
    }
    glVertexPointer(3,GL_DOUBLE,3*sizeof(double),Surface);
    DrTriangles(NPoint);
  }
  {
    //ShowTexture();
  }
  glPushMatrix();
  glTranslatef(.5,.3,.4);
  glColor4f(.0,.0,1.,1.);
  glutSolidSphere(.2,Detail,Detail);
  glPopMatrix();
  glPushMatrix();
  glTranslatef(.5,.7,.6);
  glColor4f(.0,1.0,.0,1.);
  glutSolidSphere(.2,Detail,Detail);
  glPopMatrix();
  glPopMatrix();//Rotation
  glPopMatrix();//Translation
  glutSwapBuffers();
}
void Draw::DrTriangles(int NPoint){
  vector <GLuint> VecIndices; 
  for(int i=0;i<NPoint-1;i++){
    for(int j=0;j<NPoint-1;j++){
      VecIndices.push_back(i*NPoint + j);
      VecIndices.push_back((i+1)*NPoint + j);
      VecIndices.push_back((i+1)*NPoint + (j+1));
      VecIndices.push_back(i*NPoint + j);
      VecIndices.push_back((i+1)*NPoint + (j+1));
      VecIndices.push_back((i)*NPoint + (j+1));
    }
  }
  int NIndex = VecIndices.size();
  GLuint *Indices = new GLuint [NIndex];
  for(int i=0;i<NIndex;i++)Indices[i] = VecIndices[i];
  //  GLuint *Indices; Indices = &VecIndices[0];
  VecIndices.clear();
  glDrawElements(GL_TRIANGLES,NIndex,GL_UNSIGNED_INT,Indices);
  free(Indices);
}
void Draw::ShowImage(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glDrawPixels(ImWidth,ImHeight,GL_RGBA,GL_UNSIGNED_BYTE,pixel);
  glutSwapBuffers();
}

int Draw::ApplyTexture(){
  glEnable(GL_TEXTURE_2D);
  printf("Texture %d\n",glIsEnabled(GL_TEXTURE_2D));
  if(OpenImage("Texture.tif")==0) {/*printf("Ciccia!\n")*/;return 0;}
  glGenTextures(1,&Texture);
  glBindTexture(1,Texture);
  glTexImage2D(GL_TEXTURE_2D, 0, 4, ImWidth, ImHeight , 0, GL_RGBA, GL_UNSIGNED_BYTE, pixel);
  //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  //glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
  //gluBuild2DMipmaps(GL_TEXTURE_2D, 4, infoheader.biWidth, infoheader.biHeight, GL_RGBA, GL_UNSIGNED_BYTE, l_texture);  
  return 1;
}
int Draw::ShowTexture(){
  glBindTexture(GL_TEXTURE_2D,Texture);
  glBegin(GL_QUADS);
  glTexCoord2f(0.,0.);
  glTexCoord2f(0.,0.);glVertex3f(0.,0.,1.0);
  glTexCoord2f(1.,0.);glVertex3f(1.,0.,1.0);
  glTexCoord2f(1.,1.);glVertex3f(1.,1.,1.0);
  glTexCoord2f(0.,1.);glVertex3f(0.,1.,1.0);
  glEnd();
  glRasterPos3f(0.,0.,0.);
  glDrawPixels(ImWidth,ImHeight,GL_RGBA,GL_UNSIGNED_BYTE,pixel);
}
void Draw::Transform(){
  glLineWidth(2);
  glDisable(GL_LIGHTING);
  Vettore u(.3,.1,.6);
  Vettore v(.1,.5,.4);
  Vettore w(u.x[0]-v.x[0],u.x[1]-v.x[1],u.x[2]-v.x[2]);
  Vettore o(.5,.5,.5);
  //Vettore o(.0,.0,.0);
  Vettore n(3);
  Vettore Norm1(3);
  Vettore Norm2(3);
  Vettore PS(3);
  double Angle = 90.;
  Quadri qx(u.x,Angle*DEG_RAD);
  Quadri qy(v.x,Angle*DEG_RAD);
  Quadri qz(w.x,Angle*DEG_RAD);
  //glGetDoublev(GL_MODELVIEW_MATRIX,MatrOld);
  qx.Prod(qy);
  qx.Prod(qz);
  Matrice Mx(qx,4);
  Matrice My(qy,4);
  Matrice Mz(qz,4);
  //Mx.Mult(My);
  //Mx.Mult(Mz);
  la = 0;
  IfScript = 0;
  n.Copy(&v);
  n.ProjOnAxis(&u);
  //Mx.Print();
  //glMatrixMode(GL_MODELVIEW);
  glMatrixMode(GL_PROJECTION);
  glDeleteLists(Particles,1);
  Particles = glGenLists(1);
  glNewList(Particles,GL_COMPILE);
  glColor4f( .1,.1,.1,1.);
  glPushMatrix();//x
  glBegin(GL_LINES);
  glVertex3d(0.,0.,0.);
  glVertex3d(1.,0.,0.);
  glEnd();
  glPopMatrix();//Nano
  glPushMatrix();//y
  glBegin(GL_LINES);
  glVertex3d(0.,0.,0.);
  glVertex3d(0.,1.,0.);
  glEnd();
  glPopMatrix();//Nano
  glPushMatrix();//z
  glBegin(GL_LINES);
  glVertex3d(0.,0.,0.);
  glVertex3d(0.,0.,1.);
  glEnd();
  glPopMatrix();//Nano
  //origin
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(0.,0.,0.);
  glVertex3d(o[0],o[1],o[2]);
  glEnd();
  glPopMatrix();//Nano
  glColor4f( 1.,.1,.1,1.);
  //glMultMatrixd(Mx.data);
  PutString(o[0]-.05,o[1]-.05,o[2]-.05,"o");
  //axis
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(o[0],o[1],o[2]);
  glVertex3d(o[0]+u[0],o[1]+u[1],o[2]+u[2]);
  glEnd();
  glPopMatrix();//Nano
  PutString(o[0]+u[0],o[1]+u[1],o[2]+u[2],"axis");
  glPushMatrix();//Nano
  //v
  glColor4f( .1,.1,1.,1.);
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(o[0],o[1],o[2]);
  glVertex3d(o[0]+v[0],o[1]+v[1],o[2]+v[2]);
  glEnd();
  glPopMatrix();//Nano
  PutString(o[0]+v[0],o[1]+v[1],o[2]+v[2],"v");
  glPushMatrix();//Nano
  //n
  glColor4f( 1.,1.,.1,1.);
  glPopMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(o[0]+.1,o[1],o[2]);
  glVertex3d(o[0]+n[0]+.1,o[1]+n[1],o[2]+n[2]);
  glEnd();
  glPopMatrix();//Nano
  PutString(o[0]+n[0]+.1,o[1]+n[1],o[2]+n[2],"n");
  //glCallList(Cylinder);
  n.PerpTo(&v,&u);
  //glPopMatrix();//Nano
  glPushMatrix();//Nano
  //v-n
  glColor4f( .1,1.,.1,1.);
  //glPopMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(o[0]+v[0],o[1]+v[1],o[2]+v[2]);
  glVertex3d(o[0]+v[0]+n[0],o[1]+v[1]+n[1],o[2]+v[2]+n[2]);
  glEnd();
  glPopMatrix();//Nano
  PutString(o[0]+v[0]+n[0],o[1]+v[1]+n[1],o[2]+v[2]+n[2],"v-n");
  //v'
  glColor4f(.5,.5,1.,1.);
  //glPopMatrix();//Nano
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(o[0],o[1],o[2]);
  glVertex3d(o[0]+v[0]+2.*n[0],o[1]+v[1]+2.*n[1],o[2]+v[2]+2.*n[2]);
  glEnd();
  glPopMatrix();//Nano
  PutString(o[0]+v[0]+2.*n[0],o[1]+v[1]+2.*n[1],o[2]+v[2]+2.*n[2],"v'");
  //normal
  Norm1.Normal(&v,&u);
  glColor4f(1.,.1,1.,1.);
  //glPopMatrix();//Nano
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(o[0],o[1],o[2]);
  glVertex3d(o[0]+Norm1[0],o[1]+Norm1[1],o[2]+Norm1[2]);
  glEnd();
  glPopMatrix();//Nano
  PutString(o[0]+Norm1[0],o[1]+Norm1[1],o[2]+Norm1[2],"normal");
  //normal
  Norm2.Normal(&u,&Norm1);
  glColor4f(1.,.1,1.,1.);
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(o[0],o[1],o[2]);
  glVertex3d(o[0]+Norm2[0],o[1]+Norm2[1],o[2]+Norm2[2]);
  glEnd();
  glPopMatrix();//Nano
  PutString(o[0]+Norm2[0],o[1]+Norm2[1],o[2]+Norm2[2],"normal");
  //plane
  Vettore v1(.3,.6,.8);
  //Vettore v2(.2,.7,.3);//out
  Vettore v2(.4,.2,.7);//in
  Vettore v3(.5,.7,.9);
  Piano p1(&v1,&v2,&v3);
  glColor4f(.5,.1,1.,1.);
  Vettore PSurf = p1.ProjOnSurf(&o);
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(PSurf[0],PSurf[1],PSurf[2]);
  glVertex3d(o[0],o[1],o[2]);
  glEnd();
  glColor4f(.8,.1,1.,1.);
  glBegin(GL_POLYGON);
  glVertex3d(p1.P1[0],p1.P1[1],p1.P1[2]);
  glVertex3d(p1.P2[0],p1.P2[1],p1.P2[2]);
  glVertex3d(p1.P3[0],p1.P3[1],p1.P3[2]);
  glEnd();
  glPopMatrix();//Nano  
  glPushMatrix();//Nano
  if(p1.IsOnSurf(&PSurf)){
    PutString(PSurf[0],PSurf[1],PSurf[2],"in");
  }
  else{
    PutString(PSurf[0],PSurf[1],PSurf[2],"out");
  }
  glPopMatrix();//Nano  
  Vettore Vel(.1,.2,.2);
  //Vettore Vel1 = p1.Reflect(&Vel);
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(PSurf[0],PSurf[1],PSurf[2]);
  glVertex3d(PSurf[0]+Vel[0],PSurf[1]+Vel[1],PSurf[2]+Vel[2]);
  glEnd();  
  glPopMatrix();//Nano  
  glPushMatrix();//Nano
  PutString(PSurf[0]+Vel[0],PSurf[1]+Vel[1],PSurf[2]+Vel[2],"vel");
  glPopMatrix();//Nano  
  p1.Impact(&PSurf,&Vel);
  glPushMatrix();//Nano
  glBegin(GL_LINES);
  glVertex3d(PSurf[0],PSurf[1],PSurf[2]);
  glVertex3d(PSurf[0]+Vel[0],PSurf[1]+Vel[1],PSurf[2]+Vel[2]);
  glEnd();  
  glPopMatrix();//Nano  
  glPushMatrix();//Nano
  PutString(PSurf[0]+Vel[0],PSurf[1]+Vel[1],PSurf[2]+Vel[2],"vel'");
  glPopMatrix();//Nano  
  //
  glEndList();
}
void Draw::DrCube(){
  glLineWidth(3);
  glDeleteLists(Particles,1);
  Particles = glGenLists(1);
  glNewList(Particles,GL_COMPILE);
  glColor4f( 1.,.1,.1,1.);
  glEndList();
}
void Draw::PutString(double *Pos,char *String){
  glRasterPos3f((Pos[0]*InvScaleUn),
		(Pos[1]*InvScaleUn),
		(Pos[2]*InvScaleUn));
  for(int l=0; l < strlen(String); l++) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, String[l]);
    //glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,String[l]);
  }
  glCallList(Point);
}
void Draw::PutString(double Posx,double Posy,double Posz,char *String){
  glRasterPos3f((Posx*InvScaleUn),
		(Posy*InvScaleUn),
		(Posz*InvScaleUn));
  for(int l=0; l < strlen(String); l++) {
    //glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, String[l]);
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,String[l]);
  }
  glCallList(Point);
}
#endif// __glut_h__
