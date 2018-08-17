#include "Forces.h"
#ifdef __glut_h__
#include "../include/Draw.h"
Draw *Dr;
//OpenGL
int Forces::Graphics(int argc,char **argv){
  Dr = new Draw();
  for(int d=0;d<3;d++) Dr->Edge[d] = pEdge(d);
  Dr->InvScaleUn = pInvEdge(0);
  Dr->Window(argc,argv);
  // Dr->xp = 0.;Dr->zp = 0.;
  Dr->xa = -60.;Dr->ya = -00.;Dr->za = -38.;
  Dr->xf = 160.;Dr->yf = -25.;Dr->zf = 0.;
  Dr->xi = -.4;Dr->yi = -.4;Dr->zi = 0.;
  for(int n=0;n<pNNano();n++){
    Cylinder[n] = Dr->DefCylinder(Nano[n].Rad,Nano[n].Height);
  }
  glutPostRedisplay();
  glutMainLoop();
  //Interp();
}
void Figure(){
  //Dr->Draw1();
  Dr->DFigure();
  //Dr->DMinimal();
}
void ParticleRealTime(){
  return ;
  Dr->OpenImage("RadDistrNanoR1_0H0_4h32_0.tif");
  glPopMatrix();
  glRasterPos3d(-.35,-.4,-.32);
  glPixelZoom(.5,.5);
  glDrawPixels(Dr->ImWidth,Dr->ImHeight,GL_RGBA,GL_UNSIGNED_BYTE,Dr->pixel);
  glPushMatrix();
}
void reshape(int w,int h){
  Dr->Dreshape(w,h);
}
void Timer(int v){
  Dr->DTimer(v);
}
void MouseMove(int x,int y){
  Dr->DMouseMove(x,y);
}
void mouse(int button, int state,int x,int y){
  Dr->Dmouse(button,state,x,y);
};
void special(int k, int x, int y){
  Dr->Dspecial(k,x,y);
}
void Forces::DrawScene(){
  if(VAR_IF_TYPE(SysShape,SYS_2D) && !IfLine)
    DrawCarpet();
  else 
    DrawParticles();
}
void Forces::DynamicsView(void){
  double Num=0.;
  Dr->Diap++;
  Dr->tDiap=glutGet(GLUT_ELAPSED_TIME);
  double Speed=25.;
  double fps = Dr->Diap*1000.0/(Dr->tDiap-Dr->tDiapBase);
  double Ratio = 0.;
  double NRatio = 0.;
  if (Dr->tDiap - Dr->tDiapBase > 1000) {
    fprintf(stderr,"Frame per second %lf\r",fps);
    Dr->tDiapBase = Dr->tDiap;		
    Dr->Diap = 0;
  }
  if(Dr->tDiap - Dr->tDiapBase > 1000/Speed){
    for(int i=0;i<NUpdate;i++){
      Dynamics();
      //Interp();
      if(VAR_IF_TYPE(CalcMode,CALC_NcVT))
	Dr->Step += pNPCh()-1;
    }
    if(IfMovie && (pStep()%5)==0){
      Dr->Step++;
      Dr->Picture();
    }
    Dr->za += 1.;
    CurrTime = time(NULL);
    Ratio += Dr->Step/(1000.*(CurrTime-InitTime));
    NRatio += 1.;
    double v2 = 0.;
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
	v2 += SQR(pVel(p,d));
      }
    }
    sprintf(Dr->info,"NPart %d loop/s %lf acc/step %lf in/out %lf T %lf Pot %lf",pNPart(),Ratio/NRatio,(NRemoval+NInsertion)/(double)pStep(),NInsertion/(double)NRemoval,v2/(double)(3*pNPart()),OldNrgSys);
    DrawScene();
    glutPostRedisplay();
    glutTimerFunc(1000, Timer, 1);
  }
}
void Forces::DrawSoil(){
  if(!VAR_IF_TYPE(SysShape,SYS_ELECTRO)){
    return;
  }
  glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);  
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHTING);
  double InvScaleUn = Dr->InvScaleUn;
  int NBin = 30;
  double InvNBin = 1./(double)NBin;
  Vettore v1(3);
  Vettore v2(3);
  Vettore v3(3);
  Vettore v4(3);
  Vettore vN(3);
  int Typ = 1;
  for(int sx=0;sx<NBin-1;sx++){
    for(int sy=0;sy<NBin-1;sy++){
      v1.x[0] = sx*pEdge(0)*InvNBin;
      v2.x[0] = (sx+1)*pEdge(0)*InvNBin;
      v3.x[0] = (sx+1)*pEdge(0)*InvNBin;
      v4.x[0] = (sx)*pEdge(0)*InvNBin;
      v1.x[1] = sy*pEdge(1)*InvNBin;
      v2.x[1] = (sy)*pEdge(1)*InvNBin;
      v3.x[1] = (sy+1)*pEdge(1)*InvNBin;
      v4.x[1] = (sy+1)*pEdge(1)*InvNBin;
      v1.x[2] =  .2*pEdge(2)*sin(6*v1.x[0]*pInvEdge(0))*sin(4*v1.x[1]*pInvEdge(1)) + .4*pEdge(2);
      v2.x[2] =  .2*pEdge(2)*sin(6*v2.x[0]*pInvEdge(0))*sin(4*v2.x[1]*pInvEdge(1)) + .4*pEdge(2);
      v3.x[2] =  .2*pEdge(2)*sin(6*v3.x[0]*pInvEdge(0))*sin(4*v3.x[1]*pInvEdge(1)) + .4*pEdge(2);
      v4.x[2] =  .2*pEdge(2)*sin(6*v4.x[0]*pInvEdge(0))*sin(4*v4.x[1]*pInvEdge(1)) + .4*pEdge(2);
      for(int d=0;d<3;d++){
      	v1.x[d] *= InvScaleUn;
      	v2.x[d] *= InvScaleUn;
      	v3.x[d] *= InvScaleUn;
      	v4.x[d] *= InvScaleUn;
      }
      glColor4fv(ColorType[Typ]);
      glLightfv(GL_LIGHT0, GL_AMBIENT,ColorType[Typ]);
      glLightfv(GL_LIGHT1, GL_DIFFUSE,ColorType[Typ]);
      glLightfv(GL_LIGHT1, GL_POSITION,ColorType[Typ]);
      vN.Normalize();
      glPushMatrix();//Particle
      glBegin(GL_POLYGON);
      vN = v1^v2;
      glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
      glVertex3d(v1.x[0],v1.x[1],v1.x[2]);
      vN = v2^v3;
      glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
      glVertex3d(v2.x[0],v2.x[1],v2.x[2]); 
      vN = v3^v4;
      glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
      glVertex3d(v3.x[0],v3.x[1],v3.x[2]);
      vN = v4^v1;
      glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
      glVertex3d(v4.x[0],v4.x[1],v4.x[2]);
      glEnd();
      glPopMatrix();//Particle
    }
  }
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHT1);
  glDisable(GL_LIGHTING);
}
void Forces::DrawCarpet(){
  if(!VAR_IF_TYPE(SysShape,SYS_2D)){
    return;
  }
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  // glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);  
  // glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  // glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHTING);
  //glFrontFace(GL_CCW);
  //glShadeModel(GL_SMOOTH);
  double InvScaleUn = Dr->InvScaleUn;
  Vettore v1(3);
  Vettore v2(3);
  Vettore v3(3);
  Vettore vN(3);
  for(int sx=0;sx<nEdge[0]-1;sx++){
    for(int sy=0;sy<nEdge[1]-1;sy++){
      int pRight[3] = {sx*nEdge[1]+sy,(sx+1)*nEdge[1]+sy,(sx+1)*nEdge[1]+sy+1};
      for(int d=0;d<3;d++){
	v1.x[d] = pPos(pRight[0],d)*InvScaleUn;
	v2.x[d] = pPos(pRight[1],d)*InvScaleUn;
	v3.x[d] = pPos(pRight[2],d)*InvScaleUn;
      }
      int Typ = pType(sx*pNChain()+sy);
      glColor4fv(ColorType[Typ]);
      glLightfv(GL_LIGHT0, GL_AMBIENT,ColorType[Typ]);
      glLightfv(GL_LIGHT1, GL_DIFFUSE,ColorType[Typ]);
      //glLightfv(GL_LIGHT1, GL_POSITION,ColorType[Typ]);
      //vN = v1^v2;
      vN.x[0] = v1.x[1]*v2.x[2] - v1.x[2]*v2.x[1];
      vN.x[1] = v1.x[2]*v2.x[0] - v1.x[0]*v2.x[2];
      vN.x[2] = v1.x[0]*v2.x[1] - v1.x[1]*v2.x[0];
      vN.Normalize();
      glPushMatrix();//Particle
      //glColor4f(.2,1.,.2,1.);
      glBegin(GL_TRIANGLES);
      glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
      glVertex3d(v1.x[0],v1.x[1],v1.x[2]);
      glVertex3d(v2.x[0],v2.x[1],v2.x[2]);
      glVertex3d(v3.x[0],v3.x[1],v3.x[2]);
      glEnd();
      glPopMatrix();//Particle
      int pLeft[3] = {sx*pNChain()+sy,(sx+1)*pNChain()+sy+1,sx*pNChain()+sy+1};
      for(int d=0;d<3;d++){
	v1.x[d] = pPos(pLeft[0],d)*InvScaleUn;
	v2.x[d] = pPos(pLeft[1],d)*InvScaleUn;
	v3.x[d] = pPos(pLeft[2],d)*InvScaleUn;
      }
      vN.x[0] = v1.x[1]*v2.x[2] - v1.x[2]*v2.x[1];
      vN.x[1] = v1.x[2]*v2.x[0] - v1.x[0]*v2.x[2];
      vN.x[2] = v1.x[0]*v2.x[1] - v1.x[1]*v2.x[0];
      //      vN = v1^v2;
      vN.Normalize();
      glPushMatrix();//Particle
      //glColor4f(.2,.5,.2,1.);
      glBegin(GL_TRIANGLES);
      glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
      glVertex3d(v1.x[0],v1.x[1],v1.x[2]);
      glVertex3d(v2.x[0],v2.x[1],v2.x[2]);
      glVertex3d(v3.x[0],v3.x[1],v3.x[2]);
      glEnd();
      glPopMatrix();//Particle
    }
  }
  DrawNano();
  glEndList();
}
void Forces::DrawParticles(){
  double InvScaleUn = Dr->InvScaleUn;
  double Diameter = .01;
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  glDisable(GL_LIGHTING);
  glLineWidth(3);
  glPointSize(5);
  for(int p=0,link=0;p<pNPart();p++){
    float Red = Pm[p].Pos[2]*.5*pInvEdge(2);
    int Typ = pType(p) < 6 ? pType(p) : 5;
    glColor4f(Red*ColorType[Typ][0],ColorType[Typ][1],ColorType[Typ][2],1.);
    //glColor4fv(ColorType[Typ]);
    //Dr->Numera(Pm[p].Pos,p);
    if(IfLine && Ln[p].NLink >0) DrBondLine(p);
    if(Dr->IfPoint){
      glPushMatrix();//Particle
      //glBegin(GL_POINTS);	  
      glTranslatef((GLfloat)(pPos(p,0)*InvScaleUn),
		   (GLfloat)(pPos(p,1)*InvScaleUn),
		   (GLfloat)(pPos(p,2)*InvScaleUn));
      //else if(IfColour) glCallList(Quad);
      int Typ = pType(p) < 6 ? pType(p) : 5;
      glColor4fv(ColorType[Typ]);
      if(!IfSphere) glCallList(Dr->Point);
      else glutSolidSphere(Diameter,20,20);
      glPopMatrix();//Particle
    }
  }
  if(VAR_IF_TYPE(SysShape,SYS_LEAVES)){
    glBegin(GL_LINES);
    glColor4f(.7,.7,.0,1.);
    glPushMatrix();//Line
    glVertex3f(0.,.5,.6);
    glVertex3f(1.,.5,.6);
    glPopMatrix();//Line    
    glPushMatrix();//Line
    glVertex3f(0.,.5,.4);
    glVertex3f(1.,.5,.4);
    glPopMatrix();//Line    
    glPushMatrix();//Line
    glVertex3f((GLfloat)(pPos(0,0)*InvScaleUn),
	       (GLfloat)(pPos(0,1)*InvScaleUn),
	       (GLfloat)(pPos(0,2)*InvScaleUn));
    glVertex3f((GLfloat)(pPos(NEdge-1,0)*InvScaleUn),
	       (GLfloat)(pPos(NEdge-1,1)*InvScaleUn),
	       (GLfloat)(pPos(NEdge-1,2)*InvScaleUn));
    glPopMatrix();//Line    
    glPushMatrix();//Line
    glVertex3f((GLfloat)(Pm[NEdge].Pos[0]*InvScaleUn),
	       (GLfloat)(Pm[NEdge].Pos[1]*InvScaleUn),
	       (GLfloat)(Pm[NEdge].Pos[2]*InvScaleUn));
    glVertex3f((GLfloat)(Pm[2*NEdge-1].Pos[0]*InvScaleUn),
	       (GLfloat)(Pm[2*NEdge-1].Pos[1]*InvScaleUn),
	       (GLfloat)(Pm[2*NEdge-1].Pos[2]*InvScaleUn));
    glPopMatrix();//Line    
    glEnd();
  }
  if(IfLine){
    for(int p=0;p<NShow-1;p++){
      //printf("%d (%lf %lf %lf)\n",p,Pl[p].Pos[0],Pl[p].Pos[1],Pl[p].Pos[2]);
      if(Pl[p].Idx == NEdge-2) continue;
      glPushMatrix();//Particle
      glColor4f(1.,.0,.0,1.);
      //glBegin(GL_POINTS);	  
      glBegin(GL_LINES);
      glVertex3f((GLfloat)(Pl[p].Pos[0]*InvScaleUn),
		 (GLfloat)(Pl[p].Pos[1]*InvScaleUn),
		 (GLfloat)(Pl[p].Pos[2]*InvScaleUn));
      glVertex3f((GLfloat)(Pl[p+1].Pos[0]*InvScaleUn),
		 (GLfloat)(Pl[p+1].Pos[1]*InvScaleUn),
		 (GLfloat)(Pl[p+1].Pos[2]*InvScaleUn));
      glEnd();
////           glutSolidSphere(Diameter,20,20);
//           glTranslatef((GLfloat)(Pl[p].Pos[0]*InvScaleUn),
//       		 (GLfloat)(Pl[p].Pos[1]*InvScaleUn),
//       		 (GLfloat)(Pl[p].Pos[2]*InvScaleUn));
//           glCallList(Point);
      glPopMatrix();//Particle
    }
  }
  DrawNano();
  DrawSoil();
  glEndList();
}
void Forces::DrBondLine(int p){
  double InvScaleUn = Dr->InvScaleUn;
  if(pNLink()==0)return;
  glDisable(GL_LIGHTING);
  for(int l=0;l<Ln[p].NLink;l++){
    int link = Ln[p].Link[l];
    double Pos[3] = {pPos(link,0),pPos(link,1),pPos(link,2)};
    double Mid[3];
    int IfContinue = 1;
    for(int d=0;d<3;d++){
      double Bf = floor((pPos(p,d) - pPos(link,d))*pInvEdge(d)+.5)*pEdge(d);
      if(fabs(Bf) > 0. & !PeriodicImage[d]) IfContinue = 0;
      else Pos[d] += Bf;
      Mid[d] = .5*(pPos(p,d)+Pos[d])*InvScaleUn;
    }
    //if(!IfContinue) continue;
    glPushMatrix();//Line
    glBegin(GL_LINES);
    glColor4fv(ColorType[pType(p)]);
    glVertex3d((pPos(p,0)*InvScaleUn),
	       (pPos(p,1)*InvScaleUn),
	       (pPos(p,2)*InvScaleUn));
    glVertex3d(Mid[0],Mid[1],Mid[2]);
    glColor4fv(ColorType[pType(link)]);
    glVertex3d(Mid[0],Mid[1],Mid[2]);
    glVertex3d(((Pos[0])*InvScaleUn),
	       ((Pos[1])*InvScaleUn),
	       ((Pos[2])*InvScaleUn));
    glEnd();
    glPopMatrix();//Line
   }
}
void Forces::DrawNano(){
  glDisable(GL_LIGHTING);
  for(int n=0;n<pNNano();n++){
    if(Nano[n].Shape == 0) continue;
    Vettore Ax(Nano[n].Axis[0],Nano[n].Axis[1],Nano[n].Axis[2]);
    Vettore Zed(0.,0.,1.);
    double Angle = Zed.Angle(&Ax,&Zed);
    Quadri q(Ax.x,Angle);
    Matrice M(q,4);
    glPushMatrix();//Nano
    glTranslated((pNanoPos(n,0)*Dr->InvScaleUn),
		 (pNanoPos(n,1)*Dr->InvScaleUn),
		 (pNanoPos(n,2)*Dr->InvScaleUn));
    glColor4f( 1.,.1,.1,1.);
    if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_SPH)){
      glutSolidSphere(Nano[n].Rad,20,20);
    }
    else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CYL)){
      glCallList(Cylinder[n]);
    }
    glPopMatrix();//Nano
  }
}
void Forces::Menu(){
  //submenu = glutCreateMenu(processEvent);
  glutAddMenuEntry("SottoPrima",3);

  //menu = glutCreateMenu(processEvent);
  glutAddMenuEntry("Move",1);
  glutAddMenuEntry("Stop",2);
  glutAddSubMenu("Altro",submenu);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}
void Forces::keyboard(unsigned char key,int x, int y){
  double StepDiameter = 0.01;
  switch (key){
  case 'a':
    CreateInitial();
    DrawScene();
    glutPostRedisplay();
    break;
  case 'b':
    IfLine += 1;
    if(IfLine == 2) IfLine =0;
    sprintf(Dr->info,"Bonding visualisation");
    DrawScene();
    glutPostRedisplay();
    break;
  case 'c':
    IfExt++;
    if(IfExt==0)
      sprintf(Dr->info,"Increasing radius");
    else if(IfExt==1)
      sprintf(Dr->info,"Increasing height");
    else if(IfExt==2)
      sprintf(Dr->info,"Increasing angle");
    else if(IfExt==3)
      sprintf(Dr->info,"D^4 term");
    else if(IfExt==4)
      sprintf(Dr->info,"D^2 term");
    else if(IfExt==5)
      sprintf(Dr->info,"Elastic term");
    else if(IfExt==6)
      sprintf(Dr->info,"Moving center ->x");
    else if(IfExt==7)
      sprintf(Dr->info,"Moving center ->z");
    else 
      IfExt = -1;
    break;
  case 'd':
    if(IfExt == 0){
      for(int n=0;n<pNNano();n++)
	Nano[n].Rad += StepDiameter;
      sprintf(Dr->info,"Nano->Rad %lf",Nano->Rad);
    }
    else if(IfExt == 1){
      for(int n=0;n<pNNano();n++)
	Nano[n].Height += StepDiameter;
      sprintf(Dr->info,"Nano->Height %lf",Nano->Height);
    }
    else if(IfExt == 2){
      for(int n=0;n<pNNano();n++)
	Nano[n].Hamaker += 5.;
      if(Nano->Hamaker >= 90.) Nano->Hamaker = 0.;
      sprintf(Dr->info,"ExtAngle %lf",Nano->Hamaker);
    }
    else if(IfExt == 3){
      Kf.SLap *= 10.;
      sprintf(Dr->info,"D^4 %lf ratio %lf",Kf.SLap,pow(Kf.SLap/Kf.El[2],.25));
    }
    else if(IfExt == 4){
      Kf.Lap += 1.;
      sprintf(Dr->info,"D^2 %lf",Kf.Lap);
    }
    else if(IfExt == 5){
      Kf.El[2] += 10.;
     sprintf(Dr->info,"elastic %lf ratio %lf",Kf.El[2],pow(Kf.SLap/Kf.El[2],.25));
    }
    else if(IfExt == 6){
      //for(int n=0;n<pNNano();n++)
	Nano[0].Pos[0] += StepDiameter;
      sprintf(Dr->info,"Nano->Pos x %lf",Nano->Pos[0]);
    }
    else if(IfExt == 7){
      for(int n=0;n<pNNano();n++)
	Nano[n].Pos[2] += StepDiameter;
      sprintf(Dr->info,"Nano->Pos z %lf",Nano->Pos[2]);
    }
    for(int n=0;n<pNNano();n++)
      Cylinder[n] = Dr->DefCylinder(Nano[n].Rad,Nano[n].Height);
    AddRigid();
    IfFillMatrix = 1;
    glutPostRedisplay();
    break;
  case 'D':
    if(IfExt == 0){
      for(int n=0;n<pNNano();n++)
      Nano[n].Rad -= StepDiameter;
      sprintf(Dr->info,"Rad %lf",Nano->Rad);
    }
    else if(IfExt == 1){
      for(int n=0;n<pNNano();n++)
	Nano[n].Height -= StepDiameter;
      sprintf(Dr->info,"Height %lf",Nano->Height);
    }
    else if(IfExt == 2){
      if(Nano->Hamaker < 0.) Nano->Hamaker = 0.;
      Nano->Hamaker -= 5.;
      sprintf(Dr->info,"ExtAngle %lf",Nano->Hamaker);
    }
    else if(IfExt == 3){
      Kf.SLap /= 10.;
      sprintf(Dr->info,"D^4 %lf ratio %lf",Kf.SLap,pow(Kf.SLap/Kf.El[2],.25));
    }
    else if(IfExt == 4){
      Kf.Lap -= 1.;
      sprintf(Dr->info,"D^2 %lf",Kf.Lap);
    }
    else if(IfExt == 5){
      Kf.El[2] -= 10.;
      sprintf(Dr->info,"elastic %lf ratio %lf",Kf.El[2],pow(Kf.SLap/Kf.El[2],.25));
    }
    else if(IfExt == 6){
      //for(int n=0;n<pNNano();n++)
	Nano[0].Pos[0] -= StepDiameter;
      sprintf(Dr->info,"Nano->Pos x %lf",Nano->Pos[0]);
    }
    else if(IfExt == 7){
      for(int n=0;n<pNNano();n++)
	Nano[n].Pos[2] -= StepDiameter;
      sprintf(Dr->info,"Nano->Pos z %lf",Nano->Pos[2]);
    }
    for(int n=0;n<pNNano();n++)
      Cylinder[n] = Dr->DefCylinder(Nano[n].Rad,Nano[n].Height);
    AddRigid();
    IfFillMatrix = 1;
    glutPostRedisplay();
    break;
  case 'e':
    IfSpline=1;
    IfInterp++;
    if(IfInterp == 7)
      IfInterp = 0;
    Interp();
    //Dynamics();
    //printf("IfInterp %d\n",IfInterp);
    glutPostRedisplay();    
    break;
  case 'I':
    Info();
    break;
  case 'm':
    //Menu();
    Bead2Move++;
    if(Bead2Move > pNPart()) Bead2Move = 0;
    SelectBead(Bead2Move);
    sprintf(Dr->info,"Part2Move %d",Bead2Move);
    DrawScene();
    glutPostRedisplay();
    break;
  case 'M':
    Bead2Move--;
    if(Bead2Move < 0) Bead2Move = pNPart();
    SelectBead(Bead2Move);
    sprintf(Dr->info,"Part2Move %d",Bead2Move);
    DrawScene();
    glutPostRedisplay();
    break;
  case 'o':
    Dynamics();
    DrawScene();
    glutPostRedisplay();
    break;
  case 'O':
    InitTime = time(NULL);
    glutIdleFunc(DynamicsMotion);
    break;
  case 'r':
    glutIdleFunc(NULL);
    Dr->InitConstant();
    sprintf(Dr->info,"initial configuration");
    glutPostRedisplay();
    break;
  case 'R':
    ReadConfDinamica(ConfFile);
    IfFillMatrix = 1;
    FillMatrix();
    sprintf(Dr->info,"reload configuration");
    glutPostRedisplay();
    break;
  case 'S':
    Solve();
    DrawScene();
    glutPostRedisplay();
    break;
  case 't':
    PullBead();
    DrawScene();
    glutPostRedisplay();
    break;
  case 'T':
    PushBead();
    //sprintf(info,"Perspective view");
    DrawScene();
    glutPostRedisplay();
    break;
  case 'u':
    Dynamics();
    DrawScene();
    glutPostRedisplay();
    break;
  case 'v':
    BeadType++;
    Pm[Bead2Move].Typ = BeadType;
    if(BeadType >= 4)
      BeadType = 0;
    glutPostRedisplay();
    break;
  case 'V':
    Pm[Bead2Move].Typ = 0;
    glutPostRedisplay();
    break;
  case 'w':
    VAR_REM_TYPE(SysType,VAR_SYS_TXVL);
    VAR_ADD_TYPE(SysType,VAR_SYS_XVT);
    SysFormat = VAR_SYS_TXVL;
    SetNBlock(1);
    Block[0].NChain = pNChain();
    Block[0].InitIdx = 0;
    Block[0].NPCh = pNPCh();
    Block[0].EndIdx = pNPart();
    sprintf(Block[0].Name,"GAS");
    char FileName[60];
    sprintf(FileName,"Trajectory%09d.dat",pStep());
    Write(FileName);
    break;
  case 27:
    exit(0);
    break;
  case 40:
    break;
  default:
    break;
  }
  Dr->keyboardDraw(key);
}
int Forces::Interp(){
  if(!IfSpline) return 0;
  if(IfInterp==FIT_SPLINE3){
    NShow = InterSpline3(Pm,Pl,pNPart(),NSpline);
    sprintf(Dr->info,"Interpolating via Spline3 %d",NShow);
  }
  else if(IfInterp==FIT_SPLINE4){
    NShow = InterSpline4(Pm,Pl,pNPart(),NSpline);
    sprintf(Dr->info,"Interpolating via Spline4 %d",NShow);
  }
  else if(IfInterp==FIT_PARAB2){
    NShow = InterParab2(Pm,Pl,pNPart(),NSpline);
    sprintf(Dr->info,"Interpolating via Parabola2 %d",NShow);
  }
  else if(IfInterp==FIT_CUBIC){
    NShow = InterCubica(Pm,Pl,pNPart(),NSpline);
    sprintf(Dr->info,"Interpolating via Cubic %d",NShow);
  }
  else if(IfInterp==FIT_FORTH){
    NShow = InterForth(Pm,Pl,pNPart(),NSpline);
    sprintf(Dr->info,"Interpolating via Forth %d",NShow);
  }
  else if(IfInterp==FIT_BSPLINE){
    NShow = InterBSpline(Pm,Pl,pNPart(),NSpline);
    sprintf(Dr->info,"Interpolating via BSpline %d",NShow);
  }
  else if(IfInterp==FIT_POLY){
    NShow = InterPoly(Pm,Pl,pNPart(),NSpline);
    sprintf(Dr->info,"Interpolating via %dth grade poly %d",pNPart()-1,NShow);
  }
  return 0;
}
#endif //__glut_h__
