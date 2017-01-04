/***********************************************************************
ElPoly:This progam provide a graphical visualisation of the data 
opend by VarData using openGL glut. The most important option are 
the possibility of changing the backfold of the polymers with 'c', 
see the subsequent file in the list with '>', see the bond with 'b'. 
Copyright (C) 2008 by Giovanni Marelli <sabeiro@virgilio.it>


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
***********************************************************************/
#include "ElPoly.h"

#ifdef __glut_h__
Draw *Dr;
void Figure(){
  //Dr->Draw1();
  Dr->DFigure();
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
void special(int k, int x, int y){
  Dr->Dspecial(k,x,y);
}
int ElPoly::Graphics(int argc,char **argv){
  InvScaleUn = 1./pEdge(0);
  Dr = new Draw();
  for(int d=0;d<3;d++) Dr->Edge[d] = pEdge(d);
  Dr->InvScaleUn = InvScaleUn;
  Dr->xi = -.4;Dr->yi = -.4;Dr->zi = .0;
  Dr->xp = 0.;Dr->yp = 0.;Dr->zp = 0.;
  Dr->xa = 0.;Dr->ya = 0.;Dr->za = 0.;  
  Vicinanze = 0.;
  Saturation = 1.;
  // IfColour = 0;
  // IfChains = 0;
  Quad = 0;
  Point = 0;
  Cylinder = 0;
  MetalCylinder = 0;
  Hexagon = 0;
  Cube = 0;
  DrawOutput = EL_OPENGL;
  Cylinder = (GLuint *)calloc(pNNano(),sizeof(GLuint));
  sprintf(Block2Draw,"ALL");
  Dr->Window(argc,argv);
  CreateMenu();
  ChooseDraw(cWhat2Draw);
  ElMenuChoise(What2Draw);
  ElMenuVisual(What2Draw);
  for(int n=0;n<pNNano();n++){
    Cylinder[n] = Dr->DefCylinder(Nano[n].Rad,Nano[n].Height);
  }
  //Dr->Picture();exit(0);
  glutPostRedisplay();
  glutMainLoop();
}
void ElPoly::DrRunTime(){  
  //DrSurface();
}
int ElPoly::DrIntorno(int p,double Blue){
  double Near = 0.;
  double PreChain = 0.;
  if(!IfIntorno||pType(p)==2)
    Near = 0.;
  else if(IfIntorno == 1){
    Near = Ch[pChain(p)].Pos[3]*InvScaleUn;
    if(PreChain != pChain(p)){
    }
  }
  if(Near > Vicinanze)
    return 1;
  return 0;
}
void ElPoly::DrBondOpenGl(double *Pos1,double *Pos2,float *Color){
  glColor4fv(Color);
  glPushMatrix();//Line
  glBegin(GL_LINES);
  glVertex3d(Pos1[0],Pos1[1],Pos1[2]);
  glVertex3d(Pos2[0],Pos2[1],Pos2[2]);
  glEnd();
  glPopMatrix();//Line
}
void ElPoly::DrBond(int p){
  if(!pNLink())return;
  double Pos1[3];
  double Pos2[3];
  for(int l=0;l<Ln[p].NLink;l++){
    int link = Ln[p].Link[l];
    for(int d=0;d<3;d++){
      Pos1[d] = pPos(p,d)*InvScaleUn;
      Pos2[d] = pPos(link,d)*InvScaleUn;
    }
    for(int d=0;d<3;d++){
      if(Pos1[d] - Pos2[d] > .5*pEdge(d)*InvScaleUn){
	Pos2[d] += pEdge(d)*InvScaleUn;
      }
      else if(Pos1[d] - Pos2[d] < -.5*pEdge(d)*InvScaleUn){
	Pos2[d] -= pEdge(d)*InvScaleUn;
      }
    }
    DrawBond(Pos1,Pos2,ColorType[pType(p)]);
  }
}
void ElPoly::DrawFuncHeader(){
#ifdef __glut_h__
  if(DrawOutput == EL_OPENGL){
    glDisable(GL_LIGHTING);
    //glEnable(GL_LIGHTING);
    glDeleteLists(Dr->Particles,1);
    Dr->Particles = glGenLists(1);
    if(pNPart() > 50000) NVisSkip = 2;
    else if (pNPart() > 100000) NVisSkip = 4;
    for(int n=0;n<pNNano();n++){
      glDeleteLists(Cylinder[n],1);
      Cylinder[n] = Dr->DefCylinder(Nano[n].Rad,Nano[n].Height);
      if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_WALL)){
	glDeleteLists(GlWall,1);
	GlWall = Dr->DefWall();
      }
    }
    glNewList(Dr->Particles,GL_COMPILE);
    glDisable(GL_LIGHTING);
    Draw_Part = &ElPoly::DrPartOpenGl;
    Draw_Nano = &ElPoly::DrNanoOpenGl;
    if(IfLine) Draw_Bond = &ElPoly::DrBondOpenGl;
    else Draw_Bond = &ElPoly::DrBondNo;
  }
  else if(DrawOutput == EL_POVRAY){
#endif // __glut_h_
    Draw_Part = &ElPoly::DrPartPovRay;
    Draw_Nano = &ElPoly::DrNanoPovRay;
    Draw_Bond = &ElPoly::DrBondPovRay;
    HeaderPovRay();
    SigErr(DrawOutFile != NULL,"DrawOutFile already allocated, can't use the file");
    char FName[60];
    sprintf(FName,"PovSnap%05d.pov",pStep());
    DrawOutFile = fopen(FName,"w");
    int ImSize[2] = {1000,1000};
    fprintf(DrawOutFile,"// POV 3.x input script : plot.pov\n// command: povray +W%d +H%d -I%s -O%s.tga +P +X +A +FT +C\n",ImSize[0],ImSize[1],FName,FName);
    fprintf(DrawOutFile,"#if (version < 3.5)\n#error \"POV3DisplayDevice has been compiled for POV-Ray 3.5 or above.\\nPlease upgrade POV-Ray.\"\n#end\n");
    fprintf(DrawOutFile,"#include \"PovHeader.inc\"\n");
    fprintf(DrawOutFile,"// System \n");
#ifdef __glut_h__
  }
#endif // __glut_h_
}
void ElPoly::DrawFuncFooter(){
#ifdef __glut_h__
  if(DrawOutput == EL_OPENGL){
    glEndList();
  }
  else if(DrawOutput == EL_POVRAY){
#endif // __glut_h_
    fprintf(DrawOutFile,"// End of POV-Ray 3.x generation\n");
    fclose(DrawOutFile);
#ifdef __glut_h__
  }
#endif // __glut_h_
}
void ElPoly::DrPartOpenGl(int p){
  int Typ = pType(p) < 6 ? pType(p) : 5;
  //ColorType[Typ][0] = pVel(p,3);
  glColor4fv(ColorType[Typ]);
  glPushMatrix();//Particle
  glBegin(GL_POINTS);
  glVertex3f(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn*ScaleFact);
  glEnd();
  glPopMatrix();//Particle
}
void ElPoly::DrCrossLinks(char *FileName){
  FILE *FRead = fopen(FileName,"r");
  if(FRead == NULL) return;
  char cLine[STRSIZE];
  int c1=0;
  int c2=0;
  glDisable(GL_LIGHTING);
  glColor4f(.5,0.0,.4,1.0);
  for(int k=0;!(fgets(cLine,STRSIZE,FRead)==NULL);k++){
    fscanf(FRead,"%d %d\n",&c1,&c2);
    glPushMatrix();//Line
    glBegin(GL_LINES);
    glVertex3d(pPos(c1,0)*InvScaleUn,
	       pPos(c1,1)*InvScaleUn,
	       pPos(c1,2)*InvScaleUn*ScaleFact);
    double Pos[3] = {pPos(c2,0),pPos(c2,1),pPos(c2,2)};
    for(int d=0;d<3;d++){
      if(pPos(c1,d) - Pos[d] > .5*pEdge(d))
	Pos[d] += pEdge(d);
      else if(pPos(c1,d) - Pos[d] < -.5*pEdge(d))
	Pos[d] -= pEdge(d);
    }
    glVertex3d(((Pos[0])*InvScaleUn),
	       ((Pos[1])*InvScaleUn),
	       ((Pos[2])*InvScaleUn*ScaleFact));
    glEnd();
    glPopMatrix();//Line
  }
}
void ElPoly::DrPartList(){
  double Pos1[3];
  double Pos2[3];
  DrawFuncHeader();
  for(int n=0;n<pNNano();n++) DrawNano(n);
  glDisable(GL_LIGHTING);
  DrCrossLinks("CrossLinks.dat");
  for(int b=0,NPep=0,cOff=0,pOff=0;b<pNBlock();cOff+=pNChain(b++)){
    if(strcmp(Block2Draw,"ALL")){
      if(strcmp(Block[b].Name,Block2Draw) ){
	pOff += pNChain(b)*pNPCh(b);
	continue;
      }
    }
    if(!strncmp(Block[b].Name,"PEP",3)){
      DrProtein(Nano[NPep].ArchFile,b);
      NPep++;
      continue;
    }
    for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
      int Chc = Ch[pChain(p)].Type;
      if(!CHAIN_IF_TYPE(Chc,NChType)) continue;
      if(NPType != BEAD_EVERY  && pType(p) != NPType) continue;
      // if(Pm[p].CId != 3) continue;
      // Dr->Numera(Pm[p].Pos,p);
      DrawPart(p);
      DrBond(p);
    }
    // for(int c=cOff;c<cOff+pNChain(b);c+=NVisSkip,pOff+=pNPCh(b)*NVisSkip){
    //   for(int p=pOff,link=0;p<MIN(pOff+pNPCh(b),pNPart());p++){
    // 	int Chc = Ch[pChain(p)].Type;
    // 	if(!CHAIN_IF_TYPE(Chc,NChType)) continue;
    // 	DrawPart(p);
    // 	DrBond(p);
    //   }
    // }
  }
  DrawFuncFooter();
}
void ElPoly::DrNanoOpenGl(int n){
  if(Dr->IfMaterial){
    glEnable(GL_LIGHTING);
    glEnable( GL_LIGHT0 );
  }
  else 
    glDisable(GL_LIGHTING);
  Point2Shape(Nano[n].Shape);
  double Red = fabs(Nano[n].Hamaker)/3.;
  double Alpha = 1.;//Nano[n].Hamaker/10.;
  GLfloat Color[4] = {Red,.1,.1,Alpha};
  glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,Color); 
  glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Color); 
  glColor4fv(Color);
  if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_NONE)) return;
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CLUSTER)){
    Alpha = .2;
    DrProtein(Nano[n].ArchFile,Nano[n].nBlock);
    return;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_JANUS)){
    // Point2Shape(Nano[n].Shape);
    // DrField(64,Nano[n].Rad,n);
    // continue;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_PORE)){
    PorePos();
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CYL)){
    // Point2Shape(Nano[n].Shape);
    // DrField(64,SQR(Nano[n].Rad),n);
    // continue;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_TIP)){
    DrField(32,SQR(Nano[n].Rad),n);
    return;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_BOUND)){
    DrField(32,SQR(Nano[n].Rad),n);
    return;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_TORUS)){
    double OldPos[5] = {Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],Nano->Rad,Nano->Height};
    //StalkPos(OldPos);
    DrField(64,SQR(Nano[n].Rad),n);
    return;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_STALK)){
    DrField(64,SQR(Nano[n].Rad),n);
    return;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_HARM)){
    DrField(64,SQR(Nano[n].Rad),n);
    return;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_WALL)){
    DrField(16,SQR(Nano[n].Rad),n);
    return;
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_UMBR)){
    glPushMatrix();//Nano
    glBegin(GL_QUADS);
    glNormal3f(0.,0.,-1.);
    glVertex3f(0.,0.,Nano[n].Pos[CNorm]*InvScaleUn);
    glVertex3f(0.,pEdge(1)*InvScaleUn,Nano[n].Pos[CNorm]*InvScaleUn);
    glVertex3f(pEdge(0)*InvScaleUn,pEdge(1)*InvScaleUn,Nano[n].Pos[CNorm]*InvScaleUn);
    glVertex3f(pEdge(0)*InvScaleUn,0.,Nano[n].Pos[CNorm]*InvScaleUn);
    glEnd();//GL_QUADS
    glPopMatrix();//Nano
    return;
  }
  Vettore Ax(Nano[n].Axis[0],Nano[n].Axis[1],Nano[n].Axis[2]);
  Vettore Zed(0.,0.,1.);
  Vettore Axis(3);
  Ax.Normalize();
  Axis.VetV(&Zed,&Ax);
  double Angle = Zed.Angle(&Ax,&Zed);
  //Quadri q(Axis.x,Angle);
  //double data[16];
  //q.RotMatrix(data,4);
  //Matrice M(q,4);
  Matrice M(Axis.x,Angle,4);
  glPushMatrix();//Nano
  glTranslated((pNanoPos(n,0)*InvScaleUn),
	       (pNanoPos(n,1)*InvScaleUn),
	       (pNanoPos(n,2)*InvScaleUn*ScaleFact));
  glMultMatrixd(M.data);
  //glRotatef(Angle*DEG_RAD,Axis.x[0],Axis.x[1],Axis.x[2]);
  //glCallList(Cylinder[n]);
  if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_SPH)){
    glutSolidSphere(Nano[n].Rad*InvScaleUn,20,20);
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_EXT)){
    glColor4f(Red,.1,.1,.5);
    glCallList(Cylinder[n]);
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_PORE)){
    glColor4f(Red,.1,.1,.5);
    glCallList(Cylinder[n]);
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_JANUS)){
    glCallList(Cylinder[n]);
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CYL)){
    glCallList(Cylinder[n]);
  }
  else{
    glCallList(Cylinder[n]);
  }
  glPopMatrix();//Nano
  return;
  Red = fabs(Nano[n].Hamaker)/3.;
  if(!VAR_IF_TYPE(Nano[n].Shape,SHAPE_WALL))return;
  glDisable(GL_LIGHTING);
  // double x = (pNanoPos(n,0) - .5*pEdge(0))*InvScaleUn;
  // double y = (pNanoPos(n,1) - .5*pEdge(1))*InvScaleUn;
  // double z = (pNanoPos(n,2))*InvScaleUn*ScaleFact;
  double x = pNanoPos(n,0)*InvScaleUn;
  double y = pNanoPos(n,1)*InvScaleUn;
  double z = pNanoPos(n,2)*InvScaleUn*ScaleFact;
  //Quadri q(Axis.x,Angle);
  //double data[16];
  //q.RotMatrix(data,4);
  //Matrice M(q,4);
  glPushMatrix();//Nano
  glTranslated(x,y,z);
  glMultMatrixd(M.data);
  glCallList(GlWall);
  glPopMatrix();//Nano
  // glPushMatrix();//Nano
  // glTranslated(x,y,z);
  // //glCallList(GlWall);
  // glPopMatrix();//Nano
}
void ElPoly::DrProtein(const char *ArchFile,int b){
  glDisable(GL_LIGHTING);
  FILE *Arch = fopen(ArchFile,"r");
  if(Arch == NULL) {printf("File %s not present \n",ArchFile);return;}
  char cLine[STRSIZE];
  double Pos1[3];
  double Pos2[3];
  double Mid[3];
  int l1 = 0;
  int l2 = 0;
  double Elong = 0.;
  double Alpha = 1.;
  for(int p=Block[b].InitIdx;fgets(cLine,STRSIZE,Arch);p++){
    DrawPart(p);
    sscanf(cLine,"%d %d %lf %lf\n",&l1,&l2,&Elong,&Alpha);
    l1 += Block[b].InitIdx;
    l2 += Block[b].InitIdx;
    if(l1 == l2) continue;
    for(int d=0;d<3;d++){
      Pos1[d] = pPos(l1,d)*InvScaleUn;
      Pos2[d] = pPos(l2,d)*InvScaleUn;
      Pos2[d] += floor((pPos(l1,d) - pPos(l2,d))*pInvEdge(d)+.5)*pEdge(d)*InvScaleUn;
      Mid[d] = .5*(Pos1[d]+Pos2[d]);
    }
    DrawBond(Pos1,Mid,ColorType[pType(l1)]);
    DrawBond(Mid,Pos2,ColorType[pType(l2)]);
  }
  fclose(Arch);
}
void ElPoly::DrColor(){
  glDisable(GL_LIGHTING);
  GLuint Cube = Dr->DefCube(NEdge);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  if(1==0){
    glPushMatrix();//RefPoint
    double RefPos[3] = {(ShiftPos[0]+.5)*pEdge(0),(ShiftPos[1]+.5)*pEdge(1),(ShiftPos[2]+.5)*pEdge(2)};
    glTranslated((RefPos[0]*InvScaleUn),
		 (RefPos[1]*InvScaleUn),
		 (RefPos[2]*InvScaleUn*ScaleFact));
    glColor4f( .7,.1,.1,1.);
    glutSolidSphere(0.1,20,20);
    glPopMatrix();//RefPos
  }
  int NDrawn = 0;
  for(int p=0,link=0;p<pNPart();p++){
    double Red = pVel(p,0)/pVelMax(0)*Saturation;
    if(Red < 0.) Red = -Red;
    double Green = pVel(p,1)/pVelMax(1);
    double Blue = pVel(p,2)/pVelMax(2);
    //if(Red < .3) continue;
    double Alpha = 1.;//Pm[p].Vel[3]/pVel(3);
    double Add = 1./(double)NEdge/3.;
    if(1==0){
      if(Green > .1){
	glColor4f(0.,Green,0.,Alpha);
	glPushMatrix();//Cube
	glTranslated(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn*ScaleFact);
	glCallList(Cube);
	glPopMatrix();//Cube
      }
      if(Blue > .1){
	glColor4f(0.,0.,Blue,Alpha);
	glPushMatrix();//Cube
	glTranslated(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,(pPos(p,2)*InvScaleUn+Add)*ScaleFact);
	glCallList(Cube);
	glPopMatrix();//Cube
      }
      if(Red > .1){
	glColor4f(Red,0.,0.,Alpha);
	glPushMatrix();//Cube
	glTranslated(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,(pPos(p,2)*InvScaleUn+Add*2.)*ScaleFact);
	glCallList(Cube);
	glPopMatrix();//Cube
      }
    }
    else{
      if(Green  + Blue + Red> .01)
	{
	glColor4f(Red,Green,Blue,Alpha);
	glPushMatrix();//Cube
	glTranslated(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn*ScaleFact);
	glCallList(Cube);
	glPopMatrix();//Cube
      }
    }
  }
  glEndList();
}
void ElPoly::DrVectors(){
  glDisable(GL_LIGHTING);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  ShiftRef(BF_CHAIN);
  Vettore Freccia(pEdge(0)*InvScaleUn,0.,0.);
  Vettore Axis(3);
  Vettore ChDir(3);
  Vettore ChPos(3);
  Vettore Ax0(1,0,0);
  Vettore Ax1(0,1,0);
  Vettore Ax2(0,0,1);
  for(int b=0,c1=0,NPep=0;b<pNBlock();c1+=Block[b++].NChain){
   if(strcmp(Block2Draw,"ALL")){
      if(strcmp(Block[b].Name,Block2Draw) ) continue;
    }
    if(!strncmp(Block[b].Name,"PEP",3)){
      //DrProtein(Nano[++NPep].ArchFile,b);
      continue;
    }
    //    if(strncmp(Block[b].Name,"TT",2))continue;
    for(int c=c1;c<c1+pNChain(b);c++){
      if(!CHAIN_IF_TYPE(c,NChType)) continue;
      for(int d=0;d<3;d++){
	ChDir.Set(Ch[c].Dir[d]*InvScaleUn,d);
	ChPos.Set((Ch[c].Pos[d]-Ch[c].Dir[d]*.5)*InvScaleUn,d);
      }
      double Red   = Ch[c].Angle/M_PI;
      double Green = 1.-Red;//Ax2.Angle(&Ax1,&ChDir)/M_PI;
      double Blue  = 0.;//sin(Ax0.Angle(&Ax2,&ChDir))-Green;
      double Alpha = 1.;
      if(Red > 0.35 && Red < .65);
      // else Blue = 0.;
      else continue;
      glColor4f(Red,Green,Blue,Alpha);
      glPushMatrix();//Arrow
      double Length1 = ChDir.Norm();
      double Length2 = Freccia.Norm();
      for(int d=0;d<3;d++){
	Axis.Set(.5*(ChDir.Val(d)/Length1+Freccia.Val(d)/Length2),d);
      }
      glTranslated(ChPos.Val(0),ChPos.Val(1),ChPos.Val(2));
      glRotatef(180.,Axis.Val(0),Axis.Val(1),Axis.Val(2));
      glScaled(Length1/Length2,1.0,1.0);
      glNormal3f(ChDir.Val(0),ChDir.Val(1),ChDir.Val(2));
      glCallList(Arrow);
      glPopMatrix();//Arrow
      int p1 = c*pNPCh();
      for(int ppc=0;ppc<pNPCh();ppc++){
	int p = ppc+p1;
	DrawPart(p);
	if(IfLine) DrBond(p);
      }
    }
  }
  glEndList();
}
void ElPoly::DrVector(Vettore v,Vettore Origin){
  glDisable(GL_LIGHTING);
  glPushMatrix();//Arrow
  Vettore Freccia(pEdge(0)*InvScaleUn,0.,0.);
  Vettore Axis(3);
  double Length1 = v.Norm();
  double Length2 = Freccia.Norm();
  for(int d=0;d<3;d++){
    Axis.Set(d,.5*(v[d]/Length1 + Freccia[d]/Length2));
  }
  glTranslated(Origin[0],Origin[1],Origin[2]);
  glRotatef(180.,Axis[0],Axis[1],Axis[2]);
  glScaled(Length1/Length2,1.0,1.0);
  glNormal3f(v[0],v[1],v[2]);
  glCallList(Arrow);
  glPopMatrix();//Arrow
}
void ElPoly::DrChains(){
  glDisable(GL_LIGHTING);
  double Values=200.;
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  Folding();
  for(int c=0;c<pNChain();c++){
    //Numera(Ch[c].Pos,c);
    int Chc = Ch[c].Type;
    if(!CHAIN_IF_TYPE(Chc,NChType)) continue;
    double Red = 0.;
    if(CHAIN_IF_TYPE(Chc,CHAIN_UP))
      Red = 1.;
    double Blue=0.;
    if(CHAIN_IF_TYPE(Chc,CHAIN_FLABBY))
      Blue =1.;
    double Green=.7;
    if(CHAIN_IF_TYPE(Chc,CHAIN_TILTED))
      Green =.5;
    glPushMatrix();
    glTranslatef((GLfloat)(Ch[c].Pos[0]*InvScaleUn),
		 (GLfloat)(Ch[c].Pos[1]*InvScaleUn),
		 (GLfloat)(Ch[c].Pos[2]*InvScaleUn*ScaleFact));
    glColor4f(Red,Green,Blue,1.);
    glCallList(Hexagon);
    glPopMatrix();
  }
  glEndList();
  return ;
  double Unity= 1./(double)Values*pEdge(0)*InvScaleUn;
  double *dDensity;
  dDensity = (double *)calloc(Values,sizeof(double));
  PairCorrelation(dDensity,Values,1,CHAIN_UP);
  double Max=0.;
  for(int v=0;v<Values;v++){
    if(Max < dDensity[v]){
      Max = dDensity[v];
    }
  }
  for(int v=0;v<Values;v++){
    dDensity[v] /= Max;
  }
  for(int v=0;v<Values-1;v++){
    glPushMatrix();
    glBegin(GL_LINES);
    glColor4f(1.,0.,1.,1.);
    glVertex3d((Values-v)*Unity,pEdge(1)*InvScaleUn,dDensity[v]);
    glVertex3d((Values-(v+1))*Unity,pEdge(1)*InvScaleUn,dDensity[v+1]);
    glEnd();
    glPopMatrix();
  }
  glEndList();
  free(dDensity);
}
//obsolete?
void ElPoly::DrPosCol(){
  glDisable(GL_LIGHTING);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  glPointSize(3);
  glLineWidth(3);
  for(int p=0,link=0;p<pNPart();p++){
    glColor4f(pVel(p,0)/pVelMax(0),pVel(p,1)/pVelMax(1),pVel(p,2)/pVelMax(2),pVel(p,3)/pVelMax(3));
    if(pType(p) == 2){
      glPushMatrix();//Particle
      glTranslated((pPos(p,0)*InvScaleUn),
		   (pPos(p,1)*InvScaleUn),
		   (pPos(p,2)*InvScaleUn*ScaleFact));
      glColor4f( 1.,.1,.1,1.);
      if(Nano->Shape == SHAPE_CYL)
	glCallList(Cylinder[0]);
      else 
	glutSolidSphere(Nano->Rad*InvScaleUn,20,20);
      glPopMatrix();//Particle
    }
    else{
      glPushMatrix();//Particle
      glTranslated(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn*ScaleFact);
      glCallList(Cube);
      glPopMatrix();//Particle
    }
  }
  glEndList();
}
void ElPoly::DrCrossLinks(){
  glDisable(GL_LIGHTING);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  glPointSize(3);
  glLineWidth(3);
  for(int p=0,link=0;p<pNPart();p++){
    int Chc = Ch[pChain(p)].Type;
    if(!CHAIN_IF_TYPE(Chc,NChType) && pType(p)!=2) continue;
    int Typ = pType(p) < 6 ? pType(p) : 5;
    glColor4fv(ColorType[Typ]);
    if(IfLine){
     glColor4f(0.0,1.0,0.,1.0);
      for(int l=0;l<Ln[p].NLink;l++){
	link = Ln[p].Link[l];
	double Shift[3] = {0.,0.,0.};
	if(link != p+1){ 
	  glColor4f(1.0,.0,0.,1.0);
	  for(int d=0;d<3;d++)
	    Shift[d] = -floor((pPos(p,d)-pPos(link,d))/pEdge(d)+.5)*pEdge(d);
	}
	else continue;
	glPushMatrix();//Line
	glBegin(GL_LINES);
	glVertex3d((pPos(p,0)*InvScaleUn),
		   (pPos(p,1)*InvScaleUn),
		   (pPos(p,2)*InvScaleUn*ScaleFact));
	glVertex3d(((pPos(link,0)-Shift[0])*InvScaleUn),
		   ((pPos(link,1)-Shift[1])*InvScaleUn),
		   ((pPos(link,2)-Shift[2])*InvScaleUn*ScaleFact));
	glEnd();
	glPopMatrix();//Line
      }
    }
    glPushMatrix();//Particle
    glTranslated(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn);
    //glutSolidSphere(.007,20,20);
    glCallList(Point);
    glPopMatrix();//Particle
  }
  glEndList();
}
void ElPoly::Menu(){
  //   submenu = glutCreateMenu(processEvent);
  //   glutAddMenuEntry("SottoPrima",3);

  //   menu = glutCreateMenu(processEvent);
  //   glutAddMenuEntry("Move",1);
  //   glutAddMenuEntry("Stop",2);
  //   glutAddSubMenu("Altro",submenu);
  //   glutAttachMenu(GLUT_RIGHT_BUTTON);
}
void ElPoly::ESlide1(void){
  //  if(Dr->Step > 100) glutIdleFunc(NULL);
  Dr->zw += .1;
  Dr->IfInfo = 0;
  glutPostRedisplay();
  Dr->Step++;
  Dr->Picture();
  Dr->IfInfo = 1;;
  //Dr->Picture();
  //glutTimerFunc(30, Timer, 0);
}
void ElPoly::ESlide(void){
  Dr->IfInfo = 0;
  glutPostRedisplay();
  if(quando < NFileTot)
    quando++;
    //quando += 10;
  else{
    glutIdleFunc(NULL);
  }
  OpenFile(quando);
  RenderPart();
  Dr->Step++;
  Dr->Picture();
  glutPostRedisplay();
  Dr->IfInfo = 1;
  //glutTimerFunc(30, Timer, 0);
}
#endif //__glut_h__
