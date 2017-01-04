/*
* ElPolyDrawSurf - Drawing of the surfaces defined by the chain position
 opyright (C) 2008-2010 by Giovanni Marelli <sabeiro@virgilio.it>


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
#include "ElPoly.h"
#include <vector> 

#ifdef __glut_h__
extern Draw *Dr;
void ElPoly::DrPolygon(){
  glEnable(GL_LIGHTING);
  glEnable( GL_LIGHT0 );
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  for(int n=0;n<pNNano();n++) DrawNano(n);
  double Cm[3];
  for(int p=0;p<pNPart();p+=Ln[p].NLink+1){
    //if(Pm[p].Pos[2] > .6*pEdge(2) || Pm[p].Pos[2] < .3*pEdge(2)) continue;
    double Red = 0.0;
    double Green = .4 + .3*Mat->Casuale();
    double Blue  = 0.0;
    GLfloat DrMatColor [] = {Red,Green,Blue,1.};
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,DrMatColor); 
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,DrMatColor); 
    glColor4f(Red,Green,Blue,1.);
    glPushMatrix();
    glBegin(GL_POLYGON);
    //glBegin(GL_TRIANGLE);
    //glBegin(GL_TRIANGLE_STRIP);
    glNormal3d(pVel(p,0)*InvScaleUn,pVel(p,1)*InvScaleUn,pVel(p,2)*InvScaleUn);
    glVertex3d(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn*ScaleFact);
    for(int d=0;d<3;d++) Cm[d] = pPos(p,d);
    for(int l=0;l<Ln[p].NLink;l++){
      int l1 = Ln[p].Link[l];
      glNormal3d(pVel(l1,0),pVel(l1,1),pVel(l1,2));
      glVertex3d(pPos(l1,0)*InvScaleUn,pPos(l1,1)*InvScaleUn,pPos(l1,2)*InvScaleUn*ScaleFact);
      for(int d=0;d<3;d++){
	Cm[d] += pPos(p,d);
      }
    }
    glEnd();
    glPopMatrix();
    if(IfLine){
      for(int d=0;d<3;d++){
	Cm[d] /= (double)(Ln[p].NLink+1);
      }
      glColor4f(.1,.3,1.,1.);
      glBegin(GL_LINES);
      glNormal3d(0.,0.,1.);
      glVertex3d(Cm[0]*InvScaleUn,Cm[1]*InvScaleUn,Cm[2]*InvScaleUn*ScaleFact);
      glVertex3d((Cm[0]-pVel(p,0))*InvScaleUn,(Cm[1]-pVel(p,1))*InvScaleUn,(Cm[2]-pVel(p,2))*InvScaleUn*ScaleFact);
      glEnd();
    }
  }
  glEndList();
}
void ElPoly::DrNormalPoint(int p,int NEdge){
  int link[4];
  link[0] = p + NEdge;
  if(link[0] >= pNPart()) link[0] -= pNPart();
  link[1] = p + 1;
  if(link[1] >= pNPart()) link[1] -= pNPart();
  link[2] = p - NEdge;
  if(link[2] < 0 ) link[2] += pNPart();
  link[3] = p - 1;
  if(link[3] < 0 ) link[3] += pNPart();
  double Normals[5][3];
  for(int l=0;l<4;l++){
    for(int d=0;d<3;d++){
      int d1 = (d+1)%3;
      int d2 = (d+2)%3;
      Normals[l][d] = (pPos(p,d1)-pPos(link[l],d2))*(pPos(p,d2)-pPos(link[l],d));
    }
  }
  double Norm;
  for(int d=0;d<3;d++){
    Normals[4][d] = Normals[0][d] + Normals[1][d] + Normals[2][d] + Normals[3][d];
    Norm += SQR(Normals[4][d]);
  }
  Norm = 1./sqrt(Norm);
  for(int d=0;d<3;d++){
    Normals[4][d] /= Norm;
  }
  glNormal3d(Normals[4][0],Normals[4][1],Normals[4][2]);
}
void ElPoly::DrSquareMesh(){
  if(Dr->IfMaterial){
    glEnable(GL_LIGHTING);
    glEnable( GL_LIGHT0 );
  }
  else 
    glDisable(GL_LIGHTING);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  int NEdge = (int)sqrt(pNPart());
  if(NEdge*NEdge != pNPart()) NEdge += 1;
  GLfloat Color[4];
  for(int p=0;p<pNPart();p++){
    double Depth = pPos(p,2)*pInvEdge(CNorm)*Saturation+ExtParam;
    Dr->DepthMap(Depth,Color);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,Color); 
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Color); 
    glColor4fv(Color);
    int l1 = p+1;
    if((l1%NEdge) == 0)continue;
    int l2 = p + NEdge + 1;
    if(l2 >= pNPart()-1)continue;
    int l3 = p + NEdge;
    if(l3 >= pNPart()-1)continue;
    //glPushMatrix();
    glBegin(GL_POLYGON);
    DrNormalPoint(p,NEdge);
    glVertex3d(Pm[p].Pos[0]*InvScaleUn,Pm[p].Pos[1]*InvScaleUn,Pm[p].Pos[2]*InvScaleUn*ScaleFact);
    DrNormalPoint(l1,NEdge);
    glVertex3d(Pm[l1].Pos[0]*InvScaleUn,Pm[l1].Pos[1]*InvScaleUn,Pm[l1].Pos[2]*InvScaleUn*ScaleFact);
    DrNormalPoint(l2,NEdge);
    glVertex3d(Pm[l2].Pos[0]*InvScaleUn,Pm[l2].Pos[1]*InvScaleUn,Pm[l2].Pos[2]*InvScaleUn*ScaleFact);
    DrNormalPoint(l3,NEdge);
    glVertex3d(Pm[l3].Pos[0]*InvScaleUn,Pm[l3].Pos[1]*InvScaleUn,Pm[l3].Pos[2]*InvScaleUn*ScaleFact);
    glEnd();
    //glPopMatrix();
  }
  glEndList();
  glDisable(GL_LIGHTING);
}
#include "ElPolyDrawSurf.h"
void ElPoly::DrField(int NGrid,double IsoLevel,int nNano){
  glEnable(GL_LIGHTING);
  glEnable( GL_LIGHT0 );
  glMaterialfv(GL_BACK,  GL_AMBIENT,   DrAmbientRed); 
  glMaterialfv(GL_BACK,  GL_DIFFUSE,   DrDiffuseRed); 
  glMaterialfv(GL_FRONT, GL_AMBIENT,   DrAmbientBlue); 
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   DrDiffuseBlue); 
  glMaterialfv(GL_FRONT, GL_SPECULAR,  DrSpecularWhite); 
  glMaterialf( GL_FRONT, GL_SHININESS, 25.0); 
  double CubeDist[8];
  double EdgeVertex[12][3];
  double EdgeNormal[12][3];
  double InvNGrid = 1./(double)NGrid;
  double Pos[3];
  double Pos1[3];
  glBegin(GL_TRIANGLES);
  //glBegin(GL_LINES);
  //glBegin(GL_POINTS);
  for(int gx=0;gx<NGrid;gx++){
    Pos[0] = gx*InvNGrid*pEdge(0);
    for(int gy=0;gy<NGrid;gy++){
      Pos[1] = gy*InvNGrid*pEdge(1);
      for(int gz=0;gz<NGrid;gz++){
	Pos[2] = gz*InvNGrid*pEdge(2);
	for(int v=0;v<8;v++){
	  for(int d=0;d<3;d++){
	    Pos1[d] = Pos[d] + VertCube[v][d]*InvNGrid*pEdge(d);
	  }
	  CubeDist[v] = NanoDist2(Pos1,nNano);
	}
	int Flag = 0;
	for(int v=0;v<8;v++){
	  if(CubeDist[v] <= IsoLevel)
	    Flag |= 1<<v;
	}
	int CubeShape = CubeTop[Flag];
	if(CubeShape==0) continue;
	for(int e=0;e<12;e++){
	  if(CubeShape & (1<<e)){
	    double Delta = CubeDist[EdgeConn[e][1]] - CubeDist[EdgeConn[e][0]];
	    double OffSet = (IsoLevel-CubeDist[EdgeConn[e][0]])/Delta;
	    if(Delta == 0.0){
	      OffSet = .5;
	    }
	    EdgeNormal[e][0] = NanoDist2(Pos[0]-0.01,Pos[1],Pos[2],nNano) 
	      - NanoDist2(Pos[0]+0.01,Pos[1],Pos[2],nNano);
	    EdgeNormal[e][1] = NanoDist2(Pos[0],Pos[1]-0.01,Pos[2],nNano) 
	      - NanoDist2(Pos[0],Pos[1]+0.01,Pos[2],nNano);
	    EdgeNormal[e][2] = NanoDist2(Pos[0],Pos[1],Pos[2]-0.01,nNano) 
	      - NanoDist2(Pos[0],Pos[1],Pos[2]+0.01,nNano);
	    double Norm = sqrt(SQR(EdgeNormal[e][0])+SQR(EdgeNormal[e][1])+SQR(EdgeNormal[e][2]));
	    for(int d=0;d<3;d++){
	      EdgeVertex[e][d] = Pos[d] + (VertCube[EdgeConn[e][0]][d]+OffSet*EdgeDir[e][d])*InvNGrid*pEdge(d);
	      EdgeVertex[e][d] *= InvScaleUn;
	      EdgeNormal[e][d] *= 1./Norm;
	    }
	  }
	}
	for(int t=0;t<5;t++){
	  if(TrConnTable[Flag][3*t] < 0.) break;
	  for(int d=0;d<3;d++){
	    int v = TrConnTable[Flag][3*t+d];
	    //glColor4f(1.0,0.1,0.1,1.0);
	    glNormal3d(EdgeNormal[v][0],EdgeNormal[v][1],EdgeNormal[v][2]);
	    glVertex3d(EdgeVertex[v][0],EdgeVertex[v][1],EdgeVertex[v][2]);
	  }
	}
      }
    }
  }
  glEnd();
}
void ElPoly::DrQuad1(){
  glEnable(GL_LIGHTING);
  glEnable( GL_LIGHT0 );
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  //glEnable(GL_POLYGON_STIPPLE);
  double RedMax = fabs(log(fabs(pVel(0,0))*Saturation));
  for(int p=0;p<pNPart();p++){
    if(RedMax <  fabs(log(fabs(pVel(0,0))*Saturation)) )
      RedMax  =  fabs(log(fabs(pVel(0,0))*Saturation));
  }
  int NDraw = 0;
  for(int p=0;p<pNPart();p++){
    if(pPos(p,2) < .7*pEdge(2)) continue;
    //double Red = fabs(log(fabs(pPos(p,0))*Sat))/RedMax;
    double Red = fabs(pVel(p,0)*Saturation/pVelMax(0));
    double Green = pVel(p,1)/pVelMax(1)*Saturation;
    double Blue  = pVel(p,2)/pVelMax(2)*Saturation;
    GLfloat DrMatColor [] = {Red,Green,Blue,1.};
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,DrMatColor); 
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,DrMatColor); 
    //if(Red<0.1) continue;
    glColor4f(Red,Green,Blue,1.);
    glPushMatrix();
    glTranslated(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn*ScaleFact);
    glCallList(Quad);
    glPopMatrix();
    NDraw++;
  }
  printf("%d\n",NDraw);
  glEndList();
}
void ElPoly::DrQuad(){
  glEnable(GL_LIGHTING);
  glEnable( GL_LIGHT0 );
  //CompileList();
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  int Incr = 1;//4;
  double Offset = 1./(double)NEdge;
  for(int p=0;p<pNPart();p+=Incr){
    double Red = pVel(p,0)*Saturation;// + pPos(p,2)*Sat;
    double Green = pVel(p,1)*Saturation;
    double Blue = pVel(p,2)*Saturation;
    GLfloat DrMatColor [] = {Red,Green,Blue,1.};
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,DrMatColor); 
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,DrMatColor); 
    //if(Red < 0.0001  && Green < 0.2 && Blue < 0.2) continue;
    glColor4f(Red,Green,Blue,1.0);
    glPushMatrix();
    glTranslated(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,0.);
    glCallList(Quad);
    glPopMatrix();
    // glPushMatrix();
    // glTranslated(-pPos(p,0)*InvScaleUn-Offset,pPos(p,1)*InvScaleUn,0.);
    // glCallList(Quad);
    // glPopMatrix();
  }
  glEndList();return;
}
void ElPoly::DrPotential(){
  glEnable(GL_LIGHTING);
  glEnable( GL_LIGHT0 );
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  int Nx = 0;
  int Ny = 0;
  for(int p=0;p<pNPart();p++){
    if(pPos(p,CLat2) < pPos(p+1,CLat2))
      Ny++;
    else break;
  }
  Ny++;
  Nx = (int)(pNPart()/(double)Ny);
  const int NIso = 20;
  double Min = pPos(0,CNorm);
  double Max = pPos(0,CNorm);
  for(int p=0;p<pNPart();p++){
    if(Min > pPos(p,CNorm)) Min = pPos(p,CNorm);
    if(Max < pPos(p,CNorm)) Max = pPos(p,CNorm);
  }
  Min = 1./MIN(fabs(Min),1.);
  Max = 1./MIN(1.,Max);
  for(int p=0;p<pNPart();p++){
    int l1 = p+1;
    if( !(l1%Ny) ) continue;
    int l2 = p + Ny;
    if(l2 >= pNPart() - 1) continue;
    int l3 = p + Ny + 1;
    double Red = 1.;
    double Blue = 1.;
    double Green = 1.;
    double Alpha = 1.;
    if(pPos(p,CNorm) < 0.){
      Red = 0.;
      Green = fabs(pPos(p,CNorm)*Min)*Saturation;
      //Green = fabs(pPos(p,2)*1.)*Sat;
      Blue = 0.;
      Alpha = MAX(Green-.2,0.);
    }
    else if(pPos(p,CNorm) > 0.){
      Red = 0.;
      Green = 0.;
      Blue = (pPos(p,CNorm)*pInvEdge(CNorm))*Saturation;
      //Blue = pPos(p,2)*1.*Sat;
      Alpha = MAX(Blue,0.);
    }
    else
      Alpha = 0.;
    GLfloat DrMatColor [] = {Red,Green,Blue,Alpha};
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,DrMatColor); 
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,DrMatColor); 
    glColor4d(Red,Green,Blue,Alpha);
    // glPushMatrix();//Particle
    glBegin(GL_POINTS);
    glVertex3f(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn*ScaleFact);
    glEnd();
    // glPopMatrix();//Particle
    glBegin(GL_POLYGON);
    glNormal3f(0.,0.,1.);
    glVertex3d(pPos(p,0)*InvScaleUn,pPos(p,1)*InvScaleUn,pPos(p,2)*InvScaleUn);
    glVertex3d(pPos(l1,0)*InvScaleUn,pPos(l1,1)*InvScaleUn,pPos(l1,2)*InvScaleUn);
    glVertex3d(pPos(l3,0)*InvScaleUn,pPos(l3,1)*InvScaleUn,pPos(l3,2)*InvScaleUn);
    glVertex3d(pPos(l2,0)*InvScaleUn,pPos(l2,1)*InvScaleUn,pPos(l2,2)*InvScaleUn);
    glEnd();
  }
  glEndList();
  // int Values = NEdge/4;
  // double *Border = (double *)calloc(2*Values,sizeof(double));
  // double *Norm = (double *)calloc(2*Values,sizeof(double));
  // for(int p=0;p<pNPart();p++){
  //   int v = (int)(pPos(p,0)/pEdge(0)*Values);
  //   if( v >= Values) continue;
  //   if(pPos(p,1) > .5*pEdge(1)){
  //     Border[2*v] += pPos(p,1)*pPos(p,2);
  //     Norm[2*v] += pPos(p,2);
  //   }
  //   else{
  //     Border[2*v+1] += pPos(p,1)*pPos(p,2);
  //     Norm[2*v+1] += pPos(p,2];
  //   }
  // }
  // glColor4f(.4,0.0,.6,1.0);
  // glBegin(GL_LINE_STRIP);
  // glNormal3f(0.,0.,1.);
  // for(int v=0,vv=0;v<2*Values-1;v+=2,vv++){
  //   if(!(Norm[v] > 0. ) )continue;
  //   glVertex3f(vv/(double)Values*pEdge(0)*InvScaleUn,Border[v]*InvScaleUn/Norm[v],0.);
  // }
  // glEnd();
  // glBegin(GL_LINE_STRIP);
  // for(int v=1,vv=0;v<2*Values-1;v+=2,vv++){
  //   if(!(Norm[v] > 0. ) )continue;
  //   glVertex3f(vv/(double)Values*pEdge(0)*InvScaleUn,Border[v]*InvScaleUn/Norm[v],0.);
  // }
  // glEnd();
  // glBegin(GL_LINE_STRIP);
  // for(int v=0,vv=0;v<2*Values-1;v+=2,vv++){
  //   if(!(Norm[v] > 0. ) )continue;
  //   glVertex3f(-vv/(double)Values*pEdge(0)*InvScaleUn,Border[v]*InvScaleUn/Norm[v],0.);
  // }
  // glEnd();
  // glBegin(GL_LINE_STRIP);
  // for(int v=1,vv=0;v<2*Values-1;v+=2,vv++){
  //   if(!(Norm[v] > 0. ) )continue;
  //   glVertex3f(-vv/(double)Values*pEdge(0)*InvScaleUn,Border[v]*InvScaleUn/Norm[v],0.);
  // }
  // glEnd();
  // glEndList();
}
void ElPoly::DrDensity(){
  if(Dr->IfMaterial){
    glEnable(GL_LIGHTING);
    glEnable( GL_LIGHT0 );
  }
  else 
    glDisable(GL_LIGHTING);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  double Offset = 1./(double)NEdge;
  for(int p=0,c=0,t=0;p<pNPart();p++){
    if( Ln[p].NLink != 3 ) continue;
    if(c == pChain(p)) continue;
    if(NPType != BEAD_EVERY  && pType(p) != NPType) continue;
    c = pChain(p);
    double Red = pVel(p,0)*Saturation/pVelMax(0);// + pPos(p,2)*Sat;
    double Green = pVel(p,1)*Saturation/pVelMax(1);
    double Blue = pVel(p,2)*Saturation/pVelMax(2);
    double Alpha = 1.0;//7.+(Green+Blue);
    if(Blue + Red + Green < .08){Blue = 1.; Red = 1.;Green = 1.;Alpha=0.;}
    GLfloat DrMatColor [] = {Red,Green,Blue,Alpha};
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,DrMatColor);
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,DrMatColor); 
    int l1 = Ln[p].Link[0];
    //if(l1 > Gen->NPart) continue;
    int l2 = Ln[p].Link[2];
    //if(l2 > Gen->NPart) continue;
    int l3 = Ln[p].Link[1];
    //if(Green > 0.) continue;
    //     Green = 0.; Blue = 0.;
    glColor4f(Red,Green,Blue,Alpha);
    // if(pPos(p,2)>12.) pPos(p,2) = 12.;
    // if(pPos(l1,2)>12.) pPos(l1,2) = 12.;
    // if(pPos(l2,2)>12.) pPos(l2,2) = 12.;
    // if(pPos(l3,2)>12.) pPos(l3,2) = 12.;
    Vettore v1(pPos(p,0),pPos(p,1),pPos(p,2)*ScaleFact);
    Vettore v2(pPos(l1,0),pPos(l1,1),pPos(l1,2)*ScaleFact);
    Vettore v3(pPos(l2,0),pPos(l2,1),pPos(l2,2)*ScaleFact);
    Vettore v4(pPos(l3,0),pPos(l3,1),pPos(l3,2)*ScaleFact);
    v1.Mult(InvScaleUn);
    v2.Mult(InvScaleUn);
    v3.Mult(InvScaleUn);
    v4.Mult(InvScaleUn);
    Vettore u(3);
    Vettore v(3);
    Vettore vN(3);
    vN.NormalSurf(&v1,&v2,&v3);
    glPushMatrix();//Particle
    glBegin(GL_QUADS);
    glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
    u = v2-v1;v = v3-v1;vN = u ^ v;//vN.Normalize();
    glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
    glVertex3d(v4.x[0],v4.x[1],v4.x[2]);
    u = v1-v2;v = v3-v2;vN = u ^ v;//vN.Normalize();
    glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
    glVertex3d(v3.x[0],v3.x[1],v3.x[2]);
    u = v2-v3;v = v4-v3;vN = u ^ v;//vN.Normalize();
    glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
    glVertex3d(v2.x[0],v2.x[1],v2.x[2]);
    u = v1-v4;v = v3-v4;vN = u ^ v;//vN.Normalize();
    glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
    glVertex3d(v1.x[0],v1.x[1],v1.x[2]);
    glEnd();
    glPopMatrix();//Particle
    //continue;
    glPushMatrix();//Particle
    glBegin(GL_QUADS);
    glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
    glVertex3d(-v1.x[0],v1.x[1],v1.x[2]);
    glVertex3d(-v2.x[0],v2.x[1],v2.x[2]);
    glVertex3d(-v3.x[0],v3.x[1],v3.x[2]);
    glVertex3d(-v4.x[0],v4.x[1],v4.x[2]);
    glEnd();
    glPopMatrix();//Particle
  }
  glEndList();
}
void ElPoly::DrInterpSurface(){
  glEnable(GL_LIGHTING);
  glEnable( GL_LIGHT0 );
  int NPoint = 32;
  GLfloat Deltax = 1./(GLfloat) NPoint;
  PART *PSample = (PART *)calloc(NPoint*NPoint,sizeof(PART));
  int NSample = (int)sqrt(pNPart()); 
  glColor4f(1.0,.0,1.,1.);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_VERTEX_ARRAY);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glFrontFace(GL_CCW);   //Tell OGL which orientation shall be the front face
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  //InterSpline3(Pm,PSample,NSample*NSample,NPoint*NPoint);
  glVertexPointer(3,GL_DOUBLE,sizeof(PART),&PSample[0].Pos[0]);
  glNormalPointer(GL_DOUBLE,sizeof(PART),&PSample[0].Vel[0]);
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
  GLuint *Indices = (GLuint *)calloc(NIndex,sizeof(GLuint));
  for(int i=0;i<NIndex;i++){
    Indices[i] = VecIndices[i];
  }
  glDrawElements(GL_TRIANGLES,NIndex,GL_UNSIGNED_INT,Indices);
  VecIndices.clear();
  free(Indices);
  glEndList();
  return;
}
void ElPoly::DrDoubleTria(Vettore *v00,Vettore *v01,Vettore *v11,Vettore *v10,Vettore *vN){
  glPushMatrix();//Particle
  glBegin(GL_TRIANGLES);
  // vN->NormalSurf(v00,v01,v11);
  // vN->Normalize();
  // glNormal3f(vN->x[0],vN->x[1],vN->x[2]);
  // glVertex3d(v00->x[0],v00->x[1],v00->x[2]);
  // vN->NormalSurf(v01,v00,v11);
  // vN->Normalize();
  // glNormal3f(vN->x[0],vN->x[1],vN->x[2]);
  // glVertex3d(v01->x[0],v01->x[1],v01->x[2]);
  // vN->NormalSurf(v11,v01,v10);
  // vN->Normalize();
  // glNormal3f(vN->x[0],vN->x[1],vN->x[2]);
  // glVertex3d(v11->x[0],v11->x[1],v11->x[2]);
  // glEnd();
  vN->NormalSurf(v00,v01,v11);
  vN->Normalize();
  glNormal3f(vN->x[0],vN->x[1],vN->x[2]);
  glVertex3d(v00->x[0],v00->x[1],v00->x[2]);
  glVertex3d(v01->x[0],v01->x[1],v01->x[2]);
  glVertex3d(v11->x[0],v11->x[1],v11->x[2]);
  glEnd();
  glPopMatrix();//Particle
  vN->NormalSurf(v00,v11,v10);
  vN->Normalize();
  glPushMatrix();//Particle
  glBegin(GL_TRIANGLES);
  glNormal3f(vN->x[0],vN->x[1],vN->x[2]);
  glVertex3d(v00->x[0],v00->x[1],v00->x[2]);
  glVertex3d(v11->x[0],v11->x[1],v11->x[2]);
  glVertex3d(v10->x[0],v10->x[1],v10->x[2]);
  glEnd();
  glPopMatrix();//Particle
}
void ElPoly::DrTria(Vettore *v00,Vettore *v01,Vettore *v11,Vettore *vN){
  vN->NormalSurf(v00,v01,v11);
  vN->Normalize();
  glPushMatrix();//Particle
  glBegin(GL_TRIANGLES);
  glNormal3f(vN->x[0],vN->x[1],vN->x[2]);
  glVertex3d(v00->x[0],v00->x[1],v00->x[2]);
  glVertex3d(v01->x[0],v01->x[1],v01->x[2]);
  glVertex3d(v11->x[0],v11->x[1],v11->x[2]);
  glEnd();
  glPopMatrix();//Particle
}
void ElPoly::DrTriaContour(Vettore *v00,Vettore *v01,Vettore *v11){
  glPushMatrix();//Particle
  glBegin(GL_LINES);
  glVertex3d(v00->x[0],v00->x[1],v00->x[2]);
  glVertex3d(v01->x[0],v01->x[1],v01->x[2]);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(v01->x[0],v01->x[1],v01->x[2]);
  glVertex3d(v11->x[0],v11->x[1],v11->x[2]);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(v11->x[0],v11->x[1],v11->x[2]);
  glVertex3d(v00->x[0],v00->x[1],v00->x[2]);
  glEnd();
  glPopMatrix();//Particle
}
void ElPoly::DrSmooth(double *Plot,int NSample,double Min,double Max){
  if(Dr->IfMaterial){
    glEnable(GL_LIGHTING);
    glEnable( GL_LIGHT0 );
  }
  else 
    glDisable(GL_LIGHTING);
  double FactC1 = 1./(double)NSample*InvScaleUn*pEdge(CLat1);
  double FactC2 = 1./(double)NSample*InvScaleUn*pEdge(CLat2);
  double FactCN = InvScaleUn*ScaleFact;
  GLfloat Color[4];
  Vettore v00(3);
  Vettore v01(3);
  Vettore v11(3);
  Vettore v10(3);
  Vettore vN(3);
  Min *= InvScaleUn;
  Max *= InvScaleUn;
  double Delta = 1./(Max-Min);
  for(int s=0;s<NSample-1;s++){
    for(int ss=0;ss<NSample-1;ss++){
      v00.x[CLat1]=(s+.5)*FactC1;
      v00.x[CLat2]=(ss+.5)*FactC2;
      v00.x[CNorm]=Plot[s*NSample+ss]*FactCN;
      v01.x[CLat1]=(s+1.5)*FactC1;
      v01.x[CLat2]=(ss+.5)*FactC2;
      v01.x[CNorm]=Plot[(s+1)*NSample+ss]*FactCN;
      v11.x[CLat1]=(s+1.5)*FactC1;
      v11.x[CLat2]=(ss+1.5)*FactC2;
      v11.x[CNorm]=Plot[(s+1)*NSample+(ss+1)]*FactCN;
      v10.x[CLat1]=(s+.5)*FactC1;
      v10.x[CLat2]=(ss+1.5)*FactC2;
      v10.x[CNorm]=Plot[s*NSample+ss+1]*FactCN;
      //glColor4f(.4,.4,.6,1.);
      //DrDoubleTria(&v00,&v01,&v11,&v10,&vN);
      glBegin(GL_POLYGON);
      vN.NormalSurf(&v00,&v01,&v11);
      vN.Normalize();
      glNormal3f(vN.x[0],vN.x[1],vN.x[2]);
      double Depth = (v00.x[CNorm]-Min)*Delta*Saturation+ExtParam;
      Dr->DepthMap(Depth,Color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,Color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Color);
      glColor4fv(Color);
      glVertex3d(v00.x[0],v00.x[1],v00.x[2]);
      Depth = (v01.x[CNorm]-Min)*Delta*Saturation+ExtParam;
      Dr->DepthMap(Depth,Color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,Color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Color);
      glColor4fv(Color);
      glVertex3d(v01.x[0],v01.x[1],v01.x[2]);
      Depth = (v11.x[CNorm]-Min)*Delta*Saturation+ExtParam;
      Dr->DepthMap(Depth,Color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,Color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Color);
      glColor4fv(Color);
      glVertex3d(v11.x[0],v11.x[1],v11.x[2]);
      Depth = (v10.x[CNorm]-Min)*Delta*Saturation+ExtParam;
      Dr->DepthMap(Depth,Color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,Color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Color);
      glColor4fv(Color);
      glVertex3d(v10.x[0],v10.x[1],v10.x[2]);
      glEnd();
    }
  }
  glDisable(GL_LIGHTING);
}
void ElPoly::DrSample(int NSample){
  NSample *= 4;
  MOMENTI m1 = SampleSurfaceMem(NSample);
  printf("plane min=%lf av=%lf max=%lf\n",m1.Min,m1.Uno,m1.Max);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  DrSmooth(PlotMem,NSample,.9*m1.Min,.9*m1.Max);
  glEndList();
}
void ElPoly::DrStalk(){
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  glDisable(GL_NORMALIZE);
  double Offset = 1./(double)NEdge;
  Vettore v1(3);
  Vettore v2(3);
  Vettore v3(3);
  Vettore v4(3);
  Vettore vN(3);
  glFrontFace(GL_CCW);
  //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  //glShadeModel(GL_SMOOTH);
  for(int p=0,c=0,t=0;p<pNPart();p++){
    if( Ln[p].NLink != 3 ) continue;
    if(c == pChain(p)) continue;
    if( NPType != BEAD_EVERY  && pType(p) != NPType && pType(p) != BEAD_NANO) continue;
    c = pChain(p);
    double Red = pVel(p,0)*Saturation;// + pPos(p,2)*Sat;
    double Green = pVel(p,1)*Saturation;
    double Blue = pVel(p,2)*Saturation;
    double Alpha = 1.0;//7.+(Green+Blue);
    int l1 = Ln[p].Link[0];
    //if(l1 > Gen->NPart) continue;
    int l2 = Ln[p].Link[2];
    //if(l2 > Gen->NPart) continue;
    int l3 = Ln[p].Link[1];
    //     Green = 0.; Blue = 0.;
    glColor4f(Red,Green,Blue,Alpha);
    v1.x[0]=pPos(p,0)*InvScaleUn;
    v1.x[1]=pPos(p,1)*InvScaleUn;
    v1.x[2]=pPos(p,2)*InvScaleUn*ScaleFact;
    v2.x[0]=pPos(l1,0)*InvScaleUn;
    v2.x[1]=pPos(l1,1)*InvScaleUn;
    v2.x[2]=pPos(l1,2)*InvScaleUn*ScaleFact;
    v3.x[0]=pPos(l3,0)*InvScaleUn;
    v3.x[1]=pPos(l3,1)*InvScaleUn;
    v3.x[2]=pPos(l3,2)*InvScaleUn*ScaleFact;
    v4.x[0]=pPos(l2,0)*InvScaleUn;
    v4.x[1]=pPos(l2,1)*InvScaleUn;
    v4.x[2]=pPos(l2,2)*InvScaleUn*ScaleFact;
    vN.NormalSurf(&v1,&v2,&v3);
    vN.Rescale(.2);
    glPushMatrix();//Particle
    glBegin(GL_QUADS);
    glNormal3f(vN.x[0],vN.x[1],-vN.x[2]);
    glVertex3d(v4.x[0],v4.x[1],v4.x[2]);
    glVertex3d(v3.x[0],v3.x[1],v3.x[2]);
    glVertex3d(v2.x[0],v2.x[1],v2.x[2]);
    glVertex3d(v1.x[0],v1.x[1],v1.x[2]);
    glEnd();
    glPopMatrix();//Particle
  }
  glEndList();
}
void ElPoly::DrCreateStalk(){
  int NSample = 24;
  double Threshold = 0.;
  int NLevel = 4;
  double **Plot = (double **) calloc(NLevel,sizeof(double));
  for(int l=0;l<NLevel;l++){
    Plot[l] = (double *)calloc(NSample*NSample,sizeof(double));
  }
  Stalk(NSample,NLevel,Plot,Threshold);
  //Drawing
  glColorMaterial(GL_FRONT_AND_BACK,GL_EMISSION);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  double FactC1 = 1./(double)NSample*InvScaleUn*pEdge(CLat1);
  double FactC2 = 1./(double)NSample*InvScaleUn*pEdge(CLat2);
  double FactCN = InvScaleUn*ScaleFact;
  Vettore v00(3);
  Vettore v01(3);
  Vettore v11(3);
  Vettore v10(3);
  Vettore vN(3);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  for(int s=0;s<NSample-1;s++){
    for(int ss=0;ss<NSample-1;ss++){
      for(int l=0;l<NLevel;l++){
	if(l==0 || l == 3) continue;
	glColor4f((double)l*.25,.5,.8,1.);
	v00.x[CLat1]=(s+.5)*FactC1;
	v00.x[CLat2]=(ss+.5)*FactC2;
	v00.x[CNorm]=Plot[l][s*NSample+ss]*FactCN;
	v01.x[CLat1]=(s+1.5)*FactC1;
	v01.x[CLat2]=(ss+.5)*FactC2;
	v01.x[CNorm]=Plot[l][(s+1)*NSample+ss]*FactCN;
	v11.x[CLat1]=(s+1.5)*FactC1;
	v11.x[CLat2]=(ss+1.5)*FactC2;
	v11.x[CNorm]=Plot[l][(s+1)*NSample+ss+1]*FactCN;
	v10.x[CLat1]=(s+.5)*FactC1;
	v10.x[CLat2]=(ss+1.5)*FactC2;
	v10.x[CNorm]=Plot[l][s*NSample+ss+1]*FactCN;
	if(Plot[l][(s)*NSample+ss] <= 0.){
	  if(l==0)
	    v00.x[CNorm]=Plot[3][(s)*NSample+ss+1]*FactCN;
	  if(l==1)
	    v00.x[CNorm]=Plot[2][(s)*NSample+ss+1]*FactCN;
	  if(l==2) 
	    v00.x[CNorm]=Plot[1][(s)*NSample+ss+1]*FactCN;
	  if(l==3) 
	    v00.x[CNorm]=Plot[0][(s)*NSample+ss+1]*FactCN;
	  if(v00.x[CNorm] <= 0.) continue;
	}
	if(Plot[l][(s+1)*NSample+ss] <= 0.){
	  if(l==0)
	    v01.x[CNorm]=Plot[3][(s)*NSample+ss]*FactCN;
	  if(l==1)
	    v01.x[CNorm]=Plot[2][(s)*NSample+ss]*FactCN;
	  if(l==2) 
	    v01.x[CNorm]=Plot[1][(s)*NSample+ss]*FactCN;
	  if(l==3) 
	    v01.x[CNorm]=Plot[0][(s)*NSample+ss]*FactCN;
	  if(v01.x[CNorm] <= 0.) continue;
	}
	if(Plot[l][s*NSample+ss+1] <= 0.){
	  if(l==0)
	    v10.x[CNorm]=Plot[3][(s)*NSample+ss]*FactCN;
	  if(l==1)
	    v10.x[CNorm]=Plot[2][(s)*NSample+ss]*FactCN;
	  if(l==2) 
	    v10.x[CNorm]=Plot[1][(s)*NSample+ss]*FactCN;
	  if(l==3) 
	    v10.x[CNorm]=Plot[0][(s)*NSample+ss]*FactCN;
	  if(v10.x[CNorm] <= 0.) continue;
	}
	if(Plot[l][(s+1)*NSample+(ss+1)] <= 0.){
	  if(l==0)
	    v11.x[CNorm]=Plot[3][(s)*NSample+ss+1]*FactCN;
	  if(l==1)
	    v11.x[CNorm]=Plot[2][(s)*NSample+ss+1]*FactCN;
	  if(l==2) 
	    v11.x[CNorm]=Plot[1][(s)*NSample+ss+1]*FactCN;
	  if(l==3) 
	    v11.x[CNorm]=Plot[0][(s)*NSample+ss+1]*FactCN;
	  if(v11.x[CNorm] <= 0.) continue;
	}
	DrDoubleTria(&v00,&v01,&v11,&v10,&vN);
      }
    }
  }
  glEndList();
  //Freeing
  for(int l=0;l<NLevel;l++){
    free(Plot[l]);
  }
  free(Plot);
}
void ElPoly::DrSpectrum(){
  int NSample = 32;
  double Unity = 1./(double)NSample*pEdge(0)*InvScaleUn;
  double *Plot = (double *) calloc(SQR(NSample),sizeof(double));
  double *PlotR = (double *) calloc(SQR(NSample),sizeof(double));
  SampleSurface(Plot,NSample,NChType);
  //InterBSpline2D(Plot,PlotR,NSample,NSample);
  Spettro2d(Plot,NSample);
  double Min = 0.;
  double Max = 1000000.;
  for(int s=0;s<SQR(NSample);s++){
    Plot[s] += .5*pEdge(CNorm)/ScaleFact;
    if(Max < Plot[s]) Max = Plot[s];
    if(Min > Plot[s]) Min = Plot[s];
  }
  //Drawing
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  DrSmooth(Plot,NSample,Min,Max);
  glEndList();
  //Freeing
  free(PlotR);
  free(Plot);
  return;
  double *Points = (double *)calloc(NSample,sizeof(double));
  Spettro2d(Points,NSample,NChType);
  double Max1=0.,Max2=0.;
  for(int s=0;s<NSample;s++){
    if(Max1 < Points[s]){
      Max2 = Max1;
      Max1 = Points[s];
      continue;
    }
    if(Max2 < Points[s]){
      Max2 = Points[s];
      continue;
    }
  }
  for(int s=0;s<NSample;s++){
    Points[s] /= Max2;
  }
  //SampleSurface(Plot,NSample,2,CHAIN_DOWN);
  for(int s=NSample/2;s<NSample-1;s++){
    glPushMatrix();
    glBegin(GL_LINES);
    glColor4f(1.,0.,1.,1.);
    glVertex3d((2*s-NSample)*Unity,0.,Points[s]);
    glVertex3d((2*(s+1)-NSample)*Unity,0.,Points[s+1]);
    glEnd();
    glPopMatrix();
  }
  glEndList();
  free(Points);
  free(Plot);
}
void ElPoly::DrDerivative(){
  int NSample = 20;
  double Unity= 1./(double)NSample*pEdge(0)*InvScaleUn;
  SPLINE Weight;
  Weight.a0 = 1.; Weight.a1 = 0.; Weight.a2=1.;Weight.a3=0.;Weight.a4=0.;
  Matrice *Surface= new Matrice(NSample,NSample);
  SampleSurface(Surface,NSample,NChType);
  Matrice *Resp = new Matrice(NSample,NSample);
  SpatialDerivative(Surface,Resp,Weight,NSample);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  for(int s=1;s<NSample-2;s++){
    for(int ss=1;ss<NSample-2;ss++){
      //printf("Plot[%d][%d]= %lf\n",s,ss,Plot[s][ss]);
      glPushMatrix();//Particle
      glColor4f(0.,1.,1.,1.);
      glBegin(GL_TRIANGLES);
      glNormal3f(0.,0.,1.);
      glVertex3d((s+.5)*Unity,(ss+.5)*Unity,Resp->Val(s,ss)*InvScaleUn);
      glVertex3d((s+1.5)*Unity,(ss+.5)*Unity,Resp->Val(s+1,ss)*InvScaleUn);
      glVertex3d((s+1.5)*Unity,(ss+1.5)*Unity,Resp->Val(s+1,ss+1)*InvScaleUn);      
      glEnd();
      glPopMatrix();//Particle
      glPushMatrix();//Particle
      glColor4f(0.5,.0,0.5,1.);
      glBegin(GL_TRIANGLES);
      glNormal3f(0.,0.,1.);
      glVertex3d((s+.5)*Unity,(ss+.5)*Unity,Resp->Val(s,ss)*InvScaleUn);
      glVertex3d((s+1.5)*Unity,(ss+1.5)*Unity,Resp->Val(s+1,ss+1)*InvScaleUn);
      glVertex3d((s+.5)*Unity,(ss+1.5)*Unity,Resp->Val(s,ss+1)*InvScaleUn);
      glEnd();
      glPopMatrix();//Particle
    }
  }
  delete Surface;
  delete Resp;
  glEndList();
}
void ElPoly::DrIsoipse(int Values,int NIsoipse,int CoordN){
  if(CoordN < 0 || CoordN > 3) return;
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  double InvIsoipse = 1./(double)NIsoipse;
  int Coord1 = (CoordN+1)%3;
  int Coord2 = (CoordN+2)%3;
  int *Plot = (int *) calloc(Values*Values*NIsoipse,sizeof(int));
  for(int p=0;p<pNPart();p++){
    int v = (int) (pPos(p,Coord1)*NIsoipse/pEdge(Coord1));
    if( v < 0 || v > Values) continue;
    int vv = (int) (pPos(p,Coord2)*NIsoipse/pEdge(Coord2));
    if( vv < 0 || vv > Values) continue;
    int i = (int) (pPos(p,CoordN)*NIsoipse/pEdge(CoordN));
    if( i < 0 || i > NIsoipse) continue;
    Plot[(v*Values + vv)*NIsoipse + i] = p;
  }
  for(int i=0;i<NIsoipse;i++){
    glBegin(GL_POLYGON);
    glNormal3f(0.,0.,1.);
    double Zed = i*pEdge(CoordN)*InvIsoipse*InvScaleUn*ScaleFact;
    for(int p=0;p<pNPart();p++){
      int v = (int) (pPos(p,CoordN)*NIsoipse/pEdge(CoordN));
      if( v != i) continue;
      for(int v=0;v<Values-1;v++){
	int pMin = 0;
	int pMax = 0;
	for(int vv=0;vv<Values;vv++){
	  if(Plot[(v*Values+vv)*NIsoipse+i] == 0 && Plot[(v*Values+vv+1)*NIsoipse+i] != 0) 
	    pMin = Plot[(v*Values+vv+1)*NIsoipse+i];
	  if(Plot[(v*Values+vv)*NIsoipse+i] != 0 && Plot[(v*Values+vv+1)*NIsoipse+i] == 0){
	    pMax = Plot[(v*Values+vv)*NIsoipse+i];
	    break;
	  }
	}
	glVertex3d(pPos(pMin,Coord1)*InvScaleUn,pPos(pMin,Coord2)*InvScaleUn,Zed);
	glVertex3d(pPos(pMax,Coord1)*InvScaleUn,pPos(pMax,Coord2)*InvScaleUn,Zed);
      }
      glEnd();
    }
  }
  free(Plot);
  glEndList();
}
void ElPoly::DrShell(){
  int NSample = 40;
  PART *Pn = (PART *)calloc(NSample*NSample,sizeof(PART));
}
void ElPoly::DrVoronoi(){
  int NDim = 2;
  int NCell = pNChain();
  int Batch = 1000;
  int DontCreate = 4;
  int SampleType = 1;//Halton
  int NSample = 1000;
  int NItMax = 4;
  int NItFix = 1;
  int Seme = 4532643;
  int NIt = 0;
  double ItDiff = 0.;
  double Energy = 0.;
  double *Pos = (double *)calloc(NCell*NDim,sizeof(double));
  for(int c=0;c<NCell;c++){
    for(int d=0;d<NDim;d++){
      Pos[c*NDim+d] = Ch[c].Pos[d];
    }
  }
  //cvt(NDim,NCell,Batch,DontCreate,SampleType,NSample,NItMax,NItFix,&Seme,Pos,&NIt,&ItDiff,&Energy);  
  Hexagon = Dr->DefHexagon();
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  for(int c=0;c<NCell;c++){
    glPushMatrix();
    glColor4f(.2,.5,.3,1.);
    glBegin(GL_POINTS);
    glNormal3f(0.,0.,1.);
    glVertex3d(Pos[c*NDim+0]*InvScaleUn,
	       Pos[c*NDim+1]*InvScaleUn,.5*pEdge(2)*InvScaleUn);
    glEnd();
    //glCallList(Hexagon);
    glPopMatrix();
    glColor4f(.8,.2,.5,1.);
    glPushMatrix();
    glTranslated(Ch[c].Pos[0]*InvScaleUn,
      		 Ch[c].Pos[1]*InvScaleUn,.5*pEdge(2)*InvScaleUn);
    glCallList(Hexagon);
    glPopMatrix();
  }
  glEndList();
}
void ElPoly::DrSurface(){
  DrVoronoi();
  int Vertex = 6;
  int **Triangle = (int **)calloc(Vertex+1,sizeof(int));
  for(int v=0;v<Vertex+1;v++){
    *(Triangle+v) = (int *)calloc(pNChain(),sizeof(int));
  }
  //Arrange(Triangle,Vertex);
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  for(int c=0;c<pNChain();c++){
    int Chc = Ch[c].Type;
    if(!CHAIN_IF_TYPE(Chc,NChType)) continue;
    // glBegin(GL_POLYGON);
    int l = Triangle[Vertex-1][c];
    int l1 = Triangle[0][c];
    glPushMatrix();//Particle
    // double PosX1 = .333*(Ch[l].Pos[0]+Ch[c].Pos[0]+Ch[l1].Pos[0])*InvScaleUn;
    // double PosY1 = .333*(Ch[l].Pos[1]+Ch[c].Pos[1]+Ch[l1].Pos[1])*InvScaleUn;
    // double PosZ1 = .333*(Ch[l].Pos[2]+Ch[c].Pos[2]+Ch[l1].Pos[2])*InvScaleUn;
    for(int v=0;v<Vertex;v++){
      l = Triangle[v][c];
      if(v<Vertex-1) l1 = Triangle[v+1][c];
      else l1 = Triangle[0][c];
      // double PosX = .333*(Ch[l].Pos[0]+Ch[c].Pos[0]+Ch[l1].Pos[0])*InvScaleUn;
      // double PosY = .333*(Ch[l].Pos[1]+Ch[c].Pos[1]+Ch[l1].Pos[1])*InvScaleUn;
      // double PosZ = .333*(Ch[l].Pos[2]+Ch[c].Pos[2]+Ch[l1].Pos[2])*InvScaleUn;
      double Shift[3] = {0.,0.,0.};
      for(int d=0;d<3;d++) Shift[d] = -floor((Ch[c].Pos[d]-Ch[l].Pos[d])/pEdge(d)+.5)*pEdge(d);
      glPushMatrix();//Line
      glBegin(GL_LINES);
      glNormal3f(0.,0.,1.);
      glVertex3f(Ch[c].Pos[0]*InvScaleUn,Ch[c].Pos[1]*InvScaleUn,Ch[c].Pos[2]*InvScaleUn*ScaleFact);
      glVertex3f((Ch[l].Pos[0]-Shift[0])*InvScaleUn,(Ch[l].Pos[1]-Shift[1])*InvScaleUn,(Ch[l].Pos[2]-Shift[2])*InvScaleUn*ScaleFact);
      glEnd();
      //PosX1 = PosX;PosY1=PosY;PosZ1=PosZ;
      glPopMatrix();//Line
    }
    //    printf("\n");
    //    printf("(%lf %lf %lf) \n",Ch[c].Pos[0],Ch[c].Pos[1],Ch[c].Pos[2]);    
    glColor4f(1.,1.,0.,1.);
    //glEnd();
    glPopMatrix();//Particle
    glPushMatrix();
    glTranslated(Ch[c].Pos[0]*InvScaleUn,
		 Ch[c].Pos[1]*InvScaleUn,
		 Ch[c].Pos[2]*InvScaleUn*ScaleFact);
    glColor4f(1.,0.,0.,1.);
    glCallList(Point);
    //    glutSolidSphere(.05,20,20);
    glPopMatrix();    
    glPushMatrix();//Info
    glColor3f(1.0,1.0,1.0);
    //    glTranslated(Ch[c].Pos[0]*InvScaleUn,Ch[c].Pos[1]*InvScaleUn,Ch[c].Pos[2]*InvScaleUn);
    glRasterPos3f(Ch[c].Pos[0]*InvScaleUn,
		  Ch[c].Pos[1]*InvScaleUn,
		  Ch[c].Pos[2]*InvScaleUn*ScaleFact);
    char *Numero = (char*)calloc(20,sizeof(char));
    sprintf(Numero,"%d",c);//CalcnPos(Ch[c].Pos));
    for(int i=0;i<strlen(Numero);i++)
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,Numero[i]);
    free(Numero);
    glPopMatrix();//Info
  }
  glEndList();
  free(Triangle);
  //   GLuint *Surface = (GLuint *)malloc(3*Gen->NChain*sizeof(GLuint));
  //   for(int c=0,cc=0;c<Gen->NChain;c++){
  //     for(int d=0;d<3;d++){
  //       Surface[3*c+d] = (GLuint) Ch[c].Pos[d];
  //     }
  //     printf("(%u %u %u) (%lf %lf %lf) \n",Surface[3*c+0],Surface[3*c+1],Surface[3*c+2],Ch[c].Pos[0],Ch[c].Pos[1],Ch[c].Pos[2]);
  //   }
  //glEnable(GL_VERTEX_ARRAY);
  //  glDrawElements(GL_TRIANGLES,Gen->NChain,GL_UNSIGNED_INT,Surface);
  //free(Surface);
}
void ElPoly::Tile(){
  double Dx = sqrt(pEdge(CLat1)*pEdge(CLat2)/(double)pNPart());
  int NLat1 = (int)(pEdge(CLat1)/Dx);
  int NLat2 = (int)(pEdge(CLat2)/Dx);
  printf("Dx %lf NLat1 %d NLat2 %d %d\n",Dx,NLat1,NLat2,SQR(NLat1) + SQR(NLat2));
  if( NLat1*NLat2 != pNPart()) return;
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  Vettore v1(3);
  Vettore v2(3);
  Vettore v3(3);
  Vettore v4(3);
  Vettore vN(3);
  printf("Ciccia\n");
  for(int nx = 0;nx<NLat1-2;nx+=2){
    for(int ny = 0;ny<NLat2-2;ny+=2){
      int l1 = nx*NLat2 + ny;
      int l2 = nx*NLat2 + ny + 1;
      int l3 = (nx+1)*NLat2 + ny + 1;
      int l4 = (nx+1)*NLat2 + ny;
      v1.x[0]=pPos(l1,0)*InvScaleUn;
      v1.x[1]=pPos(l1,1)*InvScaleUn;
      v1.x[2]=pPos(l1,2)*InvScaleUn*ScaleFact;
      v2.x[0]=pPos(l2,0)*InvScaleUn;
      v2.x[1]=pPos(l2,1)*InvScaleUn;
      v2.x[2]=pPos(l2,2)*InvScaleUn*ScaleFact;
      v3.x[0]=pPos(l3,0)*InvScaleUn;
      v3.x[1]=pPos(l3,1)*InvScaleUn;
      v3.x[2]=pPos(l3,2)*InvScaleUn*ScaleFact;
      v4.x[0]=pPos(l4,0)*InvScaleUn;
      v4.x[1]=pPos(l4,1)*InvScaleUn;
      v4.x[2]=pPos(l4,2)*InvScaleUn*ScaleFact;
      vN.NormalSurf(&v1,&v2,&v3);
      DrDoubleTria(&v1,&v2,&v3,&v4,&vN);
    }
  }
  glEndList();
}
void ElPoly::DrIsolevel(int NSample,double IsoLevel){
  NSample = 20;
  double VolEl = pVol()/(double)CUB(NSample);
  double *Plot = (double *)calloc(CUBE(NSample),sizeof(double));
  double Min = 0.;
  double Max = 0.;
  VAR_TRIANGLE *Tri = NULL;
  for(int p=0;p<pNPart();p++){
    //if(Pm[p].Typ != 1) continue;
    double Posx = pPos(p,0) - floor(pPos(p,0)*pInvEdge(0))*pEdge(0);
    double Posy = pPos(p,1) - floor(pPos(p,1)*pInvEdge(1))*pEdge(1);
    double Posz = pPos(p,2) - floor(pPos(p,2)*pInvEdge(2))*pEdge(2);
    int sx = (int)(Posx*pInvEdge(0)*NSample);
    int sy = (int)(Posy*pInvEdge(1)*NSample);
    int sz = (int)(Posz*pInvEdge(2)*NSample);
    int sTot = (sx*NSample+sy)*NSample+sz;
    Plot[sTot] += VolEl;
    if(Max < Plot[sTot]) Max = Plot[sTot];
    if(Min > Plot[sTot]) Min = Plot[sTot];
  }
  //IsoLevel = .1*Max;
  printf("Min %lf Max %lf Isolevel %lf\n",Min,Max,IsoLevel);
  int NTri = 0;
  Tri = MarchingCubes(Plot,NSample,IsoLevel,&NTri);
  free(Plot);
  int NVertex = CUBE(2*NSample-1);
  double *WeightL = (double *) calloc(NVertex,sizeof(double));
  double MaxW = 1.;///NormalWeight(Tri,WeightL,NSample,NTri);
  for(int t=0;t<NTri;t++){
    int vRef = Tri[t].v[0];
    //printf("%d %lf\n",t,WeightL[vRef]);
  }
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  for(int n=0;n<pNNano();n++) DrawNano(n);
  for(int t=0;t<NTri;t++){
    //if(Tri[t].p[0].x[2]> Nano->Pos[2] + 3.||Tri[t].p[0].x[2] < Nano->Pos[2] - 3.) continue;
    int vRef = Tri[t].v[0];
    double Green = WeightL[vRef]*MaxW + .3 + .4*Mat->Casuale();
    glColor4f(0.1,Green,0.2,1.);
    //glColor4f(HueUp[p].r,HueUp[p].g,HueUp[p].b,HueUp[p].a);
    glBegin(GL_TRIANGLES);
    glNormal3f(Tri[t].n[0].x[0]*InvScaleUn,
    	       Tri[t].n[0].x[1]*InvScaleUn,
    	       Tri[t].n[0].x[2]*InvScaleUn*ScaleFact);
    glVertex3f(Tri[t].p[0].x[0]*InvScaleUn,
    	       Tri[t].p[0].x[1]*InvScaleUn,
    	       Tri[t].p[0].x[2]*InvScaleUn*ScaleFact);
    glNormal3f(Tri[t].n[1].x[0]*InvScaleUn,
    	       Tri[t].n[1].x[1]*InvScaleUn,
    	       Tri[t].n[1].x[2]*InvScaleUn*ScaleFact);
    glVertex3f(Tri[t].p[1].x[0]*InvScaleUn,
    	       Tri[t].p[1].x[1]*InvScaleUn,
    	       Tri[t].p[1].x[2]*InvScaleUn*ScaleFact);
    glNormal3f(Tri[t].n[2].x[0]*InvScaleUn,
    	       Tri[t].n[2].x[1]*InvScaleUn,
    	       Tri[t].n[2].x[2]*InvScaleUn*ScaleFact);
    glVertex3f(Tri[t].p[2].x[0]*InvScaleUn,
    	       Tri[t].p[2].x[1]*InvScaleUn,
    	       Tri[t].p[2].x[2]*InvScaleUn*ScaleFact);
    glEnd();
    if(IfLine){
      glColor4f(0.1,.3,1.,1.);
      glBegin(GL_LINES);
      glVertex3d(Tri[t].c.x[0]*InvScaleUn,
		 Tri[t].c.x[1]*InvScaleUn,
		 Tri[t].c.x[2]*InvScaleUn*ScaleFact);
      glVertex3d((Tri[t].c.x[0]+Tri[t].n[0].x[0])*InvScaleUn,
		 (Tri[t].c.x[1]+Tri[t].n[0].x[1])*InvScaleUn,
		 (Tri[t].c.x[2]+Tri[t].n[0].x[2])*InvScaleUn*ScaleFact);
      glEnd();
    }
  }
  glEndList();
  free(WeightL);
}
#endif //USE_GL
