/***********************************************************************
ElPolyDrawControl: OpenGL menus, keyboards and mouse functions
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
#include "ElPoly.h"

#ifdef __glut_h__
extern Draw *Dr;
/// Type to visualize
enum ElType{
  /// Every type of chain
  EL_EVERY  = 0x0021 ,
  /// Only the polytype
  EL_POLY    ,
  /// Only the homotype
  EL_HOMO    ,
  /// Only the added
  EL_ADDED   ,
  /// Only the upwards
  EL_UP      ,
  /// Only the downwards
  EL_DOWN    ,
  /// Only the downwards
  EL_INNER   ,
  /// Only the downwards
  EL_OUTER   ,
  /// Only the stretched
  EL_STRETCH ,
  /// Only the flabby
  EL_FLABBY  ,
  /// Only the tilted
  EL_TILTED  ,
  /// Only the hydrophobic
  EL_PHOB    ,
  /// Only the hydrophilic
  EL_PHIL    ,
  /// Only the added
  EL_OIL     ,
  /// every
  EL_ALL     
};
enum ElExtra{
  EL_FULL =   0x0031,
  EL_BOND,
  EL_PICTURE,
  EL_RESTORE,
  EL_SCRIPT,
  EL_RENDER,
};
enum ElEffect{
  EL_NANOPOS = 0x0041,
  EL_DENSPOS,
  EL_CMPOS,
  EL_CICCIA
};
enum ElBackFold{
  EL_BF_NANO = 0x0051,
  EL_BF_CM,
  EL_BF_PARTICLE
};
enum ElBlend{
  EL_BLEND1 = 0x0061,
  EL_BLEND2,
  EL_BLEND3,
  EL_BLEND4,
  EL_BLEND5,
  EL_BLEND6,
  EL_BLEND7,
  EL_BLEND8,
  EL_BLEND9,
  EL_BLEND10
};
enum ElBlock{
  EL_BLOCK = 0x0071,
};
enum ElLineSize{
  EL_LSIZE1 = 0x0081,
  EL_LSIZE2  ,
  EL_LSIZE3  ,
  EL_LSIZE4  ,
  EL_LSIZE5  ,
  EL_LSIZE6  ,
  EL_LSIZE7  ,
  EL_LSIZE20 ,
};
void ElPoly::CreateMenu(){
  int SubMenu1 = 0, SubMenu2 = 0, SubMenu3 = 0;
  int SubMenu4 = 0, SubMenu5 = 0, SubMenu6 = 0;
  int SubMenu7 = 0, SubMenu8 = 0, SubMenu9 = 0;
  SubMenu1 = glutCreateMenu(MenuVisual);
  glutAddMenuEntry("Particles",EL_PART);
  glutAddMenuEntry("Chains",EL_CHAIN);
  glutAddMenuEntry("Vectors",EL_VECTORS);
  glutAddMenuEntry("Triangulate",EL_TRIA);
  glutAddMenuEntry("IsoLevel",EL_ISOLEVEL);
  glutAddMenuEntry("Color",EL_COLOR);
  glutAddMenuEntry("Quad",EL_QUAD);
  glutAddMenuEntry("Quad1",EL_QUAD1);
  glutAddMenuEntry("Potential",EL_POTENTIAL);
  glutAddMenuEntry("Density",EL_DENS);
  glutAddMenuEntry("Polygon",EL_POLYGON);
  glutAddMenuEntry("SquareMesh",EL_SQUAREMESH);
  glutAddMenuEntry("Sample",EL_SAMPLE);
  glutAddMenuEntry("Spectrum",EL_SPECTRUM);
  glutAddMenuEntry("Surface",EL_SURF);
  glutAddMenuEntry("Voronoi",EL_VORONOI);
  glutAddMenuEntry("Skin",EL_SKIN);
  glutAddMenuEntry("Mesh",EL_MESH);
  glutAddMenuEntry("Isoipse",EL_ISOIPSE);
  SubMenu2 = glutCreateMenu(MenuChoise);
  glutAddMenuEntry("Every",EL_EVERY);
  glutAddMenuEntry("Poly",EL_POLY);
  glutAddMenuEntry("Homo",EL_HOMO);
  glutAddMenuEntry("Added",EL_ADDED);
  glutAddMenuEntry("Up",EL_UP);
  glutAddMenuEntry("Down",EL_DOWN);
  glutAddMenuEntry("Inner",EL_INNER);
  glutAddMenuEntry("Outer",EL_OUTER);
  glutAddMenuEntry("Stretch",EL_STRETCH);
  glutAddMenuEntry("Flabby",EL_FLABBY);
  glutAddMenuEntry("Tilted",EL_TILTED);
  SubMenu3 = glutCreateMenu(MenuChoise);
  glutAddMenuEntry("All",EL_ALL);
  glutAddMenuEntry("Phob",EL_PHOB);
  glutAddMenuEntry("Phil",EL_PHIL);
  glutAddMenuEntry("Oil",EL_OIL);
  SubMenu4 = glutCreateMenu(MenuChoise);
  glutAddMenuEntry("NanoPos",EL_NANOPOS);
  glutAddMenuEntry("DensPos",EL_DENSPOS);
  glutAddMenuEntry("CmPos",EL_CMPOS);
  SubMenu5 = glutCreateMenu(MenuChoise);
  glutAddMenuEntry("Part",EL_BF_PARTICLE);
  glutAddMenuEntry("Nano",EL_BF_NANO);
  glutAddMenuEntry("CM",EL_BF_CM);
  SubMenu6 = glutCreateMenu(MenuChoise);
  glutAddMenuEntry("1",EL_BLEND1);
  glutAddMenuEntry("2",EL_BLEND2);
  glutAddMenuEntry("3",EL_BLEND3);
  glutAddMenuEntry("4",EL_BLEND4);
  glutAddMenuEntry("5",EL_BLEND5);
  glutAddMenuEntry("6",EL_BLEND6);
  glutAddMenuEntry("7",EL_BLEND7);
  glutAddMenuEntry("8",EL_BLEND8);
  glutAddMenuEntry("9",EL_BLEND9);
  glutAddMenuEntry("10",EL_BLEND10);
  SubMenu7 = glutCreateMenu(MenuChoise);
  glutAddMenuEntry("Every",EL_BLOCK+pNBlock());
  for(int b=0;b<pNBlock();b++){
    glutAddMenuEntry(Block[b].Name,b+EL_BLOCK);
  }
  SubMenu8 = glutCreateMenu(MenuChoise);
  glutAddMenuEntry("1",EL_LSIZE1);
  glutAddMenuEntry("2",EL_LSIZE2);
  glutAddMenuEntry("3",EL_LSIZE3);
  glutAddMenuEntry("4",EL_LSIZE4);
  glutAddMenuEntry("5",EL_LSIZE5);
  glutAddMenuEntry("6",EL_LSIZE6);
  glutAddMenuEntry("7",EL_LSIZE7);
  glutAddMenuEntry("20",EL_LSIZE20);
  glutCreateMenu(MenuChoise);
  glutAddSubMenu("Visualisation",SubMenu1);
  glutAddSubMenu("Chain",SubMenu2);
  glutAddSubMenu("Bead",SubMenu3);
  glutAddSubMenu("Move",SubMenu4);
  glutAddSubMenu("BackFold",SubMenu5);
  glutAddSubMenu("Blend",SubMenu6);
  glutAddSubMenu("Block",SubMenu7);
  glutAddSubMenu("Line size",SubMenu8);
  glutAddMenuEntry("Fullscreen",EL_FULL);
  glutAddMenuEntry("Bond",EL_BOND);
  glutAddMenuEntry("Picture",EL_PICTURE);
  glutAddMenuEntry("Restore",EL_RESTORE);
  glutAddMenuEntry("Script",EL_SCRIPT);
  glutAddMenuEntry("Render",EL_RENDER);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}
extern void ElPoly::ElMenuVisual(int option){
  InvScaleUn = 1./pEdge(0);
  for(int d=0;d<3;d++) Dr->Edge[d] = pEdge(d);
  Dr->InvScaleUn = InvScaleUn;
  Point = Dr->DefPoint();
  What2Draw = option;
  //option = EL_PART;
  switch(option){
  case EL_PART:
    DrPartList();
    break;
  case EL_COLOR:
    //glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_DST_COLOR);  
    DrColor();
    break;
  case EL_CHAIN:
    BfDefChain();
    Hexagon = Dr->DefHexagon();
    DrChains();
    break;
  case EL_VECTORS:
    Arrow = Dr->DefArrowThin();
    DrVectors();
    break;
  case EL_SURF:
    DrSurface();
    break;
  case EL_POLYGON:
    DrPolygon();
    break;
  case EL_SQUAREMESH:
    DrSquareMesh();
    break;
  case EL_SAMPLE:
    DrSample(NEdge);
    //ChooseDraw("sample");
    break;
  case EL_SKIN:
    //DrStalk();
    //DefineSkin(5);
    //Tile();
    break;
  case EL_TRIA:
    DrTriangulate();
    break;
  case EL_MESH:
    //DefineSurf();
    //DrGenMesh();
    break;
  case EL_ISOLEVEL:
    DrIsolevel(10,ExtParam);
    break;
  case EL_ISOIPSE:
    DrIsoipse(100,10,1);
    break;
  case EL_QUAD:
    Quad = Dr->DefQuad(NEdge-2);
    DrQuad();
    break;
  case EL_QUAD1:
    Quad = Dr->DefQuad(NEdge-2);
    DrQuad1();
    break;
  case EL_POTENTIAL:
    DrPotential();
    break;
  case EL_DENS:
    Dr->xp = .5*pEdge(0)*InvScaleUn;
    Dr->yp = 0.;
    Dr->zp = .5*pEdge(2)*InvScaleUn;
    DrDensity();
    break;
  case EL_CROSS:
    DrCrossLinks();
    break;
  }
}
extern void ElPoly::ElMenuChoise(int option){
  int OldIf = 0;
  InvScaleUn = 1./pEdge(0);
  for(int d=0;d<3;d++) Dr->Edge[d] = pEdge(d);
  Dr->InvScaleUn = InvScaleUn;
  Point = Dr->DefPoint();
  //option = EL_PART;
  switch(option){
  case EL_EVERY:
    NChType = CHAIN_EVERY;
    //sprintf(Dr->info,"Chain type %d",NChType);
    RenderPart();
    break;
  case EL_POLY:
    NChType = CHAIN_POLY;
    RenderPart();
    break;
  case EL_HOMO:
    NChType = CHAIN_HOMO;
    RenderPart();
    break;
  case EL_ADDED:
    NChType = CHAIN_ADDED;
    RenderPart();
    break;
  case EL_UP:
    NChType = CHAIN_UP;
    RenderPart();
    break;
  case EL_DOWN:
    NChType = CHAIN_DOWN;
    RenderPart();
    break;
  case EL_INNER:
    NChType = CHAIN_INNER;
    RenderPart();
    break;
  case EL_OUTER:
    NChType = CHAIN_OUTER;
    RenderPart();
    break;
  case EL_STRETCH:
    NChType = CHAIN_STRETCH;
    RenderPart();
    break;
  case EL_FLABBY:
    NChType = CHAIN_FLABBY;
    RenderPart();
    break;
  case EL_TILTED:
    NChType = CHAIN_TILTED;
    RenderPart();
    break;
  case EL_PHOB:
    NPType = BEAD_PHOB;
    RenderPart();
    break;
  case EL_PHIL:
    NPType = BEAD_PHIL;
    RenderPart();
    break;
  case EL_OIL:
    NPType = BEAD_OIL;
    RenderPart();
    break;
  case EL_ALL:
    NPType = BEAD_EVERY;
    RenderPart();
    break;
  case EL_FULL:
    glutFullScreen();
    glutPostRedisplay();
    break;
  case EL_BOND:
    if(!pNLink()) IfLine += 1;
    if(IfLine == 2) IfLine =0;
    RenderPart();
    break;
  case EL_PICTURE:
    Dr->Step++;
    Dr->IfInfo = 0;
    Dr->Picture();
    Dr->IfInfo = 1;
    break;
  case EL_RESTORE:
    glutIdleFunc(NULL);
    RenderPart();
    Dr->InitConstant();
    glutPostRedisplay();
    break;
  case EL_SCRIPT:
    Dr->IfScript = !Dr->IfScript;
    Dr->ReadScript();
    break;
  case EL_RENDER:
    DrawOutput = EL_POVRAY;
    RenderPart();
    DrawOutput = EL_OPENGL;
    //Conv2Povray();
    break;
  case EL_NANOPOS:
    Dr->xp = pNanoPos(0,0)*InvScaleUn;
    Dr->yp = pNanoPos(0,1)*InvScaleUn;
    Dr->zp = pNanoPos(0,2)*InvScaleUn;
    break;
  case EL_DENSPOS:
    Dr->xp = .5*pEdge(0)*InvScaleUn;
    Dr->yp = 0.;
    Dr->zp = .5*pEdge(2)*InvScaleUn;
    break;
  case EL_CMPOS:
    ShiftSys(SHIFT_CM);
    break;
  case EL_SAMPLE:
    DrSample(NEdge);
    break;
  case EL_SPECTRUM:
    DrSpectrum();
    break;
  case EL_VORONOI:
    DrVoronoi();
    break;
  case EL_BF_NANO:
    BackFold(BF_NANO);
    RenderPart();
    break;
  case EL_BF_PARTICLE:
    BackFold(BF_PART);
    RenderPart();
    break;
  case EL_BF_CM:
    BackFold(BF_CHAIN);
    RenderPart();
    break;
  case EL_BLEND1:
    Dr->ChooseBlend(1);
    break;
  case EL_BLEND2:
    Dr->ChooseBlend(2);
    break;
  case EL_BLEND3:
    Dr->ChooseBlend(3);
    break;
  case EL_BLEND4:
    Dr->ChooseBlend(4);
    break;
  case EL_BLEND5:
    Dr->ChooseBlend(5);
    break;
  case EL_BLEND6:
    Dr->ChooseBlend(6);
    break;
  case EL_BLEND7:
    Dr->ChooseBlend(7);
    break;
  case EL_BLEND8:
    Dr->ChooseBlend(8);
    break;
  case EL_BLEND9:
    Dr->ChooseBlend(9);
    break;
  case EL_BLEND10:
    Dr->ChooseBlend(10);
    break;
  case EL_LSIZE1:
    glLineWidth(1);
    glPointSize(1);
    break;
  case EL_LSIZE2:
    glLineWidth(2);
    glPointSize(2);
    break;
  case EL_LSIZE3:
    glLineWidth(3);
    glPointSize(3);
    break;
  case EL_LSIZE4:
    glLineWidth(4);
    glPointSize(4);
    break;
  case EL_LSIZE5:
    glLineWidth(5);
    glPointSize(5);
    break;
  case EL_LSIZE6:
    glLineWidth(6);
    glPointSize(6);
    break;
  case EL_LSIZE7:
    glLineWidth(7);
    glPointSize(7);
    break;
  case EL_LSIZE20:
    glLineWidth(20);
    glPointSize(20);
    break;
    //DrPosCol();
    //DrQuad1();
    //DrInterpSurface();
  }
  for(int b=0;b<pNBlock();b++){
    if(option == EL_BLOCK + b){
      sprintf(Block2Draw,"%s",Block[b].Name);
      RenderPart();
    }
  }
  if(option == EL_BLOCK + pNBlock()){
    sprintf(Block2Draw,"ALL");
    RenderPart();
  }
  //glEndList();
}
void ElPoly::RenderPart(){
  ElMenuVisual(What2Draw);
}
extern void ElPoly::ElDrawMouse(int button, int state,int x,int y){
  Dr->Dmouse(button,state,x,y);
  switch (button){
  case GLUT_WHEEL_UP:
    if( glutGetModifiers() == GLUT_ACTIVE_CTRL){
      ScaleFact += .2;
      Dr->zp -= .12*pEdge(CNorm)*InvScaleUn;
      ElMenuVisual(What2Draw);
      glutPostRedisplay();
    }
    break;
  case GLUT_WHEEL_DOWN:
    if( glutGetModifiers() == GLUT_ACTIVE_CTRL){
      ScaleFact -= .2;
      Dr->zp += .12*pEdge(CNorm)*InvScaleUn;
      ElMenuVisual(What2Draw);
      glutPostRedisplay();
    }
    break; 
  default:
    break;
  }
}
extern void ElPoly::keyboard(unsigned char key,int x, int y){
  Dr->keyboardDraw(key);
  char string[60];
  switch (key){
  case '+':
    quando+=10;
    if(quando >= NFileTot){
      quando = 0;
    }
    if(quando >= 0){
      OpenFile(quando);
      sprintf(Dr->info,"%s %d",cFile[quando],NFile[0]);
    }
    RenderPart();
    glutPostRedisplay();
    break;
  case '-':
    quando-=10;
    if(quando < 0){
      quando = NFileTot-1;
    }
    if(quando <=NFileTot){
      OpenFile(quando);
      sprintf(Dr->info,"%s %d",cFile[quando],NFile[0]);
    }
    RenderPart();
    glutPostRedisplay();
    break;
  case '>':
    quando++;
    if(quando >= NFileTot){
      quando = 0;
    }
    if(quando >= 0){
      OpenFile(quando);
      sprintf(Dr->info,"%s %d",cFile[quando],NFile[0]);
    }
    RenderPart();
    glutPostRedisplay();
    break;
  case '<':
    quando--;
    if(quando < 0){
      quando = NFileTot-1;
    }
    if(quando <=NFileTot){
      OpenFile(quando);
      sprintf(Dr->info,"%s %d",cFile[quando],NFile[0]);
    }
    RenderPart();
    glutPostRedisplay();
    break;
  case 'b':
    if(!pNLink()) IfLine = !IfLine;
    sprintf(Dr->info,"Bonding visualisation");
    RenderPart();
    glutPostRedisplay();
    break;
  case 'c':
    NBackFold++;
    if(NBackFold > 3)
      NBackFold = 0;
    if(NBackFold == BF_PART)
      sprintf(Dr->info,"Particle back fold");
    else if(NBackFold == BF_CHAIN)
      sprintf(Dr->info,"Center of mass back fold");
    else if(NBackFold == BF_NANO)
      sprintf(Dr->info,"Nano particle back fold");
    else if(NBackFold == BF_NO)
      sprintf(Dr->info,"No back fold");
    BackFold(NBackFold);
    RenderPart();
    glutPostRedisplay();
    break;
  case 'C':
    Conv2Povray();
    break;
  case 'i':
    IfIntorno += 1;
    if(IfIntorno == 3) IfIntorno = 0;
    if(IfIntorno == 1)
      sprintf(Dr->info,"Center of mass proximity");
    if(IfIntorno == 2)
      sprintf(Dr->info,"Particle proximity");
    RenderPart();
    glutPostRedisplay();
    break;
  case 'I':
    char cSystem[STRSIZE];
    SysInfo(cSystem);
    printf("%s %d %s\n",cSystem,NFile[0],cFile[NFile[0]]);
    sprintf(Dr->info,"%s %d %s\n",cSystem,NFile[0],cFile[NFile[0]]);
    glutPostRedisplay();
    break;
  case 'k':
    Dr->IfInfo = 0;
    glutPostRedisplay();
    Dr->Picture();
    Dr->IfInfo = 1;
    break;
  case 'K':
    sprintf(string,"rm %s",cFile[quando]);
    system(string);
    quando++;
    if(quando >= NFileTot){
      quando = 0;
    }
    if(quando >= 0){
      OpenFile(quando);
      sprintf(Dr->info,"%s %d",cFile[quando],NFile[0]);
    }
    RenderPart();
    glutPostRedisplay();
    break;
  case 'l':
    LineSize += 1;
    glLineWidth(LineSize);
    glPointSize(LineSize);
    RenderPart();
    glutPostRedisplay();
    break;
  case 'L':
    LineSize -= 1.;
    glLineWidth(LineSize);
    glPointSize(LineSize);
    RenderPart();
    glutPostRedisplay();
    break;
  case 'm':
    glutPostRedisplay();
    break;
  case 'n':
    Saturation *= 1.1;
    RenderPart();
    break;
  case 'N':
    Saturation /= 1.1;
    RenderPart();
    break;
  case 'p':
    ExtParam *= 1.5;
    RenderPart();
    break;
  case 'P':
    ExtParam /= 1.5;
    RenderPart();
    break;
   case 'r':
    glutIdleFunc(NULL);
    RenderPart();
    Dr->InitConstant();
    sprintf(Dr->info,"initial configuration");
    glutPostRedisplay();
    break;
  case 'R':
    RenderPart();
    glutPostRedisplay();
    break;
  case 'S':
    glutIdleFunc(Slide);
    break;
  case 't':
    Transform(2);
    RenderPart();
    break;
  case 27:
    exit(0);
    break;
  case 40:
    //glutIdleFunc(Expand);
    break;
  default:
    break;
  }
}
#endif

