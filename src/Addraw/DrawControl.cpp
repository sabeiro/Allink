#include "../../include/Draw.h"

#ifdef __glut_h__
void Draw::CameraQuat(){
//   Quaternione DxSx;
//   Quaternione SuGiu;
//   Quaternione Union;
//   float Matrix[16];
//   DxSx.AxisRotation(1.,0.,0.,AngleDxSx);
//   SuGiu.AxisRotation(0.,1.,0.,AngleSuGiu);
//   Union = DxSx * SuGiu;
//   Union.CreateMatrix(Matrix);
//   glMultMatrixf(Matrix);
}
void Draw::ChangeDxSx(GLfloat Movement){
  if( POS(Movement) < IncrVisDxSx )
    xa += Movement;
  else{
    if(Movement < 0)
      xa -= IncrVisDxSx;
    else
      xa += IncrVisDxSx;
  }
  if( xa > 360.)
    xa = -360.;
  if( xa < -360.)
    xa = 360.;
  //yp += (float)(sin(xa/PI*180.));
}
void Draw::ChangeSuGiu(GLfloat Movement){
  if( POS(Movement) < IncrVisSuGiu ){
    if(ya > 90 && ya < 270
       || (ya < -90 && ya > -270))
      ya -= Movement;
    else 
      ya += Movement;
  }
  else{
    if(Movement < 0){
      if(ya > 90 && ya < 270
	 || (ya < -90 && ya > -270))
	ya += IncrVisSuGiu;
      else 
	za -= IncrVisSuGiu;
      //ya -= IncrVisDxSx;
    }
    else {
      if(ya > 90 && ya < 270
	 || (ya < -90 && ya > -270))
	ya -= IncrVisSuGiu;
      else 
	za += IncrVisSuGiu;
      //ya += IncrVisDxSx;
    }
  }
  if( ya > 360.)
    ya = -360.;
  if( ya < -360.)
    ya = 360.;
  //xp += (float)(sin(ya/PI*180.));
  //zp -= (float)(cos(ya/PI*180.));
}
void Draw::DMouseMove(int x,int y){
  GLfloat Diff = 0.;
  if(y > yRem){
    Diff = 2.*(y - yRem);
    ChangeDxSx(Diff);
  }
  else if(y < yRem){
    Diff = 2.*(y - yRem);
    ChangeDxSx(Diff);
  }
  if(x > xRem){
    Diff = 2.*(x - xRem);
    ChangeSuGiu(Diff);
  }
  else if(x < xRem){
    Diff = 2.*(x - xRem);
    ChangeSuGiu(Diff);
  }
  xRem = x;
  yRem = y;
  //printf("(%f,%f,%lf)\n",xa,ya,za);
  glutPostRedisplay();
}
void Draw::Dmouse(int button, int state,int x,int y){
  switch (button){
  case GLUT_RIGHT_BUTTON:
    // xRem = x;
    // yRem = y;
    break;
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
    zw -= .02;
    glutPostRedisplay();
    break;
  case GLUT_WHEEL_UP:
    zw += .02;
    glutPostRedisplay();
    break;
  default:
    break;
  }
}
void Draw::Dspecial(int k, int x, int y)
{
  switch (k) {
  case GLUT_KEY_UP:
    xf += 5.0;
    break;
  case GLUT_KEY_DOWN:
    xf -= 5.0;
    break;
  case GLUT_KEY_LEFT:
    yf += 5.0;
    break;
  case GLUT_KEY_RIGHT:
    yf -= 5.0;
    break;
  default:
    return;
  }
  glutPostRedisplay();
}
void Draw::ChooseBlend(int Which){
  switch(Which){
  case 1:
    BlendSource = GL_SRC_ALPHA;
    BlendDest = GL_ONE_MINUS_SRC_ALPHA;// GL_ONE_MINUS_DST_COLOR;
    break;
  case 2:
    BlendSource = GL_SRC_ALPHA;
    BlendDest = GL_ONE_MINUS_DST_ALPHA;
    break;
  case 3:
    BlendSource = GL_DST_COLOR;
    BlendDest = GL_ONE_MINUS_SRC_ALPHA;
    break;
  case 4:
    BlendSource = GL_ONE_MINUS_SRC_ALPHA;
    BlendDest = GL_ONE_MINUS_DST_COLOR;
    break;
  case 5:
    BlendSource = GL_DST_COLOR;
    BlendDest = GL_SRC_COLOR;
    break;
  case 6:
    BlendSource = GL_ONE;
    BlendDest = GL_SRC_COLOR;
    break;
  case 7:
    BlendSource = GL_SRC_COLOR;
    BlendDest = GL_ONE_MINUS_DST_COLOR;
    break;
  case 8:
    BlendSource = GL_ONE_MINUS_SRC_ALPHA;
    BlendDest = GL_ONE_MINUS_DST_COLOR;
    break;
  case 9:
    BlendSource =GL_DST_COLOR ;
    BlendDest = GL_ZERO;
    break;
  case 10:
    BlendSource = GL_DST_COLOR;
    BlendDest = GL_ONE_MINUS_SRC_ALPHA;
    break;
  }
}
void Draw::keyboardDraw(unsigned char key){
  switch (key){
  case 'a':
    break;
  case 'A':
    IfImage=0;
    glutPostRedisplay();
    break;
  case 'B':
    if(IfBlend)
      glEnable(GL_BLEND); 
    else
      glDisable(GL_BLEND); 
    IfBlend != IfBlend;
    break;
  case 'C':
    printf("Coordinate (%.2f,%.2f,%.2f) Angolo (%.2f,%.2f,%.2f) Luce  (%.2f,%.2f,%.2f) Ruota %lf)\n"
	   ,xp,yp,zp,xa,ya,za,xf,yf,zf,zw);
    break;
  case 'f':
    glutFullScreen();
    glutPostRedisplay();
    break;
  case 'g':
    gr += 1;pr=0;
    if(gr == 3) gr=0;pr=1;
    Lista(Values);
    sprintf(info,"Grid visualized every %d system unit",GridStep);
    glutPostRedisplay();
    break;
  case 'G':
    GridStep += 1.;
    for(int d=0;d<3;d++)
      GridEdge[d] = (int)(GridStep*Edge[d]*InvScaleUn);
    Lista(Values);
    sprintf(info,"Grid visualized every %d system unit",GridStep);
    glutPostRedisplay();
    break;
  case 'h':
    glutIdleFunc(NULL);
    break;
  case 'l':
    la = !la;
    glutPostRedisplay();
    break;
  case 'L':
    IfMaterial = !IfMaterial;
    glutPostRedisplay();
    break;
  case 'k':
    Step++;
    IfInfo = 0;
    glutPostRedisplay();
    Picture();
    IfInfo = 1;
    //WritePngwriter();
    break;
  case 'm':
    //OpenImage("SplineAngolose.tif");
    //Menu();
    break;
  case 'M':
    //glutDestroyMenu(menu);
    break;
  case 'p':
    IfPoint++;
    if(IfPoint == 2) IfPoint =0;
    glutPostRedisplay();
    break;
  case 'P':
    pr += 1;
    if(pr == 2) pr =0;
    if(pr == 1)
      sprintf(info,"Prospective view");
    else if(pr == 0)
      sprintf(info,"Orthogonal view");
    glutPostRedisplay();
    break;
  case 'q':
    glutDestroyWindow(MainWindow);
    exit(0);
    break;
  case 's':
    IfScript = !IfScript;
    ReadScript();
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
  case 'w':
    Rback = 1.0 - Rback;
    Gback = 1.0 - Gback;
    Bback = 1.0 - Bback;
    Aback = 0.0;
    Lista(Values);
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
  case 40:
    //glutIdleFunc(Expand);
    break;
  default:
    break;
  }
}
#endif
