#include "../../include/Draw.h"
#ifdef __glut_h__

Draw *Dr;
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
void Slide(){}
void ParticleRealTime(){}
void ParticleList(){}
void keyboard(unsigned char key,int x, int y){
  Dr->keyboardDraw(key);
}
void Menu(){}
void Figure(){
  //Dr->Draw1();
  //Dr->ShowImage();
  Dr->Transform();
  //Dr->DrCube();
  Dr->DFigure();
}
void Figure1(){
  Dr->Draw1();
}
int main(int argc,char** argv){
  char nome[60];
  Dr = new Draw();
  if(argc > 1){
    sprintf(nome,"%s",argv[1]);
    Dr->OpenImage(nome);
  }
  else ;//Dr->ApplyTexture();
  Dr->Window(argc,argv);
  glutMainLoop();//Mantiene la finestra
  return 1;
}
#else
int main(int argc,char** argv){
  printf("OpenGL not supported\n");
  return 0;
}
#endif// __glut_h__
