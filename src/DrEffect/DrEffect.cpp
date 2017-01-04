#include "DrEffect.h"
#include "DrFinestra.h"

#ifdef __glut_h__

int main(int argc,char** argv){
  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a(argc,argv);
  DrFinestra f;
  a.setMainWidget( &f );
  f.setGeometry(100,200,800,355);
  f.show();
  for(int i=1;i<argc;i++){
    if(!strcmp(argv[i],"-d")){
      f.DataFile(argv[i+1]);
    }
    else if(!strcmp(argv[i],"-c")){
      f.ConfFile(argv[i+1]);
    }
    else if(!strncmp(argv[i],"-r",2)){
      f.Run();
    }
    else{
      f.DataFile(argv[i]);
    }
  }
  //  f.setStyle(QStyle::QPalstiqueStyle);  
  return a.exec();
}
#else
int main(int argc,char** argv){
  printf("OpenGL not supported\n");
  return 0;
}
#endif// __glut_h__
