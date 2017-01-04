#include "DrScript.h"

void DrScript::DrMessage(const char * s, ...)
{
#ifdef DEBUG
  va_list args;
  va_start(args, s);
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
#else
  return;
#endif
}
DrScript::DrScript(){
  Slide = 0;
}
DrScript::~DrScript(){
}
void DrScript::IncrSlide(){
  Slide++;
  repaint();
}
void DrScript::paintEvent(QPaintEvent *){
  DrMessage("DrScript.paintEvent");
  QPainter p(this);
  p.setPen(QPen(Qt::black,1,Qt::DashLine,Qt::RoundCap,Qt::RoundJoin));
  p.drawRect(0,0,width(),height());
  p.setFont(QFont("Helvetica",12 ));
  ReadSlide(&p);
}
void DrScript::ReadSlide(QPainter *p){
  DrMessage("DrScript.ReadSlide");
  if((File2Read = fopen("DrScript.dr","r"))==0){
    return ;
  }
  char *cLine = (char *)calloc(256,sizeof(char));
  for(int k=0;!(fgets(cLine,256,File2Read)==NULL);k++){
    if(strncmp(cLine,"Slide",5)==0){
      int NSlide = 0;
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%d\n",&NSlide);
      if(NSlide != Slide) continue;
      ExecScript(p);
    }
  } 
  free(cLine);
  fclose(File2Read);
}
void DrScript::ExecScript(QPainter *p){
  DrMessage("DrScript.ExecScript");
  char *cLine = (char *)calloc(256,sizeof(char));
  double *From = (double *)calloc(2,sizeof(double));
  int *iFrom = (int *)calloc(2,sizeof(int));
  double *To = (double *)calloc(2,sizeof(double));
  int *iTo = (int *)calloc(2,sizeof(int));
  double *Hue = (double *)calloc(4,sizeof(double));
  int *iHue = (int *)calloc(4,sizeof(int));
  for(int k=0;!(fgets(cLine,256,File2Read)==NULL);k++){
    if(strncmp(cLine,"Slide",5)==0) break;
    else if(strncmp(cLine,"Testo",5)==0){
      //Colore
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      for(int i=0;i<4;i++) iHue[i] = (int)(255.*Hue[i]);
      p->setBrush( QColor(iHue[0],iHue[1],iHue[2],iHue[3]) );
      //Posizione
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf\n",From,From+1);
      iFrom[0] = (int)(width()*From[0]);
      iFrom[1] = (int)(height()*(1.-From[1]));
      //Testo
      char *String = (char *)calloc(512,sizeof(char));
      fgets(cLine,512,File2Read);
      for(int c=0;c<strlen(cLine);c++)
	String[c] = cLine[c];
      p->drawText(iFrom[0],iFrom[1],QString(String),-1,Qt::AlignTop);
      //printf("Col %d %d %d %d Pos %d %d Text %s\n",iHue[0],iHue[1],iHue[2],iHue[3],iFrom[0],iFrom[1],String);
// 	while(1==1) {
// 	fpos_t FPos;
// 	fgetpos(File2Read,&FPos);
// 	if(fgets(cLine,256,File2Read)==0) break;
// 	if(sscanf(cLine,"%lf %lf\n",From,From+1)!=3){
// 	  fsetpos(File2Read,&FPos);
// 	  break;
// 	}
//       }
      free(String);
    }
    else if(strncmp(cLine,"Picture",4)==0){
      //Posizione
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf\n",From,From+1);
      iFrom[0] = (int)(width()*From[0]);
      iFrom[1] = (int)(height()*(1.-From[1]));
      //Posizione
      fgets(cLine,256,File2Read);
      sscanf(cLine,"%lf %lf\n",To,To+1);
      iTo[0] = (int)(width()*To[0]);
      iTo[1] = (int)(height()*(1.-To[1]));
      //Testo
      char *String = (char *)calloc(512,sizeof(char));
      fgets(cLine,512,File2Read);
      //      sprintf(String,"%s",cLine);
      for(int c=0;c<strlen(cLine)-1;c++) String[c] = cLine[c];
      //printf("Apro %s in %d %d -> %d %d\n",String,iFrom[0],iFrom[1],iTo[0],iTo[1]);
      //p->drawImage(QPoint(iFrom[0],iFrom[1]),QImage(QString(String),0),QRect(0,0,iTo[0]-iFrom[0],iTo[1]-iFrom[1]),0);
      p->drawImage(QRect(iFrom[0],iFrom[1],iTo[0],iTo[1]),QImage(QString(String),0));
      free(String);
    }
    else if(strncmp(cLine,"Background",10)==0){
    }
  }
  free(From);
  free(To);
  free(Hue);
  free(cLine);
}
