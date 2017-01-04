/***********************************************************************
Visualizza: This program interface every function of the ElementiGrafici class with the Qt3 widget for calling the plotting function of ElementiGrafici.
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
#include "ElementiGrafici.h"
#include <qapplication.h>


int main(int argc,char **argv){
  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a(argc,argv);
  Finestra f;
  f.setGeometry(100,200,800,355);
  a.setMainWidget( &f );
  f.show();
  int *FileList = (int *)calloc(argc,sizeof(int));
  int NFile = 0;
  for(int i=1;i<argc;i++){
    if(argv[i][0] != '-'){
      FileList[NFile++] = i;
    }
    if(!strcmp(argv[i],"-c") && i<argc-1){
      f.ConfFile(argv[i+1]);
      i++;
    }
  }
  f.DataFile(argv,FileList,NFile);
  //  f.setStyle(QStyle::QPalstiqueStyle);
  
  return a.exec();
}
