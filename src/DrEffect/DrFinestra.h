#ifndef DRFINESTRA_H
#define DRFINESTRA_H
#include "DrEffect.h"
#include "../Visualizza/ElementiGrafici.h"
#include "DrScript.h"
#include <qpushbutton.h>
#include <qlcdnumber.h>
#include <qfont.h>
#include <qlayout.h>
#include <qstatusbar.h>
#include <qlabel.h>
#include <QLineEdit>
#include <qspinbox.h>
#include <qradiobutton.h>
#include <qbuttongroup.h>
#include <qcheckbox.h>
#include <qslider.h>
#include <QStyle>
#include <QPlastiqueStyle>

class DrFinestra: public QWidget{
 private:
    //         Elementi Grafici
  Animation *a1;
  ElementiGrafici *e1;
  DrScript *s1;
 public:
  DrFinestra(QWidget *parent=0,const char *name=0);
  void Histo();
  void DataFile(char *FileName);
  void ConfFile(char *FileName);
  void Run();
  //  char *FileDaScrivere;
  //  void ScriviNomeFile(char *);
};
#endif
