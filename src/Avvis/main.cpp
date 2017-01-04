#include <qapplication.h>
#include <qpushbutton.h>
#include <qlcdnumber.h>
#include <qfont.h>
#include <qlayout.h>
#include <qstatusbar.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qspinbox.h>
#include <qradiobutton.h>
#include <qbuttongroup.h>
#include <qcheckbox.h>
#include <qslider.h>
#include "ElementiGrafici.h"
#include "VarElementiGrafici.h"

class Finestra: public QWidget{
public:
  Finestra(QWidget *parent=0,const char *name=0);
};
Finestra::Finestra(QWidget *parent,const char *name)
  :QWidget(parent,name){
  //         esci
  QPushButton *esci = new QPushButton("&Esci",this,"Esci");
  esci->setFont( QFont("Times",18,QFont::Bold) );
  connect(esci, SIGNAL(clicked()), qApp,SLOT(quit()) );
  //         Elementi Grafici
  ElementiGrafici *e1 = new ElementiGrafici(this,"ElementiGrafici");
  //         apri
  QPushButton *apri = new QPushButton("&Apri",this,"apri");
  connect(apri, SIGNAL(clicked() ), e1,SLOT(Apri()) );
  //         nome file
  QLineEdit *nomeFile = new QLineEdit(this,"nomeFile");
  connect(nomeFile,SIGNAL( textChanged(const QString &)),e1,SLOT( NomeFile(const QString &)) );
  connect(e1,SIGNAL( TestoCambiato(const QString &)),nomeFile,SLOT( setText(const QString &)) );
  //          segnali
  QPushButton *segnali = new QPushButton("&Segnale",this,"segnali");
  connect(segnali, SIGNAL(clicked() ), e1,SLOT(PuntiSegnale()) );
  //          distr
  QPushButton *distr = new QPushButton("&Distr",this,"distr");
  connect(distr, SIGNAL(clicked() ), e1,SLOT(PuntiDistribuzione()) );
  //          spettro
  QPushButton *spettro = new QPushButton("&Spettro",this,"Spettro");
  connect(spettro, SIGNAL(clicked() ), e1,SLOT(PuntiSpettro()) );
  //          radice
  QPushButton *radice = new QPushButton("&Radice",this,"Radice");
  connect(radice, SIGNAL(clicked() ), e1,SLOT(PuntiRadice()) );
  //           autocor
  QPushButton *autocor = new QPushButton("&ACcor",this,"autocor");
  connect(autocor, SIGNAL(clicked() ), e1,SLOT(PuntiAutocor()) );
  //           integrale
  QPushButton *integrale = new QPushButton("&Int",this,"integrale");
  connect(integrale, SIGNAL(clicked() ), e1,SLOT(PuntiIntegrale()) );
  //           interRett
  QPushButton *interRett = new QPushButton("&InRett",this,"interRett");
  connect(interRett, SIGNAL(clicked() ), e1,SLOT(PuntiInterRett()) );
  //           parabola
  QPushButton *parabola = new QPushButton("&Para",this,"parabola");
  connect(parabola, SIGNAL(clicked() ), e1,SLOT(PuntiParabola()) );
  //           mediaMobile
  QPushButton *mediaMobile = new QPushButton("&mMob",this,"mediaMobile");
  connect(mediaMobile, SIGNAL(clicked() ), e1,SLOT(PuntiMediaMob()) );
  //           correlaADue
  QPushButton *correlaADue = new QPushButton("&Cor2",this,"corrADue");
  connect(correlaADue, SIGNAL(clicked() ), e1,SLOT(PuntiCorrelaADue()) );
  //           salva
  QPushButton *salva = new QPushButton("&Salva",this,"autocor");
  connect(salva, SIGNAL(clicked() ), e1,SLOT(Salva()) );
  //         nomeSalva
  QLineEdit *nomeSalva = new QLineEdit(this,"nomeSalva");
  connect(nomeSalva,SIGNAL( textChanged(const QString &)),e1,SLOT( NomeSalva(const QString &)) );
  connect(e1,SIGNAL( SalvaCambiato(const QString &)),nomeSalva,SLOT( setText(const QString &)) );
  //           ridisegna
  QPushButton *ridisegna = new QPushButton("&R",this,"ridisegna");
  connect(ridisegna,SIGNAL(clicked() ),e1,SLOT(Ridisegna()) );
  //           sulSegnale
  QRadioButton *sulSegnale = new QRadioButton("Sul &segnale",this,"sulSegnale");
  //  sulSegnale->toggle();
  connect(sulSegnale,SIGNAL(toggled(bool) ),e1,SLOT(SulSegnale(bool)) );
  connect(e1,SIGNAL(SegnaleGrafico(bool) ),sulSegnale,SLOT(setChecked(bool)) );
  //           sulGrafico
  QRadioButton *sulGrafico = new QRadioButton("Sul &grafico",this,"sulGrafico");
  connect(sulGrafico,SIGNAL(toggled(bool) ),e1,SLOT(SulGrafico(bool)) );
  //           cosaElabora
  QButtonGroup *cosaElabora = new QButtonGroup("Elabora",this,"cosaElabora");
  cosaElabora->setRadioButtonExclusive( TRUE );
  cosaElabora->insert(sulSegnale,0);
  cosaElabora->insert(sulGrafico,1);
  cosaElabora->setButton(0);
  //           disegna linee
  QCheckBox *linee = new QCheckBox("&Linee",this,"linee");
  linee->toggle();
  connect(linee,SIGNAL(clicked() ),e1,SLOT(DisegnaLinee()) );
  //           disegna punti
  QCheckBox *punti = new QCheckBox("&Punti",this,"punti");
  connect(punti,SIGNAL(clicked() ),e1,SLOT(DisegnaPunti()) );
  //           disegna griglia
  QCheckBox *griglia = new QCheckBox("&Griglia",this,"griglia");
  connect(griglia,SIGNAL(clicked() ),e1,SLOT(DisegnaGriglia()) );
  //           LogLog
  QCheckBox *LogLog = new QCheckBox("&LogLog",this,"LogLog");
  connect(LogLog,SIGNAL(clicked() ),e1,SLOT(DisegnaLogLog()) );
  connect(e1,SIGNAL(LogLogCambiato(bool) ),LogLog,SLOT(setChecked(bool)) );
  //           cosaDisegna
  QButtonGroup *cosaDisegna = new QButtonGroup("Disegna",this,"cosaDisegna");
  cosaDisegna->insert(linee,-1);
  cosaDisegna->insert(punti,-1);
  cosaDisegna->insert(griglia,-1);
  cosaDisegna->insert(LogLog,-1);
  //           Barra
  QStatusBar *Barra = new QStatusBar(this,"Barra");
  connect(e1,SIGNAL(Stato(const QString &) ),Barra,SLOT( message(const QString &)) );
  //          Valori
  QSpinBox *Valori = new QSpinBox(2,100,20,this,"Valori");
  connect(Valori,SIGNAL(valueChanged(int)),e1,SLOT(ImpValori(int)) );
  connect(e1,SIGNAL(ValoriCambiati(int)),Valori,SLOT(setValue(int)) );
  QLabel *etValori = new QLabel(" ",this,"etValori");
  etValori->setAlignment( AlignCenter);
  etValori->setText("Valori");
  //          NMobile
  QSpinBox *NMobile = new QSpinBox(2,100,5,this,"NMobile");
  connect(NMobile,SIGNAL(valueChanged(int)),e1,SLOT(ImpNMobile(int)) );
  connect(e1,SIGNAL(NMobileCambiato(int)),NMobile,SLOT(setValue(int)) );
  QLabel *etNMobile = new QLabel(" ",this,"etNMobile");
  etNMobile->setAlignment( AlignCenter);
  etNMobile->setText("Media Mobile");
  //          NCorrela
  QSpinBox *NCorrela = new QSpinBox(2,100,2,this,"NCorrela");
  connect(NCorrela,SIGNAL(valueChanged(int)),e1,SLOT(ImpNCorrela(int)) );
  connect(e1,SIGNAL(NCorrelaCambiato(int)),NCorrela,SLOT(setValue(int)) );
  QLabel *etNCorrela = new QLabel(" ",this,"etNCorrela");
  etNCorrela->setAlignment( AlignCenter);
  etNCorrela->setText("Correla a due");
  //          NVisMin
  ImpostaValori *NVisMin = new ImpostaValori(this,"NVisMin");
  connect(NVisMin,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNVisMin(int)) );
  connect(e1,SIGNAL(NVisMinCambiato(int)),NVisMin,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterVisCambiato(int)),NVisMin,SLOT(ImpInter(int)) );
  //         NVisMass
  ImpostaValori *NVisMass = new ImpostaValori(this,"NVisMass");
  connect(NVisMass,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNVisMass(int)) );
  connect(e1,SIGNAL(NVisMassCambiato(int)),NVisMass,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterVisCambiato(int)),NVisMass,SLOT(ImpInter(int)) );
  //          NElMin
  ImpostaValori *NElMin = new ImpostaValori(this,"NElMin");
  connect(NElMin,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNElMin(int)) );
  connect(e1,SIGNAL(NElMinCambiato(int)),NElMin,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterVisCambiato(int,int)),NElMin,SLOT(ImpInter(int,int)) );
  //         NElMass
  ImpostaValori *NElMass = new ImpostaValori(this,"NElMass");
  connect(NElMass,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNElMass(int)) );
  connect(e1,SIGNAL(NElMassCambiato(int)),NElMass,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterVisCambiato(int,int)),NElMass,SLOT(ImpInter(int,int)) );

  QGridLayout *grid= new QGridLayout( this,4,2,10);
  grid->addWidget(esci,0,0);
  grid->addWidget(e1,1,1);
  grid->setColStretch(1,10);
  grid->setRowStretch(1,10);

  QVBoxLayout *ASinistra = new QVBoxLayout;
  grid->addLayout(ASinistra,1,0);
  grid->setRowSpacing(1,400);
  ASinistra->addWidget(Valori);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(etValori);
  ASinistra->addWidget(NMobile);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(etNMobile);
  ASinistra->addWidget(NCorrela);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(etNCorrela);
  ASinistra->addWidget(cosaElabora);
  ASinistra->addSpacing(5);
  ASinistra->addWidget(sulSegnale);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(sulGrafico);
  ASinistra->addWidget(apri);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(nomeFile);
  ASinistra->addWidget(salva);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(nomeSalva);
  ASinistra->addWidget(cosaDisegna);
  ASinistra->addSpacing(5);
  ASinistra->addWidget(linee);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(punti);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(griglia);
  ASinistra->addSpacing(-10);
  ASinistra->addWidget(LogLog);
  QVBoxLayout *ASinistraGiu = new QVBoxLayout;
  grid->addLayout(ASinistraGiu,3,0);


  QHBoxLayout *InAlto = new QHBoxLayout;
  grid->addLayout(InAlto,0,1);
  InAlto->addWidget(segnali);
  InAlto->addWidget(distr);
  InAlto->addWidget(spettro);
  InAlto->addWidget(radice);
  InAlto->addWidget(autocor);
  InAlto->addWidget(integrale);
  InAlto->addWidget(mediaMobile);
  InAlto->addWidget(correlaADue);
  InAlto->addWidget(interRett);
  InAlto->addWidget(parabola);

//   QHBoxLayout *MenoInAlto = new QHBoxLayout;
//   grid->addLayout(MenoInAlto,2,1);

  QHBoxLayout *InBasso = new QHBoxLayout;
  grid->addLayout(InBasso,2,1);
  InBasso->addWidget(NElMin);
  InBasso->addWidget(NElMass);
  InBasso->addWidget(ridisegna);
  InBasso->addWidget(NVisMin);
  InBasso->addWidget(NVisMass);

  QHBoxLayout *InBasso2 = new QHBoxLayout;
  grid->addLayout(InBasso2,3,1);
  InBasso2->addWidget(Barra);

}
int main(int argc,char **argv){
  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a(argc,argv);
  
  Finestra f;
  f.setGeometry(0,200,800,355);
  a.setMainWidget( &f );
  f.show();
  
  return a.exec();
}
