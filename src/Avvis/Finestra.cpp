#include "ElementiGrafici.h"
#include <qapplication.h>
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
#include <Q3VBox>
#include <Q3ButtonGroup>


Finestra::Finestra(QWidget *parent,const char *name)
  :QWidget(parent,name){
    //         esci
  QPushButton *esci = new QPushButton("&Esci",this,"Esci");
  esci->setFont( QFont("Times",18,QFont::Bold) );
  connect(esci, SIGNAL(clicked()), qApp,SLOT(quit()) );
  //         Elementi Grafici
  e1 = new ElementiGrafici(this,"ElementiGrafici");
  //         apri
  QPushButton *apri = new QPushButton("&Apri",this,"apri");
  connect(apri, SIGNAL(clicked() ), e1,SLOT(Apri()) );
  //         Aggiungi
  QPushButton *aggiungi = new QPushButton("&Aggiungi",this,"aggiungi");
  connect(aggiungi, SIGNAL(clicked() ), e1,SLOT(Aggiungi()) );
  //         nome file
  QLineEdit *nomeFile = new QLineEdit(this,"nomeFile");
  connect(nomeFile,SIGNAL( textChanged(const QString &)),e1,SLOT( NomeFile(const QString &)) );
  connect(e1,SIGNAL( TestoCambiato(const QString &)),nomeFile,SLOT( setText(const QString &)) );
  //         nome conf
  QLineEdit *nomeConf = new QLineEdit(this,"nomeConf");
  connect(nomeConf,SIGNAL( textChanged(const QString &)),e1,SLOT( NomeConf(const QString &)) );
  connect(e1,SIGNAL( ConfCambiato(const QString &)),nomeConf,SLOT( setText(const QString &)) );
  //         nomeSalva
  QLineEdit *nomeSalva = new QLineEdit(this,"nomeSalva");
  connect(nomeSalva,SIGNAL( textChanged(const QString &)),e1,SLOT( NomeSalva(const QString &)) );
  connect(e1,SIGNAL( SalvaCambiato(const QString &)),nomeSalva,SLOT( setText(const QString &)) );
  //         nomeTit
  QLineEdit *nomeTit = new QLineEdit(this,"nomeTit");
  connect(nomeTit,SIGNAL( textChanged(const QString &)),e1,SLOT( NomeTit(const QString &)) );
  connect(e1,SIGNAL( TitoloCambiato(const QString &)),nomeSalva,SLOT( setText(const QString &)) );
  //         nomeEtX
  QLineEdit *nomeEtX = new QLineEdit(this,"nomeEtX");
  connect(nomeEtX,SIGNAL( textChanged(const QString &)),e1,SLOT( NomeEtX(const QString &)) );
  connect(e1,SIGNAL( EtichettaXCambiato(const QString &)),nomeEtX,SLOT( setText(const QString &)) );
  nomeEtX->setMaximumWidth(50);
  nomeEtX->hide();
  //         nomeEtY
  QLineEdit *nomeEtY = new QLineEdit(this,"nomeEtY");
  connect(nomeEtY,SIGNAL( textChanged(const QString &)),e1,SLOT( NomeEtY(const QString &)) );
  connect(e1,SIGNAL( EtichettaYCambiato(const QString &)),nomeEtY,SLOT( setText(const QString &)) );
  nomeEtY->setMaximumWidth(50);
  nomeEtY->hide();
  QHBoxLayout *Etichette = new QHBoxLayout;
  Etichette->addWidget(nomeEtX);
  Etichette->addWidget(nomeEtY);
  Etichette->addStretch(1);
  //          segnali
  QPushButton *segnali = new QPushButton("&Segnale",this,"segnali");
  connect(segnali, SIGNAL(clicked() ), e1,SLOT(PuntiSegnale()) );
  //          distr
  QPushButton *distr = new QPushButton("&Distr",this,"distr");
  connect(distr, SIGNAL(clicked() ), e1,SLOT(PuntiDistribuzione()) );
  //          spettro
  QPushButton *spettro = new QPushButton("&Spettro",this,"Spettro");
  connect(spettro, SIGNAL(clicked() ), e1,SLOT(PuntiSpettro()) );
  //          autosim
  QPushButton *autosim = new QPushButton("&ASim",this,"autosim");
  connect(autosim, SIGNAL(clicked() ), e1,SLOT(PuntiAutosimilarita()) );
  //  autosim->hide();
  //           autocor
  QPushButton *autocor = new QPushButton("&ACcor",this,"autocor");
  connect(autocor, SIGNAL(clicked() ), e1,SLOT(PuntiAutocor()) );
  //           normalizza
  QPushButton *normalizza = new QPushButton("&Norm",this,"normalizza");
  connect(normalizza, SIGNAL(clicked() ), e1,SLOT(PuntiNormalizza()) );
  //           modulo
  QPushButton *modulo = new QPushButton("&Mod",this,"modulo");
  connect(modulo, SIGNAL(clicked() ), e1,SLOT(PuntiModulo()) );
  //           integrale
  QPushButton *integrale = new QPushButton("&Int",this,"integrale");
  connect(integrale, SIGNAL(clicked() ), e1,SLOT(PuntiIntegrale()) );
  //           derivata
  QPushButton *derivata = new QPushButton("&Der",this,"derivata");
  connect(derivata, SIGNAL(clicked() ), e1,SLOT(PuntiDerivata()) );
  //           varie
  QPushButton *varie = new QPushButton("&Varie",this,"varie");
  connect(varie, SIGNAL(clicked() ), e1,SLOT(PuntiVarie()) );
  varie->hide();
  //           sum
  QPushButton *sum = new QPushButton("&Sum",this,"sum");
  connect(sum, SIGNAL(clicked() ), e1,SLOT(PuntiSum()) );
  //           interRett
  QPushButton *interRett = new QPushButton("&InRett",this,"interRett");
  connect(interRett, SIGNAL(clicked() ), e1,SLOT(PuntiInterRett()) );
  //           interExp
  QPushButton *interExp = new QPushButton("&InExp",this,"interExp");
  connect(interExp, SIGNAL(clicked() ), e1,SLOT(PuntiInterExp()) );
  //           interGauss
  QPushButton *interGauss = new QPushButton("&InGauss",this,"interGauss");
  connect(interGauss, SIGNAL(clicked() ), e1,SLOT(PuntiInterGauss()) );
  //           parabola
  QPushButton *parabola = new QPushButton("&Para",this,"parabola");
  connect(parabola, SIGNAL(clicked() ), e1,SLOT(PuntiParabola()) );
  //parabola->hide();
  //           mediaMobile
  QPushButton *mediaMobile = new QPushButton("&mMob",this,"mediaMobile");
  connect(mediaMobile, SIGNAL(clicked() ), e1,SLOT(PuntiMediaMob()) );
  //           correlaADue
  QPushButton *correlaADue = new QPushButton("&Cor2",this,"corrADue");
  connect(correlaADue, SIGNAL(clicked() ), e1,SLOT(PuntiCorrelaADue()) );
  //           Stampa
  QPushButton *stampa = new QPushButton("&Stampa",this,"stampa");
  connect(stampa, SIGNAL(clicked() ), e1,SLOT(StampaFile()) ); 
  //           salva
  QPushButton *salva = new QPushButton("&Salva",this,"autocor");
  connect(salva, SIGNAL(clicked() ), e1,SLOT(Salva()) );
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
  Q3ButtonGroup *cosaElabora = new Q3ButtonGroup("Elabora",this,"cosaElabora");
  cosaElabora->setRadioButtonExclusive( TRUE );
  cosaElabora->hide();
  cosaElabora->insert(sulSegnale,0);
  cosaElabora->insert(sulGrafico,1);
  cosaElabora->setButton(0);
  //           suX
  QRadioButton *suX = new QRadioButton("X",this,"X");
  connect(suX,SIGNAL(toggled(bool) ),e1,SLOT(SuX(bool)) );
  //           suY
  QRadioButton *suY = new QRadioButton("Y",this,"Y");
  connect(suY,SIGNAL(toggled(bool) ),e1,SLOT(SuY(bool)) );
  //           suDX
  QRadioButton *suDX = new QRadioButton("DX",this,"DX");
  connect(suDX,SIGNAL(toggled(bool) ),e1,SLOT(SuDX(bool)) );
  //           suDY
  QRadioButton *suDY = new QRadioButton("DY",this,"DY");
  connect(suDY,SIGNAL(toggled(bool) ),e1,SLOT(SuDY(bool)) );
  Q3ButtonGroup *qualeCoord = new Q3ButtonGroup(this);
  qualeCoord->setRadioButtonExclusive( TRUE );
  qualeCoord->hide();
  qualeCoord->insert(suX,0);
  qualeCoord->insert(suY,1);
  qualeCoord->insert(suDX,2);
  qualeCoord->insert(suDY,3);
  qualeCoord->setButton(1);
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
  //           seRiscala
  QCheckBox *riscala = new QCheckBox("&Ris",this,"riscala");
  connect(riscala,SIGNAL(clicked() ),e1,SLOT(SeRiscala()) );
  //           seRisTutto
  QCheckBox *risTutto = new QCheckBox("&Tutti",this,"risTutto");
  connect(risTutto,SIGNAL(clicked() ),e1,SLOT(SeRiscalaTutto()) );
  //           seSet
  QCheckBox *nset = new QCheckBox("&Set",this,"nset");
  connect(nset,SIGNAL(clicked() ),e1,SLOT(NSet()) );
  //           Logx
  QCheckBox *Logx = new QCheckBox("&Logx",this,"Logx");
  connect(Logx,SIGNAL(clicked() ),e1,SLOT(DisegnaLogx()) );
  connect(e1,SIGNAL(LogxCambiato(bool) ),Logx,SLOT(setChecked(bool)) );
  //           Logy
  QCheckBox *Logy = new QCheckBox("&Logy",this,"Logy");
  connect(Logy,SIGNAL(clicked() ),e1,SLOT(DisegnaLogy()) );
  connect(e1,SIGNAL(LogyCambiato(bool) ),Logy,SLOT(setChecked(bool)) );
  //           cosaDisegna
  Q3ButtonGroup *cosaDisegna = new Q3ButtonGroup("cosaDisegna",this,"cosaDisegna");
  cosaDisegna->hide();
  cosaDisegna->insert(linee,-1);
  cosaDisegna->insert(punti,-1);
  cosaDisegna->insert(griglia,-1);
  cosaDisegna->insert(Logx,-1);
  cosaDisegna->insert(Logy,-1);
  cosaDisegna->insert(riscala,-1);
  cosaDisegna->insert(risTutto,-1);
  cosaDisegna->insert(nset,-1);
  //           Barra
  QStatusBar *Barra = new QStatusBar(this,"Barra");
  connect(e1,SIGNAL(Stato(const QString &) ),Barra,SLOT( message(const QString &)) );
  //           Barra2
  QStatusBar *Barra2 = new QStatusBar(this,"Barra");
  connect(e1,SIGNAL(StatoSequenza(const QString &) ),Barra2,SLOT( message(const QString &)) );
  //          NBin
  QSpinBox *NBin = new QSpinBox(2,1000,20,this,"NBin");
  connect(NBin,SIGNAL(valueChanged(int)),e1,SLOT(ImpNBin(int)) );
  connect(e1,SIGNAL(NBinCambiati(int)),NBin,SLOT(setValue(int)) );
  QLabel *etNBin = new QLabel(" ",this,"etNBin");
  etNBin->setAlignment( Qt::AlignHCenter);
  etNBin->setText("NBin");
  etNBin->hide();
  //          NMobile
  QSpinBox *NMobile = new QSpinBox(2,1000,5,this,"NMobile");
  connect(NMobile,SIGNAL(valueChanged(int)),e1,SLOT(ImpNMobile(int)) );
  connect(e1,SIGNAL(NMobileCambiato(int)),NMobile,SLOT(setValue(int)) );
  QLabel *etNMobile = new QLabel(" ",this,"etNMobile");
  etNMobile->setAlignment( Qt::AlignHCenter);
  etNMobile->setText("Media Mobile");
  etNMobile->hide();
  //          NCorrela
  QSpinBox *NCorrela = new QSpinBox(2,100,2,this,"NCorrela");
  connect(NCorrela,SIGNAL(valueChanged(int)),e1,SLOT(ImpNCorrela(int)) );
  connect(e1,SIGNAL(NCorrelaCambiato(int)),NCorrela,SLOT(setValue(int)) );
  QLabel *etNCorrela = new QLabel(" ",this,"etNCorrela");
  etNCorrela->setAlignment( Qt::AlignHCenter);
  etNCorrela->setText("Correla a due");
  etNCorrela->hide();
  //          NVisMin
  ImpostaNBin *NVisMin = new ImpostaNBin(this,"NVisMin");
  connect(NVisMin,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNVisMin(int)),Qt::DirectConnection);
  connect(e1,SIGNAL(NVisMinCambiato(int)),NVisMin,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterVisCambiato(int)),NVisMin,SLOT(ImpInter(int)) );
  //         NVisMax
  ImpostaNBin *NVisMax = new ImpostaNBin(this,"NVisMax");
  connect(NVisMax,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNVisMax(int)) );
  connect(e1,SIGNAL(NVisMaxCambiato(int)),NVisMax,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterVisCambiato(int)),NVisMax,SLOT(ImpInter(int)) );
  //          NElMin
  ImpostaNBin *NElMin = new ImpostaNBin(this,"NElMin");
  connect(NElMin,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNElMin(int)) );
  connect(e1,SIGNAL(NElMinCambiato(int)),NElMin,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterElCambiato(int)),NElMin,SLOT(ImpInter(int)) );
  //         NElMax
  ImpostaNBin *NElMax = new ImpostaNBin(this,"NElMax");
  connect(NElMax,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNElMax(int)) );
  connect(e1,SIGNAL(NElMaxCambiato(int)),NElMax,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterElCambiato(int)),NElMax,SLOT(ImpInter(int)) );
  //          NElMinY
  ImpostaNBin *NElMinY = new ImpostaNBin(this,"NElMinY");
  connect(NElMinY,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNElMinY(int)) );
  connect(e1,SIGNAL(NElMinYCambiato(int)),NElMinY,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterYCambiato(int,int)),NElMinY,SLOT(ImpInter(int,int)) );
  //         NElMaxY
  ImpostaNBin *NElMaxY = new ImpostaNBin(this,"NElMaxY");
  connect(NElMaxY,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpNElMaxY(int)) );
  connect(e1,SIGNAL(NElMaxYCambiato(int)),NElMaxY,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterYCambiato(int,int)),NElMaxY,SLOT(ImpInter(int,int)) );
  //         CoordX
  ImpostaNBin *CoordX = new ImpostaNBin(this,"CoordX");
  connect(CoordX,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpCoordX(int)) );
  connect(e1,SIGNAL(CoordXCambiata(int)),CoordX,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterCoordXCambiato(int,int)),CoordX,SLOT(ImpInter(int,int)) );
  //         CoordY
  ImpostaNBin *CoordY = new ImpostaNBin(this,"CoordY");
  connect(CoordY,SIGNAL(ValoreCambiato(int)),e1,SLOT(ImpCoordY(int)) );
  connect(e1,SIGNAL(CoordYCambiata(int)),CoordY,SLOT(ImpNumero(int)) );
  connect(e1,SIGNAL(InterCoordYCambiato(int,int)),CoordY,SLOT(ImpInter(int,int)) );

  QGridLayout *grid= new QGridLayout( this,4,4,10);
  grid->addWidget(esci,0,0);
  grid->addWidget(e1,1,1);
  grid->setColStretch(1,5);
  grid->setRowStretch(1,5);
  grid->setColumnMinimumWidth(1,500);
  grid->setRowMinimumHeight(1,500);


  QHBoxLayout *Coordinate = new QHBoxLayout;
  Coordinate->addWidget(qualeCoord);
  Coordinate->addSpacing(5);
  Coordinate->addWidget(suX);
  Coordinate->addSpacing(-10);
  Coordinate->addWidget(suY);
  Coordinate->addSpacing(-10);
  Coordinate->addWidget(suDX);
  Coordinate->addSpacing(-10);
  Coordinate->addWidget(suDY);

  QHBoxLayout *Disegna = new QHBoxLayout;
  QVBoxLayout *DisegnaSx = new QVBoxLayout;
  QVBoxLayout *DisegnaDx = new QVBoxLayout;
  //  Disegna->addWidget(cosaDisegna);
  DisegnaSx->addWidget(linee);
  DisegnaSx->addSpacing(-10);
  DisegnaSx->addWidget(punti);
  DisegnaSx->addSpacing(-10);
  DisegnaSx->addWidget(Logx);
  DisegnaSx->addSpacing(-10);
  DisegnaSx->addWidget(Logy);
  DisegnaSx->addSpacing(-10);
  DisegnaDx->addWidget(riscala);
  DisegnaDx->addSpacing(-10);
  DisegnaDx->addWidget(risTutto);
  DisegnaDx->addSpacing(-10);
  DisegnaDx->addWidget(nset);
  DisegnaDx->addSpacing(-10);
  DisegnaDx->addWidget(griglia);
  Disegna->addLayout(DisegnaSx);
  Disegna->addLayout(DisegnaDx);

  QVBoxLayout *ASinistra = new QVBoxLayout;
  grid->addLayout(ASinistra,1,0);
  //grid->setRowSpacing(1,400);
  ASinistra->addWidget(cosaElabora);
  ASinistra->addSpacing(5);
  ASinistra->addWidget(sulSegnale);
  ASinistra->addWidget(sulGrafico);
  ASinistra->addLayout(Coordinate);
  ASinistra->addWidget(apri);
  ASinistra->addWidget(nomeFile);
  ASinistra->addWidget(aggiungi);
  ASinistra->addWidget(nomeConf);
  //ASinistra->addWidget(nomeTit);
  nomeTit->hide();
  //  ASinistra->insertLayout(17,Etichette,-10);
  ASinistra->addWidget(salva);
  ASinistra->addWidget(nomeSalva);
  ASinistra->addWidget(stampa); 
  ASinistra->addLayout(Disegna);

  QVBoxLayout *ASinistraGiu = new QVBoxLayout;
  grid->addLayout(ASinistraGiu,3,0);
  ASinistraGiu->addWidget(Barra2);
  

  QHBoxLayout *InAlto = new QHBoxLayout;
  grid->addLayout(InAlto,0,1);
  InAlto->addWidget(NElMinY);
  InAlto->addWidget(NElMaxY);
  InAlto->addWidget(segnali);
  InAlto->addWidget(CoordX);
  InAlto->addWidget(CoordY);

  QVBoxLayout *InAlto2 = new QVBoxLayout;
  grid->addLayout(InAlto2,0,2);
  InAlto2->addWidget(interRett);
  InAlto2->addSpacing(-10);
  InAlto2->addWidget(parabola);


  QVBoxLayout *ADestra = new QVBoxLayout;
  grid->addLayout(ADestra,1,2);
  ADestra->addWidget(NBin);
  ADestra->addSpacing(-15);
  ADestra->addWidget(distr);
  ADestra->addWidget(NMobile);
  ADestra->addSpacing(-15);
  ADestra->addWidget(mediaMobile);
  ADestra->addWidget(NCorrela);
  ADestra->addSpacing(-15);
  ADestra->addWidget(correlaADue);
  ADestra->addSpacing(-10);
  ADestra->addWidget(spettro);
  ADestra->addSpacing(-10);
  ADestra->addWidget(integrale);
  ADestra->addSpacing(-10);
  ADestra->addWidget(derivata);
  ADestra->addSpacing(-10);
  ADestra->addWidget(sum);
  ADestra->addSpacing(-10);
  ADestra->addWidget(autosim);
  ADestra->addSpacing(-10);
  ADestra->addWidget(normalizza);
  ADestra->addSpacing(-10);
  ADestra->addWidget(modulo);
  ADestra->addSpacing(-10);
  ADestra->addWidget(autocor);
  ADestra->addSpacing(-10);
  ADestra->addWidget(interExp);
  ADestra->addSpacing(-10);
  ADestra->addWidget(interGauss);

  QVBoxLayout *ADestra1 = new QVBoxLayout;
  grid->addLayout(ADestra1,1,4);

//   QHBoxLayout *MenoInAlto = new QHBoxLayout;
//   grid->addLayout(MenoInAlto,2,1);

  QHBoxLayout *InBasso = new QHBoxLayout;
  grid->addLayout(InBasso,2,1);
  InBasso->addWidget(NVisMin);
  InBasso->addWidget(NVisMax);
  InBasso->addWidget(ridisegna);
  InBasso->addWidget(NElMin);
  InBasso->addWidget(NElMax);

  QHBoxLayout *InBasso2 = new QHBoxLayout;
  grid->addLayout(InBasso2,3,1);
  InBasso2->addWidget(Barra);
}

void Finestra::DataFile(char **argv,int *FileList,int NFile){
  e1->Apri(argv,FileList,NFile);
}
void Finestra::ConfFile(char *FileName){
  e1->ChooseConfFile(FileName);
}
