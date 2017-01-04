#include "DrFinestra.h"

DrFinestra::DrFinestra(QWidget *parent,const char *name)
  :QWidget(parent,name){
    //         esci
  QPushButton *esci = new QPushButton("&Esci",this,"Esci");
  esci->setFont( QFont("Times",18,QFont::Bold) );
  connect(esci, SIGNAL(clicked()), qApp,SLOT(quit()) );
  e1 = new ElementiGrafici(this);
  //           disegna linee
  QCheckBox *linee = new QCheckBox("&Linee",this,"linee");
  linee->toggle();
  connect(linee,SIGNAL(clicked() ),e1,SLOT(DisegnaLinee()) );
  //           disegna punti
  QCheckBox *punti = new QCheckBox("&Punti",this,"punti");
  connect(punti,SIGNAL(clicked() ),e1,SLOT(DisegnaPunti()) );
  //          distr
  QPushButton *distr = new QPushButton("&Distr",this,"distr");
  connect(distr, SIGNAL(clicked() ), e1,SLOT(PuntiDistribuzione()) );
  //          Valori
  QSpinBox *Valori = new QSpinBox(2,1000,20,this,"Valori");
  connect(Valori,SIGNAL(valueChanged(int)),e1,SLOT(ImpValori(int)) );
  connect(e1,SIGNAL(ValoriCambiati(int)),Valori,SLOT(setValue(int)) );
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
  //           salva
  QPushButton *salva = new QPushButton("&Salva",this,"autocor");
  connect(salva, SIGNAL(clicked() ), e1,SLOT(Salva()) );
  //         Animazione
  a1 = new Animation(this);
  connect(a1, SIGNAL(PuntaCoord(double **,int,int,int) ), e1,SLOT(ApriExt(double**,int,int,int)) );
  //connect(a1, SIGNAL( PuntaCoord(Variabili *v1) ), e1,SLOT(ApriExt(Variabili *v1)) );
//         image
  QPushButton *image = new QPushButton("&Image",this,"image");
  connect(image, SIGNAL(clicked() ), a1,SLOT(Open()) );
  //         punta
  QPushButton *punta = new QPushButton("&Punta",this,"Punta");
  connect(punta, SIGNAL(clicked() ), a1,SLOT(Punta()) );
  //         nomeFile
  QLineEdit *nomeFile = new QLineEdit(this,"nomeFile");
  connect(nomeFile,SIGNAL( textChanged(const QString &)),a1,SLOT( NomeFile(const QString &)) );
  //          filter
  QPushButton *filter = new QPushButton("&Filter",this,"filter");
  connect(filter, SIGNAL(clicked() ), a1,SLOT(Filter()) );
  //          motion
  QPushButton *motion = new QPushButton("&Motion",this,"motion");
  connect(motion, SIGNAL(clicked() ), a1,SLOT(Motion()) );
  //          run
  QPushButton *run = new QPushButton("&run",this,"run");
  connect(run, SIGNAL(clicked() ), a1,SLOT(Run()) );
  //          filterMC
  QPushButton *filterMC = new QPushButton("&MC",this,"MC");
  connect(filterMC, SIGNAL(clicked() ), a1,SLOT(FilterMC()) );
  //          filterCoarseGrain
  QPushButton *filterCoarseGrain = new QPushButton("&CG",this,"CG");
  connect(filterCoarseGrain, SIGNAL(clicked() ), a1,SLOT(FilterCoarseGrain()) );
  connect(filterCoarseGrain, SIGNAL(clicked() ), e1,SLOT(repaint()) );
  //          NGrana
  QSpinBox *NGrana = new QSpinBox(0,5,1,this,"NGrana");
  connect(NGrana,SIGNAL(valueChanged(int)),a1,SLOT(ImpGrana(int)) );
  //          filterIncrease
  QPushButton *filterIncrease = new QPushButton("&Increase",this,"Increase");
  connect(filterIncrease, SIGNAL(clicked() ), a1,SLOT(FilterIncrease()) );
  //          filterContrast
  QPushButton *filterContrast = new QPushButton("&Contrast",this,"Contrast");
  connect(filterContrast, SIGNAL(clicked() ), a1,SLOT(FilterContrast()) );
  connect(filterContrast, SIGNAL(clicked() ), e1,SLOT(repaint()) );
  //          BW
  QPushButton *blackwhite = new QPushButton("&B/W",this,"B/W");
  connect(blackwhite, SIGNAL(clicked() ), a1,SLOT(BlackWhite()) );
  connect(blackwhite, SIGNAL(clicked() ), e1,SLOT(repaint()) );
  //          binary
  QPushButton *binary = new QPushButton("&0/1",this,"0/1");
  connect(binary, SIGNAL(clicked() ), a1,SLOT(Binary()) );
  connect(binary, SIGNAL(clicked() ), e1,SLOT(repaint()) );
  //    Histo
  QPushButton *histo = new QPushButton("&Histo",this,"Histo");
  connect(histo, SIGNAL(clicked() ), a1,SLOT(Histo()) );
  //    NablaPhi
  QPushButton *nablaphi = new QPushButton("&Phi",this,"Phi");
  connect(nablaphi, SIGNAL(clicked() ), a1,SLOT(NablaPhi()) );
  connect(nablaphi, SIGNAL(clicked() ), e1,SLOT(repaint()) );
  //        Script
  s1 = new DrScript();
  //    >
  QPushButton *forward = new QPushButton("&>",this,">");
  connect(forward, SIGNAL(clicked() ), a1,SLOT(IncrSlide()) );
  connect(forward, SIGNAL(clicked() ), s1,SLOT(IncrSlide()) );
  //      Scena
  QGridLayout *grid= new QGridLayout(this,1,3,10);
  //grid->addWidget(a1,1,1);
  grid->setColStretch(1,15);
  grid->setColStretch(0,2);
  grid->setColStretch(2,1);
  grid->setRowStretch(1,5);
  grid->setRowStretch(0,0);
  grid->setRowStretch(2,2);
  grid->setRowMinimumHeight(1,200);
  grid->setRowStretch(2,0);
  grid->setRowStretch(3,0);
  QVBoxLayout *PicHisto = new QVBoxLayout;
  grid->addLayout(PicHisto,1,1);
  PicHisto->addWidget(e1);
  PicHisto->addWidget(a1);
  a1->setMinimumHeight(height());

  QHBoxLayout *Disegna = new QHBoxLayout;
  QVBoxLayout *DisegnaSx = new QVBoxLayout;
  QVBoxLayout *DisegnaDx = new QVBoxLayout;
  //  Disegna->addWidget(cosaDisegna);
  DisegnaSx->addWidget(linee);
  DisegnaSx->addSpacing(-10);
  DisegnaSx->addWidget(punti);
  DisegnaSx->addSpacing(-10);
  DisegnaSx->addWidget(griglia);
  DisegnaSx->addSpacing(-10);
  DisegnaDx->addWidget(riscala);
  DisegnaDx->addSpacing(-10);
  DisegnaDx->addWidget(risTutto);
  DisegnaDx->addSpacing(-10);
  DisegnaDx->addWidget(nset);
  Disegna->addLayout(DisegnaSx);
  Disegna->addLayout(DisegnaDx);

  //QHBoxLayout *Coordinate = new QHBoxLayout;


  QVBoxLayout *ASinistra = new QVBoxLayout;
  grid->addLayout(ASinistra,1,0);
  ASinistra->addLayout(Disegna);
  ASinistra->addWidget(distr);
  ASinistra->addWidget(Valori);
  ASinistra->addWidget(s1);
  ASinistra->addWidget(salva);
  //grid->addWidget(s1,2,0);

  QVBoxLayout *ASinistraGiu = new QVBoxLayout;
  grid->addLayout(ASinistraGiu,2,0);
  
  QHBoxLayout *InAlto = new QHBoxLayout;
  grid->addLayout(InAlto,0,1);
  InAlto->addWidget(histo);
  InAlto->addWidget(blackwhite);
  InAlto->addWidget(binary);
  InAlto->addWidget(nablaphi);
  InAlto->addWidget(forward);
  InAlto->addWidget(punta);

  QVBoxLayout *InAlto2 = new QVBoxLayout;
  grid->addLayout(InAlto2,0,2);

  QVBoxLayout *ADestra = new QVBoxLayout;
  grid->addLayout(ADestra,1,2);
  ADestra->addWidget(esci);
  ADestra->addWidget(image);
  ADestra->addWidget(nomeFile);
  ADestra->addWidget(filter);
  ADestra->addWidget(motion);
  ADestra->addWidget(run);
  ADestra->addWidget(filterMC);
  ADestra->addWidget(NGrana);
  ADestra->addWidget(filterCoarseGrain);
  ADestra->addWidget(filterIncrease);
  ADestra->addWidget(filterContrast);
  ADestra->setSpacing(0);
  
  QVBoxLayout *ADestra1 = new QVBoxLayout;
  grid->addLayout(ADestra1,2,2);

//   QHBoxLayout *MenoInAlto = new QHBoxLayout;
//   grid->addLayout(MenoInAlto,2,1);

  QHBoxLayout *InBasso = new QHBoxLayout;
  grid->addLayout(InBasso,2,1);
}
void DrFinestra::DataFile(char *FileName){
  a1->NomeFile(FileName);
}
void DrFinestra::ConfFile(char *FileName){
  printf("Ciccia %s\n",FileName);
}
void DrFinestra::Run(){
  a1->Run();
}
void DrFinestra::Histo(){
  //a1->Histo();
  //e1->ApriExt(a1->st,256,NLevel,100);
}
