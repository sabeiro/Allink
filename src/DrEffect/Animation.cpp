#include "DrEffect.h"
#ifdef __glut_h__
DrEffect *DrE;
void Animation::DrMessage(const char * s, ...)
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
void ParticleList(){
  //Pol->RenderPart();
}
void ParticleRealTime(){
  return;
}
void reshape(int w,int h){
  DrE->Dreshape(w,h);
}
void Timer(int v){
  DrE->DTimer(v);
}
void MouseMove(int x,int y){
  DrE->DMouseMove(x,y);
}
void mouse(int button, int state,int x,int y){
  DrE->Dmouse(button,state,x,y);
};
void special(int k, int x, int y){
  DrE->Dspecial(k,x,y);
}
void Slide(){}
void DrawParticles(){}
void Particle(){}
void keyboard(unsigned char key,int x, int y){
  DrE->keyboardDraw(key);
}
void Menu(){}
void Figure(){
  //Dr->Draw1();
  DrE->ShowImage();
}
void Figure1(){
  DrE->Draw1();
}
Animation::Animation(QWidget *parent): QGLWidget(parent){
  FileName = new char[160];
  DrE = new DrEffect();
  st = DrE->Hist;
  sprintf(FileName,"RadDistrNanoR3.8H0.5.tif");
  xRot = 0;
  yRot = 0;
  zRot = 0;
  Grana = 0;
  NVar=0;
  NMass=0;
  Valori = 20;
  v1 = new VarDatFile(NMass,NVar,Valori);
}
void Animation::initializeGL(){
  DrMessage("Animation.InitialiezeGL");
  DrE->Window();return;  //Dr->OpenImage(FileName);
  //qglClearColor(trolltechPurple.dark());
  //object = makeObject();
  glClearColor(.0,.0,.0,.0);//colore sfondo
  glShadeModel(GL_SMOOTH);// Enables Smooth Shading
  glClearDepth(1.0f);// Depth Buffer Setup
  glEnable(GL_DEPTH_TEST);//Controllo di profondit`a
  //  glDepthMask(GL_FALSE);
  //  glDepthRange(0,-10.);
  glDepthFunc(GL_LEQUAL);// The Type Of Depth Test To Do
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glEnable(GL_CULL_FACE);
}
void Animation::NomeFile(const QString &S){
  DrMessage("Animation.NomeFile");
  sprintf(FileName,"%s",S.ascii());
}
void Animation::NomeFile(char *ExtName){
  DrMessage("Animation.NomeFile");
  sprintf(FileName,"%s",ExtName);
  Open();
}
void Animation::Open(){
  DrMessage("Animation.Open");
  DrE->OpenImage(FileName);
  paintGL();
  repaint();
}
void Animation::Filter(){
  DrMessage("Animation.Filter");
  //DrE->Edges();
  DrE->EffectFilter();
  paintGL();
  repaint();
}
void Animation::Motion(){
  DrMessage("Animation.Motion");
  DrE->EffectMotion();
  paintGL();
  repaint();
}
void Animation::FilterContrast(){
  DrMessage("Animation.FilterContrast");
  DrE->Contrast();
  paintGL();
  repaint();
}
void Animation::FilterIncrease(){
  DrMessage("Animation.FilterIncrease");
  DrE->EffectIncrease();
  paintGL();
  repaint();
}
void Animation::FilterMC(){
  DrMessage("Animation.FilterMC");
  DrE->EffectMC();
  paintGL();
  repaint();
}
void Animation::IncrSlide(){
  DrMessage("Animation.IncrSlide");
  Slide++;
  paintGL();
  repaint();
}
void Animation::FilterCoarseGrain(){
  DrMessage("Animation.FilterCoarseGrain");
  DrE->EffectCoarseGrain(Grana);
  paintGL();
  repaint();
}
void Animation::Run(){
  DrMessage("Animation.Run");
  DrE->Run();
}
void Animation::ImpGrana(int ExtGrana){
  if(Grana < 0 || Grana > 32) return;
  Grana = ExtGrana;
}
void Animation::resizeGL(int width,int height){
  DrMessage("Animation.resizeGL");
  DrE->Dreshape(width,height);return;
  int side = qMin(width, height);
  glViewport((width - side) / 2, (height - side) / 2, side, side);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60,(GLfloat) width/(GLfloat) height,.1,200.);//ang,rapp,zmin,zmax  
  //glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0);
  glMatrixMode(GL_MODELVIEW);
}
void Animation::paintGL(){
  DrMessage("Animation.paintGL");
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  //glPushMatrix();
  glTranslatef(0.,0.,-.6);
  glRotated(xRot / 16.0, 1.0, 0.0, 0.0);
  glRotated(yRot / 16.0, 0.0, 1.0, 0.0);
  glRotated(zRot / 16.0, 0.0, 0.0, 1.0);
  glColor4f(.0,1.0,.0,1.);
  DrE->ShowImage();
  //  glutWireCube((GLfloat).5);
  //  Dr->Draw1();
  //glutSolidSphere(.1,20,20);
  //glPopMatrix();
}
void Animation::paintEvent(QPaintEvent (event)){
  DrMessage("Animation.paintEvent");
  QPainter p(this);
  //paintGL();
} 
void Animation::mousePressEvent(QMouseEvent *event)
{
  lastPos = event->pos();
}
void Animation::mouseMoveEvent(QMouseEvent *event){
  //  Dr->DMouseMove(event->x(),event->y());return;
  int dx = event->x() - lastPos.x();
  int dy = event->y() - lastPos.y();
  
  if (event->buttons() & Qt::LeftButton) {
    setXRotation(xRot + 8 * dy);
    setYRotation(yRot + 8 * dx);
  } else if (event->buttons() & Qt::RightButton) {
    setXRotation(xRot + 8 * dy);
    setZRotation(zRot + 8 * dx);
  }
  lastPos = event->pos();
}
void Animation::setXRotation(int angle)
{
  //normalizeAngle(&angle);
  if (angle != xRot) {
    xRot = angle;
    emit xRotationChanged(angle);
    updateGL();
  }
}
void Animation::setYRotation(int angle)
{
  //normalizeAngle(&angle);
  if (angle != yRot) {
    yRot = angle;
    emit yRotationChanged(angle);
    updateGL();
  }
}
void Animation::setZRotation(int angle)
{
  //normalizeAngle(&angle);
  if (angle != zRot) {
    zRot = angle;
    emit zRotationChanged(angle);
    updateGL();
  }
}
void Animation::Histo(){
  DrMessage("Animation.Histo");
  DrE->Histo();
  st = DrE->Hist;
  NMass = DrE->NChar;
  NVar = DrE->NLevel;
  Valori = 20;
  return ;
  delete v1;
  v1 = new VarDatFile(st,NMass,NVar,Valori);
  v1->ScriviTutto("Histo.dat",0,1.,1.);
  //Punta();
}
void Animation::NablaPhi(){
  DrMessage("Animation.NablaPhi");
  DrE->NablaPhi();
  st = DrE->Phi;
  NMass = DrE->BuffSize();
  NVar = 2;
  Valori = 100;
  delete v1;
  v1 = new VarDatFile(st,NMass,NVar,Valori);
  v1->ScriviTutto("NablaPhi.dat",0,1.,1.);
  v1->ScriviTutto("NablaPhiLogLog.dat",0,1.,1.);
  //Punta();
}
void Animation::Punta(){
  delete v1;
  v1 = new VarDatFile(st,NMass,NVar,Valori);
  DrMessage("Animation.Punta");
  emit PuntaCoord(st,NMass,NVar,Valori);
  //emit PuntaCoord(v1);
}
void Animation::BlackWhite(){
  DrMessage("Animation.BlackWhite");
  DrE->BlackWhite();
  repaint();
}
void Animation::Binary(){
  DrMessage("Animation.Binary");
  DrE->Binary(0);
  repaint();
}

#endif// __glut_h__
