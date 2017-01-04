#include "../include/Matematica.h"
Vettore::Vettore(int N){
  x = NULL;
  if(N<0)return;
  NDim = N;
  //x = (double *)calloc(NDim,sizeof(double));
  x = new double[NDim];
}
Vettore::Vettore(){
  x = NULL;
  NDim = 3;
  //x = (double *)calloc(NDim,sizeof(double));
  x = new double[NDim];
}
Vettore::Vettore(double *Pos,int N){
  x = NULL;
  NDim = N;
  x = new double[NDim];
    //  x = (double *)calloc(NDim,sizeof(double));
  for(int d=0;d<3;d++)
    x[d] = Pos[d];
};
Vettore::Vettore(double xx,double yy){
  x = NULL;
  NDim = 2;
  x = new double[NDim];
  //  x = (double *)calloc(NDim,sizeof(double));
  x[0] = xx;
  x[1] = yy;
}
Vettore::Vettore(double xx,double yy,double zz){
  x = NULL;
  NDim = 3;
  x = new double[NDim];
  //  x = (double *)calloc(NDim,sizeof(double));
  x[0] = xx;
  x[1] = yy;
  x[2] = zz;
}
Vettore::~Vettore(){
  if(x) delete [] x;
  //if(x) free(x);
}
Vettore Vettore::operator+(const Vettore &u){
  assert(NDim == u.NDim);
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    Resp.x[d] = this->x[d] + u.x[d];
  }
  return Resp;
}
Vettore Vettore::operator+(const Vettore &u) const{
  assert(NDim == u.NDim);
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    Resp.x[d] = this->x[d] + u.x[d];
  }
  return Resp;
}
Vettore& Vettore::operator+=(const Vettore &u){
  for(int d=0;d<this->NDim;d++){
    this->x[d] += u.x[d];
  }
  return *this;
}
Vettore Vettore::operator-(const Vettore &u){
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    Resp.x[d] = this->x[d] - u.x[d];
  }
  return Resp;
}
Vettore Vettore::operator-(const Vettore &u) const{
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    Resp.x[d] = this->x[d] - u.x[d];
  }
  return Resp;
}
Vettore& Vettore::operator-=(const Vettore &u){
  for(int d=0;d<this->NDim;d++){
    this->x[d] -= u.x[d];
  }
  return *this;
}
Vettore Vettore::operator*(const Vettore &u){
  assert(NDim == u.NDim);
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    Resp.x[d] = this->x[d] * u.x[d];
  }
  return Resp;
}
Vettore& Vettore::operator*=(const Vettore &u){
  for(int d=0;d<this->NDim;d++){
    this->x[d] *= u.x[d];
  }
  return *this;
}
Vettore Vettore::operator*(const double Fact){
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    Resp.x[d] = this->x[d]*Fact;
  }
  return Resp;
}
Vettore& Vettore::operator*=(const double Fact){
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    this->x[d] *= Fact;
  }
  return *this;
}
//inline 
Vettore operator*(double Fact, const Vettore& vec) {
  Vettore Resp(vec.NDim);
  for(int d=0;d<vec.NDim;d++){
    Resp.x[d] = vec.x[d]*Fact;
  }
  return Resp;
}
Vettore Vettore::operator/(const Vettore &u){
  assert(NDim == u.NDim);
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    Resp.x[d] = this->x[d] / u.x[d];
  }
  return Resp;
}
Vettore& Vettore::operator/=(const Vettore &u){
  for(int d=0;d<this->NDim;d++){
    this->x[d] /= u.x[d];
  }
  return *this;
}
double Vettore::operator%(const Vettore &u){
  assert(NDim == u.NDim);
  double Resp = 0.;
  for(int d=0;d<this->NDim;d++){
    Resp += this->x[d] * u.x[d];
  }
  return Resp;
}
Vettore Vettore::operator=(const Vettore &u){
  Vettore Resp(this->NDim);
  for(int d=0;d<this->NDim;d++){
    Resp.x[d] = u.x[d];
  }
  return Resp;
}
Vettore Vettore::operator^(const Vettore& u){
  assert(NDim == u.NDim);
  Vettore Resp(NDim);
  for(int d=0;d<NDim;d++){
    int NUno = (d+1)%NDim;
    int NDue = (d+2)%NDim;
    Resp.x[d] = u.x[NUno] * x[NDue] - u.x[NDue]*x[NUno];
  }
  return Resp;
}
Vettore& Vettore::operator^=(const Vettore& u){
  for(int d=0;d<NDim;d++){
    int NUno = (d+1)%NDim;
    int NDue = (d+2)%NDim;
    this->x[d] = this->x[NUno] * u.x[NDue] - this->x[NDue]*u.x[NUno];
  }
  return *this;
}
void Vettore::Print(){
  for(int d=0;d<NDim;d++)
    printf("%d) %lf\n",d,x[d]);
}
double Vettore::Abs(){
  double Resp = 0.;
  for(int d=0;d<NDim;d++)
    Resp += x[d];
  return ABS(Resp);
}
double Vettore::Norm(){
  double Resp=0.;
  for(int d=0;d<NDim;d++)
    Resp += QUAD((x[d]));
  return sqrt(Resp);
}
double Vettore::Normalize(){
  double Div = Norm();
  if(Div > 0.0)
    for(int d=0;d<NDim;d++)
      x[d] /= Div;
  return Div;
}
void Vettore::Subs(const Vettore *u,const Vettore *v){
  for(int d=0;d<NDim;d++)
    x[d] = u->x[d] - v->x[d];
}
void Vettore::Mult(double Fact){
  for(int d=0;d<NDim;d++)
    x[d] = Fact*x[d];
}
double Vettore::ScalS(const Vettore *u,const Vettore *v){
#ifdef VECT_DEBUG
  if(NDim != u->NDim || NDim != v->NDim){
    printf("Incompatible vectors. Dim: %d %d\n",NDim,u->NDim); 
    return 0.;
  }
#endif
  double Resp=0.;
  for(int d=0;d<NDim;d++)
    Resp += u->x[d] * v->x[d];
  return Resp;
}
double Vettore::Col(int N){
#ifdef VECT_DEBUG
  if( N >= NDim || N < 0) return 0.;
#endif
  return x[N];
}
void Vettore::Set(double Val,int N){
#ifdef VECT_DEBUG
  if( N >= NDim || N < 0) return 0.;
#endif
  x[N] = Val;
}
double Vettore::CosAngle(Vettore *u,Vettore *v){
  double Resp = ScalS(u,v);
  Resp /= v->Norm()*u->Norm();
  return Resp;
}
double Vettore::CosAngle(Vettore *u){
#ifdef VECT_DEBUG
#endif
  double Resp = 0.;
  for(int d=0;d<NDim;d++)
    Resp += x[d] * u->x[d];
  Resp /= Norm();
  Resp /= u->Norm();
  return Resp;
}
double Vettore::SinAngle(Vettore *u,Vettore *v){
#ifdef VECT_DEBUG
#endif
  Vettore w(u->NDim);
  w.VetV(u,v);
  double Resp = w.Norm()/(u->Norm()*v->Norm());
  return Resp;
}
double Vettore::SinAngle(Vettore *u){
#ifdef VECT_DEBUG
#endif
  double Resp = 0.;
  for(int d=0;d<NDim;d++){
    int NUno = (d+1)%NDim;
    int NDue = (d+2)%NDim;
    Resp += SQR(x[NUno] * u->x[NDue] - x[NDue]*u->x[NUno]);
  }
  Resp = sqrt(Resp)/(Norm()*u->Norm());
  return Resp;
}
double Vettore::Angle(Vettore *u,Vettore *v){
  return acos(CosAngle(u,v));
}
double Vettore::Angle(Vettore *u){
  return acos(CosAngle(u));
}
void Vettore::Axis(Vettore *u,Vettore *v){
#ifdef VECT_DEBUG
  if(NDim != u->NDim || NDim != v->NDim){
    printf("Incompatible vectors. Dim: %d %d\n",NDim,u->NDim); 
    return ;
  }
#endif
  for(int d=0;d<NDim;d++) x[d] = .5*(u->x[d] + v->x[d]);
}
void Vettore::Normal(const Vettore *u,const Vettore *v){
#ifdef VECT_DEBUG
  if(v->NDim != u->NDim){
    printf("Incompatible vectors. Dim: %d %d\n",u->NDim,v->NDim); 
    return ;
  }
#endif
  VetV(u,v);
  Normalize();
}
void Vettore::NormalSurf(const Vettore *u,const Vettore *v,const Vettore *w){
#ifdef VECT_DEBUG
  if(v->NDim != u->NDim){
    printf("Incompatible vectors. Dim: %d %d\n",u->NDim,v->NDim); 
    return ;
  }
#endif
  for(int d=0;d<NDim;d++){
    int NUno = (d+1)%NDim;
    int NDue = (d+2)%NDim;
    double dUno = (u->x[NUno] - w->x[NUno])*(v->x[NDue] - w->x[NDue]);
    double dDue = (u->x[NDue] - w->x[NDue])*(v->x[NUno] - w->x[NUno]);
    x[d] = dUno - dDue;
    //printf("%d -(%d %d)/%d   %lf\n",d,NUno,NDue,NDim,x[d]);
  }
}
void Vettore::ScalV(const Vettore *u,const Vettore *v){
#ifdef VECT_DEBUG
  if(NDim != u->NDim || NDim != v->NDim){
    printf("Incompatible vectors. Dim: %d %d\n",NDim,u->NDim); 
    return ;
  }
#endif
  for(int d=0;d<NDim;d++)
    x[d] = u->x[d] * v->x[d];
}
double Vettore::VetV(const Vettore *u,const Vettore *v){
  assert(NDim == u->NDim);
  double Area = 0.;
  for(int d=0;d<NDim;d++){
    int NUno = (d+1)%NDim;
    int NDue = (d+2)%NDim;
    x[d] = u->x[NUno] * v->x[NDue] - u->x[NDue]*v->x[NUno];
    Area += x[d];
  }
  return ABS(Area);
}
double Vettore::VetV(const Vettore *u){
  assert(NDim == u->NDim);
  double Area = 0.;
  Vettore w(NDim);
  for(int d=0;d<NDim;d++){
    int NUno = (d+1)%NDim;
    int NDue = (d+2)%NDim;
    w.x[d] = x[NUno] * u->x[NDue] - x[NDue]*u->x[NUno];
    Area += w.x[d];
  }
  for(int d=0;d<NDim;d++){
    x[d] = w.x[d];
  }
  return ABS(Area);
}
double Vettore::VetV3(const Vettore *u,const Vettore *v){
  double Area = 0.;
  x[0] = u->x[1]*v->x[2] - u->x[2]*v->x[1];
  x[1] = u->x[2]*v->x[0] - u->x[0]*v->x[2];
  x[2] = u->x[0]*v->x[1] - u->x[1]*v->x[0];
  for(int d=0;d<3;d++){
    Area += x[d];
  }
  return ABS(Area);
}
void Vettore::Rescale(double Length){
  double Fact = Length/Norm();
  for(int d=0;d<NDim;d++)
    x[d] *= Fact;
}
double Vettore::ProjOnAxis(Vettore *Axis){
  double Pre = 0.;
  double Norma = 0.;
  for(int d=0;d<NDim;d++){
    Pre += x[d]*Axis->x[d];
    Norma += SQR(Axis->x[d]);
  }
  double Length = 0.;
  for(int d=0;d<NDim;d++)
    Length += x[d] = Axis->x[d]*Pre/Norma;
  return Length;
}
double Vettore::ProjOnAxis(Vettore *Pos,Vettore *Axis){
  double Cos = CosAngle(Pos,Axis);
  double InvNorma = 1./Axis->Norm();
  return Cos*InvNorma;
}
void Vettore::ApplyOn(Vettore *o){
  for(int d=0;d<NDim;d++)
    x[d] = o->x[d] - x[d];
}
void Vettore::Copy(Vettore *c){
  for(int d=0;d<NDim;d++)
    x[d] = c->x[d];
}
void Vettore::Export(double *xx){
  for(int d=0;d<NDim;d++)
    xx[d] = x[d];
}
void Vettore::PerpTo(Vettore *a){
  double Sin = SinAngle(a);
  double Norma = Norm();
  Vettore w(NDim);
  for(int d=0;d<NDim;d++){
    int NUno = (d+1)%NDim;
    int NDue = (d+2)%NDim;
    w.x[d] = x[NUno]*a->x[NDue] - x[NDue]*a->x[NUno];
  }
  VetV(&w,a);
  Normalize();
  Rescale(Sin*Norma);
}
double Vettore::PerpTo(Vettore *Pos,Vettore *Axis){
  Vettore Perp1(Axis->NDim);
  Vettore Perp2(Axis->NDim);
  Perp1.VetV(Pos,Axis);
  Perp2.VetV(&Perp1,Axis);
  double Distance = Perp1.Norm()/Axis->Norm();
  double NormaInv = 1./Perp2.Norm();
  for(int d=0;d<Axis->NDim;d++)
    x[d] = Perp2[d]*Distance*NormaInv;
  return Distance;
  // double Sin = SinAngle(Axis,Pos);
  // double Norma = Norm();
  // Vettore w(NDim);
  // for(int d=0;d<NDim;d++){
  //   int NUno = (d+1)%NDim;
  //   int NDue = (d+2)%NDim;
  //   w.x[d] = Axis->x[NUno] * a.x[NDue] - x[NDue]*a->x[NUno];
  // }
  // VetV(w,a);
  // Normalize();
  // Rescale(Sin*Norma);
}
double Vettore::PerpTo3(Vettore *Pos,Vettore *Axis){
  Vettore Perp1(3);
  Vettore Perp2(3);
  Perp1.VetV3(Pos,Axis);
  Perp2.VetV3(&Perp1,Axis);
  double Distance = Perp1.Norm()/Axis->Norm();
  double NormaInv = 1./Perp2.Norm();
  for(int d=0;d<3;d++)
    x[d] = Perp2[d]*Distance*NormaInv;
  return Distance;
}
double Vettore::ProjOnSurf(Vettore *S1,Vettore *S2,Vettore *S3,Vettore *P){
  double Known[3];
  double UnKnown[3];
  Vettore Normal(3);
  Normal.NormalSurf(S1,S2,S3);
  Normal.Normalize();
  Matrice Mat(3,3);
  for(int d=0;d<3;d++){
    Known[d] = S2->Val(d) + P->Val(d);
    Mat.Set(0,d,S1->Val(d)-S2->Val(d));
    Mat.Set(1,d,S2->Val(d)-S3->Val(d));
    Mat.Set(2,d,Normal.Val(d));
  }
  Mat.Solve(Known,UnKnown);
  //double Norma = 0.;
  for(int d=0;d<3;d++){
    //Norma += SQR(UnKnown[d]);
    x[d] = UnKnown[d];
  }
  return Norm();
}
