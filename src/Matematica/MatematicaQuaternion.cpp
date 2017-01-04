#include <Matematica.h>


Quaternione::Quaternione(){
  m_x = m_y = m_z = 0.;
  m_w = 1.0;
}
void Quaternione::AxisRotation(double x,double y,double z, double degrees){
  double Angle = double((degrees / 180.0) * PI);
    // Here we calculate the sin( theta / 2) once for optimization
  double Result = (double)sin( Angle / 2.0 );
  // Calcualte the w value by cos( theta / 2 )
  m_w = (double)cos( Angle / 2.0 );
  // Calculate the x, y and z of the quaternion
  m_x = double(x * Result);
  m_y = double(y * Result);
  m_z = double(z * Result);
  //printf("%lf %lf %lf %lf\n",m_x,m_y,m_z,m_w);
}
void Quaternione::CreateMatrix(double *pMatrix)
{
    // Make sure the matrix has allocated memory to store the rotation data
    if(!pMatrix) return;
    // First row
    pMatrix[ 0] = 1.0 - 2.0 * ( m_y * m_y + m_z * m_z ); 
    pMatrix[ 1] = 2.0 * (m_x * m_y + m_z * m_w);
    pMatrix[ 2] = 2.0 * (m_x * m_z - m_y * m_w);
    pMatrix[ 3] = 0.0;  
    // Second row
    pMatrix[ 4] = 2.0 * ( m_x * m_y - m_z * m_w );  
    pMatrix[ 5] = 1.0 - 2.0 * ( m_x * m_x + m_z * m_z ); 
    pMatrix[ 6] = 2.0 * (m_z * m_y + m_x * m_w );  
    pMatrix[ 7] = 0.0;  
    // Third row
    pMatrix[ 8] = 2.0 * ( m_x * m_z + m_y * m_w );
    pMatrix[ 9] = 2.0 * ( m_y * m_z - m_x * m_w );
    pMatrix[10] = 1.0 - 2.0 * ( m_x * m_x + m_y * m_y );  
    pMatrix[11] = 0.0;  
    // Fourth row
    pMatrix[12] = 0;  
    pMatrix[13] = 0;  
    pMatrix[14] = 0;  
    pMatrix[15] = 1.0f;
    // Now pMatrix[] is a 4x4 homogeneous matrix that can be applied to an OpenGL Matrix
}

Quaternione Quaternione::operator *(Quaternione q)
{
    Quaternione r;
    r.m_w = m_w*q.m_w - m_x*q.m_x - m_y*q.m_y - m_z*q.m_z;
    r.m_x = m_w*q.m_x + m_x*q.m_w + m_y*q.m_z - m_z*q.m_y;
    r.m_y = m_w*q.m_y + m_y*q.m_w + m_z*q.m_x - m_x*q.m_z;
    r.m_z = m_w*q.m_z + m_z*q.m_w + m_x*q.m_y - m_y*q.m_x;
	
    return(r);
}
Quadri::Quadri(){
}
Quadri::Quadri(double *Vett,double Angle){
  w = cos(Angle*.5);
  double NormaInv = NormInv(Vett);
  double Sin = sin(Angle*.5);
  x = Sin*Vett[0]*NormaInv;
  y = Sin*Vett[1]*NormaInv;
  z = Sin*Vett[2]*NormaInv;
  //printf("%lf %lf %lf %lf\n",x,y,z,w);
}
Quadri::Quadri(double Pitch,double Yaw,double Roll){
  double cp = cos(Pitch*.5);
  double sp = sin(Pitch*.5);
  double cy = cos(Yaw*.5);
  double sy = sin(Yaw*.5);
  double cr = cos(Roll*.5);
  double sr = sin(Roll*.5);
  x = sr*cp*cy - cr*sp*sy;
  y = cr*sp*cy + sr*cp*sy;
  z = cr*cp*sy - sr*sp*cy;
  w = cr*cp*cy + sr*sp*sy;
  Normalize();
}
Quadri::Quadri(double ww,double xx,double yy,double zz){
  x = xx;
  y = yy;
  z = zz;
  w = ww;
}
void Quadri::RotMatrix(double *data,int dim){
  //FIXME: the determinant is not zero!
  if(dim == 4){
    int NRow = 4;
    data[NRow*0+0]  = w*w + x*x - y*y - z*z;
    data[NRow*0+1]  = 2.*x*y + 2.*w*z;
    data[NRow*0+2]  = 2.*x*z - 2.*w*y;
    data[NRow*0+3]  = 0.;
    
    data[NRow*1+0]  = 2.*x*y - 2.*w*z;
    data[NRow*1+1]  = w*w - x*x + y*y - z*z;
    data[NRow*1+2]  = 2.*y*z + 2.*w*x;
    data[NRow*1+3]  = 0.;
    
    data[NRow*2+0]  = 2.*x*z + 2.*w*y;
    data[NRow*2+1]  = 2.*y*z - 2.*w*x;
    data[NRow*2+2]  = w*w - x*x - y*y + z*z;
    data[NRow*2+3]  = 0.;
    
    data[NRow*3+0]  = 0.;
    data[NRow*3+1]  = 0.;
    data[NRow*3+2]  = 0.;
    data[NRow*3+3]  = w*w + x*x + y*y + z*z;
  }
  else{
    int NRow = 3;
    data[NRow*0+0]  = 1. - 2.*SQR(y) - 2.*SQR(z);
    data[NRow*0+1]  = 2.*x*y + 2.*w*z;
    data[NRow*0+2]  = 2.*x*z - 2.*w*y;
    
    data[NRow*1+0]  = 2.*x*y - 2.*w*z;
    data[NRow*1+1]  = 1. - 2.*SQR(x) - 2.*SQR(z);
    data[NRow*1+2]  = 2.*y*z + 2.*w*x;
    
    data[NRow*2+0]  = 2.*x*z + 2.*w*y;
    data[NRow*2+1]  = 2.*y*z - 2.*w*x;
    data[NRow*2+2]  = 1. - 2.*SQR(x) - 2.*SQR(y);
  }  
}
void Quadri::Basis(double a,double b,double c,double d,double *M){
  // first row
  M[0] = a;
  M[1] = -b;
  M[2] = d;
  M[3] = c;
  //
  M[4] = b;
  M[5] = a;
  M[6] = c;
  M[7] = -d;
  //
  M[8] = -d;
  M[9] = -c;
  M[10]= a;
  M[11]= -b;
  //
  M[12] = -c;
  M[13] = d;
  M[14] = b;
  M[15] = a;
}
void Quadri::PrintMatrix(double *M){
  printf("|%lf %lf %lf %lf|\n",M[0],M[4],M[8],M[12]);
  printf("|%lf %lf %lf %lf|\n",M[1],M[5],M[9],M[13]);
  printf("|%lf %lf %lf %lf|\n",M[2],M[6],M[10],M[14]);
  printf("|%lf %lf %lf %lf|\n",M[3],M[7],M[11],M[15]);
}
void Quadri::Conj(){
  x = -x;
  y = -y;
  z = -z;
  w = w;
}
Quadri Quadri::GetConj(){
  Quadri Resp;
  Resp.x = -x;
  Resp.y = -y;
  Resp.z = -z;
  Resp.w = w;
  return Resp;
}
double Quadri::Norm(){
  double Resp = 0.;
  Resp = QUAD(x) + QUAD(y) + QUAD(z) + QUAD(w);
  return sqrt(Resp);
}
double Quadri::Normalize(){
  double Den = 1./Norm();
  x = x*Den;
  y = y*Den;
  z = z*Den;
  w = w*Den;
  return Den;
}
double Quadri::Normalize(double *Vett){
  double Norm = 0.;
  for(int d=0;d<3;d++)
    Norm += SQR(Vett[d]);
  Norm = Norm > 0. ? 1./sqrt(Norm) : 1.;
  for(int d=0;d<3;d++)
    Vett[d] *= Norm;
  return 1./Norm;
}
double Quadri::NormInv(double *Vett){
  double Norm = 0.;
  for(int d=0;d<3;d++)
    Norm += SQR(Vett[d]);
  Norm = Norm > 0. ? 1./sqrt(Norm) : 1.;
  return Norm;
}
Quadri Quadri::operator* (const Quadri &rq) const
{
  // the constructor takes its arguments as (x, y, z, w)
  return Quadri(w * rq.x + x * rq.w + y * rq.z - z * rq.y,
		    w * rq.y + y * rq.w + z * rq.x - x * rq.z,
		    w * rq.z + z * rq.w + x * rq.y - y * rq.x,
		    w * rq.w - x * rq.x - y * rq.y - z * rq.z);
}
Quadri Quadri::operator= (const Quadri &rq) const
{
  // the constructor takes its arguments as (x, y, z, w)
  return Quadri(rq.w,rq.x,rq.y,rq.z);
}
// Multiplying a quaternion q with a vector v applies the q-rotation to v
double *Quadri::operator* (const double *Vet) const
{
  // product v' = q v q*
  double Resp[3];
  // for(int d=0;d<3;d++)
  //   Resp[d] = Vet[d];
  // Normalize(Resp);
  // Quadri vecQuat, resQuat;
  // vecQuat.x = Resp[0];
  // vecQuat.y = Resp[1];
  // vecQuat.z = Resp[2];
  // vecQuat.w = 0.0;
  
  // resQuat = vecQuat * GetConj();
  // resQuat = *this * resQuat;
  // Resp[0] = resQuat.x;
  // Resp[1] = resQuat.y;
  // Resp[2] = resQuat.z;

  return Resp;
}
void Quadri::Prod(Quadri q){
  double Uno = w*q.w - x*q.x - y*q.y - z*q.z;
  double Due = w*q.x + q.w*x + y*q.z - z*q.y;
  double Tre = w*q.y + q.w*y - x*q.z + z*q.x;
  double Qua = w*q.z + q.w*z + x*q.y - y*q.x;
  w = Uno;
  x = Due;
  y = Tre;
  z = Qua;
}
Quadri Quadri::Prod(Quadri p,Quadri q){
  Quadri Resp;
  Resp.w = p.w*q.w - p.x*q.x - p.y*q.y - p.z*q.z;
  Resp.x = p.w*q.x + q.w*p.x + p.y*q.z - p.z*q.y;
  Resp.y = p.w*q.y + q.w*p.y - p.x*q.z + p.z*q.x;
  Resp.z = p.w*q.z + q.w*p.z + p.x*q.y - p.y*q.x;
  //return Quadri;
}
double Quadri::Sqr(){
  return w*w + x*x + y*y * z*z; 
}
void Quadri::Inv(){
  Conj();
  double Num = 1./Sqr();
  w = w*Num;
  x = x*Num;
  y = y*Num;
  z = z*Num;
}
void Quadri::Matrix3x3(double *M){
  // first row
  M[0] = 1. - 2.*x*x - 2.*z*z;
  M[1] = 2.*x*y + 2.*w*z;
  M[2] = 2.*x*z - 2.*w*y;
  //
  M[3] = 2.*x*y - 2.*w*z;
  M[4] = 1.-2.*x*x-2.*z*z;
  M[5] = 2.*x*y + 2.*w*x;
  //
  M[6] = 2.*x*z + 2.*w*y;
  M[7] = 2.*y*z - 2.*w*x;
  M[8] = 1.-2.*x*x-2.*y*y;
}
void Quadri::Matrix4x4(double *M){
  // first row
  M[0] = 1. - 2.*x*x - 2.*z*z;
  M[1] = 2.*x*y + 2.*w*z;
  M[2] = 2.*x*z - 2.*w*y;
  M[3] = 0.;
  //
  M[4] = 2.*x*y - 2.*w*z;
  M[5] = 1.-2.*x*x-2.*z*z;
  M[6] = 2.*x*y + 2.*w*x;
  M[7] = 0.;
  //
  M[8] = 2.*x*z + 2.*w*y;
  M[9] = 2.*y*z - 2.*w*x;
  M[10] = 1.-2.*x*x-2.*y*y;
  M[11] = 0.;
  //
  M[12] = 0.;
  M[13] = 0.;
  M[14] = 0.;
  M[15] = 1.;
}
double *Quadri::Axis(){
  double v[3];
  double Den = Norm();
  v[0] = x/Den;
  v[1] = y/Den;
  v[2] = z/Den;
  return v;
}
double Quadri::Angle(){
  return 2.*acos(w);
}
