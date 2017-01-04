/*  colored_face.C   */     
/*  ---------------- */
#include <VarData.h>
//#ifndef __CGAL_h__
#ifdef USE_CGAL
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace CGAL;
using namespace std;

/* A facet with a color member variable. */
template < class Gt >
class My_face_base : public Triangulation_face_base_2<Gt>
{
public:
  Color color;
  My_face_base() :
    Triangulation_face_base_2<Gt>() {}
  My_face_base(void* v0, void* v1, void* v2) : 
    Triangulation_face_base_2<Gt>(v0,v1,v2) {}
  My_face_base(void* v0, void* v1, void* v2, void* n0, void* n1, void* n2) : 
    Triangulation_face_base_2<Gt>(v0,v1,v2,n0,n1,n2) {}
};

typedef Cartesian<double> Rp;
typedef Triangulation_euclidean_traits_2<Rp> Gt;
typedef Triangulation_vertex_base_2<Gt> Vb;
typedef My_face_base<Gt> Fb;
typedef Triangulation_default_data_structure_2<Gt,Vb,Fb > Tds;
typedef Triangulation_2<Gt,Tds> Triangulation;
//typedef Delaunay_graph DG;
typedef Delaunay_triangulation_2<Gt,Tds> DT;
typedef Delaunay_triangulation_adaptation_traits_2<DT> AT;
typedef Delaunay_triangulation_caching_degeneracy_removal_policy_2 <DT> AP;
typedef Voronoi_diagram_2<DT,AT,AP> VD;
typedef Point_2<Rp>  Point;
typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Face_iterator Face_iterator;
typedef Triangulation::Vertex Vertex;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Vertex_iterator Vertex_iterator;
// typedef VD::Face_handle Face_handle;
// typedef VD::Face_iterator Face_iterator;
// typedef VD::Vertex Vertex;
// typedef VD::Vertex_handle Vertex_handle;
// typedef VD::Vertex_iterator Vertex_iterator;
typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;

void VarData::AreaDistr(double *Distr,double *RadDistr,int Values) {  
  double BoxRadInv = pInvEdge(3);
  //  VD vd;
  double AreaMean = 3.*.5*Gen->NChain/(Gen->Edge[CLat1]*Gen->Edge[CLat2]);
  DT dtUp;
  DT dtDown;
  //ofstream File2Write("AreaSnapshot.dat");
  //Triangulation dt;
  //Point p;
  for(int c=0;c<Gen->NChain;c++){
    if(CHAIN_IF_TYPE(Ch[c].Type,CHAIN_UP) ){
      Point ChPos(Ch[c].Pos[CLat1],Ch[c].Pos[CLat2]);
      dtUp.insert(ChPos);
    }
    else{
      Point ChPos(Ch[c].Pos[CLat1],Ch[c].Pos[CLat2]);
      dtDown.insert(ChPos);      
    }
  }
  //  Face_iterator fc = t.faces_begin();
  Face_iterator fcTrUp = dtUp.finite_faces_begin();
  Face_iterator fcTrDown = dtDown.finite_faces_begin();
  //VD::Face_iterator fcVor = vd.faces_begin();
  //VD::Face_iterator fcVor = vd.faces_begin();
  double *Norma = new double[Values];
  memset(RadDistr,0.,Values*sizeof(double));
  vector <double> Areas;
  for(;fcTrUp != dtUp.faces_end(); ++fcTrUp){
    Vertex_handle vf1 = fcTrUp->vertex(0),
      vf2 = fcTrUp->vertex(1),vf3 = fcTrUp->vertex(2);
    Point pf1 = vf1->point();
    Point pf2 = vf2->point();
    Point pf3 = vf3->point();
    Vettore v1(pf1.x() - pf2.x(),pf1.y() - pf2.y(),0.);
    Vettore v2(pf3.x() - pf2.x(),pf3.y() - pf2.y(),0.);
    double Area = v1.VetV(&v2);
    double Distx = Nano->Pos[CLat1] - (pf1.x() + pf2.x() + pf3.x())/3.;
    double Disty = Nano->Pos[CLat2] - (pf1.y() + pf2.y() + pf3.y())/3.;
    double Rad = sqrt( SQR(Distx) + SQR(Disty) );
    if(Rad < Nano->Rad){continue;}
    int v = (int)(Rad*BoxRadInv*Values);
    //assert(v >= 0 && v < Values);
    if(v < 0 || v >= Values){continue;}
    RadDistr[v] += Area;
    Areas.push_back(Area);
    Norma[v] += 1.;
    //File2Write << vf1->point() << endl <<  vf2->point()<< endl << vf3->point() << endl << endl << endl ;
  }
  for(;fcTrDown != dtDown.faces_end(); ++fcTrDown){
    Vertex_handle vf1 = fcTrDown->vertex(0),
      vf2 = fcTrDown->vertex(1),vf3 = fcTrDown->vertex(2);
    Point pf1 = vf1->point();
    Point pf2 = vf2->point();
    Point pf3 = vf3->point();
    Vettore v1(pf1.x() - pf2.x(),pf1.y() - pf2.y(),0.);
    Vettore v2(pf3.x() - pf2.x(),pf3.y() - pf2.y(),0.);
    double Area = v1.VetV(&v2);
    double Distx = Nano->Pos[CLat1] - (pf1.x() + pf2.x() + pf3.x())/3.;
    double Disty = Nano->Pos[CLat2] - (pf1.y() + pf2.y() + pf3.y())/3.;
    double Rad = sqrt( SQR(Distx) + SQR(Disty) );
    //if(Rad < Nano->Rad){continue;}
    int v = (int)(Rad*BoxRadInv*Values);
    if(v < 0 || v >= Values){continue;}
    RadDistr[v] += Area;
    Areas.push_back(Area);
    Norma[v] += 1.;
    //File2Write << vf1->point() << endl <<  vf2->point()<< endl << vf3->point() << endl << endl << endl ;
  }
  // for(int v=0;v<Values;v++){
  //   RadDistr[v] /= Norma[v] > 0 ? Norma[v] : 1.;
  // }
  vector <double>::iterator iArea;
  for(iArea=Areas.begin();iArea!=Areas.end();iArea++){
    int v = (int)(*iArea/AreaMean*Values);
    if(v < 0 || v >= Values){continue;}
    Distr[v] += 1.;
  }
  delete [] Norma;
  //File2Write.close();
  return;
  //Finite_vertices_iterator vc = dt.finite_vertices_begin();
  // while(vc != dtUp.finite_vertices_end()){
  //   Face_handle fv = vc->face();
  //   Vertex_handle vf1 = fv->vertex(0),
  //     vf2 = fv->vertex(1),vf3 = fv->vertex(2);
  //   //File2Write << vf1->point() << endl << vf2->point() << endl << vf3->point() << endl << endl << endl;
  //   ++vc;
  // }
  // while(fcVor != vd.faces_end()){
  //   //    VD::Face fc1 = *fcVor;
  //   //VD::Vertex_handle vf1 = fcVor->vertex(0),//fcVor->vertex(0),
  //   //vf2 = fcVor->vertex(1),vf3 = fcVor->vertex(2);    
  //   //File2Write << vf1->point() << endl << vf2->point() << endl << vf3->point() << endl << endl << endl;
  //   ++fcVor;
  // }

  // Finite_vertices_iterator it1 = dt.finite_vertices_begin(),it2(it1), it3(it1);
  // ++it2;
  // ++it3; ++it3;
  // //    Vertex_handle vf0 = fc.vertex(0);
  // while( it3 != dt.finite_vertices_end()) {
  //   //File2Write << it1->point() << endl << it2->point() << endl << it3->point() << endl << endl << endl;
  //   ++it1; ++it2; ++it3; 
  // }
 
  // return 0;
}
#else
void VarData::AreaDistr(double *Distr,double *RadDistr,int Values){  
  printf("AreaDistr without CGAL not implemented\n");
}

#endif //USE_CGAL
