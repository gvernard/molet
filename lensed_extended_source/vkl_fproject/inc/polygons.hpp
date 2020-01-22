#ifndef TRIANGLES_HPP
#define TRIANGLES_HPP

#include <vector>
#include "massModels.hpp"

struct point {
  double x;
  double y;
};

struct triangle {
  point a;
  point b;
  point c;
};

int pnpoly(int nvert,double* vertx,double* verty,double testx,double testy);
void findBarycenter(int len,double* x,double* y,double& xc,double& yc);
std::vector<triangle> srcTriangles(int Ni,int Nj,double w,double h,double* x,double* y,CollectionMassModels* mycollection);
std::vector<triangle> createTriangles(int Ni,int Nj,double w,double h);
std::vector<triangle> selectTriangles(std::vector<triangle> triangles,double rlim);
void deflectTriangles(std::vector<triangle> &triangles,CollectionMassModels* mycollection);

#endif /* TRIANGLES_HPP */
