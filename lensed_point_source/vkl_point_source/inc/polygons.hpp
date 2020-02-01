#ifndef TRIANGLES_HPP
#define TRIANGLES_HPP

#include <vector>

class ImagePlane;
class CollectionMassModels;

struct point {
  double x;
  double y;
};

struct triangle {
  point a;
  point b;
  point c;
};

struct itriangle {
  int ia;
  int ib;
  int ic;
};


std::vector<triangle> imagePlaneToTriangles(ImagePlane* image);
std::vector<itriangle> imagePlaneToTriangleIndices(ImagePlane* image);
void deflectTriangles(const std::vector<triangle>& triangles_in,std::vector<triangle>& triangles_out,CollectionMassModels* mycollection);
bool pointInTriangle(point p0,point p1,point p2,point p3);
double determinant3x3(std::vector<double> row1,std::vector<double> row2,std::vector<double> row3);
void circumcircle(point A,point B,point C,double& xc,double& xy,double& r);

#endif /* TRIANGLES_HPP */
