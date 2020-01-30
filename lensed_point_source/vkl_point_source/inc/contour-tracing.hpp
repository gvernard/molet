#ifndef CONTOURTRACING_HPP
#define CONTOURTRACING_HPP

#include <cstdlib>
#include <string>
#include <vector>

class ImagePlane;
class point;

class Contour {
public:
  std::vector<double> x;  
  std::vector<double> y;

  Contour(){};
  ~Contour(){};
};

// Contour tracing: Moore Neighbor Tracing
void mooreNeighborTracing(ImagePlane* image,std::vector<Contour*>& contours);
void padImage(ImagePlane* image,ImagePlane* paddedImage,double paddingColor);

// Image multiplicity: Winding number
double isLeft(point P0,point P1,point P2);
int windingNumber(point P,Contour* C);

#endif /* CONTOURTRACING_HPP */
