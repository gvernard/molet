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

void mooreNeighborTracing(ImagePlane* image,std::vector<Contour*>& contours);
void padImage(ImagePlane* image,ImagePlane* paddedImage,double paddingColor);

#endif /* CONTOURTRACING_HPP */
