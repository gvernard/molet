#ifndef CONTOURTRACING_HPP
#define CONTOURTRACING_HPP

#include <cstdlib>
#include <string>
#include <vector>

class RectGrid;
class point;

class Contour {
public:
  std::vector<double> x;  
  std::vector<double> y;

  Contour(){};
  ~Contour(){};
};

void mooreNeighborTracing(RectGrid* image,std::vector<Contour*>& contours);
void padImage(RectGrid* image,RectGrid* paddedImage,double paddingColor);

#endif /* CONTOURTRACING_HPP */
