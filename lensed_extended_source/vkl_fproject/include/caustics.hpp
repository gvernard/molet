#ifndef CAUSTICS_HPP
#define CAUSTICS_HPP

#include <cstdlib>
#include <string>
#include <vector>

class RectGrid;

class Contour {
public:
  std::vector<double> x;  
  std::vector<double> y;

  Contour(){};
  Contour(const Contour& other);
  ~Contour(){};
};

void outputContours(std::vector<Contour> contours,std::string filepath);
std::vector<Contour> mooreNeighborTracing(RectGrid* image);
void padImage(RectGrid* image,RectGrid* paddedImage,double paddingColor);


#endif /* CAUSTICS_HPP */
