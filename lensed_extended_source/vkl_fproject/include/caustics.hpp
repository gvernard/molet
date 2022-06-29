#ifndef CAUSTICS_HPP
#define CAUSTICS_HPP

#include <cstdlib>
#include <string>
#include <vector>
#include "vkllib.hpp"

class Contour {
public:
  std::vector<double> x;  
  std::vector<double> y;

  Contour(){};
  Contour(const Contour& other);
  ~Contour(){};
};

void outputContours(std::vector<Contour> contours,std::string filepath);
std::vector<Contour> mooreNeighborTracing(vkl::RectGrid* image);
void padImage(vkl::RectGrid* image,vkl::RectGrid* paddedImage,double paddingColor);


#endif /* CAUSTICS_HPP */
