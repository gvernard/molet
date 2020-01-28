/**
 * This algorithm is called Moore Neighbor Tracing
 * An explanation of the algorithm can be found here:
 * http://www.thebigblob.com/moore-neighbor-tracing-algorithm-in-c/
 * http://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/moore.html
 *
 * @author Erik Smistad <smistad@idi.ntnu.no>
 */

#ifndef CONTOURTRACING_HPP
#define CONTOURTRACING_HPP

#include <cstdlib>
#include <string>
#include <vector>

class ImagePlane;

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
