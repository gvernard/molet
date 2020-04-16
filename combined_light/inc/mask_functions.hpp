#ifndef MASK_HPP
#define MASK_HPP

#include <string>

class ImagePlane;

void createMask(ImagePlane* image,double smear,double threshold,std::string outfile);

#endif /* MASK_HPP */
