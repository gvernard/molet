#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <string>

class ImagePlane;

class PSF {
public:
  ImagePlane* original_psf = NULL;
  ImagePlane* scaled_psf = NULL;

  PSF(std::string fname,int pix_x,int pix_y,double width,double height,ImagePlane* mydata){
    this->original_psf = new ImagePlane(fname,pix_x,pix_y,width,height);
    interpolatePSF(mydata);
  }

  ~PSF(){
    delete(original_psf);
    delete(scaled_psf);
  }

  void interpolatePSF(ImagePlane* mydata);
};

#endif /* AUXILIARY_HPP */
