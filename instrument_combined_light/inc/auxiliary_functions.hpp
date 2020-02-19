#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <string>

class ImagePlane;

class offsetPSF {
public:
  int offset_image;
  int offset_cropped;
  int nj;
  int ni;
};

class PSF {
public:
  ImagePlane* original_psf = NULL;
  ImagePlane* scaled_psf   = NULL;
  ImagePlane* cropped_psf  = NULL;
  double* kernel           = NULL;
  
  PSF(std::string fname,int pix_x,int pix_y,double width,double height,ImagePlane* mydata){
    this->original_psf = new ImagePlane(fname,pix_x,pix_y,width,height);
    interpolatePSF(mydata);
  }

  ~PSF(){
    delete(original_psf);
    delete(scaled_psf);
    delete(cropped_psf);
    free(kernel);
  }

  void interpolatePSF(ImagePlane* image);
  void cropPSF(double threshold);
  void createKernel(int Ni,int Nj);
  void convolve(ImagePlane* image);
  offsetPSF offsetPSFtoPosition(double x,double y,ImagePlane* image);
};

#endif /* AUXILIARY_HPP */
