#include <fftw3.h>

#include "vkllib.hpp"
#include "auxiliary_functions.hpp"

void PSF::interpolatePSF(ImagePlane* mydata){
  //    double newPixSize  = (mydata->xmax - mydata->xmin)/mydata->Nj;
  double newPixSize  = (mydata->width)/(mydata->Nj);
  double origPixSize = this->original_psf->width/this->original_psf->Nj;
  
  // Decide on the profile width and height in pixels based on the input profile
  int newNj,newNi;
  if( mydata->width < this->original_psf->width ){
    newNj = mydata->Nj;
    newNi = mydata->Ni;
  } else {
    newNj = floor( mydata->Nj*(0.901*this->original_psf->width)/mydata->width );
    newNi = floor( mydata->Ni*(0.901*this->original_psf->height)/mydata->height );
  }
  
  double neww    = newNj*newPixSize;
  double xoffset = (this->original_psf->width - neww)/2.0;
  double newh    = newNi*newPixSize;
  double yoffset = (this->original_psf->height - newh)/2.0;
  
  this->scaled_psf = new ImagePlane(newNi,newNj,newh,neww);
  
  
  double sum = 0.0;
  double x,y,xp,yp,dx,dy,ddx,ddy,w00,w10,w01,w11,f00,f10,f01,f11;
  int ii,jj;
  
  for(int i=0;i<this->scaled_psf->Ni;i++){
    y  = yoffset+i*newPixSize;
    ii = floor( y/origPixSize );
    yp = ii*origPixSize;
    dy = (y - yp)/origPixSize;
    ddy = (1.0 - dy);
    
    for(int j=0;j<this->scaled_psf->Nj;j++){
      x  = xoffset+j*newPixSize;
      jj = floor( x/origPixSize );
      xp = jj*origPixSize;
      dx = (x - xp)/origPixSize;
      ddx = (1.0 - dx);
      
      // first index: i (y direction) second index: j (x direction)
      w00 = ddx*ddy;
      w01 = dx*ddy;
      w10 = dy*ddx;
      w11 = dx*dy;
      
      f00 = this->original_psf->img[ii*this->original_psf->Nj+jj];
      f01 = this->original_psf->img[ii*this->original_psf->Nj+jj+1];
      f10 = this->original_psf->img[(ii+1)*this->original_psf->Nj+jj];
      f11 = this->original_psf->img[(ii+1)*this->original_psf->Nj+jj+1];
      
      this->scaled_psf->img[i*this->scaled_psf->Nj+j] = f00*w00 + f10*w10 + f01*w01 + f11*w11;
      sum += this->scaled_psf->img[i*this->scaled_psf->Nj+j];
    }
  }
  
  for(int i=0;i<this->scaled_psf->Nm;i++){
    this->scaled_psf->img[i] /= sum;
  }
  
}


void PSF::cropPSF(double threshold){
  int Ncropx = 50;
  int Ncropy = 50;

  double* blur = NULL;
  while( true ){
    int loffx = floor(Ncropx/2.0);
    int toffy = floor(Ncropy/2.0);
    double sum = 0.0;
    blur = (double*) calloc(Ncropx*Ncropy,sizeof(double));
    int offset = (floor(this->scaled_psf->Ni/2.0)-toffy)*this->scaled_psf->Ni + (floor(this->scaled_psf->Nj/2.0)-loffx);
    for(int i=0;i<Ncropy;i++){
      for(int j=0;j<Ncropx;j++){
	blur[i*Ncropx+j] = this->scaled_psf->img[offset+i*this->scaled_psf->Ni+j];
	sum += blur[i*Ncropx+j];
      }
    }

    //std::cout << sum << std::endl;
    if( sum < threshold ){
      Ncropx += 2;
      Ncropy += 2;
      free(blur);
    } else {
      break;
    }
  }


  double psf_pix_size_x = this->scaled_psf->width/this->scaled_psf->Nj;
  double psf_pix_size_y = this->scaled_psf->height/this->scaled_psf->Ni;
  this->cropped_psf = new ImagePlane(Ncropy,Ncropx,Ncropy*psf_pix_size_y,Ncropx*psf_pix_size_x);
  for(int i=0;i<this->cropped_psf->Nm;i++){
    this->cropped_psf->img[i] = blur[i];
  }
  free(blur);
}

void PSF::createKernel(int Ni,int Nj){
  int bNx = this->cropped_psf->Nj/2.0;
  int bNy = this->cropped_psf->Ni/2.0;
  this->kernel = (double*) calloc(Ni*Nj,sizeof(double));
  for(int j=0;j<bNy;j++){
    for(int i=0;i<bNx;i++){
      this->kernel[j*Ni+i]                    = this->cropped_psf->img[bNy*2*bNx+bNx+j*2*bNx+i];
      this->kernel[Ni-bNx+j*Ni+i]             = this->cropped_psf->img[bNy*2*bNx+j*2*bNx+i];
      this->kernel[Ni*(Nj-bNy)+j*Ni+i]        = this->cropped_psf->img[bNx+2*bNx*j+i];
      this->kernel[Ni*(Nj-bNy)+Ni-bNx+j*Ni+i] = this->cropped_psf->img[2*bNx*j+i];
    }
  }
}


void PSF::convolve(ImagePlane* image){
  int Ni = image->Ni;
  int Nj = image->Nj;

  fftw_complex* f_image  = (fftw_complex*) fftw_malloc(Ni*Nj*sizeof(fftw_complex));
  fftw_complex* f_kernel = (fftw_complex*) fftw_malloc(Ni*Nj*sizeof(fftw_complex));
  
  fftw_plan p1;
  p1 = fftw_plan_dft_r2c_2d(Ni,Nj,this->kernel,f_kernel,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  p1 = fftw_plan_dft_r2c_2d(Ni,Nj,image->img,f_image,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  double dum1,dum2;
  for(int i=0;i<Ni;i++) {
    for(int j=0;j<Nj;j++) {
      dum1 = f_image[i*Nj+j][0]*f_kernel[i*Nj+j][0] - f_image[i*Nj+j][1]*f_kernel[i*Nj+j][1];
      dum2 = f_image[i*Nj+j][0]*f_kernel[i*Nj+j][1] + f_image[i*Nj+j][1]*f_kernel[i*Nj+j][0];
      f_image[i*Nj+j][0] = dum1;
      f_image[i*Nj+j][1] = dum2;
    }
  }
  
  p1 = fftw_plan_dft_c2r_2d(Ni,Nj,f_image,image->img,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  fftw_free(f_image);
  fftw_free(f_kernel);
  
  // Normalize output
  for(int i=0;i<Ni*Nj;i++){
    image->img[i] /= (Ni*Nj);
  }
}

offsetPSF PSF::offsetPSFtoPosition(double x,double y,ImagePlane* image){
  double pix_size_x = image->width/image->Nj;
  double pix_size_y = image->height/image->Ni;

  double dum_x = x - image->xmin;
  double dum_y = -(y - image->ymax);
  int ix = static_cast<int>(floor( dum_x/pix_size_x ));
  int iy = static_cast<int>(floor( dum_y/pix_size_y ));
  
  int j_img = ix - this->cropped_psf->Nj/2;
  int i_img = iy - this->cropped_psf->Ni/2;
  int offset_img = i_img*image->Nj + j_img;
  
  int offset_psf = 0;
  if( offset_img < 0 ){
    int j_psf = 0;
    int i_psf = 0;
    if( j_img < 0 ){
      j_psf = abs(j_img);
    }
    if( i_img < 0 ){
      i_psf = abs(i_img);
    }
    offset_img = 0;
    offset_psf = i_psf*this->cropped_psf->Nj + j_psf;
  }
  
  int nj,ni;
  int dum_j = image->Nj - (j_img+this->cropped_psf->Nj);
  if( dum_j < 0 ){
    nj = this->cropped_psf->Nj + dum_j; // dum_j is already negative
  } else {
    nj = this->cropped_psf->Nj;
  }
  int dum_i = image->Ni - (i_img+this->cropped_psf->Ni);
  if( dum_i < 0 ){
    ni = this->cropped_psf->Ni + dum_i; // dum_i is already negative
  } else {
    ni = this->cropped_psf->Ni;
  }
  
  offsetPSF PSFoffset;
  PSFoffset.offset_image = offset_img;
  PSFoffset.offset_cropped = offset_psf;
  PSFoffset.nj = nj;
  PSFoffset.ni = ni;
  return PSFoffset;
}
