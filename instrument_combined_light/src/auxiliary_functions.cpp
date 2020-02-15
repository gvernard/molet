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

    
    //    double sum = 0.0;
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
	//	sum += this->scaled_psf->img[i*this->scaled_psf->Nj+j];
      }
    }

    //    for(int i=0;i<this->scaled_psf->Nm;i++){
    //      this->scaled_psf->img[i] /= sum;
    //    }

  }
