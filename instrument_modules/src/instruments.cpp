#include <fftw3.h>
#include <algorithm>
#include <string>
#include <fstream>

#include "vkllib.hpp"
#include "json/json.h"

#include "instruments.hpp"
#include "noise.hpp"

std::string Instrument::path = INSTRUMENT_PATH;

// START:OFFSETPSF =================================================================================================
void offsetPSF::print(){
  printf("  %20s %d\n","Offset image:",this->offset_image);
  printf("  %20s %d\n","Offset cropped:",this->offset_cropped);
  printf("  %20s %d\n","nj:",this->nj);
  printf("  %20s %d\n","ni:",this->ni);
}


void offsetPSF::printFrame(FILE* fh,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  int* x = (int*) malloc(5*sizeof(int));
  int* y = (int*) malloc(5*sizeof(int));

  x[0] = this->offset_image % Nx;
  y[0] = (int)this->offset_image/Nx;
  x[1] = x[0] + this->nj;
  y[1] = y[0];
  x[2] = x[0] + this->nj;
  y[2] = y[0] + this->ni;
  x[3] = x[0];
  y[3] = y[0] + this->ni;
  x[4] = x[0];
  y[4] = y[0];

  double w = xmax - xmin;
  double h = ymax - ymin;
  for(int i=0;i<5;i++){
    fprintf(fh,"%10.4f %10.4f\n",xmin + w*x[i]/Nx,ymin + h*y[i]/Ny);
  }
  free(x);
  free(y);
}
// END:OFFSETPSF ===================================================================================================


// START:INSTRUMENT =================================================================================================
Instrument::Instrument(std::string name,double ZP,Json::Value noise_pars):name(name){
  this->ZP = ZP;
  this->common_constructor(noise_pars);
}


Instrument::Instrument(std::string name,Json::Value noise_pars):name(name){
  this->ZP = 0.0;
  this->common_constructor(noise_pars);
}


void Instrument::common_constructor(Json::Value noise_pars){
  std::string full_path = this->path + this->name + "/";
  
  Json::Value specs;
  std::ifstream fin(full_path+"specs.json",std::ifstream::in);
  fin >> specs;
  fin.close();

  this->lambda_min = specs["lambda_min"].asDouble();
  this->lambda_max = specs["lambda_max"].asDouble();
  this->resolution = specs["resolution"].asDouble();
  this->readout    = specs["readout"].asDouble();
  
  int pix_x  = specs["psf"]["pix_x"].asInt();
  int pix_y  = specs["psf"]["pix_y"].asInt();
  int width  = specs["psf"]["width"].asDouble();
  int height = specs["psf"]["height"].asDouble();
  this->original_psf = new RectGrid(pix_x,pix_y,0,width,0,height,full_path+"psf.fits");

  this->noise = FactoryNoiseModel::getInstance()->createNoiseModel(noise_pars,this);
}


std::string Instrument::createNewInstrument(Json::Value pars,std::string path_to_psf){
  std::cout << INSTRUMENT_PATH << std::endl;

  std::string name = "mapa";
  return name;
}


Instrument::~Instrument(){
  delete(original_psf);
  delete(scaled_psf);
  delete(cropped_psf);
  free(kernel);
  delete(noise);
}


double Instrument::getResolution(std::string name){
  Json::Value specs;
  
  std::ifstream fin(path + name + "/specs.json",std::ifstream::in);
  fin >> specs;
  fin.close();

  double res = specs["resolution"].asDouble();
  return res;
}


std::string Instrument::getName(){
  return this->name;
}


void Instrument::interpolatePSF(RectGrid* grid){
  //    double newPixSize  = (mydata->xmax - mydata->xmin)/mydata->Nj;
  double newPixSize  = (grid->width)/(grid->Nx);
  double origPixSize = this->original_psf->width/this->original_psf->Nx;
  
  // Decide on the profile width and height in pixels based on the input profile
  int newNx,newNy;
  if( grid->width < this->original_psf->width ){
    newNx = grid->Nx;
    newNy = grid->Ny;
  } else {
    //newNx = floor( grid->Nx*(0.901*this->original_psf->width)/grid->width );
    //newNy = floor( grid->Ny*(0.901*this->original_psf->height)/grid->height );
    //newNx = 10*this->original_psf->Nx;
    //newNy = 10*this->original_psf->Ny;
    newNx = floor(this->original_psf->width/newPixSize);
    newNy = floor(this->original_psf->height/newPixSize);
  }


  
  double neww    = newNx*newPixSize;
  //double neww    = this->original_psf->width;
  double xoffset = (this->original_psf->width - neww)/2.0;
  double newh    = newNy*newPixSize;
  //double newh    = this->original_psf->height;
  double yoffset = (this->original_psf->height - newh)/2.0;
  
  this->scaled_psf = new RectGrid(newNx,newNy,0,neww,0,newh);
  std::cout << grid->width << std::endl;
  std::cout << this->original_psf->width << std::endl;
  std::cout << neww << std::endl;

  std::cout << grid->width/grid->Nx << std::endl;
  std::cout << this->original_psf->width/this->original_psf->Nx << std::endl;
  std::cout << this->scaled_psf->width/newNx << std::endl;

  
  double sum = 0.0;
  for(int i=0;i<this->scaled_psf->Ny;i++){
    double y = this->scaled_psf->center_y[i];
    for(int j=0;j<this->scaled_psf->Nx;j++){
      double x = this->scaled_psf->center_x[j];
      double val = this->original_psf->interp2d_bilinear(x,y,this->original_psf->z);
      //this->scaled_psf->z[i*this->scaled_psf->Nx+j] = (this->original_psf->RectGrid::interp2d)(x,y,this->original_psf->z);
      this->scaled_psf->z[i*this->scaled_psf->Nx+j] = val;
      sum += val;
    }
  }

  /*  
  double sum = 0.0;
  double x,y,xp,yp,dx,dy,ddx,ddy,w00,w10,w01,w11,f00,f10,f01,f11;
  int ii,jj;

  for(int i=0;i<this->scaled_psf->Ny;i++){
    y  = yoffset+i*newPixSize;
    ii = floor( y/origPixSize );
    yp = ii*origPixSize;
    dy = (y - yp)/origPixSize;
    ddy = (1.0 - dy);
    
    for(int j=0;j<this->scaled_psf->Nx;j++){
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
      
      f00 = this->original_psf->z[ii*this->original_psf->Nx+jj];
      f01 = this->original_psf->z[ii*this->original_psf->Nx+jj+1];
      f10 = this->original_psf->z[(ii+1)*this->original_psf->Nx+jj];
      f11 = this->original_psf->z[(ii+1)*this->original_psf->Nx+jj+1];
      
      this->scaled_psf->z[i*this->scaled_psf->Nx+j] = f00*w00 + f10*w10 + f01*w01 + f11*w11;
      sum += this->scaled_psf->z[i*this->scaled_psf->Nx+j];
    }
  }
  */
  
  for(int i=0;i<this->scaled_psf->Nz;i++){
    this->scaled_psf->z[i] /= sum;
  }
}


void Instrument::cropPSF(double threshold){
  int Ncropx = 50;
  int Ncropy = 50;

  double* blur = NULL;
  while( true ){
    int loffx = floor(Ncropx/2.0);
    int toffy = floor(Ncropy/2.0);
    double sum = 0.0;
    blur = (double*) calloc(Ncropx*Ncropy,sizeof(double));
    int offset = (floor(this->scaled_psf->Ny/2.0)-toffy)*this->scaled_psf->Ny + (floor(this->scaled_psf->Nx/2.0)-loffx);
    for(int i=0;i<Ncropy;i++){
      for(int j=0;j<Ncropx;j++){
	blur[i*Ncropx+j] = this->scaled_psf->z[offset+i*this->scaled_psf->Ny+j];
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

  double psf_pix_size_x = this->scaled_psf->width/this->scaled_psf->Nx;
  double psf_pix_size_y = this->scaled_psf->height/this->scaled_psf->Ny;
  this->cropped_psf = new RectGrid(Ncropx,Ncropy,0,Ncropx*psf_pix_size_x,0,Ncropy*psf_pix_size_y);
  for(int i=0;i<this->cropped_psf->Nz;i++){
    this->cropped_psf->z[i] = blur[i];
  }
  free(blur);
}


void Instrument::createKernel(int Nx,int Ny){
  int bNx = this->cropped_psf->Nx/2.0;
  int bNy = this->cropped_psf->Ny/2.0;
  this->kernel = (double*) calloc(Nx*Ny,sizeof(double));
  for(int j=0;j<bNy;j++){
    for(int i=0;i<bNx;i++){
      this->kernel[j*Ny+i]                    = this->cropped_psf->z[bNy*2*bNx+bNx+j*2*bNx+i];
      this->kernel[Ny-bNx+j*Ny+i]             = this->cropped_psf->z[bNy*2*bNx+j*2*bNx+i];
      this->kernel[Ny*(Nx-bNy)+j*Ny+i]        = this->cropped_psf->z[bNx+2*bNx*j+i];
      this->kernel[Ny*(Nx-bNy)+Ny-bNx+j*Ny+i] = this->cropped_psf->z[2*bNx*j+i];
    }
  }
}


void Instrument::convolve(RectGrid* grid){
  int Nx = grid->Nx;
  int Ny = grid->Ny;

  fftw_complex* f_image  = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  fftw_complex* f_kernel = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  
  fftw_plan p1;
  p1 = fftw_plan_dft_r2c_2d(Nx,Ny,this->kernel,f_kernel,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  p1 = fftw_plan_dft_r2c_2d(Nx,Ny,grid->z,f_image,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  double dum1,dum2;
  for(int i=0;i<Ny;i++) {
    for(int j=0;j<Nx;j++) {
      dum1 = f_image[i*Nx+j][0]*f_kernel[i*Nx+j][0] - f_image[i*Nx+j][1]*f_kernel[i*Nx+j][1];
      dum2 = f_image[i*Nx+j][0]*f_kernel[i*Nx+j][1] + f_image[i*Nx+j][1]*f_kernel[i*Nx+j][0];
      f_image[i*Nx+j][0] = dum1;
      f_image[i*Nx+j][1] = dum2;
    }
  }
  
  p1 = fftw_plan_dft_c2r_2d(Nx,Ny,f_image,grid->z,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  fftw_free(f_image);
  fftw_free(f_kernel);
  
  // Normalize output
  for(int i=0;i<Nx*Ny;i++){
    grid->z[i] /= (Nx*Ny);
  }
}


offsetPSF Instrument::offsetPSFtoPosition(double x,double y,RectGrid* grid){
  int Nx_img   = grid->Nx;
  int Ny_img   = grid->Ny;
  int Nx_psf   = this->cropped_psf->Nx;
  int Ny_psf   = this->cropped_psf->Ny;
  double w_psf = this->cropped_psf->width;  // xmin of the cropped psf is always 0
  double h_psf = this->cropped_psf->height; // ymin of the cropped psf is always 0

  double dx = grid->step_x;
  double dy = grid->step_y;

  // Everything below is calculated in the reference frame centered on the multiple image position
  // The bottom left corner of the image
  double xz = grid->xmin - x;
  double yz = grid->ymin - y;
  // x0,y0 is the bottom left corner of the PSF
  double x0 = -w_psf/2.0;
  double y0 = -h_psf/2.0;

  int offset_img,offset_psf,ni,nj,ii,jj;
  if( x0<xz and xz<(x0+w_psf) and y0<yz and yz<(y0+h_psf) ){
    // case 1: the image bottom left corner is inside the PSF
    offset_img = 0;
    ii = static_cast<int>(floor( (yz - y0)/dy ));
    jj = static_cast<int>(floor( (xz - x0)/dx ));
    offset_psf = ii*Nx_psf + jj;
    nj = static_cast<int>(floor( (x0 + w_psf - xz)/dx ));
    ni = static_cast<int>(floor( (y0 + h_psf - yz)/dy ));
  } else if( xz<x0 and yz<y0 ){
    // case 2: the PSF is entirely inside the image
    offset_psf = 0;
    ii = static_cast<int>(floor( (y0 - yz)/dy ));
    jj = static_cast<int>(floor( (x0 - xz)/dx ));
    offset_img = ii*Nx_img + jj;
    if( (jj+Nx_psf) > Nx_img ){
      nj = Nx_img - jj;
    } else {
      nj = Nx_psf;
    }
    if( (ii+Ny_psf) > Ny_img ){
      ni = Ny_img - ii;
    } else {
      ni = Ny_psf;
    }
  } else if( x0<xz and xz<(x0+w_psf) and yz<y0 ){
    // case 3: a part of the PSF is outside the left side of the image
    ii = static_cast<int>(floor( (y0 - yz)/dy ));
    jj = static_cast<int>(floor( (xz - x0)/dx ));
    offset_psf = jj;
    offset_img = ii*Nx_img;
    nj = Nx_psf - jj;
    if( (ii+Ny_psf) > Ny_img ){
      ni = Ny_img - ii;
    } else {
      ni = Ny_psf;
    }
  } else if( xz<x0 and y0<yz and yz<(y0+h_psf) ){
    // case 4: a part of the PSF is outside the bottom side of the image
    ii = static_cast<int>(floor( (yz - y0)/dy ));
    jj = static_cast<int>(floor( (x0 - xz)/dx ));
    offset_img = jj;
    offset_psf = ii*Nx_psf;
    ni = Ny_psf - ii;
    if( (jj+Nx_psf) > Nx_img ){
      nj = Nx_img - jj;
    } else {
      nj = Nx_psf;
    }
  } else {
    // the PSF is outside the image (the multiple image is entirely outside the simulated field of view)
    offset_img = 0; //irrelevant
    offset_psf = 0; //irrelevant
    ni = 0;
    nj = 0;
  }

  offsetPSF PSFoffset;
  PSFoffset.offset_image = offset_img;
  PSFoffset.offset_cropped = offset_psf;
  PSFoffset.nj = nj;
  PSFoffset.ni = ni;
  return PSFoffset;  
}


void Instrument::replacePSF(std::string path_to_file){
  // the new PSF must have exactly the same dimensions and number of pixels as the one it is replacing
  RectGrid* temp = this->original_psf;
  this->original_psf = new RectGrid(this->original_psf->Nx,this->original_psf->Ny,0,this->original_psf->width,0,this->original_psf->height,path_to_file);
  delete(temp);
}


void Instrument::preparePSF(RectGrid* grid,double ratio){
  delete(this->scaled_psf);
  delete(this->cropped_psf);
  delete(this->kernel);
  this->interpolatePSF(grid);
  this->cropPSF(ratio);
  this->createKernel(grid->Nx,grid->Ny);
}


double Instrument::sumPSF(offsetPSF* psf_offset){
  double sum = 0.0;
  for(int i=0;i<psf_offset->ni;i++){
    for(int j=0;j<psf_offset->nj;j++){
      int index_psf = i*this->cropped_psf->Nx + j;
      sum += this->cropped_psf->z[index_psf];
    }
  }
  return sum;
}
// END:INSTRUMENT ===================================================================================================

