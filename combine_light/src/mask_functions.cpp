#define _USE_MATH_DEFINES

#include <algorithm>
#include <string>
#include <cmath>
#include <fstream>
#include <memory>
#include <iostream>

#include <fftw3.h>

#include "vkllib.hpp"

void createMask(RectGrid* mydata,double smear,double threshold,std::string outfile){
  //=============== BEGIN:INITIALIZATION =======================
  double img_max = 0.0;
  for(int i=0;i<mydata->Nz;i++){
    if( mydata->z[i] > img_max ){
      img_max = mydata->z[i];
    }
  }
  RectGrid* image = new RectGrid(mydata->Nx,mydata->Ny,mydata->xmin,mydata->xmax,mydata->ymin,mydata->ymax);
  double threshold_brightness = img_max*threshold;
  for(int i=0;i<mydata->Nz;i++){  
    if( mydata->z[i] > threshold_brightness ){
      image->z[i] = 1;
    } else {
      image->z[i] = 0;
    }
  }
  //================= END:INITIALIZATION =======================
  



  //=============== BEGIN:FIND RING RADII =======================
  // INNER RADIUS
  int i = 2; // starting with a 2 pixel radius
  int imax = (int) floor(image->Ny/2.0);
  bool clause = true;
  while( clause && i<imax ){
    double R = i*image->height/image->Ny;
    int Ntheta = (int) ceil( M_PI*R/(2.0*image->height/image->Ny) );
    double dtheta = M_PI/(2.0*Ntheta);
    for(int j=0;j<Ntheta;j++){
      double theta = j*dtheta;
      int iy = (int) floor(sin(theta)*R/(image->height/image->Ny));
      int ix = (int) floor(cos(theta)*R/(image->height/image->Ny));
      
      int index1 = (imax - iy)*image->Nx + image->Nx/2 + ix;
      int index2 = (imax - iy)*image->Nx + image->Nx/2 - ix;
      int index3 = (imax + iy)*image->Nx + image->Nx/2 + ix;
      int index4 = (imax + iy)*image->Nx + image->Nx/2 - ix;
      
      if( image->z[index1] == 1 || image->z[index2] == 1 || image->z[index3] == 1 || image->z[index4] == 1 ){
	clause = false;
	imax = i;
	break;
      }
    }
    
    i++;
  }
  double inner_radius = imax*image->height/image->Ny;
  //std::cout << "Inner radius is: " << imax << std::endl;



  // OUTER RADIUS
  //  FILE* fhh = fopen("indices.dat","w");
  i =  imax+2; // starting with the inner radius
  imax = (int) floor(image->Ny/2.0);
  clause = true;
  while( clause && i<imax ){
    clause = false;
    double R = i*image->height/image->Ny;
    int Ntheta = (int) ceil( M_PI*R/(2.0*image->height/image->Ny) );
    double dtheta = M_PI/(2.0*Ntheta);
    for(int j=0;j<Ntheta;j++){
      double theta = j*dtheta;
      int iy = (int) floor(sin(theta)*R/(image->height/image->Ny));
      int ix = (int) floor(cos(theta)*R/(image->height/image->Ny));
      
      int index1 = (imax - iy)*image->Nx + image->Nx/2 + ix;
      int index2 = (imax - iy)*image->Nx + image->Nx/2 - ix;
      int index3 = (imax + iy)*image->Nx + image->Nx/2 + ix;
      int index4 = (imax + iy)*image->Nx + image->Nx/2 - ix;

      //fprintf(fhh,"%d %d %d %d %f %f %f %f\n",index1,index2,index3,index4,image->img[index1],image->img[index2],image->img[index3],image->img[index4]);
      //fprintf(fhh,"%d %d\n",ix,iy);
      
      if( image->z[index1] == 1 || image->z[index2] == 1 || image->z[index3] == 1 || image->z[index4] == 1 ){
	clause = true;
	break;
      }
    }

    i++;
  }
  double outer_radius = (i-1)*image->height/image->Ny;
  //  std::cout << "Outer radius is: " << i*image->height/image->Ni << std::endl;
  //fclose(fhh);
  //=============== END:FIND RING RADII =======================


  


  //=============== BEGIN:CREATE KERNEL =======================
  RectGrid blur = *image;
  int Nx = image->Nx;
  int Ny = image->Ny;

  smear = smear*(outer_radius - inner_radius)/3.0; // set the 3 sigma of the Guassian
  
  double dx = image->step_x;
  double dy = image->step_y;
  double factor1 = 1.0/(2.0*M_PI*pow(smear,2));
  double factor2 = 1.0/(2.0*pow(smear,2));
  for(int i=0;i<Ny;i++){
    for(int j=0;j<Nx;j++){
      double x = (j-Nx/2)*dx;
      double y = (i-Ny/2)*dy;
      double e = -factor2*(pow(x,2) + pow(y,2));
      blur.z[i*Nx+j] = factor1*exp(e);
    }
  }

  double blur_sum = 0.0;
  for(int i=0;i<Nx*Ny;i++){
    blur_sum += blur.z[i];
  }
  for(int i=0;i<Nx*Ny;i++){
    blur.z[i] /= blur_sum;
  }

  int bNx = Nx/2.0;
  int bNy = Ny/2.0;
  double* kernel = (double*) calloc(image->Nx*image->Ny,sizeof(double));
  for(int j=0;j<bNy;j++){
    for(int i=0;i<bNx;i++){
      kernel[j*Ny+i]                    = blur.z[bNy*2*bNx+bNx+j*2*bNx+i];
      kernel[Ny-bNx+j*Ny+i]             = blur.z[bNy*2*bNx+j*2*bNx+i];
      kernel[Ny*(Nx-bNy)+j*Ny+i]        = blur.z[bNx+2*bNx*j+i];
      kernel[Ny*(Nx-bNy)+Ny-bNx+j*Ny+i] = blur.z[2*bNx*j+i];
    }
  }
  //=============== END:CREATE KERNEL =========================

    


  //=============== BEGIN:CONVOLUTION =========================
  fftw_complex* f_image  = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  fftw_complex* f_kernel = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  
  fftw_plan p1;
  p1 = fftw_plan_dft_r2c_2d(Nx,Ny,kernel,f_kernel,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  p1 = fftw_plan_dft_r2c_2d(Nx,Ny,image->z,f_image,FFTW_ESTIMATE);
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
  
  p1 = fftw_plan_dft_c2r_2d(Nx,Ny,f_image,image->z,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  fftw_free(f_image);
  fftw_free(f_kernel);
  free(kernel);
  
  //Normalize output
  for(int i=0;i<Nx*Ny;i++){
    image->z[i] /= (Nx*Ny);
  }
  //=============== END:CONVOLUTION ===========================

    

  //=============== BEGIN:OUTPUT =======================
  //  FILE* fh = fopen("radii.dat","w");
  //  fprintf(fh,"%f\n",inner_radius);
  //  fprintf(fh,"%f\n",outer_radius);
  //  fclose(fh);


  double convolved_img_max = 0.0;
  for(int i=0;i<image->Nz;i++){
    if( image->z[i] > convolved_img_max ){
      convolved_img_max = image->z[i];
    }
  }

  threshold_brightness = convolved_img_max*threshold;
  for(int i=0;i<image->Nz;i++){  
    if( image->z[i] > threshold_brightness ){
      //if( image->img[i] > 0.00001 ){
      image->z[i] = 1;
    } else {
      image->z[i] = 0;
    }
  }

  FitsInterface::writeFits(image->Nx,image->Ny,image->z,outfile);
  delete(image);
  //================= END:OUTPUT =======================
}
