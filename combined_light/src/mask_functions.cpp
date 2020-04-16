#define _USE_MATH_DEFINES

#include <algorithm>
#include <string>
#include <cmath>
#include <fstream>
#include <memory>
#include <iostream>

#include <fftw3.h>

#include "vkllib.hpp"

void createMask(ImagePlane* mydata,double smear,double threshold,std::string outfile){
  //=============== BEGIN:INITIALIZATION =======================
  double img_max = 0.0;
  for(int i=0;i<mydata->Nm;i++){
    if( mydata->img[i] > img_max ){
      img_max = mydata->img[i];
    }
  }
  ImagePlane* image = new ImagePlane(mydata->Ni,mydata->Nj,mydata->width,mydata->height);
  double threshold_brightness = img_max*threshold;
  for(int i=0;i<mydata->Nm;i++){  
    if( mydata->img[i] > threshold_brightness ){
      image->img[i] = 1;
    } else {
      image->img[i] = 0;
    }
  }
  //================= END:INITIALIZATION =======================
  



  //=============== BEGIN:FIND RING RADII =======================
  // INNER RADIUS
  int i = 2; // starting with a 2 pixel radius
  int imax = (int) floor(image->Ni/2.0);
  bool clause = true;
  while( clause && i<imax ){
    double R = i*image->height/image->Ni;
    int Ntheta = (int) ceil( M_PI*R/(2.0*image->height/image->Ni) );
    double dtheta = M_PI/(2.0*Ntheta);
    for(int j=0;j<Ntheta;j++){
      double theta = j*dtheta;
      int iy = (int) floor(sin(theta)*R/(image->height/image->Ni));
      int ix = (int) floor(cos(theta)*R/(image->height/image->Ni));
      
      int index1 = (imax - iy)*image->Nj + image->Nj/2 + ix;
      int index2 = (imax - iy)*image->Nj + image->Nj/2 - ix;
      int index3 = (imax + iy)*image->Nj + image->Nj/2 + ix;
      int index4 = (imax + iy)*image->Nj + image->Nj/2 - ix;
      
      if( image->img[index1] == 1 || image->img[index2] == 1 || image->img[index3] == 1 || image->img[index4] == 1 ){
	clause = false;
	imax = i;
	break;
      }
    }
    
    i++;
  }
  double inner_radius = imax*image->height/image->Ni;
  //std::cout << "Inner radius is: " << imax << std::endl;



  // OUTER RADIUS
  //  FILE* fhh = fopen("indices.dat","w");
  i =  imax+2; // starting with the inner radius
  imax = (int) floor(image->Ni/2.0);
  clause = true;
  while( clause && i<imax ){
    clause = false;
    double R = i*image->height/image->Ni;
    int Ntheta = (int) ceil( M_PI*R/(2.0*image->height/image->Ni) );
    double dtheta = M_PI/(2.0*Ntheta);
    for(int j=0;j<Ntheta;j++){
      double theta = j*dtheta;
      int iy = (int) floor(sin(theta)*R/(image->height/image->Ni));
      int ix = (int) floor(cos(theta)*R/(image->height/image->Ni));
      
      int index1 = (imax - iy)*image->Nj + image->Nj/2 + ix;
      int index2 = (imax - iy)*image->Nj + image->Nj/2 - ix;
      int index3 = (imax + iy)*image->Nj + image->Nj/2 + ix;
      int index4 = (imax + iy)*image->Nj + image->Nj/2 - ix;

      //fprintf(fhh,"%d %d %d %d %f %f %f %f\n",index1,index2,index3,index4,image->img[index1],image->img[index2],image->img[index3],image->img[index4]);
      //fprintf(fhh,"%d %d\n",ix,iy);
      
      if( image->img[index1] == 1 || image->img[index2] == 1 || image->img[index3] == 1 || image->img[index4] == 1 ){
	clause = true;
	break;
      }
    }

    i++;
  }
  double outer_radius = (i-1)*image->height/image->Ni;
  //  std::cout << "Outer radius is: " << i*image->height/image->Ni << std::endl;
  //fclose(fhh);
  //=============== END:FIND RING RADII =======================


  


  //=============== BEGIN:CREATE KERNEL =======================
  ImagePlane blur = *image;
  int Ni = image->Ni;
  int Nj = image->Nj;

  smear = smear*(outer_radius - inner_radius)/3.0; // set the 3 sigma of the Guassian
  
  double dx = image->width/image->Nj;
  double dy = image->height/image->Ni;
  double factor1 = 1.0/(2.0*M_PI*pow(smear,2));
  double factor2 = 1.0/(2.0*pow(smear,2));
  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      double x = (j-Nj/2)*dx;
      double y = (i-Ni/2)*dy;
      double e = -factor2*(pow(x,2) + pow(y,2));
      blur.img[i*Nj+j] = factor1*exp(e);
    }
  }

  double blur_sum = 0.0;
  for(int i=0;i<Ni*Nj;i++){
    blur_sum += blur.img[i];
  }
  for(int i=0;i<Ni*Nj;i++){
    blur.img[i] /= blur_sum;
  }

  int bNx = Nj/2.0;
  int bNy = Ni/2.0;
  double* kernel = (double*) calloc(image->Ni*image->Nj,sizeof(double));
  for(int j=0;j<bNy;j++){
    for(int i=0;i<bNx;i++){
      kernel[j*Ni+i]                    = blur.img[bNy*2*bNx+bNx+j*2*bNx+i];
      kernel[Ni-bNx+j*Ni+i]             = blur.img[bNy*2*bNx+j*2*bNx+i];
      kernel[Ni*(Nj-bNy)+j*Ni+i]        = blur.img[bNx+2*bNx*j+i];
      kernel[Ni*(Nj-bNy)+Ni-bNx+j*Ni+i] = blur.img[2*bNx*j+i];
    }
  }
  //=============== END:CREATE KERNEL =========================

    


  //=============== BEGIN:CONVOLUTION =========================
  fftw_complex* f_image  = (fftw_complex*) fftw_malloc(Ni*Nj*sizeof(fftw_complex));
  fftw_complex* f_kernel = (fftw_complex*) fftw_malloc(Ni*Nj*sizeof(fftw_complex));
  
  fftw_plan p1;
  p1 = fftw_plan_dft_r2c_2d(Ni,Nj,kernel,f_kernel,FFTW_ESTIMATE);
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
  free(kernel);
  
  //Normalize output
  for(int i=0;i<Ni*Nj;i++){
    image->img[i] /= (Ni*Nj);
  }
  //=============== END:CONVOLUTION ===========================

    

  //=============== BEGIN:OUTPUT =======================
  //  FILE* fh = fopen("radii.dat","w");
  //  fprintf(fh,"%f\n",inner_radius);
  //  fprintf(fh,"%f\n",outer_radius);
  //  fclose(fh);


  double convolved_img_max = 0.0;
  for(int i=0;i<image->Nm;i++){
    if( image->img[i] > convolved_img_max ){
      convolved_img_max = image->img[i];
    }
  }

  threshold_brightness = convolved_img_max*threshold;
  for(int i=0;i<image->Nm;i++){  
    if( image->img[i] > threshold_brightness ){
      //if( image->img[i] > 0.00001 ){
      image->img[i] = 1;
    } else {
      image->img[i] = 0;
    }
  }

  image->writeImage(outfile);
  delete(image);
  //================= END:OUTPUT =======================
}
