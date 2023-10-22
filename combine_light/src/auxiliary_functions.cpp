#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <dirent.h>
#include <algorithm>

#include "auxiliary_functions.hpp"

#include "vkllib.hpp"
#include "instruments.hpp"
#include "noise.hpp"



// START:IMAGE MANIPULATION FUNCTIONS ========================================================================================
vkl::RectGrid createObsStatic(Instrument* mycam,vkl::RectGrid* super_extended,vkl::RectGrid* super_lens_light,int res_x,int res_y,double& F_conv_extended,double& F_conv_lens,bool convolve_lens){  
  // supersim is just carrying the resolution and extent of the image.
  int super_res_x = super_extended->Nx;
  int super_res_y = super_extended->Ny;
  double xmin = super_extended->xmin;
  double xmax = super_extended->xmax;
  double ymin = super_extended->ymin;
  double ymax = super_extended->ymax;
  
  // Get the fixed extended lensed light. This is in units of flux!
  vkl::RectGrid convolved_extended = *super_extended;
  mycam->convolve(super_extended,&convolved_extended);
  //writeCutout(&convolved_extended,out_path+"output/psf_lensed_image_super.fits");
  
  // Get the fixed lens galaxy light. This is in units of flux!
  vkl::RectGrid convolved_lens_light = *super_lens_light;
  if( convolve_lens ){
    mycam->convolve(super_lens_light,&convolved_lens_light);
    //writeCutout(&convolved_lens_light,out_path+"output/psf_lens_light_super.fits");  
  }
  
  // Calculate convolved fluxes
  F_conv_lens = convolved_lens_light.integrate();
  F_conv_extended = convolved_extended.integrate();
  
  // Combined light of the observed base image (integrating flux density from 'super' to observed resolution)
  vkl::RectGrid* base = new vkl::RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax); 
  for(int i=0;i<base->Nz;i++){
    base->z[i] = convolved_lens_light.z[i] + convolved_extended.z[i];
  }
  vkl::RectGrid obs_static = base->embeddedNewGrid(res_x,res_y,"average"); // This is in units of flux!
  delete(base);

  return obs_static;
}

vkl::RectGrid createObsPS(vkl::RectGrid* supersim,std::vector<double> image_signal,std::vector<offsetPSF>& PSFoffsets,std::vector<Instrument*>& instrument_list,std::vector<double>& psf_partial_sums,int res_x,int res_y){
  // Loop over the truncated PSF (through PSF_offsets) for each image, and add their light to the ps_light image that contains all the point source light.
  // All vectors must have the same size, equal to the number of multiple images.

  // supersim is just carrying the resolution and extent of the image.
  vkl::RectGrid super_ps_light(supersim->Nx,supersim->Ny,supersim->xmin,supersim->xmax,supersim->ymin,supersim->ymax);
  
  for(int q=0;q<PSFoffsets.size();q++){
    for(int i=0;i<PSFoffsets[q].ni;i++){
      for(int j=0;j<PSFoffsets[q].nj;j++){
	int index_img = PSFoffsets[q].offset_image + i*super_ps_light.Nx + j;
	int index_psf = PSFoffsets[q].offset_cropped + i*instrument_list[q]->cropped_psf->Nx + j;
	//super_ps_light.z[index_img] += 1.0;
	super_ps_light.z[index_img] += image_signal[q]*instrument_list[q]->cropped_psf->z[index_psf]/psf_partial_sums[q];
      }
    }
  }

  // Bin image from 'super' to observed resolution
  vkl::RectGrid obs_ps_light = super_ps_light.embeddedNewGrid(res_x,res_y,"average"); // This is (and should be) still in flux units!!!
  
  return obs_ps_light;  
}

void writeCutout(vkl::RectGrid* obs,std::string fname){
  std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
  std::vector<std::string> values{std::to_string(obs->xmin),std::to_string(obs->xmax),std::to_string(obs->ymin),std::to_string(obs->ymax)};
  std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
  vkl::FitsInterface::writeFits(obs->Nx,obs->Ny,obs->z,keys,values,descriptions,fname);
}


void writeAllCutouts(std::vector<double> tobs,Json::Value images,Json::Value samp_LC,vkl::RectGrid* super_extended,vkl::RectGrid* super_lens_light,Instrument* mycam,std::vector<offsetPSF>& PSFoffsets,std::vector<Instrument*>& instrument_list,std::vector<double>& psf_partial_sums,int res_x,int res_y,std::string mock,bool convolve_lens,std::string out_path){

  // Create list of PSF file names +++++++++++++++++++	    
  //std::vector<std::string> psf_fnames = getFileNames(tobs,argv[4]); // user-provided function to get the file names of the PSFs as a function of timestep t
  //++++++++++++++++++++++++++++++++++++++++++++++++++	    
  
  for(int t=0;t<tobs.size();t++){
    
    // Time-dependent PSF goes here ++++++++++++++++++++
    // //std::cout << psf_fnames[t] << std::endl;
    // mycam.replacePSF(psf_fnames[t]);
    // mycam.preparePSF(&supersim,0.999);
    // for(int q=0;q<instrument_list.size();q++){
    // 	PSFoffsets[q] = instrument_list[q]->offsetPSFtoPosition(images[q]["x"].asDouble(),images[q]["y"].asDouble(),&supersim);
    // }
    // for(int q=0;q<instrument_list.size();q++){
    // 	psf_partial_sums[q] = instrument_list[q]->sumPSF(&PSFoffsets[q]);
    // }
    // vkl::RectGrid obs_base_t = createObsBase(&mycam,&supersim,res_x,res_y,out_path);
    // ptr_obs_base = &obs_base_t;
    //++++++++++++++++++++++++++++++++++++++++++++++++++

    // The following line is only for a non-time-varying PSF
    double F_conv_extended,F_conv_lens;
    vkl::RectGrid obs_static = createObsStatic(mycam,super_extended,super_lens_light,res_x,res_y,F_conv_extended,F_conv_lens,convolve_lens); // This is in flux units!
    
    // construct a vector with the point source brightness in each image at the given time step
    std::vector<double> image_signal(images.size());
    for(int q=0;q<images.size();q++){
      image_signal[q] = samp_LC[q]["signal"][t].asDouble();
    }
    vkl::RectGrid obs_ps_light = createObsPS(super_extended,image_signal,PSFoffsets,instrument_list,psf_partial_sums,res_x,res_y);  // This is in flux units!

    // Calculate total flux of the PS light only
    //double F_conv_ps,F_conv_ps_mag;
    //obs_ps_light.integrate(F_conv_ps,F_conv_ps_mag,mycam->ZP);
    
    // The two grids, obs_ps_light and obs_static, have the same resolution, so I can simply add their fluxes
    for(int i=0;i<obs_ps_light.Nz;i++){
      obs_ps_light.z[i] += obs_static.z[i]; // <------ BE CAREFUL: 'obs_ps_light' now also contains the static light
    }
    
    // Adding a simple time-dependent noise here, no need to save the noise realizations
    convertGridFromFlux(&obs_ps_light,mycam->ZP);
    mycam->noise->initializeFromData(&obs_ps_light);
    mycam->noise->calculateNoise();
    mycam->noise->addNoise(&obs_ps_light);
    
    char buffer[100];
    sprintf(buffer,"%s%s/OBS_%s_%03d.fits",out_path.c_str(),mock.c_str(),mycam->name.c_str(),t);
    std::string fname(buffer);
    writeCutout(&obs_ps_light,fname);
  }
}

void convertGridFromFlux(vkl::RectGrid* obs,double ZP){
  for(int i=0;i<obs->Nz;i++){
    obs->z[i] = convertFromFlux(obs->z[i],ZP); // Convert to electrons/(s arcsec^2)
  }
}

double convertFromFlux(double flux,double ZP){
  if( flux == 0.0 ){
    return 0.0; // Flux is zero (cannot convert to electrons/(s arcsec^2))
  } else {
    return -2.5*log10(flux) + ZP; // Convert to electrons/(s arcsec^2)
  }
}
// END:IMAGE MANIPULATION FUNCTIONS ==========================================================================================




// START:TRANSFORM PSF =====================================================================================
TransformPSF::TransformPSF(double a,double b,double c,bool d,bool e){
  this->x0 = a;
  this->y0 = b;
  this->rot = c * 0.01745329251;// converting to rad
  this->flip_x = d;
  this->flip_y = e;

  this->cosrot = cos(this->rot);
  this->sinrot = sin(this->rot);
}

void TransformPSF::applyTransform(double xin,double yin,double& xout,double& yout){
  xout =  (xin-this->x0)*this->cosrot + (yin-this->y0)*this->sinrot;
  yout = -(xin-this->x0)*this->sinrot + (yin-this->y0)*this->cosrot;
  if( this->flip_x ){
    xout *= -1.0;
  }
  if( this->flip_y ){
    yout *= -1.0;
  }
}

double TransformPSF::interpolateValue(double x,double y,PSF* mypsf){
  return 1.0;
}
// END:TRANSFORM PSF =====================================================================================







// user provided function for getting the file names of a time dependent PSF
std::vector<std::string> getFileNames(std::vector<double> tobs,std::string path){
  DIR* dir;
  struct dirent* diread;
  std::vector<std::string> files;
  if( (dir=opendir(path.c_str())) != nullptr ){
    while( (diread=readdir(dir)) != nullptr ){
      files.push_back(diread->d_name);
    }
    closedir(dir);
  } else {
    // directory not found
  }
  std::sort(files.begin(),files.end(),std::greater<std::string>()); // alphanumerically sorted in descending order
  files.pop_back(); // removing '.' and '..'
  files.pop_back();
  
  // for(auto file:files){
  //   std::cout << file << std::endl;
  // }
  
  std::vector<std::string> names(tobs.size());
  for(int t=0;t<tobs.size();t++){
    int index = t % files.size();
    names[t] = path+files[index];
    //names[t] = "/home/george/Desktop/My Papers/unpublished/2021_molet/extract_with_cosmouline/psfs/3_ECAM.2019-03-01T04:54:57.000_psf_abcdef.fits"; // reference image
  }
  return names;
}

