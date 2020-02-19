#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "vkllib.hpp"

#include "auxiliary_functions.hpp"

int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  std::ifstream fin;
  Json::Value::Members jmembers;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string input_path = argv[2];

  std::string output = argv[3];

  Json::Value images;
  if( root.isMember("point_source") ){
    // Read the image parameters
    fin.open(argv[4],std::ifstream::in);
    fin >> images;
    fin.close();
  }
  //================= END:PARSE INPUT =======================



  // Loop over the bands
  const Json::Value bands = root["instrument"]["bands"];
  for(int b=0;b<bands.size();b++){

    // Set output image plane in super-resolution
    double width  = bands[b]["field-of-view_x"].asDouble();
    double height = bands[b]["field-of-view_y"].asDouble();
    int super_res_x = 10*( static_cast<int>(ceil(width/bands[b]["resolution"].asDouble())) );
    int super_res_y = 10*( static_cast<int>(ceil(height/bands[b]["resolution"].asDouble())) );
    ImagePlane mysim(super_res_x,super_res_y,width,height);

    
    // Get the psf in super-resolution, crop it, and create convolution kernel
    std::string psf_path = input_path + "psf_" + bands[b]["name"].asString() + ".fits";
    double psf_width  = bands[b]["psf"]["width"].asDouble();
    double psf_height = bands[b]["psf"]["height"].asDouble();
    int psf_Nx = bands[b]["psf"]["pix_x"].asInt();
    int psf_Ny = bands[b]["psf"]["pix_y"].asInt();
    PSF mypsf(psf_path,psf_Nx,psf_Ny,psf_width,psf_height,&mysim);
    mypsf.cropPSF(0.99);
    mypsf.createKernel(mysim.Ni,mysim.Nj);

    
    //=============== BEGIN:CREATE THE FIXED EXTENDED LENSED LIGHT ====================
    ImagePlane* extended = new ImagePlane(output+"lensed_image_super.fits",super_res_x,super_res_y,width,height);
    mypsf.convolve(extended);
    //extended->writeImage(output+"psf_lensed_image_super.fits");
    //================= END:CREATE THE FIXED EXTENDED LENSED LIGHT ====================

    //=============== BEGIN:CREATE THE FIXED LENS GALAXY LIGHT ====================
    ImagePlane* lens_light = new ImagePlane(output+"lens_light_super.fits",super_res_x,super_res_y,width,height);
    mypsf.convolve(lens_light);
    //lens_light->writeImage(output+"psf_lens_light_super.fits");
    //================= END:CREATE THE FIXED LENS GALAXY LIGHT ====================


    //=============== BEGIN:CREATE FIXED MASKS ====================
    //================= END:CREATE FIXED MASKS ====================


    // Create combined extended lensed features and lens light
    ImagePlane base(super_res_x,super_res_y,width,height); 
    for(int i=0;i<base.Nm;i++){
      base.img[i] = lens_light->img[i] + extended->img[i];
    }
    delete(extended);
    delete(lens_light);


    
    if( root.isMember("point_source") ){
	
      //=============== BEGIN:CREATE THE TIME VARYING LIGHT ====================
      
      // Set the PSF related offsets for each image
      std::vector<offsetPSF> PSFoffsets(images.size());
      for(int q=0;q<images.size();q++){
	PSFoffsets[q] = mypsf.offsetPSFtoPosition(images[q]["x"].asDouble(),images[q]["y"].asDouble(),&mysim);
	//printf("%d %d %d %d\n",PSFoffsets[q].offset_image,PSFoffsets[q].offset_cropped,PSFoffsets[q].nj,PSFoffsets[q].ni);
      }

      // Read the intrinsic light curve
      
      // Read the microlensing light curves


      
      // Loop over time starts here
      ImagePlane pp_light(super_res_x,super_res_y,width,height);

      double micro_mag = 500.0;
      double intrinsic = 1.0;
      
      for(int q=0;q<images.size();q++){
	double macro_mag = abs(images[q]["mag"].asDouble());
	double factor = intrinsic*micro_mag*macro_mag;
	for(int i=0;i<PSFoffsets[q].ni;i++){
	  for(int j=0;j<PSFoffsets[q].nj;j++){
	    int index_img = pp_light.Nj*i + j;
	    int index_psf = i*mypsf.cropped_psf->Nj + j;
	    //	    pp_light.img[PSFoffsets[q].offset_image + pp_light.Nj*i + j] = 1.0;
	    pp_light.img[PSFoffsets[q].offset_image + index_img] += factor*mypsf.cropped_psf->img[PSFoffsets[q].offset_cropped + index_psf];
	  }
	}
      }

      for(int i=0;i<pp_light.Nm;i++){
      	pp_light.img[i] = pp_light.img[i] + base.img[i];
      }
      pp_light.writeImage(output + "OBS" + bands[b]["name"].asString() + ".fits");
      // Loop over time ends here, pp_light is destroyed

      //================= END:CREATE THE TIME VARYING LIGHT ====================

    } else {
      base.writeImage(output + "OBS" + bands[b]["name"].asString() + ".fits");
    }
    
  }
  
  return 0;
}
