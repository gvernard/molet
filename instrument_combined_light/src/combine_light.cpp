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

    
    // Get the psf in super-resolution
    std::string psf_path = input_path + "psf_" + bands[b]["name"].asString() + ".fits";
    double psf_width  = bands[b]["psf"]["width"].asDouble();
    double psf_height = bands[b]["psf"]["height"].asDouble();
    int psf_Nx = bands[b]["psf"]["pix_x"].asInt();
    int psf_Ny = bands[b]["psf"]["pix_y"].asInt();
    PSF mypsf(psf_path,psf_Nx,psf_Ny,psf_width,psf_height,&mysim);

    std::cout << "mpip" << std::endl;
    
    //=============== BEGIN:CREATE THE FIXED EXTENDED LENSED LIGHT ====================
    //================= END:CREATE THE FIXED EXTENDED LENSED LIGHT ====================

    
    //=============== BEGIN:CREATE THE FIXED LENS GALAXY LIGHT ====================
    //================= END:CREATE THE FIXED LENS GALAXY LIGHT ====================


    //=============== BEGIN:CREATE FIXED MASKS ====================
    //================= END:CREATE FIXED MASKS ====================

    
    if( root.isMember("point_source") ){
	
	//=============== BEGIN:CREATE THE TIME VARYING LIGHT ====================
	
	// Create output super-resolved image that has only the variable light

	// Loop over time
	//   Loop over images
	//   Combine light
	
	//================= END:CREATE THE TIME VARYING LIGHT ====================

    } else {
      // Combine light
    }
    
  }

 

  return 0;
}
