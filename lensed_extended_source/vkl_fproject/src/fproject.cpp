#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <map>

#include "json/json.h"

#include "vkllib.hpp"
#include "instruments.hpp"

int main(int argc,char* argv[]){
  /*
    Requires:
    - angular_diameter_distances.json
  */

  //=============== BEGIN:PARSE INPUT =======================
  std::ifstream fin;
  Json::Value::Members jmembers;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string in_path = argv[2];
  std::string input   = in_path+"input_files/";
  
  std::string out_path = argv[3];
  std::string output   = out_path+"output/";
  
  // Read the cosmological parameters
  Json::Value cosmo;
  fin.open(output+"angular_diameter_distances.json",std::ifstream::in);
  fin >> cosmo;
  fin.close();

  // Initialize image plane
  double xmin = root["instruments"][0]["field-of-view_xmin"].asDouble();
  double xmax = root["instruments"][0]["field-of-view_xmax"].asDouble();
  double ymin = root["instruments"][0]["field-of-view_ymin"].asDouble();
  double ymax = root["instruments"][0]["field-of-view_ymax"].asDouble();
  double resolution = Instrument::getResolution(root["instruments"][0]["name"].asString());
  int super_res_x = 10*( static_cast<int>(ceil((xmax-xmin)/resolution)) );
  int super_res_y = 10*( static_cast<int>(ceil((ymax-ymin)/resolution)) );
  RectGrid mysim(super_res_x,super_res_y,xmin,xmax,ymin,ymax);
  //================= END:PARSE INPUT =======================


  //=============== BEGIN:CREATE THE LENSES ====================
  Json::Value jlens = root["lenses"][0];


  // Initialize mass model collection
  CollectionMassModels mass_collection = JsonParsers::parse_mass_model(jlens["mass_model"]);
  
  // Scale dpsi mass models if necessary
  for(int k=0;k<jlens["mass_model"].size();k++){
    if( jlens["mass_model"][k]["pars"].isMember("scale_factor") ){
      double scale_factor = jlens["mass_model"][k]["pars"]["scale_factor"].asDouble();
      Pert* pert = static_cast<Pert*> (mass_collection.models[k]);
      for(int m=0;m<pert->Sm;m++){
	pert->z[m] *= scale_factor;
      }
      pert->updateDerivatives();
    }
  }
  //================= END:CREATE THE LENSES ====================


  //=============== BEGIN:CREATE THE SOURCES =======================
  CollectionProfiles profile_collection = JsonParsers::parse_profile(root["source"]["light_profile"]);
  //================= END:CREATE THE SOURCES =======================

  
  //=============== BEGIN:PRODUCE IMAGE USING RAY-SHOOTING =======================
  double xdefl,ydefl;
  for(int i=0;i<mysim.Ny;i++){
    for(int j=0;j<mysim.Nx;j++){
      mass_collection.all_defl(mysim.center_x[j],mysim.center_y[i],xdefl,ydefl);
      mysim.z[i*mysim.Nx+j] = profile_collection.all_values(xdefl,ydefl);
    }
  }
  //================= END:PRODUCE IMAGE USING RAY-SHOOTING =======================
  
  
  //=============== BEGIN:OUTPUT =======================
  // Super-resolved lensed image
  FitsInterface::writeFits(mysim.Nx,mysim.Ny,mysim.z,output + "lensed_image_super.fits");
  // Super-resolved source image
  //mysource->outputProfile(output + "source_super.fits");
  //================= END:OUTPUT =======================


  return 0;
}
