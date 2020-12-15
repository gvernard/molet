#include <cstdlib>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "json/json.h"

#include "vkllib.hpp"
#include "instruments.hpp"
#include "caustics.hpp"

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
  double xdefl,ydefl;
  //================= END:PARSE INPUT =======================


  
  //=============== BEGIN:CREATE THE LENSES ====================
  Json::Value jlens = root["lenses"][0];

  // Initialize mass model collection
  CollectionMassModels mass_collection = JsonParsers::parse_mass_model(jlens["mass_model"],input);
  
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



  //=============== BEGIN:GET CRITICAL LINES AND CAUSTICS =======================
  RectGrid detA(super_res_x,super_res_y,xmin,xmax,xmin,xmax);

  for(int i=0;i<detA.Ny;i++){
    for(int j=0;j<detA.Nx;j++){
      double mydet = mass_collection.detJacobian(detA.center_x[j],detA.center_y[i]);
      if( mydet > 0 ){
	detA.z[i*detA.Nx+j] = 0;
      } else {
	detA.z[i*detA.Nx+j] = 1;
      }
    }
  }

  // Get detA contours on the image plane
  std::vector<Contour> contours = mooreNeighborTracing(&detA);

  // Add first point as last and close the polygon
  for(int i=0;i<contours.size();i++){
    contours[i].x.push_back( contours[i].x[0] );
    contours[i].y.push_back( contours[i].y[0] );
  }

  // Create caustic contours, and deflect the contours to fill them
  std::vector<Contour> caustics(contours.size());
  for(int k=0;k<contours.size();k++){
    caustics[k] = contours[k];
    for(int i=0;i<contours[k].x.size();i++){
      mass_collection.all_defl(contours[k].x[i],contours[k].y[i],xdefl,ydefl);
      caustics[k].x[i] = xdefl;
      caustics[k].y[i] = ydefl;
    }
  } 
  //================= END:GET CRITICAL LINES AND CAUSTICS =======================

  

  //=============== BEGIN:CREATE THE SOURCES =======================
  CollectionProfiles profile_collection = JsonParsers::parse_profile(root["source"]["light_profile"],input);
  //================= END:CREATE THE SOURCES =======================



  //=============== BEGIN:PRODUCE IMAGE USING RAY-SHOOTING =======================
  for(int i=0;i<mysim.Ny;i++){
    for(int j=0;j<mysim.Nx;j++){
      mass_collection.all_defl(mysim.center_x[j],mysim.center_y[i],xdefl,ydefl);
      mysim.z[i*mysim.Nx+j] = profile_collection.all_values(xdefl,ydefl);
    }
  }
  //================= END:PRODUCE IMAGE USING RAY-SHOOTING =======================


  
  //=============== BEGIN:OUTPUT =======================
  // Super-resolved lensed image
  std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
  std::vector<std::string> values{std::to_string(mysim.xmin),std::to_string(mysim.xmax),std::to_string(mysim.ymin),std::to_string(mysim.ymax)};
  std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
  FitsInterface::writeFits(mysim.Nx,mysim.Ny,mysim.z,keys,values,descriptions,output + "lensed_image_super.fits");
  
  // Super-resolved source image
  profile_collection.write_all_profiles(output + "source_super.fits");

  // Image plane magnification (0:positive, 1:negative)
  FitsInterface::writeFits(detA.Nx,detA.Ny,detA.z,output + "detA.fits");

  // Caustics and critical curves
  outputContours(contours,output+"criticals.json");
  outputContours(caustics,output+"caustics.json");
  //================= END:OUTPUT =======================
  

  
  return 0;
}
