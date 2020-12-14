#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "vkllib.hpp"
#include "instruments.hpp"

int main(int argc,char* argv[]){

  /*
    Requires:
    - angular_diameter_distances.json
    - multiple_images.json
  */
  
  //=============== BEGIN:PARSE INPUT =======================
  std::ifstream fin;
  Json::Value::Members jmembers;

  // Read the main parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string in_path = argv[2];
  std::string input   = in_path+"input_files/";

  std::string out_path = argv[3];
  std::string output = out_path+"output/";
  
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
  double res  = Instrument::getResolution(root["instruments"][0]["name"].asString());
  int super_res_x = 10*( static_cast<int>(ceil((xmax-xmin)/res)) );
  int super_res_y = 10*( static_cast<int>(ceil((ymax-ymin)/res)) );
  //================= END:PARSE INPUT =======================







  //=============== BEGIN:CREATE LENS LIGHT =======================
  Json::Value all_lenses;
  for(int k=0;k<root["lenses"].size();k++){
    for(int m=0;m<root["lenses"][k]["light_profile"].size();m++){
      all_lenses.append(root["lenses"][k]["light_profile"][m]);
    }
  }
  CollectionProfiles light_collection = JsonParsers::parse_profile(all_lenses,input);

  RectGrid mylight(super_res_x,super_res_y,xmin,xmax,ymin,ymax);
  for(int i=0;i<mylight.Ny;i++){
    for(int j=0;j<mylight.Nx;j++){
      mylight.z[i*mylight.Nx+j] = light_collection.all_values(mylight.center_x[j],mylight.center_y[i]);
    }
  }

  // Super-resolved lens light profile image
  FitsInterface::writeFits(mylight.Nx,mylight.Ny,mylight.z,output + "lens_light_super.fits");


  // Confirm that the total brightness is conserved (by numerical integration)
  /*
  double sum = 0.0;
  for(int i=0;i<mylight.Nm;i++){
    sum += mylight.img[i];
  }
  double fac = (width/super_res_x)*(height/super_res_y);
  sum *= fac;
  printf("Itot = %15.10f  Mtot = %15.10f\n",sum,-2.5*log10(sum));
  */
  //================= END:CREATE LENS LIGHT ================






  //=============== BEGIN:CREATE LENS COMPACT MASS ================
  if( root.isMember("point_source") ){
    // Factor to convert surface mass density to kappa
    double Dl  = cosmo[0]["Dl"].asDouble();
    double Ds  = cosmo[0]["Ds"].asDouble();
    double Dls = cosmo[0]["Dls"].asDouble();
    double sigma_crit = 3472.8*Ds/(Dl*Dls); // the critical density: c^2/(4pi G)  Ds/(Dl*Dls), in units of kg/m^2
    
    Json::Value all_compact;
    for(int k=0;k<root["lenses"].size();k++){
      for(int m=0;m<root["lenses"][k]["compact_mass_model"].size();m++){
	all_compact.append(root["lenses"][k]["compact_mass_model"][m]);
      }
    }
    CollectionProfiles compact_collection = JsonParsers::parse_profile(all_compact);

    // Write overall kappa_star field
    RectGrid kappa_star(super_res_x,super_res_y,xmin,xmax,ymin,ymax);
    for(int i=0;i<kappa_star.Ny;i++){
      for(int j=0;j<kappa_star.Nx;j++){
	kappa_star.z[i*kappa_star.Nx+j] = compact_collection.all_values(kappa_star.center_x[j],kappa_star.center_y[i])/sigma_crit;
      }
    }
    // Super-resolved lens compact mass profile image
    FitsInterface::writeFits(kappa_star.Nx,kappa_star.Ny,kappa_star.z,output + "lens_kappa_star_super.fits");

    
    // Read the image parameters
    Json::Value images;
    fin.open(output+"multiple_images.json",std::ifstream::in);
    fin >> images;
    fin.close();
      
    // Find the kappa_star at the multiple images
    for(int j=0;j<images.size();j++){
      double x = images[j]["x"].asDouble();
      double y = images[j]["y"].asDouble();
      double kappa_star = compact_collection.all_values(x,y)/sigma_crit;
      images[j]["s"] = 1.0 - kappa_star/images[j]["k"].asDouble();
    }

    // Overwrite multiple images file with values of s
    std::ofstream file_images(output+"multiple_images.json");
    file_images << images;
    file_images.close();
  }  
  //================= END:CREATE LENS COMPACT MASS ================


  return 0;
}
