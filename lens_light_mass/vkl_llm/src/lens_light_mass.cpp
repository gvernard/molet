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
  //================= END:PARSE INPUT =======================







  //=============== BEGIN:CREATE LENS LIGHT =======================
  for(int k=0;k<root["instruments"].size();k++){
    std::string name = root["instruments"][k]["name"].asString();

    double xmin = root["instruments"][k]["field-of-view_xmin"].asDouble();
    double xmax = root["instruments"][k]["field-of-view_xmax"].asDouble();
    double ymin = root["instruments"][k]["field-of-view_ymin"].asDouble();
    double ymax = root["instruments"][k]["field-of-view_ymax"].asDouble();
    double res  = Instrument::getResolution(name);
    int super_factor = 10;
    if( root["output_options"].isMember("super_factor") ){
      super_factor = root["output_options"]["super_factor"].asInt();
    }
    int super_res_x = super_factor*( static_cast<int>(ceil((xmax-xmin)/res)) );
    int super_res_y = super_factor*( static_cast<int>(ceil((ymax-ymin)/res)) );
    vkl::RectGrid mylight(super_res_x,super_res_y,xmin,xmax,ymin,ymax);
    
    Json::Value all_lenses;
    for(int i=0;i<root["lenses"].size();i++){
      for(int j=0;j<root["lenses"][i]["light_profile"][name].size();j++){
	all_lenses.append(root["lenses"][i]["light_profile"][name][j]);
      }
    }
    vkl::CollectionProfiles light_collection = vkl::JsonParsers::parse_profile(all_lenses,root["instruments"][k]["ZP"].asDouble(),input);

    // Confirm that the total brightness is conserved (by numerical integration)
    double total_flux = 0.0;
    for(int i=0;i<light_collection.profiles.size();i++){
      double integral = light_collection.profiles[i]->integrate(1000);
      double integral_mag = -2.5*log10(integral) + root["instruments"][k]["ZP"].asDouble();
      std::cout << "Total flux from lens profile " << i << ": " << integral << " " << integral_mag << std::endl;
      total_flux += integral;
    }
    double total_flux_mag = -2.5*log10(total_flux) + root["instruments"][k]["ZP"].asDouble();
    std::cout << "Total flux all lens profiles: " << total_flux << " " << total_flux_mag << std::endl;

    
    for(int i=0;i<mylight.Ny;i++){
      for(int j=0;j<mylight.Nx;j++){
	mylight.z[i*mylight.Nx+j] = light_collection.all_values(mylight.center_x[j],mylight.center_y[i]);
      }
    }

    // Super-resolved lens light profile image
    std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
    std::vector<std::string> values{std::to_string(mylight.xmin),std::to_string(mylight.xmax),std::to_string(mylight.ymin),std::to_string(mylight.ymax)};
    std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
    vkl::FitsInterface::writeFits(mylight.Nx,mylight.Ny,mylight.z,keys,values,descriptions,output + name + "_lens_light_super.fits");
  }
  //================= END:CREATE LENS LIGHT ================






  //=============== BEGIN:CREATE LENS COMPACT MASS ================
  if( root.isMember("point_source") ){
    // Factor to convert surface mass density to kappa
    double Dl  = cosmo[0]["Dl"].asDouble();
    double Ds  = cosmo[0]["Ds"].asDouble();
    double Dls = cosmo[0]["Dls"].asDouble();
    double sigma_crit = 3472.8*Ds/(Dl*Dls); // the critical density: c^2/(4pi G)  Ds/(Dl*Dls), in units of kg/m^2
    
    Json::Value all_compact; // Is not necessarily the same size as 'lenses' because one lens can have several 'compact_mass_model' profiles
    std::vector<std::string> flag;
    std::vector<double> zps;
    std::vector<double> areas;
    for(int i=0;i<root["lenses"].size();i++){
      if( root["lenses"][i].isMember("compact_mass_model") ){
	for(int m=0;m<root["lenses"][i]["compact_mass_model"].size();m++){
	  all_compact.append( root["lenses"][i]["compact_mass_model"][m] );
	  flag.push_back("direct_mass");
	  zps.push_back(0.0); // dummy, not needed
	  areas.push_back(0.0); // dummy, not needed
	}
      } else {
	for(int k=0;k<root["instruments"].size();k++){
	  std::string name = root["instruments"][k]["name"].asString();	
	  for(int m=0;m<root["lenses"][i]["light_profile"].size();m++){
	    if( root["lenses"][i]["light_profile"][name][m].isMember("mass-to-light") ){
	      all_compact.append( root["lenses"][i]["light_profile"][name][m] );
	      flag.push_back("mass-to-light");
	      zps.push_back(root["instruments"][k]["ZP"].asDouble());
	      double dum = Instrument::getResolution(name);
	      areas.push_back(dum*dum);
	    }
	  }
	}
      }
    }
    vkl::CollectionProfiles compact_collection = vkl::JsonParsers::parse_profile(all_compact,zps,input); // These are all LightProfile types (Sersic, Gauss, Custom)

    // Create overall kappa_star field
    double xmin,xmax,ymin,ymax;
    compact_collection.getExtent(xmin,xmax,ymin,ymax);
    vkl::RectGrid kappa_star(300,300,xmin,xmax,ymin,ymax);
    for(int i=0;i<kappa_star.Ny;i++){
      for(int j=0;j<kappa_star.Nx;j++){
	double value = 0.0;
	for(int m=0;m<compact_collection.profiles.size();m++){
	  if( flag[m] == "direct_mass" ){
	    value += compact_collection.profiles[m]->value(kappa_star.center_x[j],kappa_star.center_y[i]);
	  } else {
	    value += areas[m]*compact_collection.profiles[m]->value_to_mass(kappa_star.center_x[j],kappa_star.center_y[i]);
	  }
	}
	kappa_star.z[i*kappa_star.Nx+j] = value/sigma_crit;
      }
    }
    // Super-resolved lens compact mass profile image
    std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
    std::vector<std::string> values{std::to_string(kappa_star.xmin),std::to_string(kappa_star.xmax),std::to_string(kappa_star.ymin),std::to_string(kappa_star.ymax)};
    std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
    vkl::FitsInterface::writeFits(kappa_star.Nx,kappa_star.Ny,kappa_star.z,keys,values,descriptions,output + "lens_kappa_star_super.fits");      


    // Read the image parameters
    Json::Value images;
    fin.open(output+"multiple_images.json",std::ifstream::in);
    fin >> images;
    fin.close();

    // Find the kappa_star at the multiple images
    for(int j=0;j<images.size();j++){
      double x = images[j]["x"].asDouble();
      double y = images[j]["y"].asDouble();
      double value = 0.0;
      for(int m=0;m<compact_collection.profiles.size();m++){
	if( flag[m] == "direct_mass" ){
	  value += compact_collection.profiles[m]->value(x,y);
	} else {
	  value += areas[m]*compact_collection.profiles[m]->value_to_mass(x,y);
	}
      }
      double k_star = value/sigma_crit;
      double s =  1.0 - k_star/images[j]["k"].asDouble();
      std::cout << images[j]["k"].asDouble() << " " << k_star << " " << s << std::endl;
      if( s < 0 ){
	images[j]["s"] = 0.0;
      } else {
	images[j]["s"] = s;	
      }
    }
        
    // Overwrite multiple images file with values of s
    std::ofstream file_images(output+"multiple_images.json");
    file_images << images;
    file_images.close();
  }  
  //================= END:CREATE LENS COMPACT MASS ================


  return 0;
}
