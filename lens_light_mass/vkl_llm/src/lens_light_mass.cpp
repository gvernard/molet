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
  //================= END:PARSE INPUT =======================







  //=============== BEGIN:CREATE LENS LIGHT =======================
  Json::Value lens_fluxes;
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
    std::vector<std::string> lens_profile_ind; // Needed to identify each profile with a lens when storing its total flux
    for(int i=0;i<root["lenses"].size();i++){
      for(int j=0;j<root["lenses"][i]["light_profile"][name].size();j++){
	all_lenses.append(root["lenses"][i]["light_profile"][name][j]);
	lens_profile_ind.push_back( std::to_string(i) + "-" + std::to_string(j) );
      }
    }
    vkl::CollectionProfiles light_collection = vkl::JsonParsers::parse_profile(all_lenses,root["instruments"][k]["ZP"].asDouble(),input);

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


    // Calculate total flux per profile and in total
    double total_flux = 0.0;
    double ZP = root["instruments"][k]["ZP"].asDouble();
    Json::Value profiles = Json::Value(Json::arrayValue);
    for(int i=0;i<light_collection.profiles.size();i++){      
      double prof_flux = light_collection.profiles[i]->integrate(xmin,xmax,ymin,ymax,500);
      total_flux += prof_flux;

      Json::Value profile;
      profile["flux"]  = prof_flux;
      profile["mag"]   = -2.5*log10(prof_flux) + ZP;
      profile["type"]  = light_collection.profiles[i]->profile_type;
      profile["index"] = lens_profile_ind[i];
      profiles.append(profile);
    }
    double total_flux_mag = 0.0;
    if( total_flux != 0.0 ){
      total_flux_mag = -2.5*log10(total_flux) + ZP;
    }

    lens_fluxes[name]["profiles"]      = profiles;
    lens_fluxes[name]["total"]["flux"] = total_flux;
    lens_fluxes[name]["total"]["mag"]  = total_flux_mag;
  }


  // Write fluxes
  Json::Value fluxes;
  fin.open(output+"fluxes.json",std::ifstream::in);
  if( fin.fail() ){
    fluxes["lens_flux"] = lens_fluxes;
    std::ofstream file_fluxes(output+"fluxes.json");
    file_fluxes << fluxes;
    file_fluxes.close();
  } else {
    // Overwrite fluxes after adding the lens fluxes
    Json::Value fluxes;
    fin >> fluxes;
    fluxes["lens_flux"] = lens_fluxes;
    std::ofstream file_fluxes(output+"fluxes.json");
    file_fluxes << fluxes;
    file_fluxes.close();
  }
  fin.close();
  //================= END:CREATE LENS LIGHT ================






  //=============== BEGIN:CREATE LENS COMPACT MASS ================
  if( root.isMember("point_source") ){
    Json::Value all_compact; // Is not necessarily the same size as 'lenses' because one lens can have several 'compact_mass_model' profiles
    std::vector<std::string> flag;
    std::vector<double> zps; // zero-points
    for(int i=0;i<root["lenses"].size();i++){
      if( root["lenses"][i].isMember("compact_mass_model") ){
	for(int m=0;m<root["lenses"][i]["compact_mass_model"].size();m++){
	  all_compact.append( root["lenses"][i]["compact_mass_model"][m] );
	  flag.push_back("direct_mass");
	  zps.push_back(0.0); // dummy
	}
      } else {
	for(int k=0;k<root["instruments"].size();k++){
	  std::string name = root["instruments"][k]["name"].asString();	
	  for(int m=0;m<root["lenses"][i]["light_profile"].size();m++){
	    if( root["lenses"][i]["light_profile"][name][m].isMember("mass-to-light") ){
	      all_compact.append( root["lenses"][i]["light_profile"][name][m] );
	      flag.push_back("mass-to-light");
	      zps.push_back(root["instruments"][k]["ZP"].asDouble());
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
    // double xmin = root["instruments"][0]["field-of-view_xmin"].asDouble();
    // double xmax = root["instruments"][0]["field-of-view_xmax"].asDouble();
    // double ymin = root["instruments"][0]["field-of-view_ymin"].asDouble();
    // double ymax = root["instruments"][0]["field-of-view_ymax"].asDouble();
    // vkl::RectGrid kappa_star(500,500,xmin,xmax,ymin,ymax);
    for(int i=0;i<kappa_star.Ny;i++){
      for(int j=0;j<kappa_star.Nx;j++){
	double value = 0.0;
	for(int m=0;m<compact_collection.profiles.size();m++){
	  if( flag[m] == "direct_mass" ){
	    value += compact_collection.profiles[m]->value(kappa_star.center_x[j],kappa_star.center_y[i]);
	  } else {
	    value += compact_collection.profiles[m]->value_to_mass(kappa_star.center_x[j],kappa_star.center_y[i]);
	  }
	}
	kappa_star.z[i*kappa_star.Nx+j] = value;
      }
    }
    // Super-resolved lens compact mass profile image
    std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
    std::vector<std::string> values{std::to_string(kappa_star.xmin),std::to_string(kappa_star.xmax),std::to_string(kappa_star.ymin),std::to_string(kappa_star.ymax)};
    std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
    vkl::FitsInterface::writeFits(kappa_star.Nx,kappa_star.Ny,kappa_star.z,keys,values,descriptions,output + "lens_kappa_star_super.fits");      

    // double total_kappa_star,dummy;
    // kappa_star.integrate(total_kappa_star,dummy,0.0);
    // std::cout << "Total k_star: " << total_kappa_star << std::endl;
    
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
	  value += compact_collection.profiles[m]->value_to_mass(x,y);
	}
      }
      double k_star = value;
      double s =  1.0 - k_star/images[j]["k"].asDouble();
      double kappa = images[j]["k"].asDouble();
      images[j]["k_star"] = k_star;
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
