#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "vkllib.hpp"


int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  std::ifstream fin;
  Json::Value::Members jmembers;

  // Read the main parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();
  
  // Read the cosmological parameters
  Json::Value cosmo;
  fin.open(argv[2],std::ifstream::in);
  fin >> cosmo;
  fin.close();

  // Read the image parameters
  Json::Value images;
  fin.open(argv[3],std::ifstream::in);
  fin >> images;
  fin.close();
  
  std::string output = argv[4];
  
  // Initialize image plane
  const Json::Value jiplane = root["instrument"]["bands"][0];
  double width  = jiplane["field-of-view_x"].asDouble();
  double height = jiplane["field-of-view_y"].asDouble();
  int super_res_x = 10*( static_cast<int>(ceil(width/jiplane["resolution"].asDouble())) );
  int super_res_y = 10*( static_cast<int>(ceil(height/jiplane["resolution"].asDouble())) );
  //================= END:PARSE INPUT =======================







  //=============== BEGIN:CREATE LENS LIGHT =======================
  const Json::Value jll = root["lenses"][0]["light_profile"];
  BaseProfile* lens_light = NULL;

  std::string light_model = jll["type"].asString();
  if( light_model == "analytic" ){

    std::vector<std::string> names;
    std::vector<std::map<std::string,double> > all_pars;
    for(int i=0;i<jll["pars"].size();i++){
      std::string name = jll["pars"][i]["type"].asString();
      names.push_back(name);

      std::map<std::string,double> pars;
      const Json::Value::Members jpars = jll["pars"][i].getMemberNames();
      for(int j=0;j<jpars.size();j++){
	if( jpars[j] != "type" ){
	  pars[jpars[j]] = jll["pars"][i][jpars[j]].asDouble();
	}
	if( jpars[j] == "pa" ){
	  pars[jpars[j]] = jll["pars"][i][jpars[j]].asDouble() + 90.0;
	}
      }
      all_pars.push_back(pars);
    }
    lens_light = new Analytic(names,all_pars);

  } else if( light_model == "custom" ){

    std::string filename = jll["pars"]["filename"].asString();
    int Ni               = jll["pars"]["Ni"].asInt();
    int Nj               = jll["pars"]["Nj"].asInt();
    double height        = jll["pars"]["height"].asDouble();
    double width         = jll["pars"]["width"].asDouble();
    double x0            = jll["pars"]["x0"].asDouble();
    double y0            = jll["pars"]["y0"].asDouble();
    lens_light = new fromFITS(filename,Ni,Nj,height,width,x0,y0);

  } else {

    std::cout << "Unknown light profile type" << std::endl;
    return 1;

  }

  ImagePlane mylight(super_res_x,super_res_y,width,height);
  for(int i=0;i<mylight.Nm;i++){
    mylight.img[i] = lens_light->value(mylight.x[i],mylight.y[i]);
  }

  // Super-resolved lens light profile image
  mylight.writeImage(output + "lens_light_super.fits");
  //================= END:CREATE LENS LIGHT ================







  //=============== BEGIN:CREATE LENS COMPACT MASS ================
  if( root["lenses"][0].isMember("compact_mass_model") ){
    // Factor to convert surface mass density to kappa
    double Dl  = cosmo[0]["Dl"].asDouble();
    double Ds  = cosmo[0]["Ds"].asDouble();
    double Dls = cosmo[0]["Dls"].asDouble();
    double sigma_crit = 3472.8*Ds/(Dl*Dls); // the critical density: c^2/(4pi G)  Ds/(Dl*Dls), in units of kg/m^2

    
    BaseProfile* lens_compact = NULL;
    const Json::Value jlm = root["lenses"][0]["compact_mass_model"];
    std::string compact_model = jlm["type"].asString();
    if( compact_model == "mass_to_light" ){

      lens_compact = lens_light;
      double ratio = jlm["pars"]["ratio"].asDouble(); // in units of Y_solar
      ratio *= 5133; // in units of kg/m^2 per W
      sigma_crit /= ratio; // in units of W kg/m^2, because I will use it to divide a flux (mass-to-light ratio) as opposed to the other purely mass profiles below
      
    } else if( compact_model == "analytic" ){
      
      std::vector<std::string> names;
      std::vector<std::map<std::string,double> > all_pars;
      for(int i=0;i<jlm["pars"].size();i++){
	std::string name = jlm["pars"][i]["type"].asString();
	names.push_back(name);
	
	std::map<std::string,double> pars;
	const Json::Value::Members jpars = jlm["pars"][i].getMemberNames();
	for(int j=0;j<jpars.size();j++){
	  if( jpars[j] != "type" ){
	    pars[jpars[j]] = jlm["pars"][i][jpars[j]].asDouble();
	  }
	}
	all_pars.push_back(pars);
      }
      lens_compact = new Analytic(names,all_pars);
      
    } else if( compact_model == "custom" ){
      
      std::string filename = jlm["pars"]["filename"].asString();
      int Ni               = jlm["pars"]["Ni"].asInt();
      int Nj               = jlm["pars"]["Nj"].asInt();
      double height        = jlm["pars"]["height"].asDouble();
      double width         = jlm["pars"]["width"].asDouble();
      double x0            = jlm["pars"]["x0"].asDouble();
      double y0            = jlm["pars"]["y0"].asDouble();
      lens_compact = new fromFITS(filename,Ni,Nj,height,width,x0,y0);
      
    } else {
      
      std::cout << "Unknown compact mass profile type" << std::endl;
      return 1;
      
    }


    // Write overall kappa_star field
    ImagePlane kappa_star(super_res_x,super_res_y,width,height);
    for(int i=0;i<kappa_star.Nm;i++){
      kappa_star.img[i] = lens_compact->value(kappa_star.x[i],kappa_star.y[i])/sigma_crit;
    }
    // Super-resolved lens compact mass profile image
    kappa_star.writeImage(output + "lens_kappa_star_super.fits");


    // Find the kappa_star at the multiple images
    for(int j=0;j<images.size();j++){
      double x = images[j]["x"].asDouble();
      double y = images[j]["y"].asDouble();
      double kappa_star = lens_compact->value(x,y)/sigma_crit;
      images[j]["s"] = 1.0 - kappa_star/images[j]["k"].asDouble();
    }

    // Overwrite multiple images file with values of s
    std::ofstream file_images(output+"multiple_images.json");
    file_images << images;
    file_images.close();

    if( compact_model != "mass_to_light" ){
      delete(lens_compact);
    }
  }  
  //================= END:CREATE LENS COMPACT MASS ================



  









  


  delete(lens_light); // May be needed by the compact mass profile if set to 'mass_to_light'


  return 0;
}
