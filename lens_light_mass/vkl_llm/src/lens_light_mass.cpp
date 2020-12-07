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

  std::string out_path = argv[2];
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
  double res    = Instrument::getResolution(root["instruments"][0]["name"].asString());
  int super_res_x = 10*( static_cast<int>(ceil((xmax-xmin)/res)) );
  int super_res_y = 10*( static_cast<int>(ceil((ymax-ymin)/res)) );
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
    int Nx               = jll["pars"]["Nx"].asInt();
    int Ny               = jll["pars"]["Ny"].asInt();
    double xmin          = jll["pars"]["xmin"].asDouble();
    double xmax          = jll["pars"]["xmax"].asDouble();
    double ymin          = jll["pars"]["ymin"].asDouble();
    double ymax          = jll["pars"]["ymax"].asDouble();
    double Mtot          = jll["pars"]["M_tot"].asDouble();
    lens_light = new Custom(filename,Nx,Ny,xmin,xmax,ymin,ymax,Mtot,"bilinear");

  } else {

    std::cout << "Unknown light profile type" << std::endl;
    return 1;

  }

  RectGrid mylight(super_res_x,super_res_y,xmin,xmax,ymin,ymax);
  for(int i=0;i<mylight.Ny;i++){
    for(int j=0;j<mylight.Nx;j++){
      mylight.z[i*mylight.Nx+j] = lens_light->value(mylight.center_x[j],mylight.center_y[i]);
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

    
    BaseProfile* lens_compact = NULL;
    const Json::Value jlm = root["lenses"][0]["compact_mass_model"];
    std::string compact_model = jlm["type"].asString();
    if( compact_model == "mass_to_light" ){

      lens_compact = lens_light;
      double ratio = jlm["pars"]["ratio"].asDouble(); // in units of Y_solar
      ratio *= 5133; // in units of kg/W
      sigma_crit /= ratio; // now Sigma_crit is in W/m^2

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
      int Nx               = jlm["pars"]["Nx"].asInt();
      int Ny               = jlm["pars"]["Ny"].asInt();
      double xmin          = jlm["pars"]["xmin"].asDouble();
      double xmax          = jlm["pars"]["xmax"].asDouble();
      double ymin          = jlm["pars"]["ymin"].asDouble();
      double ymax          = jlm["pars"]["ymax"].asDouble();
      double Mtot          = jlm["pars"]["M_tot"].asDouble();
      lens_compact = new Custom(filename,Nx,Ny,xmin,xmax,ymin,ymax,Mtot,"bilinear");
      
    } else {
      
      std::cout << "Unknown compact mass profile type" << std::endl;
      return 1;
      
    }


    // Write overall kappa_star field
    RectGrid kappa_star(super_res_x,super_res_y,xmin,xmax,ymin,ymax);
    for(int i=0;i<kappa_star.Ny;i++){
      for(int j=0;j<kappa_star.Nx;j++){
	kappa_star.z[i*kappa_star.Nx+j] = lens_compact->value(kappa_star.center_x[j],kappa_star.center_y[i])/sigma_crit;
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
