#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "vkllib.hpp"


using std::cout;
using std::endl;




double myRandomNumber(double min,double max){
  double r = (double) rand() / (double)RAND_MAX;
  return min + r * (max - min);
}


int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  //Read in the json parameters
  std::ifstream fin;
  Json::Value::Members jmembers;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();
  
  // Read the cosmological parameters
  Json::Value cosmo;
  fin.open(argv[2],std::ifstream::in);
  fin >> cosmo;
  fin.close();

  std::string output = argv[3];
  
  // Initialize image plane
  const Json::Value jiplane = root["instrument"]["bands"][0];
  double width  = jiplane["field-of-view_x"].asDouble();
  double height = jiplane["field-of-view_y"].asDouble();
  int super_res_x = 10*( static_cast<int>(ceil(width/jiplane["resolution"].asDouble())) );
  int super_res_y = 10*( static_cast<int>(ceil(height/jiplane["resolution"].asDouble())) );
  ImagePlane mysim(super_res_x,super_res_y,width,height);
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
  //================= END:CREATE LENS LIGHT =======================





  //=============== BEGIN:OUTPUT =======================
  for(int i=0;i<mysim.Nm;i++){
    mysim.img[i] = lens_light->value(mysim.x[i],mysim.y[i]);
  }


  // Super-resolved lens light profile image
  mysim.writeImage(output + "lens_light_super.fits");
  //================= END:OUTPUT =======================





  //=============== BEGIN:CLEAN UP =======================
  delete(lens_light);
  //================= END:CLEAN UP =======================


  return 0;
}
