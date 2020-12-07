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
  const Json::Value jlens = root["lenses"][0];

  // Initialize mass model physical parameters
  CollectionMassModels* mycollection = new CollectionMassModels();

  // Initialize main mass model
  mycollection->models.resize(jlens["mass_model"].size());
  for(int k=0;k<jlens["mass_model"].size();k++){
    std::string mmodel = jlens["mass_model"][k]["type"].asString();

    if( mmodel == "custom" ){

      std::string filename = input + jlens["mass_model"][k]["pars"]["filename"].asString();
      int dpsi_Ny = jlens["mass_model"][k]["pars"]["Ny"].asInt();
      int dpsi_Nx = jlens["mass_model"][k]["pars"]["Nx"].asInt();
      Pert* custom = new Pert(dpsi_Nx,dpsi_Ny,mysim.xmin,mysim.xmax,mysim.ymin,mysim.ymax,filename);

      if( jlens["mass_model"][k]["pars"].isMember("scale_factor") ){
	double scale_factor = jlens["mass_model"][k]["pars"]["scale_factor"].asDouble();
	for(int m=0;m<custom->Sm;m++){
	  custom->z[m] *= scale_factor;
	}
	custom->updateDerivatives();
      }
      
      mycollection->models[k] = custom;
      
    } else if ( mmodel == "eagle" ){
      
    } else {
      
      jmembers = jlens["mass_model"][k]["pars"].getMemberNames();
      std::map<std::string,double> pars;
      for(int i=0;i<jmembers.size();i++){
	pars.insert( std::pair<std::string,double>(jmembers[i],jlens["mass_model"][k]["pars"][jmembers[i]].asDouble()) );
      }
      mycollection->models[k] = FactoryParametricMassModel::getInstance()->createParametricMassModel(mmodel,pars);
    }
  }

  

  //  for(int i=0;i<mycollection->models.size();i++){
  //    mycollection->models[i]->printMassPars();
  //  }
  //  mycollection->printPhysPars();
  //================= END:CREATE THE LENSES ====================






  //=============== BEGIN:CREATE THE SOURCES =======================
  const Json::Value jsource = root["source"]["light_profile"];
  BaseProfile* mysource = NULL;

  std::string smodel = jsource["type"].asString();
  if( smodel == "analytic" ){

    std::vector<std::string> names;
    std::vector<std::map<std::string,double> > all_pars;
    for(int i=0;i<jsource["pars"].size();i++){
      std::string name = jsource["pars"][i]["type"].asString();
      names.push_back(name);

      std::map<std::string,double> pars;
      const Json::Value::Members jpars = jsource["pars"][i].getMemberNames();
      for(int j=0;j<jpars.size();j++){
	if( jpars[j] != "type" ){
	  pars[jpars[j]] = jsource["pars"][i][jpars[j]].asDouble();
	}
      }
      all_pars.push_back(pars);
    }
    mysource = new Analytic(names,all_pars);

  } else if( smodel == "delaunay" ){

    std::string filename = jsource["pars"]["filename"].asString();
    mysource = new myDelaunay(filename);

  } else if( smodel == "custom" ){

    std::string filename = input + jsource["pars"]["filename"].asString();
    int Nx               = jsource["pars"]["Nx"].asInt();
    int Ny               = jsource["pars"]["Ny"].asInt();
    double xmin          = jsource["pars"]["xmin"].asDouble();
    double xmax          = jsource["pars"]["xmax"].asDouble();
    double ymin          = jsource["pars"]["ymin"].asDouble();
    double ymax          = jsource["pars"]["ymax"].asDouble();
    double Mtot          = jsource["pars"]["M_tot"].asDouble();
    mysource = new Custom(filename,Nx,Ny,xmin,xmax,ymin,ymax,Mtot,"bilinear");
    
  } else {

    std::cout << "Unknown source profile type" << std::endl;
    return 1;

  }
  //================= END:CREATE THE SOURCES =======================





  
  //=============== BEGIN:PRODUCE IMAGE USING RAY-SHOOTING =======================
  double xdefl,ydefl;
  for(int i=0;i<mysim.Ny;i++){
    for(int j=0;j<mysim.Nx;j++){
      mycollection->all_defl(mysim.center_x[j],mysim.center_y[i],xdefl,ydefl);
      mysim.z[i*mysim.Nx+j] = mysource->value(xdefl,ydefl);
    }
  }
  //================= END:PRODUCE IMAGE USING RAY-SHOOTING =======================
  
  

  
  
  
  //=============== BEGIN:OUTPUT =======================
  // Super-resolved lensed image
  FitsInterface::writeFits(mysim.Nx,mysim.Ny,mysim.z,output + "lensed_image_super.fits");
  // Super-resolved source image
  mysource->outputProfile(output + "source_super.fits");
  //================= END:OUTPUT =======================






  //=============== BEGIN:CLEAN UP =======================
  delete(mycollection);
  delete(mysource);
  //================= END:CLEAN UP =======================


  return 0;
}
