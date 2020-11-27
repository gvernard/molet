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
  double width  = root["instruments"][0]["field-of-view_x"].asDouble();
  double height = root["instruments"][0]["field-of-view_y"].asDouble();
  double resolution = Instrument::getResolution(root["instruments"][0]["name"].asString());
  int super_res_x = 10*( static_cast<int>(ceil(width/resolution)) );
  int super_res_y = 10*( static_cast<int>(ceil(height/resolution)) );
  ImagePlane mysim(super_res_x,super_res_y,width,height);
  //================= END:PARSE INPUT =======================




  //=============== BEGIN:CREATE THE LENSES ====================
  const Json::Value jlens = root["lenses"][0];

  // Initialize mass model physical parameters
  jmembers = jlens["external_shear"].getMemberNames();
  
  std::vector<Nlpar*> ext_pars;
  for(int i=0;i<jmembers.size();i++){
    double value = jlens["external_shear"][jmembers[i]].asDouble();
    ext_pars.push_back( new Nlpar(jmembers[i],0,0,value,0,0,0) );
  }
  CollectionMassModels* mycollection = new CollectionMassModels(ext_pars);
  for(int i=0;i<ext_pars.size();i++){ delete(ext_pars[i]); }
  ext_pars.clear();

  // Initialize main mass model
  mycollection->models.resize(jlens["mass_model"].size());
  for(int k=0;k<jlens["mass_model"].size();k++){
    std::string mmodel = jlens["mass_model"][k]["type"].asString();

    if( mmodel == "custom" ){

      std::string filename = input + jlens["mass_model"][k]["pars"]["filename"].asString();
      int dpsi_Ni = jlens["mass_model"][k]["pars"]["Ni"].asInt();
      int dpsi_Nj = jlens["mass_model"][k]["pars"]["Nj"].asInt();
      double dpsi_width,dpsi_height;
      if( jlens["mass_model"][k]["pars"].isMember("width") ){
	dpsi_width = jlens["mass_model"][k]["pars"]["width"].asDouble();
      } else {
	dpsi_width = width;
      }
      if( jlens["mass_model"][k]["pars"].isMember("height") ){
	dpsi_height = jlens["mass_model"][k]["pars"]["height"].asDouble();
      } else {
	dpsi_height = height;
      }
      std::string reg = "identity"; // dummy argument
      
      Pert* custom = new Pert(filename,dpsi_Ni,dpsi_Nj,dpsi_width,dpsi_height,reg);

      if( jlens["mass_model"][k]["pars"].isMember("scale_factor") ){
	double scale_factor = jlens["mass_model"][k]["pars"]["scale_factor"].asDouble();
	for(int m=0;m<custom->dpsi->Sm;m++){
	  custom->dpsi->src[m] *= scale_factor;
	}
	custom->updatePert();
      }
      
      mycollection->models[k] = custom;
      
    } else if ( mmodel == "eagle" ){
      
    } else {
      
      jmembers = jlens["mass_model"][k]["pars"].getMemberNames();
      std::vector<Nlpar*> pars;
      for(int i=0;i<jmembers.size();i++){
	pars.push_back( new Nlpar(jmembers[i],0,0,jlens["mass_model"][k]["pars"][jmembers[i]].asDouble(),0,0,0) ); // only nam and val have meaning in this call
      }
      mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(mmodel,pars,cosmo[0]["Dls"].asDouble(),cosmo[0]["Ds"].asDouble());
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
    int Ni               = jsource["pars"]["Ni"].asInt();
    int Nj               = jsource["pars"]["Nj"].asInt();
    double height        = jsource["pars"]["height"].asDouble();
    double width         = jsource["pars"]["width"].asDouble();
    double x0            = jsource["pars"]["x0"].asDouble();
    double y0            = jsource["pars"]["y0"].asDouble();
    double Mtot          = jsource["pars"]["M_tot"].asDouble();
    mysource = new fromFITS(filename,Ni,Nj,height,width,x0,y0,Mtot,"bilinear");
    
  } else {

    std::cout << "Unknown source profile type" << std::endl;
    return 1;

  }
  //================= END:CREATE THE SOURCES =======================




  
  //=============== BEGIN:PRODUCE IMAGE USING RAY-SHOOTING =======================
  double xdefl,ydefl;
  for(int i=0;i<mysim.Nm;i++){
    mycollection->all_defl(mysim.x[i],mysim.y[i],xdefl,ydefl);
    mysim.img[i] = mysource->value(xdefl,ydefl);
  }
  //================= END:PRODUCE IMAGE USING RAY-SHOOTING =======================



  
  
  //=============== BEGIN:OUTPUT =======================
  // Super-resolved lensed image
  mysim.writeImage(output + "lensed_image_super.fits");
  // Super-resolved source image
  mysource->outputProfile(output + "source_super.fits");
  //================= END:OUTPUT =======================





  //=============== BEGIN:CLEAN UP =======================
  delete(mycollection);
  delete(mysource);
  //================= END:CLEAN UP =======================


  return 0;
}
