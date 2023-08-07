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
  
  double xdefl,ydefl;
  //================= END:PARSE INPUT =======================
  


  //=============== BEGIN:CREATE THE LENSES ====================
  Json::Value jlens = root["lenses"][0];

  // Initialize mass model collection
  vkl::CollectionMassModels mass_collection = vkl::JsonParsers::parse_mass_model(jlens["mass_model"],input);
  
  // Scale dpsi mass models if necessary
  for(int k=0;k<jlens["mass_model"].size();k++){
    if( jlens["mass_model"][k]["pars"].isMember("scale_factor") ){
      double scale_factor = jlens["mass_model"][k]["pars"]["scale_factor"].asDouble();
      vkl::Pert* pert = static_cast<vkl::Pert*> (mass_collection.models[k]);
      for(int m=0;m<pert->Sm;m++){
	pert->z[m] *= scale_factor;
      }
      pert->updateDerivatives();
    }
  }
  //================= END:CREATE THE LENSES ====================
  


  //=============== BEGIN:GET CRITICAL LINES AND CAUSTICS =======================
  // Calculate detA
  double xmin,xmax,ymin,ymax;
  mass_collection.getExtent(xmin,xmax,ymin,ymax);
  vkl::RectGrid detA(300,300,xmin,xmax,xmin,xmax);
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
  
  // Check nonzero contours and add first point as last to close the polygon
  std::vector<Contour> clean_contours;
  for(int i=0;i<contours.size();i++){
    if( contours[i].x.size() > 4 ){
      contours[i].x.push_back( contours[i].x[0] );
      contours[i].y.push_back( contours[i].y[0] );
      clean_contours.push_back( contours[i]);
    }
  }
  contours.clear();
  contours = clean_contours;
  clean_contours.clear();

  // Create caustic contours, and deflect the critical contours to fill them
  std::vector<Contour> caustics(contours.size());
  for(int k=0;k<contours.size();k++){
    caustics[k] = contours[k];
    for(int i=0;i<contours[k].x.size();i++){
      mass_collection.all_defl(contours[k].x[i],contours[k].y[i],xdefl,ydefl);
      caustics[k].x[i] = xdefl;
      caustics[k].y[i] = ydefl;
    }
  } 

  // Image plane magnification (0:positive, 1:negative), caustics, and critical curves
  std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
  std::vector<std::string> values{std::to_string(detA.xmin),std::to_string(detA.xmax),std::to_string(detA.ymin),std::to_string(detA.ymax)};
  std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
  vkl::FitsInterface::writeFits(detA.Nx,detA.Ny,detA.z,keys,values,descriptions,output + "detA.fits");
  outputContours(contours,output+"criticals.json");
  outputContours(caustics,output+"caustics.json");
  //================= END:GET CRITICAL LINES AND CAUSTICS =======================

  

  //=============== BEGIN:LENS THE SOURCES =======================
  Json::Value lensed_flux;
  Json::Value unlensed_flux;
  for(int k=0;k<root["instruments"].size();k++){
    std::string name = root["instruments"][k]["name"].asString();
    
    // Initialize image plane
    double xmin = root["instruments"][k]["field-of-view_xmin"].asDouble();
    double xmax = root["instruments"][k]["field-of-view_xmax"].asDouble();
    double ymin = root["instruments"][k]["field-of-view_ymin"].asDouble();
    double ymax = root["instruments"][k]["field-of-view_ymax"].asDouble();
    double resolution = Instrument::getResolution(name);
    int super_factor = 10;
    if( root["output_options"].isMember("super_factor") ){
      super_factor = root["output_options"]["super_factor"].asInt();
    }
    int super_res_x = super_factor*( static_cast<int>(ceil((xmax-xmin)/resolution)) );
    int super_res_y = super_factor*( static_cast<int>(ceil((ymax-ymin)/resolution)) );
    vkl::RectGrid mysim(super_res_x,super_res_y,xmin,xmax,ymin,ymax);


    
    // Create the source and output super-resolved image
    vkl::CollectionProfiles profile_collection = vkl::JsonParsers::parse_profile(root["source"]["light_profile"][name],root["instruments"][k]["ZP"].asDouble(),input);
    profile_collection.write_all_profiles(output + name + "_source_super.fits");


    
    // Produce image using ray-shooting
    for(int i=0;i<mysim.Ny;i++){
      for(int j=0;j<mysim.Nx;j++){
	mass_collection.all_defl(mysim.center_x[j],mysim.center_y[i],xdefl,ydefl);
	mysim.z[i*mysim.Nx+j] = profile_collection.all_values(xdefl,ydefl);
      }
    }
          
    // Output super-resolved lensed image
    std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
    std::vector<std::string> values{std::to_string(mysim.xmin),std::to_string(mysim.xmax),std::to_string(mysim.ymin),std::to_string(mysim.ymax)};
    std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
    vkl::FitsInterface::writeFits(mysim.Nx,mysim.Ny,mysim.z,keys,values,descriptions,output + name + "_lensed_image_super.fits");



    // Calculate total unlensed flux per profile and in total.
    double total_flux = 0.0;
    double prof_flux,prof_flux_mag;
    Json::Value profiles = Json::Value(Json::arrayValue);
    for(int i=0;i<profile_collection.profiles.size();i++){      
      prof_flux = profile_collection.profiles[i]->integrate(xmin,xmax,ymin,ymax,500);
      prof_flux_mag = -2.5*log10(prof_flux) + root["instruments"][k]["ZP"].asDouble();
      total_flux += prof_flux;

      Json::Value profile;
      profile["flux"]  = prof_flux;
      profile["mag"]   = prof_flux_mag;
      profile["type"]  = profile_collection.profiles[i]->profile_type;
      profile["index"] = i;
      profiles.append(profile);
    }
    double total_flux_mag = -2.5*log10(total_flux) + root["instruments"][k]["ZP"].asDouble();

    // Calculate total lensed flux
    double total_lensed_flux = 0.0;
    double total_lensed_flux_mag;
    mysim.integrate(total_lensed_flux,total_lensed_flux_mag,root["instruments"][k]["ZP"].asDouble());

    unlensed_flux[name]["profiles"]      = profiles;
    unlensed_flux[name]["total"]["flux"] = total_flux;
    unlensed_flux[name]["total"]["mag"]  = total_flux_mag;    
    lensed_flux[name]["total"]["flux"]   = total_lensed_flux;
    lensed_flux[name]["total"]["mag"]    = total_lensed_flux_mag;

    
    
    /*
    // Output deflections
    vkl::RectGrid defl_x(super_res_x/super_factor,super_res_y/super_factor,xmin,xmax,ymin,ymax);
    vkl::RectGrid defl_y(super_res_x/super_factor,super_res_y/super_factor,xmin,xmax,ymin,ymax);
    for(int i=0;i<defl_x.Ny;i++){
      for(int j=0;j<defl_x.Nx;j++){
	mass_collection.all_defl(defl_x.center_x[j],defl_x.c    double sum = 0.0;
    double area = (mysim.center_x[1]-mysim.center_x[0])*(mysim.center_y[1] - mysim.center_y[0]);
    for(int i=0;i<mysim.Nz;i++){
      sum += mysim.z[i];
    }
    sum *= area;
    double sum_mag = -2.5*log10(sum) + root["instruments"][k]["ZP"].asDouble();enter_y[i],xdefl,ydefl);
	defl_x.z[i*defl_x.Nx+j] = xdefl;
	defl_y.z[i*defl_x.Nx+j] = ydefl;
      }
    }
    vkl::FitsInterface::writeFits(defl_x.Nx,defl_x.Ny,defl_x.z,keys,values,descriptions,output + name + "_defl_x.fits");
    vkl::FitsInterface::writeFits(defl_y.Nx,defl_y.Ny,defl_y.z,keys,values,descriptions,output + name + "_defl_y.fits");

    double xxx[4] = {-1.75,1.715,-1.75,1.715};
    double yyy[4] = {1.75,1.75,-1.715,-1.715};
    for(int k=0;k<4;k++){
      mass_collection.all_defl(xxx[k],yyy[k],xdefl,ydefl);
      printf("%10.4f %10.4f - %10.4f %10.4f\n",xxx[k],yyy[k],xdefl,ydefl);
    }
    mass_collection.models[0]->printMassPars();
    */
    
    
  }
  //================= END:LENS THE SOURCES =======================



  // Write fluxes
  Json::Value fluxes;
  fluxes["lensed_source_flux"] = lensed_flux;
  fluxes["unlensed_source_flux"] = unlensed_flux;
  std::ofstream file_fluxes(output+"fluxes.json");
  file_fluxes << fluxes;
  file_fluxes.close();

  
  return 0;
}
