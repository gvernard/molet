// Checks that the duration of the following time vectors is consistent (everything is in the observer's time frame):
// - intrinsic light curve time
// - observing time vector for each instrument
// To compare these, we also need the maximum time delay.

// This code outputs the final minimum and maximum observational time encompassing all instruments

#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <algorithm>

#include "json/json.h"
#include "CCfits/CCfits"

#include "vkllib.hpp"
#include "gerlumph.hpp"
#include "instruments.hpp"
#include "noise.hpp"


int main(int argc,char* argv[]){
  std::ifstream fin;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string in_path = argv[2];
  std::string out_path = argv[3];




  // Read the multiple images' parameters from JSON to get the maximum time delay
  Json::Value images;
  fin.open(out_path+"output/multiple_images.json",std::ifstream::in);
  fin >> images;
  fin.close();
  // Get maximum image time delay
  double td_max = 0.0;
  for(int q=0;q<images.size();q++){
    double td = images[q]["dt"].asDouble();
    if( td > td_max ){
      td_max = td;
    }
  }


  // Find the tmin and tmax across all instrument times.
  std::vector<double> tmins;
  std::vector<double> tmaxs;
  std::vector<double> lrest;
  std::vector<std::string> names;
  double zs = root["source"]["redshift"].asDouble();
  for(int b=0;b<root["instruments"].size();b++){
    int Ntime = root["instruments"][b]["time"].size();
    tmins.push_back( root["instruments"][b]["time"][0].asDouble() );
    tmaxs.push_back( root["instruments"][b]["time"][Ntime-1].asDouble() );
    names.push_back( root["instruments"][b]["name"].asString() );
    Instrument mycam(root["instruments"][b]["name"].asString(),root["instruments"][b]["noise"]);
    double lobs  = (mycam.lambda_min + mycam.lambda_max)/2.0;
    lrest.push_back( lobs/(1.0+zs) );
  }
  double tobs_max = tmaxs[0];
  double tobs_min = tmins[0];
  for(int k=0;k<tmins.size();k++){
    if( tobs_max < tmaxs[k] ){
      tobs_max = tmaxs[k];
    }
    if( tobs_min > tmins[k] ){
      tobs_min = tmins[k];
    }
  }
  

  

  bool check = false;
  


  std::string ex_type = root["point_source"]["variability"]["extrinsic"]["type"].asString();
  std::string in_type = root["point_source"]["variability"]["intrinsic"]["type"].asString();

  
  // Check intrinsic light curves
  //================================================================================================================
  if( ex_type == "moving_fixed_source" || ex_type == "moving_fixed_source_custom" ){

    if( in_type == "custom" ){
      // Check if all intrinsic light curves begin at t < tobs_min - td_max
      // and that the maximum observational time is at least equal to the maximum intrinsic time (in the observer's frame).
      for(int n=0;n<names.size();n++){
	std::string iname = names[n];      
	Json::Value json_intrinsic;
	fin.open(in_path+"input_files/"+iname+"_LC_intrinsic.json",std::ifstream::in);
	fin >> json_intrinsic;
	fin.close();
		
	for(int i=0;i<json_intrinsic.size();i++){
	  if( json_intrinsic[i]["time"][0].asDouble() > (tobs_min-td_max) ){
	    fprintf(stderr,"Instrument %s: Intrinsic light curve %d requires an earlier starting time by at least %f days!\n",iname.c_str(),i,json_intrinsic[i]["time"][0].asDouble() - (tobs_min-td_max));
	    check = true;
	  }
	  int Ntime = json_intrinsic[i]["time"].size();
	  if( json_intrinsic[i]["time"][Ntime-1].asDouble() < tobs_max ){
	    fprintf(stderr,"Instrument %s: Intrinsic light curve %d requires a later ending time by at least %f days!\n",iname.c_str(),i, tobs_max - json_intrinsic[i]["time"][Ntime-1].asDouble());
	    check = true;
	  }
	}
      }
    } else {
      fprintf(stderr,"Intrinsic variability of a moving fixed source must be given!");
      check = true;
    }

  } else if ( ex_type == "expanding_source" ){

    if( in_type == "custom" ){
      // Check that all the intrinsic light curves per instrument have the same starting and ending time.
      std::vector<double> tmp_tins_min(names.size());
      std::vector<double> tmp_tins_max(names.size());
      for(int n=0;n<names.size();n++){
	std::string iname = names[n];      
	Json::Value json_intrinsic;
	fin.open(in_path+"input_files/"+iname+"_LC_intrinsic.json",std::ifstream::in);
	fin >> json_intrinsic;
	fin.close();

	std::vector<double> tstart(json_intrinsic.size());	
	std::vector<double> tend(json_intrinsic.size());	
	for(int i=0;i<json_intrinsic.size();i++){
	  int N = json_intrinsic[i]["time"].size();
	  tstart[i] = json_intrinsic[i]["time"][0].asDouble();
	  tend[i] = json_intrinsic[i]["time"][N-1].asDouble();
	}

	tmp_tins_min[n] = *std::min_element(tstart.begin(),tstart.end());
	tmp_tins_max[n] = *std::max_element(tend.begin(),tend.end());
      }
      double tins_min =  *std::min_element(tmp_tins_min.begin(),tmp_tins_min.end());
      double tins_max =  *std::max_element(tmp_tins_max.begin(),tmp_tins_max.end());

      // Change tobs limits as needed
      if( tins_min < tobs_min ){
	fprintf(stdout,"Extending observing time DOWN to the earliest start of the intrinsic light curves!");
	tobs_min = tins_min;
      } 
      if( tobs_max < (tins_max + td_max) ){
	fprintf(stdout,"Extending observing time UP to the latest end of the intrinsic light curves + the maximum time delay!");
	tobs_max = tins_max + td_max;
      } 

      // Output tmin and tmax (in the observerational frame) of the intrinsic light curves in all instruments
      Json::Value out;
      out["tins_min"] = tins_min;
      out["tins_max"] = tins_max;
      std::ofstream file_tins(out_path+"output/intrinsic_times.json");
      file_tins << out;
      file_tins.close();
      
    } else {
      fprintf(stderr,"Intrinsic variability of an expanding source must be given!");
      check = true;
    }
    
  } else if ( ex_type == "moving_variable_source" ){

    // NO instrinsic variability here (incorporated in the extrinsic model). The checks happen in the extrinsic part below. 
    for(int n=0;n<names.size();n++){
      if( root["point_source"]["variability"]["extrinsic"].isMember(names[n]) ){

	if( !root["point_source"]["variability"]["extrinsic"][names[n]].isMember("pixSize") ){
	  fprintf(stderr,"Pixel size for instrument %s is not given.\n",names[n].c_str());
	  check = true;
	}

	int Nx,Ny;
	if( root["point_source"]["variability"]["extrinsic"][names[n]].isMember("Nx") && root["point_source"]["variability"]["extrinsic"][names[n]].isMember("Ny") ){
	  Nx = root["point_source"]["variability"]["extrinsic"][names[n]]["Nx"].asInt();
	  Ny = root["point_source"]["variability"]["extrinsic"][names[n]]["Ny"].asInt();
	  if( Nx != Ny ){
	    fprintf(stderr,"Width and height of the source snapshots for instrument %s must be equal (square profile).\n",names[n].c_str());
	    check = true;
	  }
	} else {
	  fprintf(stderr,"Width and height of the source snapshots in pixels for instrument %s is not given.\n",names[n].c_str());
	  check = true;
	}

	// Check timestep vector compared to the observed time vector of the instrument and the time delay. All times are in days in the observer's frame.
	int Nsteps = root["point_source"]["variability"]["extrinsic"][names[n]]["time"].size();
	double tmax_source = root["point_source"]["variability"]["extrinsic"][names[n]]["time"][Nsteps - 1].asDouble();
	if( tmax_source < (tmaxs[n] - tmins[n] + td_max) ){
	  fprintf(stderr,"Maximum timestep (%f) for instrument %s must be larger than the duration of the observation + the longest time delay, viz. >%f days.\n",tmax_source,names[n].c_str(),tmaxs[n]-tmins[n]+td_max);
	  check = true;
	}	

	// Check snapshot number and names
	std::vector<std::string> fnames;
	std::vector<std::string> paths;
	for(const auto& entry : std::filesystem::directory_iterator(in_path+"input_files/vs_"+names[n]) ){
	  std::filesystem::path p(entry);
	  fnames.push_back(p.stem());
	  paths.push_back(p);
	}
	int Nsnapshots = fnames.size();

	if( Nsnapshots != Nsteps ){
	  fprintf(stderr,"Number of timesteps (%d) and number of snapshots (%d) for instrument %s do not match.\n",Nsteps,Nsnapshots,names[n].c_str());
	  check = true;
	} else {

	  // Check that snapshot names are sequential
	  std::vector<std::string> seq(Nsteps);
	  char str[11];
	  for(int i=0;i<Nsteps;i++){
	    int dum = sprintf(str,"%04d",i);
	    seq[i] = str;
	  }
	  std::vector<std::string> diff;
	  std::sort(fnames.begin(),fnames.end());
	  std::set_difference(seq.begin(),seq.end(),fnames.begin(),fnames.end(),std::inserter(diff,diff.begin()));

	  if( diff.size() > 0 ){
	    // Report missing timesteps 
	    fprintf(stderr,"The snapshot names for instrument %s must be in a sequence. The following timesteps are missing:\n",names[n].c_str());
	    for(int i=0;i<diff.size();i++){
	      fprintf(stderr," %s",diff[i].c_str());
	    }
	    fprintf(stderr,"\n");
	    check = true;
	  } else {
	    // Read the fits headers and make sure the agree with the given width and height
	    for(int i=0;i<Nsteps;i++){
	      std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(paths[i],CCfits::Read,true));
	      CCfits::PHDU& image = pInfile->pHDU();
	      image.readAllKeys();
	      int fNy = image.axis(0);
	      int fNx = image.axis(1);
	      if( fNx != Nx || fNy != Ny ){
		fprintf(stderr,"The dimensions of snapshot %s (%dx%d) for instrument %s do not match the given ones (%dx%d).\n",fnames[i].c_str(),fNx,fNy,names[n].c_str(),Nx,Ny);
		check = true;
	      }
	    }
	  }
	  
	}	
	
      } else {
	fprintf(stderr,"Variability properties for instrument %s are not given.\n",names[n].c_str());
	check = true;
      }      
    }

  } else if( ex_type == "custom" ){
    // Do nothing here
  } else {
    fprintf(stderr,"Unknown variability model: %s\n",ex_type.c_str());
    fprintf(stderr,"Allowed options are:");
    fprintf(stderr," - expanding_source");
    fprintf(stderr," - moving_fixed_source");
    fprintf(stderr," - moving_variable_source");
    check = true;
  }


  
  // Perform checks related to the microlensing variability model
  //================================================================================================================
  if( ex_type == "custom" ){

    // Check custom extrinsic light curves, they must extend from below tobs_min to above tobs_max
    for(int n=0;n<names.size();n++){
      std::string iname = names[n];      
      Json::Value json_extrinsic;
      fin.open(in_path+"input_files/"+iname+"_LC_extrinsic.json",std::ifstream::in);
      fin >> json_extrinsic;
      fin.close();

      for(int q=0;q<json_extrinsic.size();q++){
	if( json_extrinsic[q].size() > 0 ){
	  for(int i=0;i<json_extrinsic[q].size();i++){
	    if( json_extrinsic[q][i]["time"][0] > tobs_min ){
	      fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must have a starting time earlier than the minimum observing time, viz. <%f days.\n",i,iname.c_str(),tobs_min);
	      check = true;
	    }
	    int Ntime = json_extrinsic[q][i]["time"].size();
	    if( json_extrinsic[q][i]["time"][Ntime-1] < tobs_max ){
	      fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must have a later ending time than the maximum observing time, viz. >%f days.\n",i,iname.c_str(),tobs_max);
	      check = true;
	    }
	  }
	}
      }
    }


    
  } else if( ex_type == "expanding_source" ){

    double cutoff_fac = root["point_source"]["variability"]["extrinsic"]["size_cutoff"].asDouble();  // in Rein
    if( cutoff_fac > 7.0 ){
      fprintf(stderr,"Cutoff size for the largest source profile too big. Consider reducing it below 7 Einstein radii.\n");
      check = true;
    }

    Json::Value cosmo;
    fin.open(out_path+"output/angular_diameter_distances.json",std::ifstream::in);
    fin >> cosmo;
    fin.close();
    double Dl  = cosmo[0]["Dl"].asDouble();
    double Ds  = cosmo[0]["Ds"].asDouble();
    double Dls = cosmo[0]["Dls"].asDouble();
    double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
    double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm
    
    for(int n=0;n<names.size();n++){
      if( root["point_source"]["variability"]["extrinsic"]["v_expand"].isMember(names[n]) ){
	// Loop over the maximum extent of the observations in each filter
	double R_max = (tmaxs[n]-tmins[n])*root["point_source"]["variability"]["extrinsic"]["v_expand"][names[n]].asDouble()*8.64; // in 10^14 cm
	if( R_max/Rein > cutoff_fac ){
	  fprintf(stderr,"Maximum physical size of the source in instrument %s (%f) is above the cutoff size of %f Einstein radii.\n",names[n].c_str(),R_max/Rein,cutoff_fac);
	  check = true;
	}
      } else {
	fprintf(stderr,"Expansion velocity for instrument %s not found!\n",names[n].c_str());
	check = true;
      }
    }

    

  } else if( ex_type == "moving_variable_source" ){

    // Check accretion disc size with respect to the map?




  } else if( ex_type == "moving_fixed_source_custom" ){

    for(int n=0;n<names.size();n++){
      if( root["point_source"]["variability"]["extrinsic"].isMember(names[n]) ){
	  
	if( !root["point_source"]["variability"]["extrinsic"][names[n]].isMember("pixSize") ){
	  fprintf(stderr,"Pixel size for instrument %s is not given.\n",names[n].c_str());
	  check = true;
	}
	  
	int Nx,Ny;
	if( root["point_source"]["variability"]["extrinsic"][names[n]].isMember("Nx") && root["point_source"]["variability"]["extrinsic"][names[n]].isMember("Ny") ){
	  Nx = root["point_source"]["variability"]["extrinsic"][names[n]]["Nx"].asInt();
	  Ny = root["point_source"]["variability"]["extrinsic"][names[n]]["Ny"].asInt();
	  if( Nx != Ny ){
	    fprintf(stderr,"Width and height of the source snapshots for instrument %s must be equal (square profile).\n",names[n].c_str());
	    check = true;
	  }
	} else {
	  fprintf(stderr,"Width and height of the source snapshots in pixels for instrument %s is not given.\n",names[n].c_str());
	  check = true;
	}

	// Read the fits headers and make sure the agree with the given width and height
	std::string filepath = in_path+"input_files/cs_"+names[n]+".fits";
	if( std::filesystem::exists(filepath) ){
	  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
	  CCfits::PHDU& image = pInfile->pHDU();
	  image.readAllKeys();
	  int fNy = image.axis(0);
	  int fNx = image.axis(1);
	  if( fNx != Nx || fNy != Ny ){
	    fprintf(stderr,"The dimensions of snapshot %s (%dx%d) for instrument %s do not match the given ones (%dx%d).\n",filepath.c_str(),fNx,fNy,names[n].c_str(),Nx,Ny);
	    check = true;
	  }
	} else {
	    fprintf(stderr,"The custom profile input file for instrument %s does not exist.\n",names[n].c_str());
	    check = true;
	}
	
      } else {
	fprintf(stderr,"Variability properties for instrument %s are not given.\n",names[n].c_str());
	check = true;
      }

    }

    if( check != true ){	
      // Check accretion disc size with respect to the map - it must be possible to create a profile
      std::map<std::string,std::string> main_map;
      main_map.insert( std::pair<std::string,std::string>("type","custom") );
      main_map.insert( std::pair<std::string,std::string>("shape","custom") );
      main_map.insert( std::pair<std::string,std::string>("pixSizePhys","1") ); // The size of the map pixels is dummy here
      main_map.insert( std::pair<std::string,std::string>("incl","0") );
      main_map.insert( std::pair<std::string,std::string>("orient","0") );
      main_map.insert( std::pair<std::string,std::string>("filename","") );
      main_map.insert( std::pair<std::string,std::string>("profPixSizePhys","") );

      std::vector<double> rhalfs(names.size());
      for(int i=0;i<names.size();i++){
	main_map["filename"] = in_path+"input_files/cs_"+names[i]+".fits";
	main_map["profPixSizePhys"] = std::to_string( root["point_source"]["variability"]["extrinsic"][names[i]]["pixSize"].asDouble() );
	gerlumph::BaseProfile* profile = gerlumph::FactoryProfile::getInstance()->createProfileFromPars(main_map);
	rhalfs[i] = profile->getHalfRadius(); // the second argument is the lrest, which is dummy in this call
	delete(profile);
      }

      Json::Value cosmo;
      fin.open(out_path+"output/angular_diameter_distances.json",std::ifstream::in);
      fin >> cosmo;
      fin.close();
      double Dl  = cosmo[0]["Dl"].asDouble();
      double Ds  = cosmo[0]["Ds"].asDouble();
      double Dls = cosmo[0]["Dls"].asDouble();
      double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
      double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm
      for(int i=0;i<names.size();i++){
	if( rhalfs[i]/Rein > 7 ){
	  fprintf(stderr,"Accretion disc size for instrument %s is too big. Consider reducing rest wavelength %f so that the disc becomes smaller than 7 Einstein radii.\n",names[i].c_str(),lrest[i]);
	  check = true;	
	}
      }
    }    

    

  } else if( ex_type == "moving_fixed_source" ){
      
    // Quick loop to check if the accretion dize does not become too large for some wavelength, e.g. above several Rein
    Json::Value::Members json_members = root["point_source"]["variability"]["extrinsic"]["profiles"].getMemberNames();
    std::map<std::string,std::string> main_map; // size should be json_members.size()+2
    for(int i=0;i<json_members.size();i++){
      main_map.insert( std::pair<std::string,std::string>(json_members[i],root["point_source"]["variability"]["extrinsic"]["profiles"][json_members[i]].asString()) );
    }
    main_map.insert( std::pair<std::string,std::string>("rhalf","") );
    
    std::vector<double> rhalfs(names.size());
    for(int i=0;i<names.size();i++){
      if( main_map["type"] == "vector" ){
	rhalfs[i] = root["point_source"]["variability"]["extrinsic"]["profiles"]["rhalf"][i].asDouble();
      } else {
	rhalfs[i] = gerlumph::BaseProfile::getSize(main_map,lrest[i]);
      }
    }
    
    Json::Value cosmo;
    fin.open(out_path+"output/angular_diameter_distances.json",std::ifstream::in);
    fin >> cosmo;
    fin.close();
    double Dl  = cosmo[0]["Dl"].asDouble();
    double Ds  = cosmo[0]["Ds"].asDouble();
    double Dls = cosmo[0]["Dls"].asDouble();
    double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
    double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm
    for(int i=0;i<names.size();i++){
      if( rhalfs[i]/Rein > 7 ){
	fprintf(stderr,"Accretion disc size for instrument %s is too big. Consider reducing rest wavelength %f so that the disc becomes smaller than 7 Einstein radii.\n",names[i].c_str(),lrest[i]);
	check = true;	
      }
    }
    

    
  } else {
    fprintf(stderr,"Unknown variability model: %s\n",ex_type.c_str());
    fprintf(stderr,"Allowed options are:");
    fprintf(stderr," - custom");
    fprintf(stderr," - expanding_source");
    fprintf(stderr," - moving_fixed_source");
    fprintf(stderr," - moving_fixed_source_custom");
    fprintf(stderr," - moving_variable_source");
    check = true;
  }



  if( root["point_source"]["variability"].isMember("unmicro") ){
    for(int n=0;n<names.size();n++){
      // Check flux ratio
      double flux_ratio = root["point_source"]["variability"][names[n]]["flux_ratio"].asDouble();
      if( flux_ratio < 0 || 1 < flux_ratio ){
	fprintf(stderr,"Flux ratio in instrument %s must be bettwen 0 and 1 (current value: %f).\n",names[n].c_str(),flux_ratio);
	check = true;
      }
    }
  }

  
      
 
  //=======================================================================================================================    

  
  if( check ){
    return 1;
  }

  
  // Output tmin and tmax for the observerational frame in all instruments
  Json::Value out;
  out["tobs_min"] = tobs_min;
  out["tobs_max"] = tobs_max;
  std::ofstream file_tobs(out_path+"output/tobs.json");
  file_tobs << out;
  file_tobs.close();
  
  return 0;
}
