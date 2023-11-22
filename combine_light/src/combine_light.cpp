#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>
#include <filesystem>

#include "json/json.h"

#include "vkllib.hpp"

#include "auxiliary_functions.hpp"
#include "instruments.hpp"
#include "noise.hpp"

int main(int argc,char* argv[]){

  /*
    Requires:
    - main input
    - fluxes
    - extended lensed source light   (1 per instrument, from fproject.cpp)
    - lens light                     (1 per instrument, from lens_light_mass.cpp)

    If point source exists:
    - multiple_images.json
    - tobs.json
  */

  
  //=============== BEGIN:PARSE INPUT ===========================================
  std::ifstream fin;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);                                           // INPUT FILE 1: main MOLET input file
  fin >> root;
  fin.close();

  std::string in_path = argv[2];
  std::string out_path = argv[3];
  double zs = root["source"]["redshift"].asDouble();
  bool convolve_lens = false;
  if( root["output_options"].isMember("convolve_lens") ){
    convolve_lens = root["output_options"]["convolve_lens"].asBool();
  }
  bool conserve_flux = false;
  if( root["output_options"].isMember("conserve_flux") ){
    conserve_flux = root["output_options"]["conserve_flux"].asBool();
  }
  
  Json::Value fluxes;
  fin.open(out_path+"output/fluxes.json",std::ifstream::in);                     // INPUT FILE 2: fluxes
  fin >> fluxes;
  fin.close();


  Json::Value images;
  if( root.isMember("point_source") ){
    // Read the multiple images' parameters from JSON
    fin.open(out_path+"output/multiple_images.json",std::ifstream::in);          // INPUT FILE 3 (only for point source): properties of the multiple images
    fin >> images;
    fin.close();
  }
  //=============================================================================
  








  
  
  // Loop over the instruments
  // ===================================================================================================================
  // ===================================================================================================================
  for(int b=0;b<root["instruments"].size();b++){
    const Json::Value instrument = root["instruments"][b];
    std::string instrument_name = root["instruments"][b]["name"].asString();
    Instrument mycam(instrument_name,root["instruments"][b]["ZP"].asDouble(),root["instruments"][b]["noise"],conserve_flux);
    double F_conv_extended,F_conv_lens,F_conv_ps;



    
    //=============== START: INITIALIZE FINAL IMAGE GRIDS ====================
    // Set output image plane in super-resolution
    double xmin = root["instruments"][b]["field-of-view_xmin"].asDouble();
    double xmax = root["instruments"][b]["field-of-view_xmax"].asDouble();
    double ymin = root["instruments"][b]["field-of-view_ymin"].asDouble();
    double ymax = root["instruments"][b]["field-of-view_ymax"].asDouble();
    int res_x = static_cast<int>(ceil((xmax-xmin)/mycam.resolution));
    int res_y = static_cast<int>(ceil((ymax-ymin)/mycam.resolution));
    int super_factor = 10;
    if( root["output_options"].isMember("super_factor") ){
      super_factor = root["output_options"]["super_factor"].asInt();
    }
    int super_res_x = super_factor*res_x;
    int super_res_y = super_factor*res_y;
    vkl::RectGrid supersim(super_res_x,super_res_y,xmin,xmax,ymin,ymax);    

    // Get the psf in super-resolution, crop it, and create convolution kernel
    mycam.preparePSF(&supersim,0.999);
    //vkl::FitsInterface::writeFits(mycam.scaled_psf->Nx,mycam.scaled_psf->Ny,mycam.scaled_psf->z,out_path + "output/supersampled_psf.fits");
    //vkl::FitsInterface::writeFits(mycam.cropped_psf->Nx,mycam.cropped_psf->Ny,mycam.cropped_psf->z,out_path + "output/cropped_psf.fits");
    //vkl::FitsInterface::writeFits(supersim.Nx,supersim.Ny,mycam.kernel,out_path + "output/kernel_psf.fits");

    double total = mycam.original_psf->integrate();
    std::cout << "Original PSF flux: " << total << " " << -2.5*log10(total)+mycam.ZP << std::endl;
    double total_scaled = mycam.scaled_psf->integrate();
    std::cout << "Scaled PSF flux:   " << total_scaled << " " << -2.5*log10(total_scaled)+mycam.ZP << std::endl;
    double total_cropped = mycam.cropped_psf->integrate();
    std::cout << "Cropped PSF flux:  " << total_cropped << " " << -2.5*log10(total_cropped)+mycam.ZP << " (" << 100*(total_cropped/total_scaled) << ")" << std::endl;

    
    
    vkl::RectGrid super_extended = vkl::RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax,out_path+"output/"+instrument_name+"_lensed_image_super.fits");        // INPUT FILE 5: super-resolved lensed source
    vkl::RectGrid super_lens_light = vkl::RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax,out_path+"output/"+instrument_name+"_lens_light_super.fits");        // INPUT FILE 6: super-resolved lens light
    
    // Create static image (lens light and extended lensed source)
    vkl::RectGrid obs_static = createObsStatic(&mycam,&super_extended,&super_lens_light,res_x,res_y,F_conv_extended,F_conv_lens,convolve_lens); // This is in units of electrons/(s arcsec^2)!

    // Calculate total flux of the static image
    double F_final_static = obs_static.integrate();
    
    // Assign noise grid
    mycam.noise->setGrid(&obs_static);


    fluxes["lensed_source_flux"][instrument_name]["total_convolved"]["flux"] = F_conv_extended;
    fluxes["lensed_source_flux"][instrument_name]["total_convolved"]["mag"]  = convertFromFlux(F_conv_extended,mycam.ZP);
    fluxes["lens_flux"][instrument_name]["total_convolved"]["flux"] = F_conv_lens;
    fluxes["lens_flux"][instrument_name]["total_convolved"]["mag"]  = convertFromFlux(F_conv_lens,mycam.ZP);
    fluxes["final_static"][instrument_name]["flux"] = F_final_static;
    fluxes["final_static"][instrument_name]["mag"]  = convertFromFlux(F_final_static,mycam.ZP);

    
    // All the static light components have been created.
    // The resulting observed image (not super-resolved) is "obs_static".
    // If there is no time dimension required, the code just outputs the obs_static image and stops.
    //===============   END: INITIALIZE FINAL IMAGE GRIDS ====================
    



    //=============== CREATE THE TIME VARYING LIGHT ====================
    if( root.isMember("point_source") ){

      
      //==================================== Configure the PSF for the point source ===========================================
      // Perturb the PSF at each image location
      std::vector<Instrument*> instrument_list(images.size());
      for(int q=0;q<images.size();q++){
	//TransformPSF* dum = new TransformPSF(images[q]["x"].asDouble(),images[q]["y"].asDouble(),0.0,false,false);
	//transPSF[q] = dum;
	instrument_list[q] = &mycam; // just copy the same instrument pointer per image
      }
      // Set the PSF related offsets for each image
      std::vector<offsetPSF> PSFoffsets(images.size());
      for(int q=0;q<images.size();q++){
	PSFoffsets[q] = instrument_list[q]->offsetPSFtoPosition(images[q]["x"].asDouble(),images[q]["y"].asDouble(),&supersim);
      }
      FILE* fh = fopen((out_path+"output/psf_locations.dat").c_str(),"w");
      for(int q=0;q<images.size();q++){
	PSFoffsets[q].printFrame(fh,supersim.Nx,supersim.Ny,supersim.xmin,supersim.xmax,supersim.ymin,supersim.ymax);
      }
      fclose(fh);
      // Calculate the appropriate PSF sums
      std::vector<double> psf_partial_sums(images.size());
      for(int q=0;q<images.size();q++){
	psf_partial_sums[q] = instrument_list[q]->sumPSF(&PSFoffsets[q]);
      }
      //=======================================================================================================================


      
      // *********************** Product: Observed sampled cut-outs (images) *******************************************
      if( root["output_options"]["output_PS_cutouts"].asBool() ){
	std::vector<double> tobs;
	for(int t=0;t<instrument["time"].size();t++){
	  tobs.push_back(instrument["time"][t].asDouble());
	}

	std::vector<std::string> mocks;
	for(auto const& entry : std::filesystem::directory_iterator(out_path)){
	  std::filesystem::path p(entry);
	  std::string dir_name = p.stem();
	  if( dir_name.compare(0,4,"mock") == 0 ){
	    mocks.push_back( p.stem() );
	  }
	}
  
	for(int m=0;m<mocks.size();m++){
	  Json::Value samp_LC;
	  fin.open(out_path+"/"+mocks[m]+"/"+instrument_name+"_LC_sampled.json",std::ifstream::in);
	  fin >> samp_LC;
	  fin.close();
	  writeAllCutouts(tobs,images,samp_LC,&super_extended,&super_lens_light,&mycam,PSFoffsets,instrument_list,psf_partial_sums,res_x,res_y,mocks[m],convolve_lens,out_path);
	}
      }
      // *********************** End of product ************************************************************************

      
      
      // *********************** Product: Just one cutout with only the macromagnification *****************************
      // For each image we use the same brightness value, which is a multiple of the total extended lensed source flux.
      // We gather the entire extended lensed source flux in one pixel and multiply by a factor
      double M_ps_unlensed = root["point_source"]["M_tot_unlensed"].asDouble();
      double F_ps_unlensed = pow(10,-0.4*(M_ps_unlensed - mycam.ZP)); // This is a total flux in [electrons/s]
      double area_super = supersim.step_x*supersim.step_y;
      double F_ps = 0.0;
      std::vector<double> image_signal(images.size());
      for(int q=0;q<images.size();q++){
	double macro_mag = fabs(images[q]["mag"].asDouble());
	image_signal[q] = macro_mag*(F_ps_unlensed/area_super); // This needs to be flux density, i.e. divided by the area of a pixel [electrons/(s arcsec^2)]
	F_ps += image_signal[q]*area_super; // This is the total magnified flux, i.e. multiplied by the area of a pixel
      }
      vkl::RectGrid obs_ps_light = createObsPS(&supersim,image_signal,PSFoffsets,instrument_list,psf_partial_sums,res_x,res_y); // This is in [electrons/(s arcsec^2)]

      // Calculate total flux of the PS light only
      F_conv_ps = obs_ps_light.integrate();

      // The two grids, obs_ps_light and obs_static, have the same resolution, so I can simply add their fluxes
      for(int i=0;i<obs_ps_light.Nz;i++){
	obs_ps_light.z[i] += obs_static.z[i]; // <------ BE CAREFUL: 'obs_ps_light' now also contains the static light
      }      

      // Calculate total flux of the 'PS macro' image
      double F_final_with_ps = obs_ps_light.integrate();

      // Output the noiseless static image with a PS with macromagnification only
      writeCutout(&obs_ps_light,out_path+"output/OBS_"+instrument_name+"_ps_macro_noiseless.fits");
      
      // Convert to units of electrons/s, must be done before adding noise
      double area = obs_ps_light.step_x*obs_ps_light.step_y;
      for(int i=0;i<obs_ps_light.Nz;i++){
	obs_ps_light.z[i] *= area;
      }
      
      // Adding time-dependent noise here
      mycam.noise->initializeFromData(&obs_ps_light);
      mycam.noise->calculateNoise();
      mycam.noise->addNoise(&obs_ps_light); // <------ BE CAREFUL: Noise is added directly to the 'obs_ps_light' pixels.
      mycam.noise->outputNoiseProperties(out_path + "output/",instrument_name+"_ps_macro");

      // Output the observed static image with a PS with macromagnification only
      writeCutout(&obs_ps_light,out_path+"output/OBS_"+instrument_name+"_ps_macro.fits");
      //vkl::FitsInterface::writeFits(obs_static.Nx,obs_static.Ny,obs_static.z,fname);

      fluxes["lensed_ps_flux"][instrument_name]["flux"]   = F_ps;
      fluxes["lensed_ps_flux"][instrument_name]["mag"]    = convertFromFlux(F_ps,mycam.ZP);
      fluxes["lensed_ps_flux"][instrument_name]["flux"]   = F_conv_ps;
      fluxes["lensed_ps_flux"][instrument_name]["mag"]    = convertFromFlux(F_conv_ps,mycam.ZP);
      fluxes["unlensed_ps_flux"][instrument_name]["flux"] = F_ps_unlensed;
      fluxes["unlensed_ps_flux"][instrument_name]["mag"]  = convertFromFlux(F_ps_unlensed,mycam.ZP);
      fluxes["final_with_ps"][instrument_name]["flux"]    = F_final_with_ps;
      fluxes["final_with_ps"][instrument_name]["mag"]     = convertFromFlux(F_final_with_ps,mycam.ZP);
      // *********************** End of product ************************************************************************

    }
    //================= END:CREATE THE TIME VARYING LIGHT ====================

    


    //=============== START: CREATE STATIC LIGHT ====================
    // Output the noiseless static image
    writeCutout(&obs_static,out_path + "output/OBS_" + instrument_name + "_noiseless.fits");
    
    // Convert to units of electrons/s, must be done before adding noise
    double area = obs_static.step_x*obs_static.step_y;
    for(int i=0;i<obs_static.Nz;i++){
      obs_static.z[i] *= area;
    }
    
    // Adding noise here and output realization
    mycam.noise->initializeFromData(&obs_static);
    mycam.noise->calculateNoise();
    mycam.noise->addNoise(&obs_static); // <------ BE CAREFUL: Noise is added directly to the obs_static pixels.
    mycam.noise->outputNoiseProperties(out_path + "output/",instrument_name);
    
    // Output the observed static image
    writeCutout(&obs_static,out_path + "output/OBS_" + instrument_name + ".fits");
    //===============   END: CREATE STATIC LIGHT ====================

    
    
  }
  // Loop over the instruments ends here
  // ===================================================================================================================
  // ===================================================================================================================


  // Add convolved and final fluxes
  std::ofstream file_fluxes(out_path+"output/fluxes.json");
  file_fluxes << fluxes;
  file_fluxes.close();


  
  return 0;
}
