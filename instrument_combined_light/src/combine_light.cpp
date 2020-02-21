#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "vkllib.hpp"

#include "auxiliary_functions.hpp"

int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  std::ifstream fin;
  Json::Value::Members jmembers;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string path = argv[2];

  std::string output     = path + "output/";
  std::string input_path = path + "input_files/";

  Json::Value images;
  Json::Value intrinsic_lc;
  Json::Value extrinsic_lc;
  if( root.isMember("point_source") ){
    // Read the image parameters
    fin.open(output+"multiple_images.json",std::ifstream::in);
    fin >> images;
    fin.close();
    
    // Read intrinsic light curve
    if( root["point_source"]["variability"]["intrinsic"]["type"].asString() == "custom" ){
      fin.open(input_path+"intrinsic_light_curves.json",std::ifstream::in);
    } else {
      fin.open(output+"intrinsic_light_curves.json",std::ifstream::in);
    }
    fin >> intrinsic_lc;
    fin.close();

    // Read extrinsic light curve
    if( root["point_source"]["variability"]["extrinsic"]["type"].asString() == "custom" ){
      fin.open(input_path+"extrinsic_light_curves.json",std::ifstream::in);
    } else {
      fin.open(output+"extrinsic_light_curves.json",std::ifstream::in);
    }
    fin >> extrinsic_lc;
    fin.close();
  }
  //================= END:PARSE INPUT =======================



  // Loop over the bands
  const Json::Value bands = root["instrument"]["bands"];
  for(int b=0;b<bands.size();b++){
    std::string band_name = bands[b]["name"].asString();
    
    // Set output image plane in super-resolution
    double width  = bands[b]["field-of-view_x"].asDouble();
    double height = bands[b]["field-of-view_y"].asDouble();
    int res_x = static_cast<int>(ceil(width/bands[b]["resolution"].asDouble()));
    int res_y = static_cast<int>(ceil(height/bands[b]["resolution"].asDouble()));
    int super_res_x = 10*res_x;
    int super_res_y = 10*res_y;
    ImagePlane mysim(super_res_x,super_res_y,width,height);

    
    // Get the psf in super-resolution, crop it, and create convolution kernel
    std::string psf_path = input_path + "psf_" + band_name + ".fits";
    double psf_width  = bands[b]["psf"]["width"].asDouble();
    double psf_height = bands[b]["psf"]["height"].asDouble();
    int psf_Nx = bands[b]["psf"]["pix_x"].asInt();
    int psf_Ny = bands[b]["psf"]["pix_y"].asInt();
    PSF mypsf(psf_path,psf_Nx,psf_Ny,psf_width,psf_height,&mysim);
    mypsf.cropPSF(0.99);
    mypsf.createKernel(mysim.Ni,mysim.Nj);

    
    //=============== BEGIN:CREATE THE FIXED EXTENDED LENSED LIGHT ====================
    ImagePlane* extended = new ImagePlane(output+"lensed_image_super.fits",super_res_x,super_res_y,width,height);
    mypsf.convolve(extended);
    //extended->writeImage(output+"psf_lensed_image_super.fits");
    //================= END:CREATE THE FIXED EXTENDED LENSED LIGHT ====================

    //=============== BEGIN:CREATE THE FIXED LENS GALAXY LIGHT ====================
    ImagePlane* lens_light = new ImagePlane(output+"lens_light_super.fits",super_res_x,super_res_y,width,height);
    mypsf.convolve(lens_light);
    //lens_light->writeImage(output+"psf_lens_light_super.fits");
    //================= END:CREATE THE FIXED LENS GALAXY LIGHT ====================


    //=============== BEGIN:CREATE FIXED MASKS ====================
    //================= END:CREATE FIXED MASKS ====================


    // Create combined extended lensed features and lens light
    ImagePlane base(super_res_x,super_res_y,width,height); 
    for(int i=0;i<base.Nm;i++){
      base.img[i] = lens_light->img[i] + extended->img[i];
    }
    delete(extended);
    delete(lens_light);


    
    if( root.isMember("point_source") ){
	
      //=============== BEGIN:CREATE THE TIME VARYING LIGHT ====================
      
      // Read observed time vector
      std::vector<double> obs_time;
      for(int t=0;t<bands[b]["time"].size();t++){
	obs_time.push_back(bands[b]["time"][t].asDouble());
      }
      
      // Get maximum image time delay
      double td_max = 0.0;
      for(int q=0;q<images.size();q++){
	double td = images[q]["dt"].asDouble();
	if( td > td_max ){
	  td_max = td;
	}
      }
	
      // Read the light curves
      LightCurve LC_intrinsic(intrinsic_lc[band_name]);
      std::vector<LightCurve> LC_extrinsic;
      for(int q=0;q<images.size();q++){
	LC_extrinsic.push_back(extrinsic_lc[q][band_name]);
      }

      // Create interpolated intrinsic + extrinsic signal for each image
      double* intrinsic_signal = (double*) malloc(obs_time.size()*sizeof(double));
      double** img_signal = (double**) malloc(images.size()*sizeof(double*));
      for(int q=0;q<images.size();q++){
	LC_intrinsic.interpolate(obs_time,td_max - images[q]["dt"].asDouble(),intrinsic_signal);
	img_signal[q] = (double*) malloc(obs_time.size()*sizeof(double));
	if( LC_extrinsic[q].time.size() > 0 ){
	  LC_extrinsic[q].interpolate(obs_time,0.0,img_signal[q]);
	  for(int t=0;t<obs_time.size();t++){
	    // Do nothing in this loop in order to keep the microlensing signal only
	    //img_signal[q][t] = intrinsic_signal[t]; // this line includes only the intrinsic signal and excludes microlensing
	    img_signal[q][t] *= intrinsic_signal[t]; // this line includes both intrinsic and microlensing signals
	  }
	} else {
	  for(int t=0;t<obs_time.size();t++){
	    img_signal[q][t] = intrinsic_signal[t]; // this line includes only the intrinsic signal and excludes microlensing
	  }
	}
      }
      free(intrinsic_signal);


      for(int t=0;t<obs_time.size();t++){
	printf("%3d: ",t);
	for(int q=0;q<images.size();q++){
	  printf("%f ",img_signal[q][t]);
	}
	printf("\n");
      }

      
      // Set the PSF related offsets for each image
      std::vector<offsetPSF> PSFoffsets(images.size());
      for(int q=0;q<images.size();q++){
	PSFoffsets[q] = mypsf.offsetPSFtoPosition(images[q]["x"].asDouble(),images[q]["y"].asDouble(),&mysim);
	//printf("%d %d %d %d\n",PSFoffsets[q].offset_image,PSFoffsets[q].offset_cropped,PSFoffsets[q].nj,PSFoffsets[q].ni);
      }

      // Loop over time starts here
      for(int t=0;t<obs_time.size();t++){
	ImagePlane pp_light(super_res_x,super_res_y,width,height);
	
	for(int q=0;q<images.size();q++){
	  double macro_mag = abs(images[q]["mag"].asDouble());
	  //double macro_mag = 10.0;
	  double factor = macro_mag * img_signal[q][t];
	  for(int i=0;i<PSFoffsets[q].ni;i++){
	    for(int j=0;j<PSFoffsets[q].nj;j++){
	      int index_img = pp_light.Nj*i + j;
	      int index_psf = i*mypsf.cropped_psf->Nj + j;
	      //	    pp_light.img[PSFoffsets[q].offset_image + pp_light.Nj*i + j] += 1.0;
	      pp_light.img[PSFoffsets[q].offset_image + index_img] += factor*mypsf.cropped_psf->img[PSFoffsets[q].offset_cropped + index_psf];
	    }
	  }
	}
	
	for(int i=0;i<pp_light.Nm;i++){
	  pp_light.img[i] = pp_light.img[i] + base.img[i];
	}

	// Bin image from 'super' to observed resolution
	ImagePlane obs_img(res_x,res_y,width,height);
	int* counts = (int*) calloc(obs_img.Nm,sizeof(int));
	double inf_dx = width/super_res_x;
	double inf_dy = height/super_res_y;
	double obs_dx = width/res_x;
	double obs_dy = height/res_y;
	for(int i=0;i<pp_light.Ni;i++){
	  int ii = (int) floor(i*inf_dy/obs_dy);
	  for(int j=0;j<pp_light.Nj;j++){
	    int jj = (int) floor(j*inf_dx/obs_dx);
	    obs_img.img[ii*obs_img.Nj + jj] += pp_light.img[i*pp_light.Nj + j];
	    counts[ii*obs_img.Nj + jj] += 1;
	  }
	}
	for(int i=0;i<obs_img.Nm;i++){
	  obs_img.img[i] = obs_img.img[i]/counts[i];
	}
	free(counts);

	char buf[4];
	sprintf(buf,"%03d",t);
	std::string timestep = buf;
	obs_img.writeImage(output + "OBS_" + band_name + "_" + timestep + ".fits");
      }

      for(int q=0;q<images.size();q++){
	free(img_signal[q]);
      }
      free(img_signal);
      //================= END:CREATE THE TIME VARYING LIGHT ====================

    } else {

      // Bin image from 'super' to observed resolution
      ImagePlane obs_img(res_x,res_y,width,height);
      int* counts = (int*) calloc(obs_img.Nm,sizeof(int));
      double inf_dx = width/super_res_x;
      double inf_dy = height/super_res_y;
      double obs_dx = width/res_x;
      double obs_dy = height/res_y;
      for(int i=0;i<base.Ni;i++){
	int ii = (int) floor(i*inf_dy/obs_dy);
	for(int j=0;j<base.Nj;j++){
	  int jj = (int) floor(j*inf_dx/obs_dx);
	  obs_img.img[ii*obs_img.Nj + jj] += base.img[i*base.Nj + j];
	  counts[ii*obs_img.Nj + jj] += 1;
	}
      }
      for(int i=0;i<obs_img.Nm;i++){
	obs_img.img[i] = obs_img.img[i]/counts[i];
      }
      free(counts);
      obs_img.writeImage(output + "OBS_" + band_name + ".fits");
      
    }
    
  }
  
  return 0;
}
