#include <algorithm>
#include <stdlib.h>

#include "vkllib.hpp"

#include "instruments.hpp"
#include "noise.hpp"

// START: BaseNoise ===============================
BaseNoise::BaseNoise(double texp){
  this->texp = texp;
}

BaseNoise::~BaseNoise(){
  delete(noise_realization);
}

void BaseNoise::addNoise(RectGrid* mydata){
  for(int i=0;i<mydata->Nz;i++){
    mydata->z[i] += this->noise_realization->z[i];
  }
}

void BaseNoise::outputNoiseRealization(std::string output,std::string instrument_name){
  FitsInterface::writeFits(this->noise_realization->Nx,this->noise_realization->Ny,this->noise_realization->z,output + "OBS_" + instrument_name + "_noise_realization.fits");
}
// END: BaseNoise =================================




// START: PoissonNoise ============================
PoissonNoise::PoissonNoise(double texp,double Msb,double zp,double readout,double res): BaseNoise(texp),Msb(Msb),ZP(zp){
  this->res2 = res*res;
  //this->Ibg = pow(10.0,-0.4*(Msb-ZP))*res2*texp + pow(readout,2); // the sky background in electrons
  //double sky_bg = pow( pow(10.0,-0.4*(Msb-ZP))*res2 ,2)/texp;
  double sky_bg = pow(10.0,-0.4*(Msb-ZP))*res2/texp;
  double readout_noise = pow(readout,2)/pow(texp,2);
  this->Ibg =  sky_bg + readout_noise; // the sky background in electrons
}

PoissonNoise::~PoissonNoise(){
  delete(sigma_map);
}

void PoissonNoise::setGrid(RectGrid* obs_grid){
  this->noise_realization = new RectGrid(obs_grid->Nx,obs_grid->Ny,obs_grid->xmin,obs_grid->xmax,obs_grid->ymin,obs_grid->ymax);
  this->sigma_map = new RectGrid(obs_grid->Nx,obs_grid->Ny,obs_grid->xmin,obs_grid->xmax,obs_grid->ymin,obs_grid->ymax);
}

void PoissonNoise::calculateNoise(RectGrid* mydata){
  this->seed += 2; // increment seed at each call
  srand48(seed);
  
  double z1,z2,u1,u2,sdev,electrons;
  //Applying the Box-Muller transformation
  for(int i=0;i<mydata->Nz;i++){
    u1 = drand48();
    u2 = drand48();
    z1 = sqrt(-2.0 * log(u1)) * cos(this->two_pi * u2);
    //    z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);    

    // The units of obs_base are electrons/(s pix), so multiply it by exposure time
    electrons = mydata->z[i];//*this->texp;
    sdev = sqrt(electrons/this->texp + this->Ibg);
    this->sigma_map->z[i] = sdev;
    this->noise_realization->z[i] = z1*sdev;
  }
}

void PoissonNoise::outputNoiseProperties(std::string output,std::string instrument_name){
  this->outputNoiseRealization(output,instrument_name);
  FitsInterface::writeFits(this->sigma_map->Nx,this->sigma_map->Ny,this->sigma_map->z,output + "OBS_" + instrument_name + "_sigma_map.fits");
  // output the background sigma
}
// END: PoissonNoise ==============================





// START: UniformGaussian ============================
UniformGaussian::UniformGaussian(double sn): BaseNoise(1){
  this->sn = sn;
}

void UniformGaussian::setGrid(RectGrid* obs_grid){
  this->noise_realization = new RectGrid(obs_grid->Nx,obs_grid->Ny,obs_grid->xmin,obs_grid->xmax,obs_grid->ymin,obs_grid->ymax);
}

void UniformGaussian::calculateNoise(RectGrid* mydata){
  this->seed += 2; // increment seed at each call

  double maxdata = *std::max_element(mydata->z,mydata->z+mydata->Nz);
  double sigma = maxdata/this->sn;
  srand48(seed);

  double min_noise = sigma; // just a starting value
  double z1,z2,u1,u2,noise;
  //Applying the Box-Muller transformation
  for(int i=0;i<this->noise_realization->Nz;i++){
    u1 = drand48();
    u2 = drand48();
    z1 = sqrt(-2.0 * log(u1)) * cos(this->two_pi * u2);
    //    z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);    

    noise = z1*sigma;
    this->noise_realization->z[i] += noise;
    if( noise < min_noise ){
      min_noise = noise;
    }
  }

  // Renormalize by adding the minimum (negative) noise value
  // This is needed because the pixels can contain nan's if mag scale is used
  for(int i=0;i<this->noise_realization->Nz;i++){
    this->noise_realization->z[i] += abs(min_noise);
  }
}

void UniformGaussian::outputNoiseProperties(std::string output,std::string instrument_name){
  this->outputNoiseRealization(output,instrument_name);
  // Output the sdev of the gaussian
}
// END: UniformGaussian ==============================
