#include <algorithm>
#include <stdlib.h>
#include <fstream>

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

std::vector<double> BaseNoise::getBoxMullerVector(int N){
  this->seed += 2; // increment seed at each call
  srand48(seed);
  std::vector<double> z(N);
  double z1,z2,u1,u2;
  //Applying the Box-Muller transformation
  for(int i=0;i<N;i++){
    u1 = drand48();
    u2 = drand48();
    z1 = sqrt(-2.0 * log(u1)) * cos(this->two_pi * u2);
    //    z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    z[i] = z1;
  }
  return z;
}

void BaseNoise::addNoise(RectGrid* mydata){
  for(int i=0;i<mydata->Nz;i++){
    mydata->z[i] += this->noise_realization->z[i];
  }
}

void BaseNoise::outputNoiseRealization(std::string output,std::string instrument_name){
  FitsInterface::writeFits(this->noise_realization->Nx,this->noise_realization->Ny,this->noise_realization->z,output + instrument_name + "_noise_realization.fits");
}
// END: BaseNoise =================================




// START: PoissonNoise ============================
PoissonNoise::PoissonNoise(double texp,double Msb,double zp,double readout,double res): BaseNoise(texp),Msb(Msb),ZP(zp){
  this->res2 = res*res;
  //this->Ibg = pow(10.0,-0.4*(Msb-ZP))*res2*texp + pow(readout,2); // the sky background in electrons
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

void PoissonNoise::initializeFromData(RectGrid* mydata){
  for(int i=0;i<mydata->Nz;i++){
    double electrons = mydata->z[i];//*this->texp;
    double sigma = sqrt(electrons/this->texp + this->Ibg);
    this->sigma_map->z[i] = sigma;
  }
}

void PoissonNoise::calculateNoise(){
  std::vector<double> z = this->getBoxMullerVector(this->noise_realization->Nz);
  for(int i=0;i<this->noise_realization->Nz;i++){
    this->noise_realization->z[i] = z[i]*this->sigma_map->z[i];
  }
}

void PoissonNoise::outputNoiseProperties(std::string output,std::string instrument_name){
  this->outputNoiseRealization(output,instrument_name);
  FitsInterface::writeFits(this->sigma_map->Nx,this->sigma_map->Ny,this->sigma_map->z,output + instrument_name + "_sigma_map.fits");

  Json::Value json_obj;
  json_obj["type"] = "PoissonNoise";
  json_obj["sigma_bg"] = sqrt(this->Ibg);  
  std::ofstream noise_properties(output + instrument_name + "_noise_properties.json");
  noise_properties << json_obj;
  noise_properties.close();
}
// END: PoissonNoise ==============================





// START: UniformGaussian ============================
UniformGaussian::UniformGaussian(double sn): BaseNoise(1){
  this->sn = sn;
}

void UniformGaussian::setGrid(RectGrid* obs_grid){
  this->noise_realization = new RectGrid(obs_grid->Nx,obs_grid->Ny,obs_grid->xmin,obs_grid->xmax,obs_grid->ymin,obs_grid->ymax);
}

void UniformGaussian::initializeFromData(RectGrid* mydata){
  double maxdata = *std::max_element(mydata->z,mydata->z+mydata->Nz);
  this->sigma = maxdata/this->sn;
}
  
void UniformGaussian::calculateNoise(){
  std::vector<double> z = this->getBoxMullerVector(this->noise_realization->Nz);
  double noise,min_noise;
  for(int i=0;i<this->noise_realization->Nz;i++){
    noise = z[i]*this->sigma;
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

  Json::Value json_obj;
  json_obj["type"] = "UniformGaussian";
  json_obj["sigma"] = this->sigma;  
  std::ofstream noise_properties(output + instrument_name + "_noise_properties.json");
  noise_properties << json_obj;
  noise_properties.close();
}
// END: UniformGaussian ==============================
