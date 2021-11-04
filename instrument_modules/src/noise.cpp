#include <algorithm>
#include <stdlib.h>

#include "vkllib.hpp"

#include "instruments.hpp"
#include "noise.hpp"

// START: NoNoise ====================================
NoNoise::NoNoise(){}
void NoNoise::addNoise(RectGrid* mydata){}
// END: NoNoise ======================================


// START: PoissonNoise ============================
PoissonNoise::PoissonNoise(double a,double b,double c,double readout,double res): BaseNoise(a),Msb(b),ZP(c){
  this->res2 = res*res;
  this->Ibg = pow(10.0,-0.4*(Msb-ZP))*res2*texp + pow(readout,2); // the sky background in electrons
}

void PoissonNoise::addNoise(RectGrid* mydata){
  this->seed += 2; // increment seed at each call
  srand48(seed);
  
  double z1,z2,u1,u2,sdev,electrons;
  //Applying the Box-Muller transformation
  for(int i=0;i<mydata->Nz;i++){
    u1 = drand48();
    u2 = drand48();
    z1 = sqrt(-2.0 * log(u1)) * cos(this->two_pi * u2);
    //    z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);    

    // The units of obs_base are electrons/(s arcsec^2), so multiply it by the area of a pixel and exposure time
    electrons = mydata->z[i]*this->res2*this->texp;
    sdev = sqrt(electrons + this->Ibg);
    mydata->z[i] += z1*sdev;
  }
} 
// END: PoissonNoise ==============================




// START: UniformGaussian ============================
UniformGaussian::UniformGaussian(double sn){
  this->sn = sn;
}

void UniformGaussian::addNoise(RectGrid* mydata){
  this->seed += 2; // increment seed at each call

  double maxdata = *std::max_element(mydata->z,mydata->z+mydata->Nz);
  double sigma = maxdata/this->sn;
  srand48(seed);

  double min_noise = sigma; // just a starting value
  double z1,z2,u1,u2,noise;
  //Applying the Box-Muller transformation
  for(int i=0;i<mydata->Nz;i++){
    u1 = drand48();
    u2 = drand48();
    z1 = sqrt(-2.0 * log(u1)) * cos(this->two_pi * u2);
    //    z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);    

    noise = z1*sigma;
    mydata->z[i] += noise;
    if( noise < min_noise ){
      min_noise = noise;
    }
  }

  // Renormalize by adding the minimum (negative) noise value
  // This is needed because the pixels can contain nan's if mag scale is used
  for(int i=0;i<mydata->Nz;i++){
    mydata->z[i] += abs(min_noise);
  }
}
// END: UniformGaussian ==============================
