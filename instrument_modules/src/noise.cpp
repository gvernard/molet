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
PoissonNoise::PoissonNoise(double a,double b,double c,double d,double readout,double res): sn(a),texp(b),Msb(c),ZP(d){
  this->Ibg = pow(10.0,-0.4*(Msb-ZP))*pow(res,2)*texp + pow(readout,2); // the sky background in electrons
}

void PoissonNoise::addNoise(RectGrid* mydata){
  
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
