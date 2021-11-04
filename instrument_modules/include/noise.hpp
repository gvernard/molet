#ifndef NOISE_HPP
#define NOISE_HPP

#define _USE_MATH_DEFINES

#include <string>

#include "json/json.h"

class RectGrid;

class BaseNoise {
public:
  int seed = 123;
  double texp; // exposure time in seconds
  BaseNoise(){};
  BaseNoise(double texp): texp(texp){};
  ~BaseNoise(){};
  virtual void addNoise(RectGrid* mydata) = 0;
};

class NoNoise: public BaseNoise {
public:
  NoNoise();
  void addNoise(RectGrid* mydata);
};

class UniformGaussian: public BaseNoise {
public:
  const double two_pi = 2.0*M_PI;
  double sn;   // signal to noise ratio
  UniformGaussian(double sn);
  void addNoise(RectGrid* mydata);
};

class PoissonNoise: public BaseNoise {
public:
  const double two_pi = 2.0*M_PI;
  double Msb;  // the sky background in mag/arcsec^2
  double ZP;   // zero-point, Vega's aparent magnitude in V
  double Ibg;  // the sky background in electrons
  double res2; // pixel area in arcsec^2
  PoissonNoise(double texp,double Msb,double ZP,double readout,double res);
  void addNoise(RectGrid* mydata);
};

class FactoryNoiseModel{//This is a singleton class.
public:
  FactoryNoiseModel(FactoryNoiseModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryNoiseModel const&) = delete;

  static FactoryNoiseModel* getInstance(){
    static FactoryNoiseModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseNoise* createNoiseModel(Json::Value noise_pars,Instrument* instrument){
    std::string type = noise_pars["type"].asString();
    if( type == "NoNoise" ){
      return new NoNoise();
    } else if( type == "UniformGaussian" ){
      double sn = noise_pars["sn"].asDouble();
      return new UniformGaussian(sn);
    } else if( type == "PoissonNoise" ){
      double texp = noise_pars["texp"].asDouble();
      double Msb = noise_pars["Msb"].asDouble();
      double ZP = noise_pars["ZP"].asDouble();
      return new PoissonNoise(texp,Msb,ZP,instrument->readout,instrument->resolution);
    } else {
      return new NoNoise();
    }
  }

private:
  FactoryNoiseModel(){};
};


#endif /* NOISE_HPP */
