#ifndef NOISE_HPP
#define NOISE_HPP

#define _USE_MATH_DEFINES

#include <string>
#include "json/json.h"
#include "vkllib.hpp"

class BaseNoise {
public:
  int seed = 123;
  const double two_pi = 2.0*M_PI;
  double texp; // exposure time in seconds
  vkl::RectGrid* noise_realization = NULL;

  BaseNoise(){};
  BaseNoise(double texp);
  ~BaseNoise();

  void addNoise(vkl::RectGrid* mydata);
  
  virtual void setGrid(vkl::RectGrid* obs_grid) = 0;
  virtual void initializeFromData(vkl::RectGrid* mydata) = 0;
  virtual void calculateNoise() = 0;
  virtual void outputNoiseProperties(std::string out_path,std::string instrument_name) = 0;

protected:
  void outputNoiseRealization(std::string output,std::string instrument_name);
  std::vector<double> getBoxMullerVector(int N);
};



class NoNoise: public BaseNoise {
public:
  NoNoise(){};
  void setGrid(vkl::RectGrid* obs_grid);
  void initializeFromData(vkl::RectGrid* mydata){};
  void calculateNoise(){};
  void outputNoiseProperties(std::string out_path,std::string instrument_name){};
};

class UniformGaussian: public BaseNoise {
public:
  double sn;   // signal to noise ratio
  double sigma;
  double min_noise; // the minimum value of the noise
  UniformGaussian(double sn);
  void setGrid(vkl::RectGrid* obs_grid);
  void initializeFromData(vkl::RectGrid* mydata);
  void calculateNoise();
  void outputNoiseProperties(std::string out_path,std::string instrument_name);
};

class PoissonNoise: public BaseNoise {
public:
  double Msb;  // the sky background in mag/arcsec^2
  double ZP;   // zero-point, Vega's aparent magnitude in V
  double Ibg;  // the sky background in electrons
  double res2; // pixel area in arcsec^2
  vkl::RectGrid* sigma_map = NULL;
  PoissonNoise(double texp,double Msb,double ZP,double readout,double res);
  ~PoissonNoise();
  void setGrid(vkl::RectGrid* obs_grid);
  void initializeFromData(vkl::RectGrid* mydata);
  void calculateNoise();
  void outputNoiseProperties(std::string out_path,std::string instrument_name);
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
      return new PoissonNoise(texp,Msb,instrument->ZP,instrument->readout,instrument->resolution);
    } else {
      return new NoNoise();
    }
  }

private:
  FactoryNoiseModel(){};
};


#endif /* NOISE_HPP */
