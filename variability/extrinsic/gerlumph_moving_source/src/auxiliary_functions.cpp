#include <iostream>
#include <string>
#include <ctime>

#include "json/json.h"
#include "gerlumph.hpp"

#include "auxiliary_functions.hpp"

double dateDifference(std::string start,std::string end){
  struct std::tm tm_start;
  strptime(start.c_str(),"%H:%M:%S %d-%m-%Y",&tm_start);
  struct std::tm tm_end;
  strptime(end.c_str(),"%H:%M:%S %d-%m-%Y",&tm_end);

  std::time_t x = std::mktime(&tm_start);
  std::time_t y = std::mktime(&tm_end);
  if( x != (std::time_t)(-1) && y != (std::time_t)(-1) ){
    double difference = std::difftime(y,x)/(60*60*24); // in days
    return difference;
  } else {
    return 0.0;
  }
}

BaseProfile* createProfileFromJson(Json::Value json_profile,double pixSizePhys){
  std::string profile_type = json_profile["type"].asString();

  BaseProfile* profile;
  if( profile_type == "uniform" ){
    double rhalf  = json_profile["rhalf"].asDouble();
    double incl   = json_profile["incl"].asDouble();
    double orient = json_profile["orient"].asDouble();
    profile = new UniformDisc(pixSizePhys,rhalf/0.707,incl,orient);
  } else if( profile_type == "gaussian" ){
    double rhalf  = json_profile["rhalf"].asDouble();
    double incl   = json_profile["incl"].asDouble();
    double orient = json_profile["orient"].asDouble();
    profile = new Gaussian(pixSizePhys,rhalf/1.18,incl,orient);
  } else {
    std::string filename   = json_profile["filename"].asString();
    double profPixSizePhys = json_profile["profPixSizePhys"].asDouble();
    double incl            = json_profile["incl"].asDouble();
    double orient          = json_profile["orient"].asDouble();
    profile = new Custom(pixSizePhys,filename,profPixSizePhys,incl,orient);
  }

  return profile;
}
