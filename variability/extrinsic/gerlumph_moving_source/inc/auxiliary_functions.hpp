#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <string>
#include "json/json.h"
#include "gerlumph.hpp"

int filterMaxVelTot(std::vector<double> vtot);
double dateDifference(std::string start,std::string end);
BaseProfile* createProfileFromJson(Json::Value profile,double pixSizePhys);

#endif /* AUXILIARY_HPP */
