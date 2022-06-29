#include <string>
#include <iostream>
#include "gerlumph.hpp"

int main(){
  gerlumph::MagnificationMap mymap;
  std::string path = mymap.printMapPath();
  std::cout << path << std::endl;  
  return 0;
}
