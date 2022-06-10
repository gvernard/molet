#include <string>
#include <fstream>
#include <filesystem>

#include "json/json.h"
#include "instruments.hpp"

int main(int argc,char* argv[]){
  if( argc != 3 ){
    fprintf(stderr,"Two command line arguments are required (%d given), the paths to the specs.json and the psf files!\n",argc-1);
    return 1;
  }
  
  std::string path_to_specs = argv[1];
  std::string path_to_psf   = argv[2];

  // Make sure files exist
  if( !(std::filesystem::exists(path_to_specs) && std::filesystem::exists(path_to_psf)) ){
    fprintf(stderr,"One of the given paths does not correspond to an existing file!\n");
    return 1;
  }

  // Read in json
  Json::Value specs;
  std::ifstream fin(path_to_specs,std::ifstream::in);
  fin >> specs;
  fin.close();
  std::string name = specs["name"].asString() + "-" + specs["band"].asString();

 
  // Check if instrument exists
  if( Instrument::checkInstrumentExists(name) ){
    fprintf(stderr,"Instrument %s already exists! Rename the one you try to add or use the existing one.\n",name.c_str());
    return 1;
  } else {
    Instrument::createNewInstrument(specs,path_to_psf);
  }
  
  return 0;
}
