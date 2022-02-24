// If light curves are given as input, then this checks to ensure that the computed multiple images are the same in number.

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include "json/json.h"




int main(int argc,char* argv[]){
  std::ifstream fin;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string in_path = argv[2];
  std::string out_path = argv[3];



  
  // Read the multiple images' parameters from JSON to get the maximum time delay
  Json::Value images;
  fin.open(out_path+"output/multiple_images.json",std::ifstream::in);
  fin >> images;
  fin.close();
  int n_mult = images.size();

  // Loop over instruments and check the number of given extrinsic light curves
  bool check = false;
  for(int i=0;i<root["instruments"].size();i++){
    std::string name = root["instruments"][i]["name"].asString();

    Json::Value lcs_json;
    fin.open(in_path+"/input_files/"+name+"_LC_extrinsic.json",std::ifstream::in);
    fin >> lcs_json;
    fin.close();
    int n_lc = lcs_json.size();

    std::cout << n_mult << " " << n_lc << std::endl;
    if( n_mult != n_lc ){
      check = true;
      fprintf(stderr,"For instrument %s: Number of given light curves (%d) does not match the computed number of multiple images (%d).\n",name.c_str(),n_lc,n_mult);
    }
  }
  
  if( check ){
    return 1;
  }

  return 0;
}
