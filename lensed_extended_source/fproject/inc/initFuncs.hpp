#ifndef INIT_FUNCS_HPP
#define INIT_FUNCS_HPP

#define DEFAULT_TARGET_DIR  ""
#define DEFAULT_X0          0.0
#define DEFAULT_Y0          0.0
#define DEFAULT_PA          30.0
#define DEFAULT_Q           0.7
#define DEFAULT_B           1.3
#define DEFAULT_W           4.
#define DEFAULT_H           4.
#define DEFAULT_PIXW        200
#define DEFAULT_PIXH        200
#define DEFAULT_XS          0.0
#define DEFAULT_YS          0.0
#define DEFAULT_S           0.5


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include <cstdlib>
#include <cmath>

typedef struct{
  std::string target_dir;
  double x0;
  double y0;
  double pa;
  double q;
  double b;
  double w;
  double h;
  int pixw;
  int pixh;
  double xs;
  double ys;
  double s;
} params;


template<typename T>
inline bool fromString(const std::string &s, T &t){
  std::istringstream ss(s);
  ss >> t;
  if( ss.fail() ){
    return false;
  } else {
    return true;
  }
}

void printUsage(char* argv[]);
int  setParams(int argc,char* argv[],params &parameters);

#endif /* INIT_FUNCS_HPP */
