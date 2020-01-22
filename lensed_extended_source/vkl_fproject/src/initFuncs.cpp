#include "initFuncs.hpp"

using std::cout;
using std::endl;

void printUsage(char* argv[]){
  cout << "Usage: " << argv[0] << " [options]" << endl;
  cout << "  -h,        --help                        Show this usage information" << endl;
  cout << "  -td        --target_dir                  Set the I/O directory to an existing directory (string)" << endl;
  cout << "  -x0                                      Set the abscissa of the lens (double)" << endl;
  cout << "  -y0                                      Set the ordinate of the lens (double)" << endl;
  cout << "  -pa                                      Set the position angle of the lens (double)" << endl;
  cout << "  -q                                       Set the ellipticity of the lens (double)" << endl;
  cout << "  -b                                       Set the Einstein radius of the lens (double)" << endl;
  cout << "  -w                                       Set the width of the image in arcsec (double)" << endl;
  cout << "  -h                                       Set the height of the image in arcsec (double)" << endl;
  cout << "  -pixw                                    Set the width of the image in pixels (int)" << endl;
  cout << "  -pixh                                    Set the height of the image in pixels (int)" << endl;
  cout << "  -xs                                      Set the abscissa of the source (double)" << endl;
  cout << "  -ys                                      Set the ordinate of the source (double)" << endl;
  cout << "  -s                                       Set the size of the source (double)" << endl;
}

int setParams(int argc,char* argv[],params &parameters){
  std::string target_dir = DEFAULT_TARGET_DIR;
  double x0              = DEFAULT_X0;
  double y0              = DEFAULT_Y0;
  double pa              = DEFAULT_PA;
  double q               = DEFAULT_Q;
  double b               = DEFAULT_B;
  double w               = DEFAULT_W;
  double h               = DEFAULT_H;
  int pixw               = DEFAULT_PIXW;
  int pixh               = DEFAULT_PIXH;
  double xs              = DEFAULT_XS;
  double ys              = DEFAULT_YS;
  double s               = DEFAULT_S;

  // Print usage information when no params are passed
  //  if( 1 == argc ) {printUsage(argv);return 1;}
  
  // 'Simple' command-line parsing
  for(int i=1;i<argc;i++){

    if( strcmp(argv[i], "--target_dir") == 0 || strcmp(argv[i], "-td") == 0 ){
      fromString(argv[++i],target_dir);
      char *temp = (char*) target_dir.c_str();
      if( opendir(temp) == NULL ){
	cout << "Error: target directory does not exist" << endl;
	printUsage(argv);
	return 1;
      }
      if( target_dir.substr(target_dir.length()-1,1).compare("/") != 0 ){
	target_dir += "/";
      }

    } else if( strcmp(argv[i], "-x0") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],x0) ) {
	  cout << "Error parsing param 'x0'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-y0") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],y0) ) {
	  cout << "Error parsing param 'y0'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-pa") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],pa) ) {
	  cout << "Error parsing param 'pa'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-q") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],q) ) {
	  cout << "Error parsing param 'q'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-b") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],b) ) {
	  cout << "Error parsing param 'b'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-w") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],w) ) {
	  cout << "Error parsing param 'w'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-h") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],h) ) {
	  cout << "Error parsing param 'h'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-pixw") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],pixw) ) {
	  cout << "Error parsing param 'pixw'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-pixh") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],pixh) ) {
	  cout << "Error parsing param 'pixh'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-xs") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],xs) ) {
	  cout << "Error parsing param 'xs'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-ys") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],ys) ) {
	  cout << "Error parsing param 'ys'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "-s") == 0 ){
      if( i+1 < argc ) {
	if( !fromString(argv[++i],s) ) {
	  cout << "Error parsing param 's'" << endl;
	  printUsage(argv);
	  return 1;
	}
      } else{
	printUsage(argv);
	return 1;
      }

    } else if( strcmp(argv[i], "--help")==0 || strcmp(argv[i], "-h")==0 ) {
      printUsage(argv);
      return 1;
      
    } else {
      cout << "Unknown option '" << argv[i] << "'" << endl;
      printUsage(argv);
      return 1;
    }
  }


  parameters.target_dir = target_dir;
  parameters.x0         = x0;
  parameters.y0         = y0;
  parameters.pa         = pa;
  parameters.q          = q;
  parameters.b          = b;
  parameters.w          = w;
  parameters.h          = h;
  parameters.pixw       = pixw;
  parameters.pixh       = pixh;
  parameters.xs         = xs;
  parameters.ys         = ys;
  parameters.s          = s;

  return 0;
}
