#include <cmath>

double transverse_comoving_distance(double Wm,double Wr,double Wk,double Wv,double z){
  int n = 1000;            // number of points in integrals
  double az = 1.0/(1.0+z); // 1/(1+z), the scale factor at given redshift
  
  // Radial comoving distance: do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
  double DCMR = 0.0;       // Radial comoving distance in units of c/H0
  for(int i=0;i<n;i++){
    double a    = az + (1.0-az)*(i+0.5)/n;
    double adot = sqrt( Wk + (Wm/a) + (Wr/(a*a)) + (Wv*a*a) );
    DCMR = DCMR + 1.0/(a*adot);
  }
  DCMR = (1.0-az)*DCMR/n;

  // Transverse comoving distance
  double ratio = 1.00;
  double x = sqrt(std::abs(Wk))*DCMR;
  if( x > 0.1 ){
    if( Wk > 0){
      ratio = 0.5*(exp(x)-exp(-x))/x;
    } else {
      ratio = sin(x)/x;
    }
  } else {
    double y = x*x;
    if( Wk < 0 ){
      y = -y;
    }
    ratio = 1.0 + y/6.0 + y*y/120.0;
  }
  double DCMT   = ratio*DCMR;

  return DCMT;
}
