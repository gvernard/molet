#ifndef POINTIMAGE_HPP
#define POINTIMAGE_HPP

class pointImage {
public:
  double x;
  double y;
  double dx;
  double dy;
  double k;
  double g;
  double phig;
  double s; // just here for later
  double mag;
  double dt;
  
  pointImage(double a,double b,double c,double d,double e,double f,double g,double h,double i,double j):x(a),y(b),dx(c),dy(d),k(e),g(f),phig(g),s(h),mag(i),dt(j){}
  ~pointImage(){};
};

#endif /* POINTIMAGE_HPP */
