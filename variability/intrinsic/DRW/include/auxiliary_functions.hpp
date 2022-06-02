#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <utility>


template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p){
  return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::map<B,A> flip_map(const std::map<A,B> &src){
  std::map<B,A> dst;
  std::transform(src.begin(),src.end(),std::inserter(dst,dst.begin()),flip_pair<A,B>);
  return dst;
}

std::vector< std::vector<double> > getDRWLightCurve(std::vector<double> time,double zs,double Mi,double lrest,double mean_mag,int N_in,int seed);
void writeLightCurves(std::vector<double> time,std::vector< std::vector<double> > all_lc_signals,std::string outfile);

#endif /* AUXILIARY_HPP */
