#ifndef TABLE_DEFINITION_HPP
#define TABLE_DEFINITION_HPP

struct mytriplet {
  int i;
  int j;
  double v;
};

struct mytable {
  int Ti;
  int Tj;
  std::vector<mytriplet> tri;
};

class InterpolationCell {
public:
  int size;
  int* ind;
  double* wei;

  InterpolationCell(int size){
    this->size = size;
    this->ind  = (int*) malloc(size*sizeof(int));
    this->wei  = (double*) malloc(size*sizeof(double));
  }
  ~InterpolationCell(){
    free(ind);
    free(wei);
  }

};


#endif /* TABLE_DEFINITION_HPP */
