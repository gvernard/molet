#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct { double x, y; } vec_t, *vec;
typedef struct { int len, alloc; vec v; } poly_t, *poly;

inline double dot(vec a,vec b);
inline double cross(vec a,vec b);
inline vec vsub(vec a,vec b,vec res);
int left_of(vec a,vec b,vec c);
int line_sect(vec x0,vec x1,vec y0,vec y1,vec res);
poly poly_new();
void poly_free(poly p);
void poly_append(poly p,vec v);
int poly_winding(poly p);
void poly_edge_clip(poly sub,vec x0,vec x1,int left,poly res);
poly poly_clip(poly sub,poly clip);

