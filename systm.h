#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "grid3D.h"

class systm {
 public:
  int N1,N2,N3;
  grid3D phi;
  grid3D mu;
  systm* coarse;
  systm* fine;
  systm(int,int,int);
  systm(grid3D*,grid3D*);
  ~systm();
  double l2_norm();
  systm* prolongate(int,int,int=1);
  systm* prolongate_CC(int,int,int=1);
  systm* restrict_FW();
  systm* restrict_CC();
  systm* injection();
  void dirichletBoundary(double);
  void neumannBoundary(double);
  void periodicBoundary();
  systm* getCoarseGrid();
  void operator=(const double);
  void operator=(systm&);
  void operator+=(systm&);
  void operator-=(systm&);
};

#endif