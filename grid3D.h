/* 
 * Copyright (c) 2006, Daniel A. Cogswell
 * Department of Materials Science and Engineering
 * Massachusetts Institute of Technology
 * 
 * This program not intended for distribution without the consent of the author.
 */

#ifndef _GRID3D_H_
#define _GRID3D_H_

#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;

//A macro for looping over a grid3D object
//-----------------------------------------------------------------------------
#define gridLoop3D(grid) for (int i=0; i<(grid).N1; ++i) for (int j=0; j<(grid).N2; ++j) for (int k=0; k<(grid).N3; ++k)
#define bndryGridLoop3D(grid) for (int i=-boundary; i<(grid).N1+boundry; ++i) for (int j=-boundary; j<(grid).N2+boundry; ++j) for (int k=-boundary; k<(grid).N3+boundry; ++k)

//Math Macros
//-----------------------------------------------------------------------------
#define sq(x) ((x)*(x))
#define cube(x) ((x)*sq(x))

//Phyiscal Constants
//-----------------------------------------------------------------------------
#define PI 3.14159

//Function prototypes
//-----------------------------------------------------------------------------
double Lint(double,double,double);
double Cint(double,double,double,double,double);
double Cint_fwd(double,double,double,double,double);

class grid3D{
 public:
  int N1,N2,N3,boundary;
  grid3D(int,int,int, int=1, double=0);
  grid3D(const char*,int=1);
  ~grid3D(void);
  void initializeRandom(double,double);
  void initializeGaussian(double);
  void initializeSphere(double);
  double mean(void);
  double sum();
  void periodicBoundary(void);
  void neumannBoundary(double);
  void dirichletBoundary(double);
  void xAxisPeriodicBoundary(int=0,int=0);
  void yAxisPeriodicBoundary(int=0,int=0);
  void zAxisPeriodicBoundary(int=0,int=0);
  void xAxisNeumannBoundary(double,int=1,int=1);
  void yAxisNeumannBoundary(double,int=1,int=1);
  void zAxisNeumannBoundary(double,int=1,int=1);
  void xAxisDirichletBoundary(double,int=1,int=1);
  void yAxisDirichletBoundary(double,int=1,int=1);
  void zAxisDirichletBoundary(double,int=1,int=1);
  grid3D* prolongate(int,int,int=1);
  grid3D* prolongate_cubic(int,int,int=1);
  grid3D* restrict_FW();
  grid3D* restrict_HW();
  grid3D* injection();
  grid3D* getCoarseGrid();
  double l2_norm();
  void writeToFile(const char*);
  void writeToFile(const char*,int);
  void writeToFileDx(const char*);
  double*** getGrid(){return(grid);}
  inline double& operator()(int,int,int);
  inline double operator()(double,int,int);
  inline double operator()(int,double,int);
  inline double operator()(int,int,double);
  inline double operator()(double,double,int);
  inline double operator()(double,double,double);
  inline double cubic(double,double,double);
  
  //Finite difference operators
  inline double DX(int i,int j,int k,double h){return((*this)(i+1,j,k)-(*this)(i-1,j,k))/(2*h);}
  inline double DY(int i,int j,int k,double h){return((*this)(i,j+1,k)-(*this)(i,j-1,k))/(2*h);}
  inline double DZ(int i,int j,int k,double h){return((*this)(i,j,k+1)-(*this)(i,j,k-1))/(2*h);}

  inline double DXX(int i,int j,int k,double h) {return(((*this)(i+1,j,k)+(*this)(i-1,j,k)-2*(*this)(i,j,k))/sq(h));}
  inline double DYY(int i,int j,int k,double h) {return(((*this)(i,j+1,k)+(*this)(i,j-1,k)-2*(*this)(i,j,k))/sq(h));}
  inline double DZZ(int i,int j,int k,double h) {return(((*this)(i,j,k+1)+(*this)(i,j,k-1)-2*(*this)(i,j,k))/sq(h));}

  inline double laplacian(int i, int j, int k,double h) {return(DXX(i,j,k,h)+DYY(i,j,k,h)+DZZ(i,j,k,h));}

  void operator=(const double);
  void operator=(grid3D&);
  void operator+=(grid3D&);
  void operator-=(grid3D&);
  grid3D *coarse, *fine;

 protected:
  double cubic_x(double,double,double);
  double cubic_y(double,double,double);
  double cubic_z(double,double,double);
  void neg_yAxisDirichletBoundary(double,int=0,int=0);
  void pos_yAxisDirichletBoundary(double,int=0,int=0);
  void allocate(int,int,int,int);
  double*** grid;
  double* data;
};

//-----------------------------------------------------------------------------
//This function provides access to elements in the grid
double& grid3D::operator()(int i, int j, int k){
  return(grid[i][j][k]);
}

double grid3D::operator()(double i, int j, int k){
    int i0=floor(i);
    return(i==i0 ? operator()(i0,j,k) : Lint(i-i0,operator()(i0,j,k),operator()(i0+1,j,k)));
}

double grid3D::operator()(int i, double j, int k){
    int j0=floor(j);
    return(j==j0 ? operator()(i,j0,k) : Lint(j-j0,operator()(i,j0,k),operator()(i,j0+1,k)));
}

double grid3D::operator()(int i, int j, double k){
    int k0=floor(k);
    return(k==k0 ? operator()(i,j,k0) : Lint(k-k0,operator()(i,j,k0),operator()(i,j,k0+1)));
}

double grid3D::operator()(double i, double j, int k){
  int i0=floor(i);
  int j0=floor(j);
  
  if (i==i0)
    return(operator()(i0,j,k));
  else if (j==j0)
    return(operator()(i,j0,k));
  
  return(Lint(j-j0,operator()(i,j0,k),operator()(i,j0+1,k)));
}

//Trilinear interpolation
double grid3D::operator()(double i, double j, double k){
  int i0=floor(i);
  int j0=floor(j);
  int k0=floor(k);

  if (k==k0)
    return(operator()(i,j,k0));

  return(Lint(k-k0,operator()(i,j,k0),operator()(i,j,k0+1)));
}

//Tricubic interpolation
double grid3D::cubic(double i, double j, double k){
  int i0=(int)i;
  int j0=(int)j;
  int k0=(int)k;

  //If the doubles are equal to integers, return the position
  if (i0==i && j0==j && k0==k)
    return(grid[i0][j0][k0]);
  else if (i0==i && j0==j)
    return(cubic_z(i,j,k));
  else if (i0==i)
    return(cubic_y(i,j,k));
  else
    return(cubic_x(i,j,k));
}

#endif
