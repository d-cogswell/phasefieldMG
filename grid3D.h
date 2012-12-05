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
#include <stdlib.h>
#include <math.h>
using namespace std;

#define X_DIM 1
#define Y_DIM 2
#define Z_DIM 3
#define BND_DIM 7
#define XY_PLANE 8
#define XZ_PLANE 9
#define YZ_PLANE 10

//A macro for looping over a grid3D object
//-----------------------------------------------------------------------------
#define gridLoop3D(grid) for (int i=0; i<(grid).getDimension(1); ++i) for (int j=0; j<(grid).getDimension(2); ++j) for (int k=0; k<(grid).getDimension(3); ++k)

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

class grid3D{
 public:
  grid3D(int,int,int, int=1, double=0);
  grid3D(char*,int=1, double=0, int=0, int=0, int=0, int=0, int=0, int=0);
  grid3D(const grid3D&);
  ~grid3D(void);
  void initializeRandom(double,double);
  void initializeGaussian(double);
  void initializeSphere(double);
  double mean(void);
  double sum();
  double squaredSum(int,int,int,int,int,int,int* = NULL);
  void periodicBoundary(void);
  void neumannBoundary(double);
  void dirichletBoundary(double);
  void xAxisPeriodicBoundary(int=0,int=0);
  void yAxisPeriodicBoundary(int=0,int=0);
  void zAxisPeriodicBoundary(int=0,int=0);
  void xAxisNeumannBoundary(double,int=0,int=0);
  void yAxisNeumannBoundary(double,int=0,int=0);
  void zAxisNeumannBoundary(double,int=0,int=0);
  void xAxisDirichletBoundary(double,int=0,int=0);
  void yAxisDirichletBoundary(double,int=0,int=0);
  void zAxisDirichletBoundary(double,int=0,int=0);
  grid3D* prolongate(int,int,int=1);
  grid3D* restrict();
  grid3D* injection();
  grid3D* getCoarseGrid();
  double l2_norm();
  void writeToFile(char*);
  void writeToFile(char*,int);
  void writeToFileDx(char*);
  double*** getGrid(){return(grid);}
  double* getPlane(int,int);
  inline int getDimension(int n){switch (n){case X_DIM: return(N1); case Y_DIM: return(N2); case Z_DIM: return(N3); case BND_DIM: return(boundary); case 4: return(N1_orig); case 5: return(N2_orig); case 6: return(N3_orig);};}
  grid3D* select(int,int,int,int,int,int);
  inline double& operator()(int,int,int);
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

  void neg_yAxisNeumannBoundary(double,int=0,int=0);
  void pos_yAxisNeumannBoundary(double,int=0,int=0);
  void neg_yAxisDirichletBoundary(double,int=0,int=0);
  void pos_yAxisDirichletBoundary(double,int=0,int=0);
  void allocate(int,int,int,int);
  void getPlane(double*,int,int,int=0,int=0);
  void setPlane(double*,int,int,int=0,int=0);
  int N1;
  int N2;
  int N3;
  int boundary;
  double*** grid;
  double* data;

  /* N_inc variables are used to keep track of the initial size of the grid 
   * that was read from the input file.  warned_range and warned_bndry are used
   * to insure that particular warning messages are only displayed once to
   * avoid filling output files with error messages.
   */
 private:
  int N1_orig;
  int N2_orig;
  int N3_orig;
  int warned_range;
  int warned_bndry;
};

//-----------------------------------------------------------------------------
//This function provides access to elements in the grid
double& grid3D::operator()(int i, int j, int k){
  return(grid[i][j][k]);
}

//Trilinear interpolation
double grid3D::operator()(double i, double j, double k){
  int i0=(int)i;
  int j0=(int)j;
  int k0=(int)k;

  //If the doubles are equal to integers, return the position
  if (i0==i && j0==j && k0==k)
    return(grid[i0][j0][k0]);

  double v1_k0=Lint(i-i0,grid[i0][j0][k0],grid[i0+1][j0][k0]);
  double v2_k0=Lint(i-i0,grid[i0][j0+1][k0],grid[i0+1][j0+1][k0]);
  double v1_k1=Lint(i-i0,grid[i0][j0][k0+1],grid[i0+1][j0][k0+1]);
  double v2_k1=Lint(i-i0,grid[i0][j0+1][k0+1],grid[i0+1][j0+1][k0+1]);
    
  double v_k0=Lint(j-j0,v1_k0,v2_k0);
  double v_k1=(j-j0,v1_k1,v2_k1);

  return(Lint(k-k0,v_k0,v_k1));
}

//Tricubic interpolation
double grid3D::cubic(double i, double j, double k){
  int i0=(int)i;
  int j0=(int)j;
  int k0=(int)k;

  //If the doubles are equal to integers, return the position
  if (i0==i && j0==j && k0==k)
    return(grid[i0][j0][k0]);
  
  double tij[4][4];
  for (int x=-1; x<3; ++x)
    for (int y=-1; y<3; ++y){
      tij[x+1][y+1]=Cint(k-k0,
              grid[i0+x][j0+y][k0-1],
              grid[i0+x][j0+y][k0],
              grid[i0+x][j0+y][k0+1],
              grid[i0+x][j0+y][k0+2]);   
    }

  double ui[4];
  for (int x=0; x<4; ++x)
      ui[x]=Cint(j-j0,tij[x][0],tij[x][1],tij[x][2],tij[x][3]);
  
  return(Cint(i-i0,ui[0],ui[1],ui[2],ui[3]));
}

#endif
