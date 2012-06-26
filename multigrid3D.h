#ifndef _MULTIGRID3D_H_
#define _MULTIGRID3D_H_

#include "phasefield3DMG.h"

//Function definitions
//-----------------------------------------------------------------------------
template <class system>
double multigrid(grid3D*,system&,system&,system&,system&,double,double,int=2,int=1);
template <class system>
double FAS_multigrid(grid3D*,system&,system&,system&,system&,double,double,int=2,int=1);
void gaussian_elimination(grid3D&,grid3D&,grid3D&);

//This is a recursive function that relaxes the error using 'max_level' grid
//levels.   max_levels=2 corresponds to the two grid scheme.
template <class system>
double multigrid(grid3D* L, system& u, system& f, system& d, system& e, double dt, double h, int max_level, int level){

  //Set number of pre and post smoothing iterations
  int v1=4;
  int v2=4;

  //Presmoothing
  for (int i=0;i<v1;++i){
    GS_LEX_CH(u,f,dt,h);
  }

  //Compute the defect and restrict it
  dfct_CH(d,u,f,dt,h);
  system& d2h=*d.restrict();
  system& e2h=*e.getCoarseGrid();

  //Direct solve on the coarsest mesh
  if (level==max_level){

     //Get the grid dimensions for the current level
    int Nx=e2h.getDimension(1);
    int Ny=e2h.getDimension(2);

    //Only allocate space for L once
    if (L==NULL)
      L = new grid3D(Nx*Ny,Nx*Ny,1);

    L_CH(*L,Nx,Ny,dt,2*h);
    gaussian_elimination(*L,e2h,d2h);
  }

  //Otherwise perform a coarse grid correction to get e2h
  else
    multigrid<system>(L,*u.getCoarseGrid(),*f.getCoarseGrid(),d2h,e2h,dt,2*h,max_level,level+1);

  //Prolongate the error to the fine mesh
  e2h.prolongate(e.getDimension(1),e.getDimension(2));

  //Compute the corrected approximation
  gridLoop3D(u)
    u(i,j,k)+=e(i,j,k);

  //Postsmoothing
  for (int i=0;i<v2;++i)
    GS_LEX_CH(u,f,dt,h);

  if (level==1){
    dfct_CH(d,u,f,dt,h);
    return(d.l2_norm());
  }
}

//This function applies the full approximation scheme for solving nonlinear problems
//-----------------------------------------------------------------------------
template <class system>
double FAS_multigrid(grid3D* L, system& u, system& f, system& d, system& e, double dt, double h, int max_level, int level){

  //Set number of iterations on the fine grid and coarse grid
  int v1=4;
  int v2=4;

  //Presmoothing
  for (int i=0;i<v1;++i)
    GS_LEX_AC(u,f,dt,h);

  //Compute the defect
  dfct_AC(d,u,f,dt,h);

  //Restrict the defect and smoothed u
  system& d2h=*d.restrict();
  system& u2h=*u.restrict();

  //Compute the RHS
  system& f2h=*f.getCoarseGrid();
  system& e2h=*e.getCoarseGrid();
  d_plus_Nu_AC(f2h,d2h,u2h,dt,2*h);

  //Direct solve on the coarsest mesh
  if (level==max_level){

    //Get the grid dimensions for the current level
    int Nx=f2h.getDimension(1);
    int Ny=f2h.getDimension(2);

    //Only allocate space for L once
    if (L==NULL)
      L = new grid3D(Nx*Ny,Nx*Ny,1);

    L_AC(*L,u2h,f2h,Nx,Ny,dt,2*h);
    gaussian_elimination(*L,e2h,f2h);

    //Compute the coarse grid correction
    gridLoop3D(e2h)
      e2h(i,j,k)-=u2h(i,j,k);
  }

  //Otherwise perform a coarse grid correction
  else
    FAS_multigrid<system>(L,u2h,f2h,d2h,e2h,dt,2*h,max_level,level+1);

  //Prolongate the error to the fine mesh
  e2h.prolongate(e.getDimension(1),e.getDimension(2));

  //Compute the corrected approximation
  gridLoop3D(u)
    u(i,j,k)+=e(i,j,k);

  //Postsmoothing
  for (int i=0;i<v2;++i)
    GS_LEX_AC(u,f,dt,h);

  if (level==1){
    dfct_AC(d,u,f,dt,h);
    return(d.l2_norm());
  }
}

#endif