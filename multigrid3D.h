#ifndef _MULTIGRID3D_H_
#define _MULTIGRID3D_H_

#include "phasefield3DMG.h"

//Function definitions
//-----------------------------------------------------------------------------
template <class system>
void multigrid(grid3D**,system&,system&,double,double,int=2,int=1);
template <class system>
void FAS_multigrid(grid3D**,system&,system&,double,double,int=2,int=1);
void gaussian_elimination(grid3D&,grid3D&,grid3D&);

//This is a recursive function that relaxes the error using 'max_level' grid
//levels.   max_levels=2 corresponds to the two grid scheme.
template <class system>
void multigrid(grid3D** L, system& u, system& f, double dt, double h, int max_level, int level){

  //Set number of pre and post smoothing iterations
  int v1=1;
  int v2=1;

  //Presmoothing
  for (int i=0;i<v1;++i)
    GS_LEX_CH(u,f,dt,h);

  //Create temporary storage space
  int Nx=u.getDimension(1);
  int Ny=u.getDimension(2);
  int Nz=u.getDimension(3);
  system d(Nx,Ny,Nz);
  system e(Nx,Ny,Nz);

  //Compute the defect and restrict it
  dfct_CH(d,u,f,dt,h);
  system& d2h=*d.restrict();
  system& e2h=*e.getCoarseGrid();

  //Direct solve on the coarsest mesh
  if (level==max_level-1){

     //Get the grid dimensions for the current level
    int Nx=e2h.getDimension(1);
    int Ny=e2h.getDimension(2);

    //Only allocate space for L once
    if (*L==NULL)
      *L = new grid3D(Nx*Ny,Nx*Ny,1);

    L_CH(**L,Nx,Ny,dt,2*h);
    gaussian_elimination(**L,e2h,d2h);
  }

  //Otherwise perform a coarse grid correction to solve for e2h
  else
    multigrid<system>(L,e2h,d2h,dt,2*h,max_level,level+1);
  
  //Prolongate the error to the fine mesh
  e2h.prolongate(e.getDimension(1),e.getDimension(2));

  //Compute the corrected approximation
  u+=e;

  //Postsmoothing
  for (int i=0;i<v2;++i)
    GS_LEX_CH(u,f,dt,h);
}

//This function applies the full approximation scheme for solving nonlinear problems
//-----------------------------------------------------------------------------
template <class system>
void FAS_multigrid(grid3D** L, system& u, system& f, double dt, double h, int max_level, int level){

  //Set number of iterations on the fine grid and coarse grid
  int v1=1;
  int v2=1;

  //Presmoothing
  for (int i=0;i<v1;++i)
    GS_LEX_CH(u,f,dt,h);

  int Nx=u.getDimension(1);
  int Ny=u.getDimension(2);
  int Nz=u.getDimension(3);
  system d(Nx,Ny,Nz);
  system e(Nx,Ny,Nz);
  
  //Compute the defect
  dfct_CH(d,u,f,dt,h);
  
  //Restrict the defect and smoothed u
  system& d2h=*d.restrict();
  system& u2h=*u.restrict();

  //Compute the RHS
  system& f2h=*f.getCoarseGrid();
  system& e2h=*e.getCoarseGrid();
  d_plus_Nu_CH(f2h,d2h,u2h,dt,2*h);

  //Direct solve on the coarsest mesh
  if (level==max_level-1){

    //Get the grid dimensions for the current level
    int Nx=f2h.getDimension(1);
    int Ny=f2h.getDimension(2);

    //Only allocate space for L once
    if (*L==NULL)
      *L = new grid3D(Nx*Ny,Nx*Ny,1);

    L_CH(**L,Nx,Ny,dt,2*h);
    gaussian_elimination(**L,e2h,f2h);
  }

  //Otherwise perform a coarse grid correction
  else{
    e2h=u2h;
    FAS_multigrid<system>(L,e2h,f2h,dt,2*h,max_level,level+1);
  }

  //Compute the coarse grid correction
  e2h-=u2h;
  
  //Prolongate the error to the fine mesh
  e2h.prolongate(e.getDimension(1),e.getDimension(2));

  //Compute the corrected approximation
  u+=e;

  //Postsmoothing
  for (int i=0;i<v2;++i)
    GS_LEX_CH(u,f,dt,h);
}

#endif
