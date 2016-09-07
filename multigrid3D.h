#ifndef _MULTIGRID3D_H_
#define _MULTIGRID3D_H_

#include "phasefield3DMG.h"

//Function definitions
//-----------------------------------------------------------------------------
template <class system>
void multigrid(system&,system&,system&,system&,double,double,int,int=2,int=1);
template <class system>
void FAS_multigrid(system&,system&,system&,system&,double,double,int,int=2,int=1);

//This is a recursive function that relaxes the error using 'max_level' grid
//levels.   max_levels=2 corresponds to the two grid scheme.
template <class system>
void multigrid(system& u, system& f, system& d, system& e, double dt, double h, int gamma, int max_level, int level){

  //Set number of pre and post smoothing iterations
  int v1=1;
  int v2=1;

  //Presmoothing
  for (int i=0;i<v1;++i)
    GS_LEX_CH(u,f,dt,h);

  //Compute the defect and restrict it
  dfct_CH(d,u,f,dt,h);
  system& d2h=*d.restrict_FW();
  system& e2h=*e.getCoarseGrid();

  //Set the initial guess for the LHS
  e2h=0;

  //Direct solve on the coarsest mesh
  if (level==max_level-1){
    GS_LEX_CH(e2h,d2h,dt,2*h);
  }

  //Otherwise perform a coarse grid correction to solve for e2h
  else{
    for (int i=gamma;i>0;--i)
      multigrid<system>(e2h,d2h,*u.getCoarseGrid(),*f.getCoarseGrid(),dt,2*h,i,max_level,level+1);
  }
  
  //Prolongate the error to the fine mesh
  e2h.prolongate(e.N1,e.N2);

  //Compute the corrected approximation
  u+=e;

  //Postsmoothing
  for (int i=0;i<v2;++i)
    GS_LEX_CH(u,f,dt,h);
}

//This function applies the full approximation scheme for solving nonlinear problems
//-----------------------------------------------------------------------------
template <class system>
void FAS_multigrid(system& u, system& f, system& d, system& w, double dt, double h, int gamma, int max_level, int level){

  //Set number of iterations on the fine grid and coarse grid
  int v1=1;
  int v2=1;

  //Presmoothing
  for (int i=0;i<v1;++i)
    GS_RB_CH(u,f,dt,h);
  
  //Compute the defect
  dfct_CH(d,u,f,dt,h);
  
  //Restrict the defect and smoothed u
  system& d2h=*d.restrict_CC();
  system& u2h=*u.restrict_CC();
  system& w2h=*w.getCoarseGrid();

  //Compute the RHS
  system& f2h=*f.getCoarseGrid();
  d_plus_Nu_CH(f2h,d2h,u2h,dt,2*h);

  //Set the initial guess for the LHS
  w2h=u2h;

  //Direct solve on the coarsest mesh
  if (level==max_level-1){
    GS_RB_CH(w2h,f2h,dt,2*h);
  }

  //Otherwise perform a coarse grid correction
  else{
    for (int i=gamma;i>0;--i)
      FAS_multigrid<system>(w2h,f2h,d2h,u2h,dt,2*h,i,max_level,level+1);
  }

  //Compute the coarse grid correction
  d2h=w2h;
  d2h-=u2h;  
  
  //Prolongate the error to the fine mesh
  d2h.prolongate_CC(d.N1,d.N2,d.N3);

  //Compute the corrected approximation
  u+=d;

  //Postsmoothing
  for (int i=0;i<v2;++i)
    GS_RB_CH(u,f,dt,h);
}

#endif
