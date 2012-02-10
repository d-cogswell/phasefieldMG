#ifndef _MULTIGRID3D_H_
#define _MULTIGRID3D_H_

#include "phasefield3DMG.h"

//Function definitions
//-----------------------------------------------------------------------------
double multigrid(grid3D*,grid3D&,grid3D&,grid3D&,grid3D&,double,double,int=2,int=1);
double FAS_multigrid(grid3D*,grid3D&,grid3D&,grid3D&,grid3D&,double,double,int=2,int=1);
void gaussian_elimination(grid3D&,grid3D&,grid3D&);

#endif