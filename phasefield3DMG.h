#ifndef _PHASEFIELD3DMG_H_
#define _PHASEFIELD3DMG_H_

#include <Magick++.h>
#include "grid3D.h"
using namespace Magick;

inline double clip(double n, double low, double high){
  return(n<low ? low : (n>high ? high : n));
}

//Function definitions
//-----------------------------------------------------------------------------
void GS_LEX_CH(grid3D&, grid3D&, double, double);
void dfct_CH(grid3D&, grid3D&, grid3D&, double, double);
void f_CH(grid3D&, grid3D&, double, double);
void L_CH(grid3D&, int, int, double, double);

void GS_LEX_heat_eqn(grid3D&, grid3D&, double, double);
void dfct_heat_eqn(grid3D&, grid3D&, grid3D&, double, double);
void L_heat_eqn(grid3D&, int, int, double, double);
void f_heat_eqn(grid3D&, grid3D&, double, double);

#endif