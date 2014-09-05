#ifndef _PHASEFIELD3DMG_H_
#define _PHASEFIELD3DMG_H_
#include "grid3D.h"

inline double clip(double n, double low, double high){
  return(n<low ? low : (n>high ? high : n));
}

//Function definitions
//-----------------------------------------------------------------------------
void GS_LEX_AC(grid3D&, grid3D&, double, double);
void dfct_AC(grid3D&, grid3D&, grid3D&, double, double);
void d_plus_Nu_AC(grid3D&, grid3D&, grid3D&, double, double);
void f_AC(grid3D&, grid3D&, double, double);
void L_AC(grid3D&, int, int, double, double);

void GS_LEX_CH(grid3D&, grid3D&, double, double);
void dfct_CH(grid3D&, grid3D&, grid3D&, double, double);
void d_plus_Nu_CH(grid3D&, grid3D&, grid3D&, double, double);
void f_CH(grid3D&, grid3D&, double, double);
void L_CH(grid3D&, int, int, double, double);

void GS_LEX_heat_eqn(grid3D&, grid3D&, double, double);
void dfct_heat_eqn(grid3D&, grid3D&, grid3D&, double, double);
void d_plus_Nu_heat_eqn(grid3D&, grid3D&, grid3D&, double, double);
void L_heat_eqn(grid3D&, int, int, double, double);
void f_heat_eqn(grid3D&, grid3D&, double, double);

//Anisotropy: 1+e4*cos(4*theta)
#define e4 .025
//-----------------------------------------------------------------------------
inline double epsilon(grid3D& xi, int i, int j){
  double XI=xi(i,j,0);
  double e=1;
  if (XI>.01 && XI<.99){
    double dx=xi.DX(i,j,0,1);
    double dy=xi.DY(i,j,0,1);
    double cos2_t=sq(dx)/(sq(dx)+sq(dy));
    double cos4t=8*sq(cos2_t)-8*cos2_t+1;
    e=1+e4*cos4t;
  }
  return(e);
}

inline double depsilon_dtheta(grid3D& xi, int i, int j){
  double XI=xi(i,j,0);
  double de_dtheta=0;
  if (XI>.01 && XI<.99){
    double dx=xi.DX(i,j,0,1);
    double dy=xi.DY(i,j,0,1);
    double sin_4t=4*dx*dy*(sq(dx)-sq(dy))/(sq(dx)+sq(dy));
    de_dtheta=-4*e4*sin_4t;
  }
  return(de_dtheta);
}

#endif