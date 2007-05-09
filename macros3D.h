#ifndef _MACROS_H_
#define _MACROS_H_

#include <math.h>

//Finite difference operations (Macros)
//-----------------------------------------------------------------------------
#define get(grid,i,j,k) (grid)(i,j,k)
#define DX3D(grid,i,j,k,h) (get(grid,i+1,j,k) - get(grid,i-1,j,k))/(2*(h))
#define DY3D(grid,i,j,k,h) (get(grid,i,j+1,k) - get(grid,i,j-1,k))/(2*(h))
#define DZ3D(grid,i,j,k,h) (get(grid,i,j,k+1) - get(grid,i,j,k-1))/(2*(h))
#define DXX3D(grd,i,j,k,h) (get(grd,i+1,j,k)+get(grd,i-1,j,k)-2*get(grd,i,j,k))/(sq(h))
#define DYY3D(grd,i,j,k,h) (get(grd,i,j+1,k)+get(grd,i,j-1,k)-2*get(grd,i,j,k))/(sq(h))
#define DZZ3D(grd,i,j,k,h) (get(grd,i,j,k+1)+get(grd,i,j,k-1)-2*get(grd,i,j,k))/(sq(h))
#define laplacian3D(grid,i,j,k,h) (DXX3D(grid,i,j,k,h) + DYY3D(grid,i,j,k,h) + DZZ3D(grid,i,j,k,h))
#define gradDot(grd1,grd2,i,j,k,h) (DX3D(grd1,i,j,k,h)*DX3D(grd2,i,j,k,h)+DY3D(grd1,i,j,k,h)*DY3D(grd2,i,j,k,h)+DZ3D(grd1,i,j,k,h)*DZ3D(grd2,i,j,k,h))
#define gradMag(grd1,i,j,h) sqrt(sq(DX3D(grd1,i,j,k,h))+sq(DY3D(grd1,i,j,k,h))+sq(DZ3D(grd1,i,j,k,h)))

//Use this within a gridLoop
#define laplacian(grid) laplacian3D(grid,i,j,k,h)

//Math Macros
//-----------------------------------------------------------------------------
#define sq(x) ((x)*(x))
#define cube(x) ((x)*sq(x))
#define abs(x) ((x)<0 ? (-(x)) : (x))

//Phyiscal Constants
//-----------------------------------------------------------------------------
#define PI 3.14159

#endif

