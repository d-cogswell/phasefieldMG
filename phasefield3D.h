#ifndef _PHASEFIELD3D_H_
#define _PHASEFIELD3D_H_

#include <iostream>
#include "grid3D.h"
#include "macros3D.h"
using namespace std;

//Boundary layers
//-----------------------------------------------------------------------------
#define bc 1

//Energy density function
//-----------------------------------------------------------------------------
/*I use the energy density function of Eggleston, McFadden, Voorhees:
   f(c)=1/4*W*c^2(1-c)^2 */
#define DfDphi(phi) (4*(cube(phi) - 1.5*sq(phi) + .5*phi))

//Mobility dependence on phi
//-----------------------------------------------------------------------------
/*If M(phi,i,j,h) is not defined, constant mobility is assumed.  
 *The constant mobility case is much faster to compute. */
//#define M(phi,i,j,h) sq(phi)*sq(phi-1)

//This model makes the mobility highest where the gradient is high
//#define M(phi,i,j,h) (sq(DX2D(phi,i,j,h))+sq(DY2D(phi,i,j,h)))

//#define M(phi,i,j,h) .5*((DX2D(phi,i,j,h)==0) ? 0 : (1+cos(4*atan(DY2D(phi,i,j,h)/DX2D(phi,i,j,h)))))
//-----------------------------------------------------------------------------

//This struct is used to define the component properties of an allow
// Tm=melting temperature
// W=energy hump between solid and liquid phases
// L=latent heat
//-----------------------------------------------------------------------------
struct component {
  double Tm;
  double W;
  double L;
};

//Function definitions
//-----------------------------------------------------------------------------
int cahn_hilliard3D(grid3D*,double,int,int);
int allen_cahn3D(grid3D*,double,int,int);
int liquid_solid3D(component,double,grid3D*,double,int,int);
int binary_alloy3D(component,component,double,grid3D*,grid3D*,double,int,int);
int ternary_alloy3D(component,component,component,double,grid3D*,grid3D*,grid3D*,double,int,int);
#endif
