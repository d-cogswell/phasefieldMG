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
//#define DfDphi(phi) (4*(cube(phi) - 1.5*sq(phi) + .5*phi))
#define DfDphi(phi) (cube(phi)-phi)

//Mobility dependence on phi
//-----------------------------------------------------------------------------
/*If M(phi,i,j,h) is not defined, constant mobility is assumed.  
 *The constant mobility case is much faster to compute. */
//#define M(phi,i,j,h) sq(phi)*sq(phi-1)

//This model makes the mobility highest where the gradient is high
//#define M(phi,i,j,h) (sq(DX2D(phi,i,j,h))+sq(DY2D(phi,i,j,h)))

//#define M(phi,i,j,h) .5*((DX2D(phi,i,j,h)==0) ? 0 : (1+cos(4*atan(DY2D(phi,i,j,h)/DX2D(phi,i,j,h)))))
//-----------------------------------------------------------------------------

//Function definitions
//-----------------------------------------------------------------------------
int cahn_hilliard3D(grid3D*,double,int,int);
int allen_cahn3D(grid3D*,double,int,int);
double multigrid(grid3D**,grid3D*,grid3D*,grid3D*,grid3D*,double,double,int=2,int=1);
void gaussian_elimination(grid3D*,grid3D*,grid3D*);


//-----------------------------------------------------------------------------
inline void GS_LEX_CH(grid3D* u, grid3D* f, double dt, double h){
  //Parameters
  double K=1.5;

  int Nx=u->getDimension(1);
  int Ny=u->getDimension(2);

  double D=1+dt*(K*20/(h*h*h*h)+8/sq(h));
  gridLoop3D(*u){
    int i1=(i+1)%Nx, i2=(i+2)%Nx, i_1=(i+Nx-1)%Nx, i_2=(i+Nx-2)%Nx;
    int j1=(j+1)%Ny, j2=(j+2)%Ny, j_1=(j+Ny-1)%Ny, j_2=(j+Ny-2)%Ny;

    (*u)(i,j,k)=(-dt*(
      K/(h*h*h*h)*(-8*(*u)(i1,j,k)-8*(*u)(i_1,j,k)-8*(*u)(i,j1,k)-8*(*u)(i,j_1,k)
                   +2*(*u)(i1,j1,k)+2*(*u)(i_1,j_1,k)+2*(*u)(i1,j_1,k)+2*(*u)(i_1,j1,k)
                   +(*u)(i2,j,k)+(*u)(i_2,j,k)+(*u)(i,j2,k)+(*u)(i,j_2,k))
      -2/sq(h)*((*u)(i1,j,k)+(*u)(i_1,j,k)+(*u)(i,j1,k)+(*u)(i,j_1,k))
      )+(*f)(i,j,k))/D;
  }
}
//-----------------------------------------------------------------------------
inline void dfct_CH(grid3D* d, grid3D* u, grid3D* f,double dt, double h){
  //Parameters
  double K=1.5;

  int Nx=u->getDimension(1);
  int Ny=u->getDimension(2);
  gridLoop3D(*d){
    int i1=(i+1)%Nx, i2=(i+2)%Nx, i_1=(i+Nx-1)%Nx, i_2=(i+Nx-2)%Nx;
    int j1=(j+1)%Ny, j2=(j+2)%Ny, j_1=(j+Ny-1)%Ny, j_2=(j+Ny-2)%Ny;
    (*d)(i,j,k)=(*f)(i,j,k)-(*u)(i,j,k)-dt*(
      K*(20*(*u)(i,j,k)
        -8*((*u)(i1,j,k)+(*u)(i_1,j,k)+(*u)(i,j1,k)+(*u)(i,j_1,k))
        +2*((*u)(i1,j1,k)+(*u)(i_1,j_1,k)+(*u)(i1,j_1,k)+(*u)(i_1,j1,k))
        +((*u)(i2,j,k)+(*u)(i_2,j,k)+(*u)(i,j2,k)+(*u)(i,j_2,k)))
      -2*laplacian(*u));
  }
}
//-----------------------------------------------------------------------------
inline void f_CH(grid3D* f, grid3D* u, double dt, double h){
  int Nx=u->getDimension(1);
  int Ny=u->getDimension(2);
  grid3D u_cubed(Nx,Ny,1);
  gridLoop3D(u_cubed){
    double u_val=(*u)(i,j,k);
    u_cubed(i,j,k)=cube(u_val);
  }
  u_cubed.periodicBoundary();

  gridLoop3D(*f){
    (*f)(i,j,k)=dt*(laplacian(u_cubed)-3*laplacian(*u))+(*u)(i,j,k);
  }
}
//-----------------------------------------------------------------------------
inline void L_CH(grid3D* L, int Nx, int Ny, double dt, double h){
  //Parameters
  double K=1.5;

  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      (*L)(row,row,0)+=1;
  }

  //4th derivative term
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      double fct=dt*K/(h*h*h*h);
      (*L)(row,row,0)+=20*fct;

      (*L)(row,j*Nx+(i+1)%Nx,0)+=-8*fct;
      (*L)(row,j*Nx+(i+Nx-1)%Nx,0)+=-8*fct;
      (*L)(row,((j+1)%Ny)*Nx+i,0)+=-8*fct;
      (*L)(row,((j+Ny-1)%Ny)*Nx+i,0)+=-8*fct;

      (*L)(row,j*Nx+(i+2)%Nx,0)+=fct;
      (*L)(row,j*Nx+(i+Nx-2)%Nx,0)+=fct;
      (*L)(row,((j+2)%Ny)*Nx+i,0)+=fct;
      (*L)(row,((j+Ny-2)%Ny)*Nx+i,0)+=fct;

      (*L)(row,((j+1)%Ny)*Nx+(i+1)%Nx,0)+=2*fct;
      (*L)(row,((j+1)%Ny)*Nx+(i+Nx-1)%Nx,0)+=2*fct;
      (*L)(row,((j+Ny-1)%Ny)*Nx+(i+1)%Nx,0)+=2*fct;
      (*L)(row,((j+Ny-1)%Ny)*Nx+(i+Nx-1)%Nx,0)+=2*fct;
  }

  //Laplacian term
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      double fct=-2*dt/sq(h);
      (*L)(row,row,0)+=-4*fct;
      (*L)(row,j*Nx+(i+1)%Nx,0)+=fct;
      (*L)(row,j*Nx+(i+Nx-1)%Nx,0)+=fct;
      (*L)(row,((j+1)%Ny)*Nx+i,0)+=fct;
      (*L)(row,((j+Ny-1)%Ny)*Nx+i,0)+=fct;
  }
}

//-----------------------------------------------------------------------------
inline void GS_LEX_heat_eqn(grid3D* u, grid3D* f, double dt, double h){
  double D=2*sq(h)+4*dt;
  gridLoop3D(*u){
    (*u)(i,j,k)=(dt*((*u)(i+1,j,k)+(*u)(i-1,j,k)+(*u)(i,j+1,k)+(*u)(i,j-1,k))+(*f)(i,j,k))/D;
  }
}
//-----------------------------------------------------------------------------
inline void dfct_heat_eqn(grid3D* d, grid3D* u, grid3D* f, double dt, double h){
  gridLoop3D(*d){
    (*d)(i,j,k)=(*f)(i,j,k)-(2*sq(h)+4*dt)*(*u)(i,j,k)+dt*((*u)(i+1,j,k)+(*u)(i-1,j,k)+(*u)(i,j+1,k)+(*u)(i,j-1,k));
  }
}
//-----------------------------------------------------------------------------
inline void L_heat_eqn(grid3D* L, int Nx, int Ny, double dt, double h){
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      (*L)(row,row,0)+=2*sq(h)+4*dt;
      (*L)(row,j*Nx+(i+1)%Nx,0)+=-dt;
      (*L)(row,j*Nx+(i+Nx-1)%Nx,0)+=-dt;
      (*L)(row,((j+1)%Ny)*Nx+i,0)+=-dt;
      (*L)(row,((j+Ny-1)%Ny)*Nx+i,0)+=-dt;
    }
}
//-----------------------------------------------------------------------------
inline void f_heat_eqn(grid3D* f, grid3D* u, double dt, double h){
  gridLoop3D(*f){
    (*f)(i,j,k)=dt*((*u)(i+1,j,k)
                   +(*u)(i-1,j,k)
                   +(*u)(i,j+1,k)
                   +(*u)(i,j-1,k))
                   +(2*sq(h)-4*dt)*(*u)(i,j,k);
  }
}


#endif
