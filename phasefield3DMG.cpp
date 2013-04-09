#include "phasefield3DMG.h"

double kappa=1.5;

/*The following functions solve the Allen-Cahn equation using Eyre's 
 *linearly stabilized splitting.*/
//-----------------------------------------------------------------------------
void GS_LEX_AC(grid3D& u, grid3D& f, double dt, double h){
  double D=dt*kappa/(h*h);
  u.periodicBoundary();
  gridLoop3D(u){
    u(i,j,k)=(f(i,j,k)+D*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)))/(1+2*dt+4*D);
  }
  u.periodicBoundary();
}
//-----------------------------------------------------------------------------
void dfct_AC(grid3D& d, grid3D& u, grid3D& f, double dt, double h){
  gridLoop3D(d){
    d(i,j,k)=f(i,j,k)-(u(i,j,k)+dt*(2*u(i,j,k)-kappa*u.laplacian(i,j,k,h)));
  }
}
//-----------------------------------------------------------------------------
void d_plus_Nu_AC(grid3D& f, grid3D& d, grid3D& u,double dt, double h){
  u.periodicBoundary();
  gridLoop3D(f){
    f(i,j,k)=d(i,j,k)+u(i,j,k)+dt*(2*u(i,j,k)-kappa*u.laplacian(i,j,k,h));
  }
}
//-----------------------------------------------------------------------------
void f_AC(grid3D& f, grid3D& u, double dt, double h){ 
  gridLoop3D(f){
    f(i,j,k)=u(i,j,k)+dt*(3*u(i,j,k)-cube(u(i,j,k)));
  }
}
//-----------------------------------------------------------------------------
void L_AC(grid3D& L, int Nx, int Ny, double dt, double h){
  double D=dt*kappa/(h*h);
  L=0;

  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      L(row,row,0)+=1+2*dt+4*D;
      L(row,j*Nx+(i+1)%Nx,0)+=-D;
      L(row,j*Nx+(i+Nx-1)%Nx,0)+=-D;
      L(row,((j+1)%Ny)*Nx+i,0)+=-D;
      L(row,((j+Ny-1)%Ny)*Nx+i,0)+=-D;
  }
}


/*The following functions solve the Cahn-Hilliard equation using Eyre's 
 *linearly stabilized splitting from "An Unconditionally Stable One-Step Scheme
 *for Gradient Systems".*/
//-----------------------------------------------------------------------------
void GS_LEX_CH(grid3D& u, grid3D& f, double dt, double h){
  int Nx=u.getDimension(1);
  int Ny=u.getDimension(2);

  double D=1+dt*(kappa*20/(h*h*h*h)+8/sq(h));
  gridLoop3D(u){
    int i1=(i+1)%Nx, i2=(i+2)%Nx, i_1=(i+Nx-1)%Nx, i_2=(i+Nx-2)%Nx;
    int j1=(j+1)%Ny, j2=(j+2)%Ny, j_1=(j+Ny-1)%Ny, j_2=(j+Ny-2)%Ny;

    u(i,j,k)=(-dt*(
      kappa/(h*h*h*h)*(-8*u(i1,j,k)-8*u(i_1,j,k)-8*u(i,j1,k)-8*u(i,j_1,k)
                   +2*u(i1,j1,k)+2*u(i_1,j_1,k)+2*u(i1,j_1,k)+2*u(i_1,j1,k)
                   +u(i2,j,k)+u(i_2,j,k)+u(i,j2,k)+u(i,j_2,k))
      -2/sq(h)*(u(i1,j,k)+u(i_1,j,k)+u(i,j1,k)+u(i,j_1,k))
      )+f(i,j,k))/D;
  }
  u.periodicBoundary();  
}
//-----------------------------------------------------------------------------
void dfct_CH(grid3D& d, grid3D& u, grid3D& f,double dt, double h){
  int Nx=u.getDimension(1);
  int Ny=u.getDimension(2);
  gridLoop3D(d){
    int i1=(i+1)%Nx, i2=(i+2)%Nx, i_1=(i+Nx-1)%Nx, i_2=(i+Nx-2)%Nx;
    int j1=(j+1)%Ny, j2=(j+2)%Ny, j_1=(j+Ny-1)%Ny, j_2=(j+Ny-2)%Ny;
    d(i,j,k)=f(i,j,k)-u(i,j,k)-dt*(
      kappa/(h*h*h*h)*(20*u(i,j,k)
        -8*(u(i1,j,k)+u(i_1,j,k)+u(i,j1,k)+u(i,j_1,k))
        +2*(u(i1,j1,k)+u(i_1,j_1,k)+u(i1,j_1,k)+u(i_1,j1,k))
        +(u(i2,j,k)+u(i_2,j,k)+u(i,j2,k)+u(i,j_2,k)))
      -2*u.laplacian(i,j,k,h));
  }
}
//-----------------------------------------------------------------------------
void d_plus_Nu_CH(grid3D& f, grid3D& d, grid3D& u, double dt, double h){
  int Nx=f.getDimension(1);
  int Ny=f.getDimension(2);
  u.periodicBoundary();
  gridLoop3D(f){
    int i1=(i+1)%Nx, i2=(i+2)%Nx, i_1=(i+Nx-1)%Nx, i_2=(i+Nx-2)%Nx;
    int j1=(j+1)%Ny, j2=(j+2)%Ny, j_1=(j+Ny-1)%Ny, j_2=(j+Ny-2)%Ny;
    f(i,j,k)=d(i,j,k)+u(i,j,k)+dt*(
      kappa/(h*h*h*h)*(20*u(i,j,k)
        -8*(u(i1,j,k)+u(i_1,j,k)+u(i,j1,k)+u(i,j_1,k))
        +2*(u(i1,j1,k)+u(i_1,j_1,k)+u(i1,j_1,k)+u(i_1,j1,k))
        +(u(i2,j,k)+u(i_2,j,k)+u(i,j2,k)+u(i,j_2,k)))
      -2*u.laplacian(i,j,k,h));
  }
}
//-----------------------------------------------------------------------------
void f_CH(grid3D& f, grid3D& u, double dt, double h){
  int Nx=u.getDimension(1);
  int Ny=u.getDimension(2);
  grid3D u_cubed(Nx,Ny,1);
  gridLoop3D(u_cubed){
    double u_val=u(i,j,k);
    u_cubed(i,j,k)=cube(u_val);
  }
  u_cubed.periodicBoundary();

  gridLoop3D(f){
    (f)(i,j,k)=dt*(u_cubed.laplacian(i,j,k,h)-3*u.laplacian(i,j,k,h))+u(i,j,k);
  }
}
//-----------------------------------------------------------------------------
void L_CH(grid3D& L, int Nx, int Ny, double dt, double h){
  L=0;

  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      L(row,row,0)+=1;
  }

  //4th derivative term
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      double fct=dt*kappa/(h*h*h*h);
      L(row,row,0)+=20*fct;

      L(row,j*Nx+(i+1)%Nx,0)+=-8*fct;
      L(row,j*Nx+(i+Nx-1)%Nx,0)+=-8*fct;
      L(row,((j+1)%Ny)*Nx+i,0)+=-8*fct;
      L(row,((j+Ny-1)%Ny)*Nx+i,0)+=-8*fct;

      L(row,j*Nx+(i+2)%Nx,0)+=fct;
      L(row,j*Nx+(i+Nx-2)%Nx,0)+=fct;
      L(row,((j+2)%Ny)*Nx+i,0)+=fct;
      L(row,((j+Ny-2)%Ny)*Nx+i,0)+=fct;

      L(row,((j+1)%Ny)*Nx+(i+1)%Nx,0)+=2*fct;
      L(row,((j+1)%Ny)*Nx+(i+Nx-1)%Nx,0)+=2*fct;
      L(row,((j+Ny-1)%Ny)*Nx+(i+1)%Nx,0)+=2*fct;
      L(row,((j+Ny-1)%Ny)*Nx+(i+Nx-1)%Nx,0)+=2*fct;
  }

  //Laplacian term
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      double fct=-2*dt/sq(h);
      L(row,row,0)+=-4*fct;
      L(row,j*Nx+(i+1)%Nx,0)+=fct;
      L(row,j*Nx+(i+Nx-1)%Nx,0)+=fct;
      L(row,((j+1)%Ny)*Nx+i,0)+=fct;
      L(row,((j+Ny-1)%Ny)*Nx+i,0)+=fct;
  }
}

/*The following functions solve the heat equation using Crank-Nicolson.*/
//-----------------------------------------------------------------------------
void GS_LEX_heat_eqn(grid3D& u, grid3D& f, double dt, double h){
  u.periodicBoundary();
  double D=dt/sq(h);
  gridLoop3D(u){
    u(i,j,k)=(f(i,j,k)+D*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)))/(2+4*D);
  }
  u.periodicBoundary();
}
//-----------------------------------------------------------------------------
void dfct_heat_eqn(grid3D& d, grid3D& u, grid3D& f, double dt, double h){
  gridLoop3D(d){
    d(i,j,k)=f(i,j,k)-(2*u(i,j,k)-dt*u.laplacian(i,j,k,h));
  }
}
//-----------------------------------------------------------------------------
void d_plus_Nu_heat_eqn(grid3D& f, grid3D& d, grid3D& u,double dt, double h){
  u.periodicBoundary();
  gridLoop3D(f){
    f(i,j,k)=d(i,j,k)+2*u(i,j,k)-dt*u.laplacian(i,j,k,h);
  }
}
//-----------------------------------------------------------------------------
void L_heat_eqn(grid3D& L, int Nx, int Ny, double dt, double h){
  L=0;
  double D=dt/sq(h);
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      L(row,row,0)+=2+4*D;
      L(row,j*Nx+(i+1)%Nx,0)+=-D;
      L(row,j*Nx+(i+Nx-1)%Nx,0)+=-D;
      L(row,((j+1)%Ny)*Nx+i,0)+=-D;
      L(row,((j+Ny-1)%Ny)*Nx+i,0)+=-D;
    }
}
//-----------------------------------------------------------------------------
void f_heat_eqn(grid3D& f, grid3D& u, double dt, double h){
  gridLoop3D(f){
    f(i,j,k)=2*u(i,j,k)+dt*u.laplacian(i,j,k,h);
  }
}