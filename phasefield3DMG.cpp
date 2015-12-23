#include "phasefield3DMG.h"

double gama=1.1547;
double delta=3.4641;
double kappa=3*delta*gama/2;
double W=12*gama/delta;

inline double dfdphi_c(double x){
	return(W*.5*cube(2*x-1));
}

inline double d2fdphi2_c(double x){
	return(W*3*sq(2*x-1));
}

inline double dfdphi_e(double x){
	return(W*(.5-x));
}

/*The following functions solve the Allen-Cahn equation using Eyre's 
 *nonlinearly stabilized splitting.
 */
//-----------------------------------------------------------------------------
inline void GS_AC_update(int i, int j, int k, grid3D& u, grid3D& f, double dt, double h){
  u(i,j,k)=(f(i,j,k)-dt*dfdphi_c(u(i,j,k))+dt*d2fdphi2_c(u(i,j,k))*u(i,j,k)
          +dt*kappa*laplacian_RHS(u,i,j,k,h))/(1+dt*d2fdphi2_c(u(i,j,k))+6*dt*kappa/sq(h));
}
//-----------------------------------------------------------------------------
void GS_RB_AC(grid3D& u, grid3D& f, double dt, double h){
  u.periodicBoundary();
  
  //Red
  #pragma omp parallel for collapse(2)
  for (int i=0; i<u.N1; ++i){
    for (int j=0; j<u.N2; ++j)
      for (int k=(i+j)%2; k<u.N3; k+=2)
        GS_AC_update(i,j,k,u,f,dt,h);
  }

  //Black
  #pragma omp parallel for collapse(2)
  for (int i=0; i<u.N1; ++i){
    for (int j=0; j<u.N2; ++j)
      for (int k=(i+j+1)%2; k<u.N3; k+=2)
        GS_AC_update(i,j,k,u,f,dt,h);
  }
}
//-----------------------------------------------------------------------------
inline double AC_LHS(grid3D& u, double dt, double h, int i, int j, int k){
  return(u(i,j,k)+dt*(dfdphi_c(u(i,j,k))-kappa*u.laplacian(i,j,k,h)));
}
//-----------------------------------------------------------------------------
void dfct_AC(grid3D& d, grid3D& u, grid3D& f, double dt, double h){
  u.periodicBoundary();
  gridLoop3D(d){
    d(i,j,k)=f(i,j,k)-AC_LHS(u,dt,h,i,j,k);
  }
  d.periodicBoundary();
}
//-----------------------------------------------------------------------------
void d_plus_Nu_AC(grid3D& f, grid3D& d, grid3D& u,double dt, double h){
  u.periodicBoundary();
  gridLoop3D(f){
    f(i,j,k)=d(i,j,k)+AC_LHS(u,dt,h,i,j,k);
  }
}
//-----------------------------------------------------------------------------
void f_AC(grid3D& f, grid3D& u, double dt, double h){ 
  u.periodicBoundary();
  gridLoop3D(f){
    f(i,j,k)=u(i,j,k)-dt*dfdphi_e(u(i,j,k));
  }
}
//-----------------------------------------------------------------------------
void L_AC(grid3D& L, grid3D& u, grid3D& f, int Nx, int Ny, double dt, double h){
  double D=dt*kappa/(h*h);
  L=0;

  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      f(i,j,0)+=-dt*dfdphi_c(u(i,j,0))+dt*d2fdphi2_c(u(i,j,0))*u(i,j,0);
      L(row,row,0)+=1+dt*d2fdphi2_c(u(i,j,0))+4*D;
      L(row,j*Nx+(i+1)%Nx,0)+=-D;
      L(row,j*Nx+(i+Nx-1)%Nx,0)+=-D;
      L(row,((j+1)%Ny)*Nx+i,0)+=-D;
      L(row,((j+Ny-1)%Ny)*Nx+i,0)+=-D;
  }
}


/*The following functions solve the Cahn-Hilliard equation using Eyre's 
 *nonlinearly stabilized splitting from "An Unconditionally Stable One-Step 
 * Scheme for Gradient Systems".
 */
//-----------------------------------------------------------------------------
inline void GS_CH_update(int i, int j, int k, systm& u, systm& f, double dt, double h){
  double A11=1;
  double A12=dt/sq(h)*6;
  double A21=-d2fdphi2_c(u.phi(i,j,k))-6*kappa/sq(h);
  double A22=1;
  
  double f1=f.phi(i,j,k)+dt*laplacian_RHS(u.mu,i,j,k,h);
  double f2=f.mu(i,j,k)+dfdphi_c(u.phi(i,j,k))-d2fdphi2_c(u.phi(i,j,k))*u.phi(i,j,k)
    -kappa*laplacian_RHS(u.phi,i,j,k,h);
  
  //solve Au=f;
  double det=A11*A22-A12*A21;
  u.phi(i,j,k)=(A22*f1-A12*f2)/det;
  u.mu(i,j,k)=(-A21*f1+A11*f2)/det;
}
//-----------------------------------------------------------------------------
void GS_RB_CH(systm& u, systm& f, double dt, double h){
  u.periodicBoundary();
  
  //Red
  #pragma omp parallel for collapse(2)
  for (int i=0; i<u.N1; ++i){
    for (int j=0; j<u.N2; ++j)
      for (int k=(i+j)%2; k<u.N3; k+=2)
        GS_CH_update(i,j,k,u,f,dt,h);
  }
    
  //Black
  #pragma omp parallel for collapse(2)
  for (int i=0; i<u.N1; ++i){
    for (int j=0; j<u.N2; ++j)
      for (int k=(i+j+1)%2; k<u.N3; k+=2)
        GS_CH_update(i,j,k,u,f,dt,h);
  }
}
//-----------------------------------------------------------------------------
inline double CH_phi_LHS(systm& u, double dt, double h, int i, int j, int k){
  return(u.phi(i,j,k)-dt*u.mu.laplacian(i,j,k,h));
}
inline double CH_mu_LHS(systm& u, double dt, double h, int i, int j, int k){
  return(u.mu(i,j,k)-dfdphi_c(u.phi(i,j,k))+kappa*u.phi.laplacian(i,j,k,h));
}
//-----------------------------------------------------------------------------
void dfct_CH(systm& d, systm& u, systm& f,double dt, double h){
  u.periodicBoundary();
  #pragma omp parallel for collapse(3)
  gridLoop3D(d){
    d.phi(i,j,k)=f.phi(i,j,k)-CH_phi_LHS(u,dt,h,i,j,k);
    d.mu(i,j,k)=f.mu(i,j,k)-CH_mu_LHS(u,dt,h,i,j,k);
  }
  d.periodicBoundary();
}
//-----------------------------------------------------------------------------
void d_plus_Nu_CH(systm& f, systm& d, systm& u, double dt, double h){
  u.periodicBoundary();
  #pragma omp parallel for collapse(3)
  gridLoop3D(f){
    f.phi(i,j,k)=d.phi(i,j,k)+CH_phi_LHS(u,dt,h,i,j,k);
    f.mu(i,j,k)=d.mu(i,j,k)+CH_mu_LHS(u,dt,h,i,j,k);
  }
}
//-----------------------------------------------------------------------------
void f_CH(systm& f, systm& u, double dt, double h){
  u.periodicBoundary();
  #pragma omp parallel for collapse(3)
  gridLoop3D(f){
    f.phi(i,j,k)=u.phi(i,j,k);
    f.mu(i,j,k)=dfdphi_e(u.phi(i,j,k));
  }
}
//-----------------------------------------------------------------------------
void L_CH(grid3D& L, grid3D& u, grid3D& f, int Nx, int Ny, double dt, double h){
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
      double fct=dt/sq(h);
      f(i,j,0)+=(dfdphi_c(u(i+1,j,0))+dfdphi_c(u(i-1,j,0))
                +dfdphi_c(u(i,j+1,0))+dfdphi_c(u(i,j-1,0))
                -4*dfdphi_c(u(i,j,0))
              -(d2fdphi2_c(u(i+1,j,0))*u(i+1,j,0)+d2fdphi2_c(u(i-1,j,0))*u(i-1,j,0)
               +d2fdphi2_c(u(i,j+1,0))*u(i,j+1,0)+d2fdphi2_c(u(i,j-1,0))*u(i,j-1,0)
               -4*d2fdphi2_c(u(i,j,0))*u(i,j,0)))*fct;
      L(row,row,0)+=4*d2fdphi2_c(u(i,j,0))*fct;
      L(row,j*Nx+(i+1)%Nx,0)+=-d2fdphi2_c(u(i+1,j,0))*fct;
      L(row,j*Nx+(i+Nx-1)%Nx,0)+=-d2fdphi2_c(u(i-1,j,0))*fct;
      L(row,((j+1)%Ny)*Nx+i,0)+=-d2fdphi2_c(u(i,j+1,0))*fct;
      L(row,((j+Ny-1)%Ny)*Nx+i,0)+=-d2fdphi2_c(u(i,j-1,0))*fct;
  }
}

