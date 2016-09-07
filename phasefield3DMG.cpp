#include "phasefield3DMG.h"

double gama=1.1547;
double delta=3.4641;
double kappa=3*delta*gama/2;
double W=12*gama/delta;

inline double dfdphi_c(double x){
  return(W*2*x);
}

inline double d2fdphi2_c(double x){
  return(W*2);
}

inline double dfdphi_e(double x){
  return(W*(-6*sq(x)+4*cube(x)));
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
  
  u.periodicBoundary();
}
//-----------------------------------------------------------------------------
inline double AC_LHS(grid3D& u, double dt, double h, int i, int j, int k){
  return(u(i,j,k)+dt*(dfdphi_c(u(i,j,k))-kappa*u.laplacian(i,j,k,h)));
}
//-----------------------------------------------------------------------------
void dfct_AC(grid3D& d, grid3D& u, grid3D& f, double dt, double h){
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


/*The following functions solve the Cahn-Hilliard equation using Eyre's 
 *nonlinearly stabilized splitting from "An Unconditionally Stable One-Step 
 * Scheme for Gradient Systems".
 */
//-----------------------------------------------------------------------------
#pragma omp declare simd uniform(u,f,dt,h,i,j) //linear(i:1,j:1,k:1)
void GS_CH_update(int i, int j, int k, systm& u, systm& f, double dt, double h){
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
  
  u.periodicBoundary();
}
//-----------------------------------------------------------------------------
#pragma omp declare simd uniform(u,dt,h,i,j) //linear(i:1,j:1,k:1)
inline double CH_phi_LHS(systm& u, double dt, double h, int i, int j, int k){
  return(u.phi(i,j,k)-dt*u.mu.laplacian(i,j,k,h));
}

#pragma omp declare simd uniform(u,dt,h,i,j) //linear(i:1,j:1,k:1)
inline double CH_mu_LHS(systm& u, double dt, double h, int i, int j, int k){
  return(u.mu(i,j,k)-dfdphi_c(u.phi(i,j,k))+kappa*u.phi.laplacian(i,j,k,h));
}
//-----------------------------------------------------------------------------
void dfct_CH(systm& d, systm& u, systm& f,double dt, double h){
  #pragma omp parallel for collapse(2)
  for (int i=0;i<d.N1;++i)
    for (int j=0;j<d.N2;++j)
      #pragma omp simd
      for (int k=0;k<d.N3;++k){
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
