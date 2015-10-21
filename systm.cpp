#include "systm.h"

systm::systm(int Nx, int Ny, int Nz) : 
N1(Nx),N2(Ny),N3(Nz),phi(Nx,Ny,Nz),mu(Nx,Ny,Nz){
  fine=NULL;
  coarse=NULL;
}

systm::systm(grid3D* PHI, grid3D* MU) :
N1(PHI->N1),N2(PHI->N2),N3(PHI->N3),phi(*PHI),mu(*MU){
  fine=NULL;
  coarse=NULL;
}

systm::~systm(){  
  fine=NULL;
  coarse=NULL;
}

double systm::l2_norm(){
  double e1,e2;
  e1=phi.l2_norm();
  e2=mu.l2_norm();
  return(sqrt(e1*e1+e2*e2));
}

systm* systm::prolongate(int Nx, int Ny, int Nz){
  grid3D *phi2h, *mu2h;
  phi2h=phi.prolongate(Nx,Ny,Nz);
  mu2h=mu.prolongate(Nx,Ny,Nz);

  if (!fine){
    fine = new systm(phi2h,mu2h);
    fine->coarse=this;
  }
  
  return(fine);
}

systm* systm::prolongate_CC(int Nx, int Ny, int Nz){
  grid3D *phi2h, *mu2h;
  phi2h=phi.prolongate_CC(Nx,Ny,Nz);
  mu2h=mu.prolongate_CC(Nx,Ny,Nz);

  if (!fine){
    fine = new systm(phi2h,mu2h);
    fine->coarse=this;
  }
  
  return(fine);
}

systm* systm::restrict_FW(){
  grid3D *phi2h, *mu2h;
  phi2h=phi.restrict_FW();
  mu2h=mu.restrict_FW();
  
  if (!coarse){
    coarse = new systm(phi2h,mu2h);
    coarse->fine=this;
  }

  return(coarse);
}

systm* systm::restrict_CC(){
  grid3D *phi2h, *mu2h;
  phi2h=phi.restrict_CC();
  mu2h=mu.restrict_CC();
  
  if (!coarse){
    coarse = new systm(phi2h,mu2h);
    coarse->fine=this;
  }

  return(coarse);
}

systm* systm::injection(){
  grid3D *phi2h,*mu2h;
  phi2h=phi.injection();
  mu2h=mu.injection();
  
  if (!coarse){
    coarse = new systm(phi2h,mu2h);
    coarse->fine=this;
  }

  return(coarse);
}

void systm::periodicBoundary(){
  phi.periodicBoundary();
  mu.periodicBoundary();
}

systm* systm::getCoarseGrid(){
  if (!coarse){
    grid3D* phi2h=phi.getCoarseGrid();
    grid3D* mu2h=mu.getCoarseGrid();
    coarse=new systm(phi2h,mu2h);
    coarse->fine=this;
  }
  return(coarse);
}

void systm::operator=(const double b){
  phi=b;
  mu=b;
}

void systm::operator=(systm& b){
  phi=b.phi;
  mu=b.mu;
}

void systm::operator+=(systm& b){
  phi+=b.phi;
  mu+=b.mu;
}

void systm::operator-=(systm& b){ 
  phi-=b.phi;
  mu-=b.mu;
}