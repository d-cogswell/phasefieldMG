#include <iostream>
#include <unistd.h>
#include "grid3D.h"
#include "phasefield3D.h"
using namespace std;

//-----------------------------------------------------------------------------
int main(int argc, char **argv){
  double h=1;
  int iterations=10;
  int outputEvery=1;

  char* filename;
  grid3D* initial_condition;
  int option_char;

  // Handle command line options
  bool inputFileSupplied=false;
  while ((option_char = getopt(argc, argv, "i:")) != -1)
    switch (option_char){
    case 'i':
      inputFileSupplied=true;
      filename=optarg;
      initial_condition = new grid3D(filename);
      break;
    }

  if (!inputFileSupplied){
    cout << "No input file supplied!" << endl;
    initial_condition = new grid3D(25,25,1);
  }

  //Cahn-hilliard and allen-cahn solvers
  initial_condition->initializeSphere(5);
  h=1./26;
  int Nx=25, Ny=25;
  double dt=.01;

  //create f
  grid3D* f = new grid3D(Nx,Ny,1);
  gridLoop3D(*f)
      (*f)(i,j,k)=(*initial_condition)(i,j,k);

  //Direct solver
  grid3D* u3 = new grid3D(Nx,Ny,1);
  gridLoop3D(*u3)
    (*u3)(i,j,k)=(*initial_condition)(i,j,k);
  grid3D L(Nx*Ny,Nx*Ny,1);
  L_heat_eqn(&L,Nx,Ny,h,dt);
  gaussian_elimination(&L,u3,f);
  u3->writeToFile("output/out.direct.phi");

  //Multigrid solver
  grid3D* u = new grid3D(Nx,Ny,1);
  gridLoop3D(*u)
    (*u)(i,j,k)=(*initial_condition)(i,j,k);
  for (int n=0; n<20; ++n)
    multigrid(u,f,h,2);
  u->writeToFile("output/out.multi.phi");

/*
  //cahn_hilliard3D(initial_condition,h,iterations,outputEvery);
  //allen_cahn3D(initial_phi,h,iterations,outputEvery);
*/
  return 0;
}
