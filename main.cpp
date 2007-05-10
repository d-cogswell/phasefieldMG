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
  initial_condition->initializeRandom(0,1);
  grid3D* u=initial_condition;
  h=1./26;
  double dt=.001;
  initial_condition->writeToFile("output/init.phi");

  //Multigrid solver
  multigrid(initial_condition,h,2);
  initial_condition->writeToFile("output/out.multi.phi");
/*
  //Construct L
  int Nx=25;
  int Ny=25;
  grid3D L(Nx*Ny,Nx*Ny,1);
  L_heat_eqn(&L,Nx,Ny,h,dt);

  //Construct f
  grid3D f(Nx,Ny,1);
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j)
      f(i,j,0)=sq(h)/dt*(*u)(i,j,0);

  gaussian_elimination(&L,u,&f);
  initial_condition->writeToFile("output/out.dir.phi");
*/
/*
  //cahn_hilliard3D(initial_condition,h,iterations,outputEvery);
  //allen_cahn3D(initial_phi,h,iterations,outputEvery);
*/
  return 0;
}
