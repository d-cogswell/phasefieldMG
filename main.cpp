#include <iostream>
#include <unistd.h>
#include "grid3D.h"
#include "phasefield3D.h"
using namespace std;

//-----------------------------------------------------------------------------
int main(int argc, char **argv){
  int Nx=33,Ny=33;
  double h=1;
  double dt=.1;
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
    initial_condition = new grid3D(Nx,Ny,1);
  }

  initial_condition->initializeRandom(0,1);
  grid3D* f = new grid3D(Nx,Ny,1);
  grid3D* u = initial_condition;

  //Performs semi-implicit timestepping
  for (int t=0; t<=iterations; ++t){

    //Write output, if necessary
    char outFile[128];
    if (!(t%outputEvery)){
      sprintf(outFile,"output/p%6.6i.phi",t);
      cout << "writing output: " << outFile << endl;
      u->writeToFile(outFile);
    }

    //Create f
    u->periodicBoundary();
    f_heat_eqn(f,u,dt,h);

/*
    //Direct solve
    grid3D L(Nx*Ny,Nx*Ny,1);
    L_heat_eqn(&L,Nx,Ny,dt,h);
    gaussian_elimination(&L,u,f);
*/

    //multigrid
    while(multigrid(u,f,dt,h,5)>1.e-7);
  }
  delete initial_condition,f;
/*
  //cahn_hilliard3D(initial_condition,h,iterations,outputEvery);
  //allen_cahn3D(initial_phi,h,iterations,outputEvery);
*/
  return 0;
}
