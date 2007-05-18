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
    initial_condition = new grid3D(129,129,1);
  }

  //Cahn-hilliard and allen-cahn solvers
  initial_condition->initializeRandom(0,1);
  initial_condition->periodicBoundary();
  int Nx=129, Ny=129;
  h=1./26;
  double dt=.001;

  //Multigrid solver
  grid3D* f = new grid3D(Nx,Ny,1);
  grid3D* u = initial_condition;

  //Performs semi-implicit timestepping
  for (int t=0; t<iterations; ++t){

    //Write output, if necessary
    char outFile[128];
    if (!(t%outputEvery)){
      sprintf(outFile,"output/p%6.6i.phi",t);
      cout << "writing output: " << outFile << endl;
      u->writeToFile(outFile);
    }

    //Create f
    gridLoop3D(*f){
      (*f)(i,j,k)=dt*((*initial_condition)(i+1,j,k)
                     +(*initial_condition)(i-1,j,k)
                     +(*initial_condition)(i,j+1,k)
                     +(*initial_condition)(i,j-1,k))
                     +(2*sq(h)-4*dt)*(*initial_condition)(i,j,k);
    }

    //Solve
    for (int n=0; n<15; ++n){
      multigrid(u,f,h,7);
    }
  }
  delete initial_condition,f;
/*
  //cahn_hilliard3D(initial_condition,h,iterations,outputEvery);
  //allen_cahn3D(initial_phi,h,iterations,outputEvery);
*/
  return 0;
}
