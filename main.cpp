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
    initial_condition = new grid3D(33,33,1);
  }

  //Cahn-hilliard and allen-cahn solvers
  initial_condition->initializeRandom(0,1);
  initial_condition->periodicBoundary();
  int Nx=33, Ny=33;
  h=1./26;
  double dt=.001;
/*
  //Iterative solver
  grid3D* u3 = new grid3D(Nx,Ny,1);
  gridLoop3D(*u3)
    (*u3)(i,j,k)=(*initial_condition)(i,j,k);
  for (int n=0; n<dt/.0001; ++n){
    u3->periodicBoundary();
    gridLoop3D(*u3){
      (*u3)(i,j,k)+=.0001/sq(h)*((*u3)(i+1,j,k)+(*u3)(i-1,j,k)+(*u3)(i,j+1,k)+(*u3)(i,j-1,k)-4*(*u3)(i,j,k));
    }
  }
  u3->writeToFile("output/out.iter.phi");

  //Direct solver
  grid3D* u2 = new grid3D(Nx,Ny,1);
  gridLoop3D(*u2)
    (*u2)(i,j,k)=(*initial_condition)(i,j,k);
  grid3D L(Nx*Ny,Nx*Ny,1);
  L_heat_eqn(&L,Nx,Ny,h,dt);
  gaussian_elimination(&L,u2,f);
  u2->writeToFile("output/out.direct.phi");
*/
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
    for (int n=0; n<10; ++n){
      multigrid(u,f,h,4);
    }
  }
/*
  //cahn_hilliard3D(initial_condition,h,iterations,outputEvery);
  //allen_cahn3D(initial_phi,h,iterations,outputEvery);
*/
  return 0;
}
