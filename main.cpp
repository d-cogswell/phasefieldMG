#include <iostream>
#include <unistd.h>
#include "grid3D.h"
#include "phasefield3D.h"
using namespace std;

//-----------------------------------------------------------------------------
//This function creates the initial block configuration for problem #2
//B should be a array of 4 integers that represent the source blocks
grid3D* createGrid(grid3D* grid, int* B, int N){

  int blockResolution=4;
  int blockWidth=(N-1)/6+1;

  //Set gridpoints corresponding to each block
  (*grid)=0;
  for (int b=0; b<4; ++b){
    int x=1+(B[b]-1)/blockResolution;
    int y=B[b]-(x-1)*blockResolution;

    for (int i=0; i<blockWidth; ++i)
      for (int j=0; j<blockWidth; ++j)
        (*grid)((blockWidth-1)*x+i,(blockWidth-1)*y+j,0)=1;
  }
}

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
  int blocks [4] = {1,7,14,16};
  createGrid(initial_condition,&blocks[0],25);
  initial_condition->initializeRandom(0,1);
  h=1./26;
  int Nx=25;
  int Ny=25;
  double dt=.5*h*h/2;
  dt=.01;

  //Construct u
  grid3D u(Nx*Ny,1,1,0);
  gridLoop3D(*initial_condition){
    u(j*Nx+i,0,0)=(*initial_condition)(i,j,k);
  }

  //Construct L
  grid3D L(Nx*Ny,Nx*Ny,0);
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j){
      int row=j*Nx+i;
      L(row,row,0)+=4+sq(h)/dt;
      L(row,j*Nx+(i+1)%Nx,0)+=-1;
      L(row,j*Nx+(i+Nx-1)%Nx,0)+=-1;
      L(row,((j+1)%Ny)*Nx+i,0)+=-1;
      L(row,((j+Ny-1)%Ny)*Nx+i,0)+=-1;
  }

  //Construct f
  grid3D f(Nx*Ny,1,0);
  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j)
      f(j*Nx+i,0,0)=sq(h)/dt*(*initial_condition)(i,j,0);

  gaussian_elimination(&L,&u,&f);

  //Convert u back to 2D grid
  gridLoop3D(*initial_condition){
    (*initial_condition)(i,j,k)=u(j*Nx+i,0,0);
  }
  initial_condition->writeToFile("output/out.phi");

/*
  multigrid(initial_condition,h,iterations,outputEvery);
  //cahn_hilliard3D(initial_condition,h,iterations,outputEvery);
  //allen_cahn3D(initial_phi,h,iterations,outputEvery);
*/
  return 0;
}
