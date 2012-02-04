#include <iostream>
#include <unistd.h>
#include "grid3D.h"
#include "phasefield3DMG.h"
using namespace std;

//-----------------------------------------------------------------------------
int main(int argc, char **argv){
  int Nx=129,Ny=129;
  double h=1;
  double dt=1;
  int iterations=100;
  int outputEvery=10;

  char* filename;
  grid3D* initial_condition;
  int option_char;
 Image img(Geometry(Nx,Ny),"black");

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

  initial_condition->initializeGaussian(.01);
  grid3D* f = new grid3D(Nx,Ny,1);
  grid3D* d = new grid3D(Nx,Ny,1);
  grid3D* e = new grid3D(Nx,Ny,1);
  grid3D* u = initial_condition;
  grid3D* L = NULL;

  //Performs semi-implicit timestepping
  for (int t=0; t<=iterations; ++t){

    //Write output, if necessary
    char outFile[128];
    if (!(t%outputEvery)){
      sprintf(outFile,"output/p%6.6i.png",t);
      cout << "writing output: " << outFile << endl;
      gridLoop3D(*u){
        double val=(1+(*u)(i,j,0))/2;
        val=MaxRGB*clip(val,0,1);
        img.pixelColor(i,j,Color(val,val,val,0));
      }
      img.write(outFile);
    }

    //Create f
    u->periodicBoundary();
    f_CH(f,u,dt,h);

    //Direct solve
    //grid3D L(Nx*Ny,Nx*Ny,1);
    //L_CH(&L,Nx,Ny,dt,h);
    //gaussian_elimination(&L,u,f);

    //multigrid
    while(multigrid(&L,u,f,d,e,dt,h,7)>1.e-5);
  }

  delete initial_condition,L,f;
  return 0;
}
