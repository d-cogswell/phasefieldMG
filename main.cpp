#include <omp.h>
#include <getopt.h>
#include <Magick++.h>
#include "multigrid3D.h"
#include "systm.h"
using namespace Magick;

//-----------------------------------------------------------------------------
int main(int argc, char **argv){
  int Nx=1,Ny=128,Nz=128;
  double h=1;
  double dt=.25;
  int iterations=100;
  int outputEvery=10;

  //Find the max dimension and use it to calculate the number of grid levels
  int max = Nx>Ny ? (Nx>Nz ? Nx : Nz) : Ny>Nz ? Ny : Nz;
  int grids=(int)(log(max)/log(2));

  printf("Max threads: %i\n",omp_get_max_threads());
  
  char* filename;
  grid3D* initial_condition;
  int option_char;
  Image img(Geometry(Ny,Nz),"black");
  const char* outDir="output";

  // Handle command line options
  bool inputFileSupplied=false;
  while ((option_char = getopt(argc, argv, "i:o:")) != -1)
    switch (option_char){
      case 'i':
        inputFileSupplied=true;
        filename=optarg;
        initial_condition = new grid3D(filename);
        break;
      case 'o':
        outDir=optarg;
        break;
    }

  if (!inputFileSupplied){
    printf("No input file supplied!\n");
    initial_condition = new grid3D(Nx,Ny,Nz);
    initial_condition->initializeGaussian(.01);
  }

  systm u(Nx,Ny,Nz), f(Nx,Ny,Nz), d(Nx,Ny,Nz), w(Nx,Ny,Nz);
  gridLoop3D(u){
    u.phi(i,j,k)=.5+.5*(*initial_condition)(i,j,k);
  }
  
  //Performs semi-implicit timestepping
  for (int t=0; t<=iterations; ++t){

    //Write output, if necessary
    char outFile[128];
    if (!(t%outputEvery)){
      sprintf(outFile,"%s/p%6.6i.jpg",outDir,t);
      printf("writing output: %s\n", outFile);
      gridLoop3D(u){
        double val=u.phi(0,j,k);
        val=QuantumRange*clip(val,0,1);
        img.pixelColor(j,k,Color(val,val,val,0));
      }
      img.write(outFile);
    }

    if (t<iterations){
      
      //Create f
      f_CH(f,u,dt,h);
      
      //Compute initial error
      dfct_CH(d,u,f,dt,h);
      double error=d.l2_norm();
      
      //multigrid
      while (error>1.e-4){
        FAS_multigrid<systm>(u,f,d,w,dt,h,1,grids);
        dfct_CH(d,u,f,dt,h);
        error=d.l2_norm();

        if (!isfinite(error)){
          printf("error is inf or nan!\n");
          exit(1);
        }
      }
    }
  }

  delete initial_condition;
  return 0;
}
