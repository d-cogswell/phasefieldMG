#include <omp.h>
#include <getopt.h>
#include <Magick++.h>
#include "multigrid3D.h"
#include "systm.h"
using namespace Magick;

//-----------------------------------------------------------------------------
int main(int argc, char **argv){
  int Nx=129,Ny=129,Nz=1;
  int grids=(int)(log(Nx>Ny ? Nx : Ny)/log(2))-1;
  double h=1;
  double dt=1;
  int iterations=100;
  int outputEvery=10;

  printf("Max threads: %i\n",omp_get_max_threads());
  
  char* filename;
  grid3D* initial_condition;
  grid3D* L = NULL;
  int option_char;
  Image img(Geometry(Nx,Ny),"black");
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

  systm u(Nx,Ny,Nz), f(Nx,Ny,Nz), d(Nx,Ny,Nz), v(Nx,Ny,Nz), w(Nx,Ny,Nz);
  u.phi=*initial_condition;
  
  //Performs semi-implicit timestepping
  for (int t=0; t<=iterations; ++t){

    //Write output, if necessary
    char outFile[128];
    if (!(t%outputEvery)){
      sprintf(outFile,"%s/p%6.6i.jpg",outDir,t);
      printf("writing output: %s\n", outFile);
      gridLoop3D(u){
        double val=(1+u.phi(i,j,0))/2;
        val=MaxRGB*clip(val,0,1);
        img.pixelColor(i,j,Color(val,val,val,0));
      }
      img.write(outFile);
    }

    //Create f
    f_CH(f,u,dt,h);

    //multigrid
    double error=1;
    if (t<iterations){
      while (error>1.e-4){
        FAS_multigrid<systm>(L,u,f,d,v,w,dt,h,1,grids);
        dfct_CH(d,u,f,dt,h);
        error=d.l2_norm();

        if (!isfinite(error)){
          printf("error is inf or nan!\n");
          exit(1);
        }
      }
    }
  }

  delete initial_condition,L;
  return 0;
}
