#include "multigrid3D.h"

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
    printf("No input file supplied!\n");
    initial_condition = new grid3D(Nx,Ny,1);
  }

  initial_condition->initializeGaussian(.01);
  grid3D f(Nx,Ny,1), d(Nx,Ny,1);
  grid3D* u = initial_condition;
  grid3D* L = NULL;
  
  //Performs semi-implicit timestepping
  for (int t=0; t<=iterations; ++t){

    //Write output, if necessary
    char outFile[128];
    if (!(t%outputEvery)){
      sprintf(outFile,"output/p%6.6i.jpg",t);
      printf("writing output: %s\n", outFile);
      gridLoop3D(*u){
        double val=(1+(*u)(i,j,0))/2;
        val=MaxRGB*clip(val,0,1);
        img.pixelColor(i,j,Color(val,val,val,0));
      }
      img.write(outFile);
    }

    //Create f
    u->periodicBoundary();
    f_CH(f,*u,dt,h);
    
    double error=1;
    double error_old=NULL;

    //multigrid
    if (t<iterations){
      while (error>1.e-3){
        FAS_multigrid<grid3D>(&L,*u,f,dt,h,6);
        dfct_CH(d,*u,f,dt,h);
        error=d.l2_norm();
	  
        if (error_old && error>=error_old){
          printf("error not decreasing!\n");
          exit(1);
        }
        error_old=error;
      }
    }
  }

  delete initial_condition,L;
  return 0;
}
