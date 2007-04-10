#include <iostream>
#include <unistd.h>
#include "grid3D.h"
#include "convexification.h"
#include "phasefield3Dmpi.h"
//#include "phasefield3D.h"
using namespace std;

//Calibrated for 10 boundary pts
double factor=.26;
double expnt = 1.0;

//-----------------------------------------------------------------------------
int main(int argc, char **argv){
  double lambda=0;
  double theta=0;
  //double h=.5;
  double h=.25;
  int iterations=200000;
  int outputEvery=1000;

  char* filename;
  grid3D* initialCond;
  int option_char;

  // Handle command line options
  bool inputFileSupplied=false;

  while ((option_char = getopt(argc, argv, "f:e:w:i:a:")) != -1)
    switch (option_char){
    case 'f':
      factor=atof(optarg);
      break;
    case 'e':
      expnt=atof(optarg);
      break;
    case 'i':
      inputFileSupplied=true;
      filename=optarg;
      //initialCond = new grid3D(filename,bc,1,0,80+50,0,0,40+50,0);
      initialCond = new grid3D(filename);
      break;
    case 'w':
      lambda=atof(optarg);
      break;
    case 'a':
      theta=atof(optarg);
      break;
    }

  if (!inputFileSupplied){
    cout << "No input file supplied!" << endl;
    initialCond = new grid3D(100,100,100);
    //initialCond->initializeSphere(25);
    //initialCond = new grid3D(80,80,80);
    initialCond->initializeSphere(20);
  }

    //phasefield3D(initialCond,h,iterations,outputEvery);
    phasefield3DMPI(argc,argv,initialCond,h,iterations,outputEvery);

  //Solve Maxwell's equations
  //maxwell3D(argc,argv,initialCond,lambda,theta,h,iterations,outputEvery);
  //maxwell3D(initialCond,lambda,theta,h,iterations,outputEvery);
  return 0;
}
