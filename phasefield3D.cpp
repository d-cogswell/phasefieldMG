#include "phasefield3D.h"

//-----------------------------------------------------------------------------
//Free energy functions
//-----------------------------------------------------------------------------

//Double-well function from Boettinger
//-----------------------------------------------------------------------------
inline double g(double phi){
  return(sq(phi*(1-phi)));
}
inline double g_prime(double phi){
  return(phi*(4*sq(phi)-6*phi+2));
}

//Interpolating function from Boettinger
//-----------------------------------------------------------------------------
inline double p(double phi){
  return(cube(phi)*(6*sq(phi)-15*phi+10));
}
inline double p_prime(double phi){
  return(30*g(phi));
}

//The free energy of a pure component
//-----------------------------------------------------------------------------
inline double f(component A, double phi, double T){
  double Tm=A.Tm;
  double W=A.W;
  double L=A.L;
  return(W*g(phi)+L*(Tm-T)/Tm*p(phi));
}

inline double dfdphi(component comp, double phi, double T){
  double Tm=comp.Tm;
  double W=comp.W;
  double L=comp.L;
  return(W*g_prime(phi)+L*(Tm-T)/Tm*p_prime(phi));
}

//The free energy of an alloy using the regular solution model
//-----------------------------------------------------------------------------
inline double dfdphi(component A, component B, double omega_l, double omega_s, double phi, double c, double T){
  return((1-c)*dfdphi(A,phi,T)+c*dfdphi(B,phi,T)+c*(1-c)*(omega_l-omega_s)*p_prime(phi));
}

inline double dfdc(component A, component B, double omega_l, double omega_s, double phi, double c, double T){
  double R=8.314/1000; //Gas constant, energy in kJ
  double p_phi=p(phi);
  return(f(B,phi,T)-f(A,phi,T)+R*T*(log(c)-log(1-c))+(1-2*c)*(omega_s*(1-p_phi)+omega_l*p_phi));
}

//Cahn-Hilliard solver
//A mobility dependence on phi can be defined in the header file
//-----------------------------------------------------------------------------
int cahn_hilliard3D(grid3D* phi, double h, int iterations, int outputEvery){

  //Global variables
  double dt=.001;
  double K=1.5;
  double M=2.0;
  char outFile[128];
  grid3D insideTerms(*phi);

  //If mobility is defined as a function of phi in the header file
  //allocate space for a mobility grid here
  #ifdef M
  grid3D mobility(Nx,Ny,Nz);
  double **M = mobility.getGrid();
  #endif

  //Begin iterating
  //---------------------------------------------------------------------------
  for (int n=0; n<iterations+1; n++){

    //Write output, if necessary
    if (!(n%outputEvery)){
      sprintf(outFile,"output/p%6.6i.phi",n);
      cout << "writing output: " << outFile << endl;
      phi->writeToFile(outFile);
    }

    //Set boundary conditions for phi
    phi->periodicBoundary();

    //Compute mobility, and the inner terms of the cahn-hilliard equation
    gridLoop3D(*phi){
      insideTerms(i,j,k)=DfDphi((*phi)(i,j,k))-K*laplacian((*phi));

      //If mobility is not constant, calculate it
      #ifdef M
      M(i,j,k)=M(grid1,i,j,h);
      #endif
    }

    /* It's necessary to impose boundary conditions on insideTerms, 
     * because we take it's laplacian in the next step, 
     * and this requires examining neighbors
     */
     insideTerms.periodicBoundary();

    //Solve the Cahn-Hilliard equation using Forward Euler time-stepping
    gridLoop3D(*phi){
      #ifdef M
      double dphidt = M(i,j,k)*laplacian(insideTerms)+gradDot(M,insideTerms,i,j,h);
      #else
      double dphidt = M*laplacian(insideTerms);
      #endif

      (*phi)(i,j,k)=(*phi)(i,j,k)+dphidt*dt;
    }
  }
  return(0);
}


//Allen-cahn solver
//-----------------------------------------------------------------------------
int allen_cahn3D(grid3D* phi, double h, int iterations, int outputEvery){

  //Global variables
  double dt=.001;
  double K=1.5;
  double M=2.0;
  char outFile[128];

  //Begin iterating
  //---------------------------------------------------------------------------
  for (int n=0; n<iterations+1; n++){

    //Write output, if necessary
    if (!(n%outputEvery)){
      sprintf(outFile,"output/p%6.6i.phi",n);
      cout << "writing output: " << outFile << endl;
      phi->writeToFile(outFile);
    }

    //Set boundary conditions for phi
    phi->periodicBoundary();

    //Solve the Allen-Cahn equation using forward euler timestepping
    gridLoop3D(*phi){
      double dphidt=-M*(DfDphi((*phi)(i,j,k))-K*laplacian((*phi)));
      (*phi)(i,j,k)=(*phi)(i,j,k)+dphidt*dt;
    }
  }
  return(0);
}


//This function solves cahn-allen for a pure material melting or solidifying
//-----------------------------------------------------------------------------
int liquid_solid3D(component A, double T, grid3D* phi, double h, int iterations, int outputEvery){

  //Global variables
  double dt=.001;
  double K=1.5;
  double M=2.0;
  char outFile[128];

  //Begin iterating
  //---------------------------------------------------------------------------
  for (int n=0; n<iterations+1; n++){

    //Write output, if necessary
    if (!(n%outputEvery)){
      sprintf(outFile,"output/p%6.6i.phi",n);
      cout << "writing output: " << outFile << endl;
      phi->writeToFile(outFile);
    }

    //Set boundary conditions for phi
    phi->periodicBoundary();

    //Solve the Allen-Cahn equation using forward euler timestepping
    gridLoop3D(*phi){
      double dphidt=-M*(dfdphi(A,(*phi)(i,j,k),T)-K*laplacian((*phi)));
      (*phi)(i,j,k)=(*phi)(i,j,k)+dphidt*dt;
    }
  }
  return(0);
}

//This function solves coupled cahn-hilliard and allen-cahn to simulate a 
//binary system
//-----------------------------------------------------------------------------
int binary_alloy3D(component A, component B, double T, grid3D* phi, grid3D* c, double h, int iterations, int outputEvery){

  //Global variables
  double dt=.001;
  double K=1.5;
  double M=2.0;
  char outFile[128];
  grid3D insideTerms(*phi);

  //Define interaction perameters, in kJ (use zero for an ideal solution)
  double omega_l=0;
  double omega_s=0;

  //Begin iterating
  //---------------------------------------------------------------------------
  for (int n=0; n<iterations+1; n++){

    //Write output, if necessary
    if (!(n%outputEvery)){
      sprintf(outFile,"output/p%6.6i.phi",n);
      cout << "writing output: " << outFile << endl;
      phi->writeToFile(outFile);

      sprintf(outFile,"output/p%6.6i.c",n);
      c->writeToFile(outFile);
    }

    //Set boundary conditions for phi and c
    phi->periodicBoundary();
    c->periodicBoundary();

    //Solve the Allen-Cahn equation for phi using forward euler timestepping
    gridLoop3D(*phi){
      double dphidt=-M*(dfdphi(A,B,omega_l,omega_s,(*phi)(i,j,k),(*c)(i,j,k),T)-K*laplacian((*phi)));
      (*phi)(i,j,k)=(*phi)(i,j,k)+dphidt*dt;
    }

    //Solve the Cahn-Hilliard equation for composition
    //Compute mobility, and the inner terms of the cahn-hilliard equation
    gridLoop3D(*c){
      insideTerms(i,j,k)=dfdc(A,B,omega_l,omega_s,(*phi)(i,j,k),(*c)(i,j,k),T)-K*laplacian((*c));
    }

    /* It's necessary to impose boundary conditions on insideTerms, 
     * because we take it's laplacian in the next step, 
     * and this requires examining neighbors
     */
     insideTerms.periodicBoundary();

    //Solve the Cahn-Hilliard equation using Forward Euler time-stepping
    gridLoop3D(*c){
      double dcdt = M*laplacian(insideTerms);
      (*c)(i,j,k)=(*c)(i,j,k)+dcdt*dt;
    }
  }
  return(0);
}
