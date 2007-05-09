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
//-----------------------------------------------------------------------------
int multigrid(grid3D* u, double h, int iterations, int outputEvery){
  int N=25;
  double dt=.5*h*h/2;

  grid3D f(N,N,1,0);
  grid3D u_smth(N,N,1,0);
  grid3D u_CGC(N,N,1,0);
  grid3D *v;

  //Set number of iterations on the fine grid and coarse grid
  int v1=2;
  int v2=0;

  grid3D r(N,N,1,0);
  grid3D e(N,N,1,0);

  for (int n=0; n<iterations; n++){

    //Write output, if necessary
    char outFile[128];
    if (!(n%outputEvery)){
      sprintf(outFile,"output/p%6.6i.phi",n);
      cout << "writing output: " << outFile << endl;
      u->writeToFile(outFile);
    }

    //Compute f
    for (int i=1; i<u->getDimension(1)-1; ++i)
      for (int j=1; j<u->getDimension(2)-1; ++j)
        f(i,j,0)=(*u)(i,j,0)+dt/sq(h)*((*u)(i+1,j,0)+(*u)(i-1,j,0)+(*u)(i,j+1,0)+(*u)(i,j-1,0)-4*(*u)(i,j,0));

    //Smooth v1 on the fine mesh
    for (int n1=0; n1<v1;++n1){
      GS_LEX(&u_smth,u,&f,h);
    }

    //Calculate the residual
    for (int i=1; i<u->getDimension(1)-1; ++i)
      for (int j=1; j<u->getDimension(2)-1; ++j)
        r(i,j,0)=f(i,j,0)-u_smth(i,j,0);

    //Control grid refinement levels here (2=two grid)
    int levels=4;
    v=multigrid_CGC(&e,&r,h,levels);

    //Add the relaxed error back to u
    gridLoop3D(*u){
      //u_CGC(i,j,k)=u_smth(i,j,k)+(*v)(i,j,k);
      (*u)(i,j,k)=u_smth(i,j,k)+(*v)(i,j,k);
      //(*u)(i,j,k)=u_CGC(i,j,k);
    }

    //Smooth v2 on the fine mesh
    for (int n2=0; n2<v2;++n2){
      GS_LEX(u,u,&f,h);
    }
  }
  return(0);
}

//-----------------------------------------------------------------------------
//This is a recursive function that relaxes the error using 'max_level' grid
//levels.   max_levels=2 corresponds to the two grid scheme.
grid3D* multigrid_CGC(grid3D* e, grid3D* r, double h, int max_level, int level){
  if (level<=max_level){

    //Create the restricted grids
    int Nx=r->getDimension(1);
    int Ny=r->getDimension(2);
    int N2x=(Nx+1)/2;
    int N2y=(Ny+1)/2;
    grid3D *e2h = new grid3D(N2x,N2y,0);
    //grid3D* e2h=e->restrict();

    //For debugging
    //cout << "Level " << level << ":(" << Nx << "," << Ny << ")" << endl;

    //Restrict the residual to a coarse mesh
    grid3D* r2h=r->restrict();

    //Solve on the coarse mesh
    L_heat_eqn(e2h,r2h,2*h);

    //Recursively call multigrid for higher levels, and add up relaxed error
    grid3D* e4h=multigrid_CGC(e2h,r2h,2*h,max_level,level+1);
    if (e4h){
      gridLoop3D(*e2h){
	(*e2h)(i,j,k)+=(*e4h)(i,j,k);
      }
    }

    //Prolongate the error to the fine mesh and return it
    return(e2h->prolongate(Nx,Ny));
  }
  return(0);
}
//-----------------------------------------------------------------------------
inline void GS_LEX(grid3D* u_smooth, grid3D* u, grid3D* f, double h){
  for (int i=1; i<u->getDimension(1)-1; ++i)
    for (int j=1; j<u->getDimension(2)-1; ++j)
      (*u_smooth)(i,j,0)=.25*(sq(h)*(*f)(i,j,0)+(*u)(i+1,j,0)+(*u_smooth)(i-1,j,0)+(*u)(i,j+1,0)+(*u_smooth)(i,j-1,0));
}
//-----------------------------------------------------------------------------
inline void L_heat_eqn(grid3D* u, grid3D* f, double h){
  for (int i=1; i<u->getDimension(1)-1; ++i)
    for (int j=1; j<u->getDimension(2)-1; ++j)
      (*u)(i,j,0)=(*f)(i,j,0);
}
//-----------------------------------------------------------------------------
void gaussian_elimination(grid3D* L, grid3D* u, grid3D* f){
  int N1=L->getDimension(1);
  int N2=L->getDimension(2);

  //Forward elimination
  for (int col=0; col<N2; ++col){
    for (int row=col+1; row<N1; ++row){
      double factor=-(*L)(row,col,0)/(*L)(col,col,0);
      for (int j=col; j<N2; ++j){
        (*L)(row,j,0)+=factor*(*L)(col,j,0);
      }
      (*f)(row,0,0)+=factor*(*f)(col,0,0);
    }
  }

  //Backward elimination
  for (int col=N2-1; col>=0; --col){
    (*f)(col,0,0)/=(*L)(col,col,0);
    for (int row=col-1; row>=0; --row){
      (*f)(row,0,0)-=(*L)(row,col,0)*(*f)(col,0,0);
    }
  }

  for (int i=0; i<N1; ++i){
    (*u)(i,0,0)=(*f)(i,0,0);
  }
}
