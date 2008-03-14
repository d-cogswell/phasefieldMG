#include "phasefield3D.h"

//Cahn-Hilliard solver
//A mobility dependence on phi can be defined in the header file
//-----------------------------------------------------------------------------
int cahn_hilliard3D(grid3D* phi, double h, int iterations, int outputEvery){

  //Global variables
  double dt=.001;
  double K=1.5;
  char outFile[128];

  int Nx=phi->getDimension(1);
  int Ny=phi->getDimension(2);
  int Nz=phi->getDimension(3);
  grid3D insideTerms(Nx,Ny,Nz);
  Image img(Geometry(Nx,Ny),"black");

  //Begin iterating
  //---------------------------------------------------------------------------
  for (int n=0; n<iterations+1; n++){

    //Write output, if necessary
    if (!(n%outputEvery)){
      sprintf(outFile,"output/p%6.6i.png",n);
      cout << "writing output: " << outFile << endl;
      gridLoop3D(*phi){
        double val=(1+(*phi)(i,j,0))/2;
        val=MaxRGB*clip(val,0,1);
        img.pixelColor(i,j,Color(val,val,val,0));
      }
      img.write(outFile);
    }

    //Set boundary conditions for phi
    phi->periodicBoundary();

    //Compute mobility, and the inner terms of the cahn-hilliard equation
    gridLoop3D(*phi){
      insideTerms(i,j,k)=DfDphi((*phi)(i,j,k))-K*laplacian((*phi));
    }

    /* It's necessary to impose boundary conditions on insideTerms, 
     * because we take it's laplacian in the next step, 
     * and this requires examining neighbors
     */
     insideTerms.periodicBoundary();

    //Solve the Cahn-Hilliard equation using Forward Euler time-stepping
    gridLoop3D(*phi){
      double dphidt = laplacian(insideTerms);
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
  char outFile[128];
  Image img(Geometry(phi->getDimension(1),phi->getDimension(2)),"black");

  //Begin iterating
  //---------------------------------------------------------------------------
  for (int n=0; n<iterations+1; n++){

    //Write output, if necessary
    if (!(n%outputEvery)){
      sprintf(outFile,"output/p%6.6i.png",n);
      cout << "writing output: " << outFile << endl;
      gridLoop3D(*phi){
        double val=(1+(*phi)(i,j,0))/2;
        val=MaxRGB*clip(val,0,1);
        img.pixelColor(i,j,Color(val,val,val,0));
      }
      img.write(outFile);
    }

    //Set boundary conditions for phi
    phi->periodicBoundary();

    //Solve the Allen-Cahn equation using forward euler timestepping
    gridLoop3D(*phi){
      double dphidt=-(DfDphi((*phi)(i,j,k))-K*laplacian((*phi)));
      (*phi)(i,j,k)=(*phi)(i,j,k)+dphidt*dt;
    }
  }
  return(0);
}
//-----------------------------------------------------------------------------
//This is a recursive function that relaxes the error using 'max_level' grid
//levels.   max_levels=2 corresponds to the two grid scheme.
double multigrid(grid3D** L, grid3D* u, grid3D* f,grid3D* d, grid3D* e, double dt, double h, int max_level, int level){

    //Get the grid dimensions for the current level
    int Nx=u->getDimension(1);
    int Ny=u->getDimension(2);

    //Set number of iterations on the fine grid and coarse grid
    int v1=1;
    int v2=1;

    //Direct solve on the coarsest mesh
    if (level==max_level){
      if (*L==NULL){
        *L = new grid3D(Nx*Ny,Nx*Ny,1);
        L_CH(*L,Nx,Ny,dt,h);
      }
      gaussian_elimination(*L,e,d);
    }

    //Otherwise perform a coarse grid correction
    else{

      //Presmoothing
      u->periodicBoundary();
      for (int i=0;i<v1;++i){
        GS_LEX_CH(u,f,dt,h);
      }

      //Compute the defect
      dfct_CH(d,u,f,dt,h);

      //solve the defect equation on a coarse mesh
      multigrid(L,e->restrict(),d->restrict(),u->getCoarseGrid(),f->getCoarseGrid(),dt,2*h,max_level,level+1);

      //Prolongate the error to the fine mesh
      e->getCoarseGrid()->prolongate(Nx,Ny);

      //Compute the corrected approximation
      gridLoop3D(*u)
        (*u)(i,j,k)+=(*e)(i,j,k);

      //Postsmoothing
      u->periodicBoundary();
      for (int i=0;i<v2;++i)
        GS_LEX_CH(u,f,dt,h);
    }

    if (level==1)
      return(d->l2_norm());
}
//-----------------------------------------------------------------------------
void gaussian_elimination(grid3D* L, grid3D* u, grid3D* f){
  int N1=L->getDimension(1);
  int N2=L->getDimension(2);
  int Nx=u->getDimension(1);

  //Convert f to a column vector
  grid3D *f_col = new grid3D(N2,1,1);
  gridLoop3D(*f){
    (*f_col)(j*Nx+i,0,0)=(*f)(i,j,0);
  }

  //Forward elimination
  for (int col=0; col<N2; ++col){
    for (int row=col+1; row<N1; ++row){
      double factor=-(*L)(row,col,0)/(*L)(col,col,0);
      if (factor!=0){
        for (int j=col; j<N2; ++j){
          (*L)(row,j,0)+=factor*(*L)(col,j,0);
        }
        (*f_col)(row,0,0)+=factor*(*f_col)(col,0,0);
      }
    }
  }

  //Backward elimination
  for (int col=N2-1; col>=0; --col){
    (*f_col)(col,0,0)/=(*L)(col,col,0);
    for (int row=col-1; row>=0; --row){
      (*f_col)(row,0,0)-=(*L)(row,col,0)*(*f_col)(col,0,0);
    }
  }

  //copy f to u
  gridLoop3D(*u){
    (*u)(i,j,0)=(*f_col)(j*Nx+i,0,0);
  }
  delete f_col;
}
