#include "phasefield3DMG.h"

//-----------------------------------------------------------------------------
//This is a recursive function that relaxes the error using 'max_level' grid
//levels.   max_levels=2 corresponds to the two grid scheme.
double multigrid(grid3D** L, grid3D* u, grid3D* f,grid3D* d, grid3D* e, double dt, double h, int max_level, int level){

  //Set number of pre and post smoothing iterations
  int v1=1;
  int v2=1;

  //Presmoothing
  u->periodicBoundary();
  for (int i=0;i<v1;++i){
    GS_LEX_CH(u,f,dt,h);
  }

  //Compute the defect and restrict it
  dfct_CH(d,u,f,dt,h);
  grid3D* d2h=d->restrict();
  grid3D* e2h=e->getCoarseGrid();

  //Direct solve on the coarsest mesh
  if (level==max_level-1){

     //Get the grid dimensions for the current level
    int Nx=e2h->getDimension(1);
    int Ny=e2h->getDimension(2);

    //Only allocate space for L once
    if (*L==NULL)
      *L = new grid3D(Nx*Ny,Nx*Ny,1);

    L_CH(*L,Nx,Ny,dt,2*h);
    gaussian_elimination(*L,e2h,d2h);
  }

  //Otherwise perform a coarse grid correction
  else
    multigrid(L,u->getCoarseGrid(),f->getCoarseGrid(),d2h,e2h,dt,2*h,max_level,level+1);

  //Prolongate the error to the fine mesh
  e2h->prolongate(e->getDimension(1),e->getDimension(2));

  //Compute the corrected approximation
  gridLoop3D(*u)
    (*u)(i,j,k)+=(*e)(i,j,k);

  //Postsmoothing
  u->periodicBoundary();
  for (int i=0;i<v2;++i)
    GS_LEX_CH(u,f,dt,h);

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
