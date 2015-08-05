#include "multigrid3D.h"

void gaussian_elimination(grid3D& L, grid3D& u, grid3D& f){
  int N1=L.N1;
  int N2=L.N2;
  int Nx=u.N1;

  //Convert f to a column vector
  grid3D f_col(N2,1,1);
  gridLoop3D(f){
    f_col(j*Nx+i,0,0)=f(i,j,0);
  }

  //Forward elimination
  for (int col=0; col<N2; ++col){
    for (int row=col+1; row<N1; ++row){
      double factor=-L(row,col,0)/L(col,col,0);
      if (factor!=0){
        for (int j=col; j<N2; ++j){
          L(row,j,0)+=factor*L(col,j,0);
        }
        f_col(row,0,0)+=factor*f_col(col,0,0);
      }
    }
  }

  //Backward elimination
  for (int col=N2-1; col>=0; --col){
    f_col(col,0,0)/=L(col,col,0);
    for (int row=col-1; row>=0; --row){
      f_col(row,0,0)-=L(row,col,0)*f_col(col,0,0);
    }
  }

  //copy f to u
  gridLoop3D(u){
    u(i,j,0)=f_col(j*Nx+i,0,0);
  }
}
