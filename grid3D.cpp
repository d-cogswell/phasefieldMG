#include "grid3D.h"

//Grid macros
//-----------------------------------------------------------------------------
#define gridLoop for(int i=0; i<N1; ++i) for (int j=0; j<N2; ++j) for (int k=0; k<N3; ++k)
#define bndryGridLoop for (int i=-boundary; i<N1+boundary; ++i) for (int j=-boundary; j<N2+boundary; ++j) for (int k=-boundary; k<N3+boundary; ++k)
#define fastBndryLoop for (int n=0;n<(N1+2*boundary)*(N2+2*boundary)*(N3+2*boundary);++n)

//-----------------------------------------------------------------------------
//Linear interpolation
double Lint(double x, double u0, double u1){
    return(u0+(u1-u0)*x);
}

//Cubic interpolation
double Cint(double x, double u_1, double u0, double u1, double u2){
    double a0=u0;
    double a1=-u_1/3-u0/2+u1-u2/6;
    double a2=u_1/2-u0+u1/2;
    double a3=-u_1/6+u0/2-u1/2+u2/6;
    return(a0+a1*x+a2*x*x+a3*x*x*x);
}

//Forward looking cubic interpolation
double Cint_fwd(double x, double u0, double u1, double u2, double u3){
    double a0=u0;
    double a1=-11./6*u0+3*u1-3./2*u2+1./3*u3;
    double a2=u0-5./2*u1+2*u2-1./2*u3;
    double a3=-1./6*u0+1./2*u1-1./2*u2+1./6*u3;
    return(a0+a1*x+a2*x*x+a3*x*x*x);
}

double grid3D::cubic_x(double i, double j, double k){
  int i0=(int)i;
  int j0=(int)j;
  int k0=(int)k;
  
  //If there aren't enough gridpoints, perform linear interpolation
  if (N1<4)
      return(Lint(i-i0,cubic_y(i0,j,k),cubic_y(i0+1,j,k)));
  
  if (i0==0)
    return(Cint_fwd(i-i0,cubic_y(i,j,k),cubic_y(i+1,j,k),cubic_y(i+2,j,k),cubic_y(i+3,j,k)));    
  else if (i0==N1-1)
    return(Cint_fwd(1-(i-i0),cubic_y(i,j,k),cubic_y(i-1,j,k),cubic_y(i-2,j,k),cubic_y(i-3,j,k)));
  else
    return(Cint(i-i0,cubic_y(i-1,j,k),cubic_y(i,j,k),cubic_y(i+1,j,k),cubic_y(i+2,j,k)));   
}


double grid3D::cubic_y(double i, double j, double k){
  int i0=(int)i;
  int j0=(int)j;
  int k0=(int)k;
  
  if (N2<4)
    return(Lint(j-j0,cubic_z(i,j,k),cubic_z(i,j+1,k)));
  
  if (j0==0)
    return(Cint_fwd(j-j0,cubic_z(i,j,k),cubic_z(i,j+1,k),cubic_z(i,j+2,k),cubic_z(i,j+3,k)));
  else if (j0==N2-1)
    return(Cint_fwd(1-(j-j0),cubic_z(i,j,k),cubic_z(i,j-1,k),cubic_z(i,j-2,k),cubic_z(i,j-3,k)));
  else
    return(Cint(j-j0,cubic_z(i,j-1,k),cubic_z(i,j,k),cubic_z(i,j+1,k),cubic_z(i,j+2,k)));  
}


double grid3D::cubic_z(double i, double j, double k){
  int i0=(int)i;
  int j0=(int)j;
  int k0=(int)k;
  
  if (N3<4)
    return(Lint(k-k0,grid[i0][j0][k0],grid[i0][j0][k0+1]));
  
  if (k0==0)
    return(Cint_fwd(k-k0,grid[i0][j0][k0],grid[i0][j0][k0+1],grid[i0][j0][k0+2],grid[i0][j0][k0+3]));
  else if (k0==N3-1)
    return(Cint_fwd(1-(k-k0),grid[i0][j0][k0],grid[i0][j0][k0-1],grid[i0][j0][k0-2],grid[i0][j0][k0-3]));
  else
    return(Cint(k-k0,grid[i0][j0][k0-1],grid[i0][j0][k0],grid[i0][j0][k0+1],grid[i0][j0][k0+2]));  
}


//-----------------------------------------------------------------------------
grid3D::grid3D(int n1, int n2, int n3, int bound, double initialVal)
:N1(n1),N2(n2),N3(n3),boundary(bound),coarse(NULL),fine(NULL){

  //Allocate space for the grid
  allocate(N1,N2,N3,boundary);
  (*this)=initialVal;
}
//This function allows a grid3D object to be read from a file
//-----------------------------------------------------------------------------
grid3D::grid3D(const char *file, int bound)
:boundary(bound),coarse(NULL),fine(NULL){
  ifstream inFile(file, ios::in);
  inFile >> N1 >> N2 >> N3;

  allocate(N1,N2,N3,boundary);

  //Read in the grid from file handle
  for (int k=0; k<N3; ++k)
    for (int j=0; j<N2; ++j)
      for (int i=0; i<N1; ++i){
        inFile >> (*this)(i,j,k);
      }
  
  inFile.close();
}

//-----------------------------------------------------------------------------
grid3D::~grid3D(void){
  grid-=boundary;
  for (int i=0; i<N1+2*boundary; ++i){
    grid[i]-=boundary;
    delete [] grid[i];
  }
  delete [] grid;
  delete [] data;

  //Delete coarse and fine grids if they've been allocated
  if (coarse){
    coarse->fine=fine;
    delete coarse;
  }

  if (fine){
    fine->coarse=coarse;
    delete fine;
  }
}

//-----------------------------------------------------------------------------
void grid3D::initializeRandom(double min, double max){
  //srand((unsigned)time(NULL));
  gridLoop{
    (*this)(i,j,k)=((double)rand()/RAND_MAX)*(max-min)+min;
  }
}
//-----------------------------------------------------------------------------
/* This method uses a polar Box-Mueller transformation to convert a uniform
 * distribution in the range 0 to 1 to a gaussian distribution */
void grid3D::initializeGaussian(double sigma){
  double x1, x2, w, y1, y2;
  //srand((unsigned)time(NULL));
  gridLoop{
    w=1;
    while (!(w < 1)){
      x1=2*(double)rand()/RAND_MAX-1;
      x2=2*(double)rand()/RAND_MAX-1;
      w=x1*x1+x2*x2;
    }
    w=sqrt(-2*log(w)/w);
    y1=sqrt(sigma)*x1*w;
    y2=sqrt(sigma)*x2*w;

    //Add the random numbers to the grid
    (*this)(i,j,k)=y1;
    k++;
    if (k<N3) {
      (*this)(i,j,k)=y2;
    }
  }
}
//-----------------------------------------------------------------------------
void grid3D::initializeSphere(double r){
  gridLoop{
    (*this)(i,j,k)=(sq(i-N1/2)+sq(j-N2/2)+sq(k-N3/2))<=sq(r) ? 1 : 0;
  }
}
//-----------------------------------------------------------------------------
double grid3D::mean(void){
  double sum=0;
  gridLoop{
    sum+=(*this)(i,j,k);
  }
  return(sum/(N1*N2*N3));
}
//-----------------------------------------------------------------------------
double grid3D::sum(){
  double sum=0;
  gridLoop{
    sum+=(*this)(i,j,k);
  }
  return(sum);
}
//-----------------------------------------------------------------------------
void grid3D::periodicBoundary(void){
  xAxisPeriodicBoundary();
  yAxisPeriodicBoundary();
  zAxisPeriodicBoundary();
}

/* The 'ext' variables are used to determine if values on other boundary should
 * be passed as part of the periodic boundary.  These variables must be either
 * one or zero.  ext1 and ext2 variables represent the directions i,j,k in that
 * order for the other two axes.
 */
//-----------------------------------------------------------------------------
void grid3D::xAxisPeriodicBoundary(int ext1, int ext2){
  for (int j=-boundary*ext1; j<N2+boundary*ext1; j++)
    for (int k=-boundary*ext2; k<N3+boundary*ext2; k++){
      (*this)(-1,j,k)=(*this)(N1-1,j,k);
      (*this)(N1,j,k)=(*this)(0,j,k);
    }
}
//-----------------------------------------------------------------------------
void grid3D::yAxisPeriodicBoundary(int ext1, int ext2){
  for (int i=-boundary*ext1; i<N1+boundary*ext1; i++)
    for (int k=-boundary*ext2; k<N3+boundary*ext2; k++){
      (*this)(i,-1,k)=(*this)(i,N2-1,k);
      (*this)(i,N2,k)=(*this)(i,0,k);
    }
}
//-----------------------------------------------------------------------------
void grid3D::zAxisPeriodicBoundary(int ext1, int ext2){
  for (int i=-boundary*ext1; i<N1+boundary*ext1; i++)
    for (int j=-boundary*ext2; j<N2+boundary*ext2; j++){
      (*this)(i,j,-1)=(*this)(i,j,N3-1);
      (*this)(i,j,N3)=(*this)(i,j,0);
    }
}
//-----------------------------------------------------------------------------
void grid3D::xAxisNeumannBoundary(double nh, int ext1, int ext2){
  for (int j=-boundary*ext1; j<N2+boundary*ext1; j++)
    for (int k=-boundary*ext2; k<N3+boundary*ext2; k++){
      (*this)(-1,j,k)=(*this)(0,j,k)-nh;
      (*this)(N1,j,k)=nh+(*this)(N1-1,j,k);
    }
}

void grid3D::yAxisNeumannBoundary(double nh, int ext1, int ext2){
  for (int i=-boundary*ext1; i<N1+boundary*ext1; i++)
    for (int k=-boundary*ext2; k<N3+boundary*ext2; k++){
      (*this)(i,-1,k)=(*this)(i,0,k)-nh;
      (*this)(i,N2,k)=nh+(*this)(i,N2-1,k);
    }
}

void grid3D::zAxisNeumannBoundary(double nh, int ext1, int ext2){
  for (int i=-boundary*ext1; i<N1+boundary*ext1; i++)
    for (int j=-boundary*ext2; j<N2+boundary*ext2; j++){
      (*this)(i,j,-1)=(*this)(i,j,0)-nh;
      (*this)(i,j,N3)=nh+(*this)(i,j,N3-1);
  }
}

void grid3D::neumannBoundary(){
  xAxisNeumannBoundary(0);
  yAxisNeumannBoundary(0);
  zAxisNeumannBoundary(0);
}
//-----------------------------------------------------------------------------
void grid3D::xAxisDirichletBoundary(double c, int ext1, int ext2){
  for (int j=-boundary*ext1; j<N2+boundary*ext1; j++)
    for (int k=-boundary*ext2; k<N3+boundary*ext2; k++){
      (*this)(-1,j,k)=c;
      (*this)(N1,j,k)=c;
    }
}

void grid3D::yAxisDirichletBoundary(double c, int ext1, int ext2){
  for (int i=-boundary*ext1; i<N1+boundary*ext1; i++)
    for (int k=-boundary*ext2; k<N3+boundary*ext2; k++){
      (*this)(i,-1,k)=c;
      (*this)(i,N2,k)=c;
    }
}

void grid3D::neg_yAxisDirichletBoundary(double c, int ext1, int ext2){
  for (int i=-boundary*ext1; i<N1+boundary*ext1; i++)
    for (int k=-boundary*ext2; k<N3+boundary*ext2; k++){
      (*this)(i,-1,k)=c;
  }
}

void grid3D::pos_yAxisDirichletBoundary(double c, int ext1, int ext2){
  for (int i=-boundary*ext1; i<N1+boundary*ext1; i++)
    for (int k=-boundary*ext2; k<N3+boundary*ext2; k++){
      (*this)(i,N2,k)=c;
    }
}

void grid3D::zAxisDirichletBoundary(double c, int ext1, int ext2){
  for (int i=-boundary*ext1; i<N1+boundary*ext1; i++)
    for (int j=-boundary*ext2; j<N2+boundary*ext2; j++){
      (*this)(i,j,-1)=c;
      (*this)(i,j,N3)=c;
    }
}

void grid3D::dirichletBoundary(double c){
  xAxisDirichletBoundary(c);
  yAxisDirichletBoundary(c);
  zAxisDirichletBoundary(c);
}
//-----------------------------------------------------------------------------
//Function to perform prolongation by interpolation
//Dimensions for the fine matrix must be provided - they can't be determined
//from the size of the coarse matrix
grid3D* grid3D::prolongate(int Nx, int Ny, int Nz){
  if (!fine){
    fine=new grid3D(Nx,Ny,Nz,boundary);
    fine->coarse=this;
  }

  #pragma omp parallel for collapse(3)
  gridLoop3D(*fine){
    (*fine)(i,j,k)=(*this)((double)i/2,(double)j/2,(double)k/2); 
  }
  return(fine);
}
grid3D* grid3D::prolongate_CC(int Nx, int Ny, int Nz){
  if (!fine){
    fine=new grid3D(Nx,Ny,Nz,boundary);
    fine->coarse=this;
  }

  #pragma omp parallel for collapse(3)
  gridLoop3D(*fine){
    //Injection
//    (*fine)(i,j,k)=(*this)((int)i/2,(int)j/2,(int)k/2);
    
    //Bilinear interpolation
    (*fine)(i,j,k)=(*this)((double)i/2-.25,(double)j/2-.25,(double)k/2-.25);
  }
  return(fine);
}
grid3D* grid3D::prolongate_cubic(int Nx, int Ny, int Nz){
  if (!fine){
    fine=new grid3D(Nx,Ny,Nz,boundary);
    fine->coarse=this;
  }

  #pragma omp parallel for collapse(3)
  gridLoop3D(*fine){
    (*fine)(i,j,k)=cubic((double)i/2,(double)j/2,(double)k/2); 
  }
  return(fine);
}
//-----------------------------------------------------------------------------
//Restriction with the FW operator
grid3D* grid3D::restrict_FW(){
  if (!coarse){
    int N2x=(N1+1)/2;
    int N2y=(N2+1)/2;
    int N2z=(N3+1)/2;

    coarse=new grid3D(N2x,N2y,N2z,boundary);
    coarse->fine=this;
  }

  #pragma omp parallel for collapse(3)
  gridLoop3D(*coarse){
    (*coarse)(i,j,k)=1./64*(8*(*this)(2*i,2*j,2*k)
     +4*((*this)(2*i+1,2*j,2*k)+(*this)(2*i-1,2*j,2*k)
        +(*this)(2*i,2*j+1,2*k)+(*this)(2*i,2*j-1,2*k)
        +(*this)(2*i,2*j,2*k+1)+(*this)(2*i,2*j,2*k-1))
     +2*((*this)(2*i+1,2*j+1,2*k)+(*this)(2*i-1,2*j-1,2*k)
        +(*this)(2*i+1,2*j-1,2*k)+(*this)(2*i-1,2*j+1,2*k)
        +(*this)(2*i+1,2*j,2*k+1)+(*this)(2*i-1,2*j,2*k-1)
        +(*this)(2*i+1,2*j,2*k-1)+(*this)(2*i-1,2*j,2*k+1)
        +(*this)(2*i,2*j+1,2*k+1)+(*this)(2*i,2*j-1,2*k-1)
        +(*this)(2*i,2*j+1,2*k-1)+(*this)(2*i,2*j-1,2*k+1))
      +(*this)(2*i+1,2*j+1,2*k+1)+(*this)(2*i-1,2*j-1,2*k-1)
      +(*this)(2*i-1,2*j+1,2*k+1)+(*this)(2*i+1,2*j-1,2*k-1)
      +(*this)(2*i+1,2*j-1,2*k+1)+(*this)(2*i-1,2*j+1,2*k-1)
      +(*this)(2*i+1,2*j+1,2*k-1)+(*this)(2*i-1,2*j-1,2*k+1));
  }
  
  return(coarse);
}
//-----------------------------------------------------------------------------
//Restriction with the HW operator
grid3D* grid3D::restrict_HW(){
  if (!coarse){
    int N2x=(N1+1)/2;
    int N2y=(N2+1)/2;
    int N2z=(N3+1)/2;

    coarse=new grid3D(N2x,N2y,N2z,boundary);
    coarse->fine=this;
  }

  #pragma omp parallel for collapse(3)
  gridLoop3D(*coarse){
    (*coarse)(i,j,k)=1./12*(6*(*this)(2*i,2*j,2*k)
     +(*this)(2*i+1,2*j,2*k)+(*this)(2*i-1,2*j,2*k)
     +(*this)(2*i,2*j+1,2*k)+(*this)(2*i,2*j-1,2*k)
     +(*this)(2*i,2*j,2*k+1)+(*this)(2*i,2*j,2*k-1));
  }
  
  return(coarse);
}
//-----------------------------------------------------------------------------
grid3D* grid3D::restrict_CC(){
  if (!coarse){
    int N2x=(N1+1)/2;
    int N2y=(N2+1)/2;
    int N2z=(N3+1)/2;

    coarse=new grid3D(N2x,N2y,N2z,boundary);
    coarse->fine=this;
  }

  #pragma omp parallel for collapse(3)
  gridLoop3D(*coarse){
    (*coarse)(i,j,k)=1./8*((*this)(2*i,2*j,2*k)
      +(*this)(2*i+1,2*j,2*k)+(*this)(2*i,2*j+1,2*k)+(*this)(2*i,2*j,2*k+1)
      +(*this)(2*i,2*j+1,2*k+1)+(*this)(2*i+1,2*j,2*k+1)+(*this)(2*i+1,2*j+1,2*k)
      +(*this)(2*i+1,2*j+1,2*k+1));
  }
  
  return(coarse);
}
//-----------------------------------------------------------------------------
//Restriction by injection
grid3D* grid3D::injection(){
  if (!coarse){
    int N2x=(N1+1)/2;
    int N2y=(N2+1)/2;
    int N2z=(N3+1)/2;

    coarse=new grid3D(N2x,N2y,N2z,boundary);
    coarse->fine=this;
  }

  #pragma omp parallel for collapse(3)
  gridLoop3D(*coarse){
    (*coarse)(i,j,k)=(*this)(2*i,2*j,2*k);
  }
  
  return(coarse);
}
//-----------------------------------------------------------------------------
grid3D* grid3D::getCoarseGrid(){
  if (!coarse){
    int N2x=(N1+1)/2;
    int N2y=(N2+1)/2;
    int N2z=(N3+1)/2;
    coarse=new grid3D(N2x,N2y,N2z,boundary);
    coarse->fine=this;
  }
  return(coarse);
}
//-----------------------------------------------------------------------------
double grid3D::l2_norm(){
  double sum=0;
  gridLoop{
      sum += sq((*this)(i,j,k));
  }
  return(sqrt(sum/(N1*N2*N3)));
}
//-----------------------------------------------------------------------------
void grid3D::writeToFile(const char* file){
  writeToFile(file, 0);
}
//-----------------------------------------------------------------------------
void grid3D::writeToFile(const char *file, int includeBndry){
  int ref_pt;
  ofstream outFile(file);
  
  if (includeBndry){
    ref_pt=boundary;
  }
  else{
    ref_pt=0;
  }
  
  outFile << N1+2*ref_pt << " " << N2+2*ref_pt << " " << N3+2*ref_pt << "\n";

  for (int k=-ref_pt; k<N3+ref_pt; k++){
    for (int j=-ref_pt; j<N2+ref_pt; j++){  
      for (int i=-ref_pt; i<N1+ref_pt; i++)
        outFile << (*this)(i,j,k) << " ";
      outFile << "\n";
    }
  }
  outFile.close();
}
//-----------------------------------------------------------------------------
void grid3D::writeToFileDx(const char* file){
  ofstream outFile(file);

  outFile << "object 1 class gridpositions counts " << N1 << " " << N2 << " " << N3 << "\n";
  outFile << "origin 0 0 0" << "\n";
  outFile << "delta 1 0 0" << "\n";
  outFile << "delta 0 1 0" << "\n";
  outFile << "delta 0 0 1" << "\n";
  outFile << "object 2 class gridconnections counts " << N1 <<  " " << N2 << " "<< N3 << "\n";
  outFile << "object 3 class array type double rank 0 items " << N1*N2*N3 << " data follows" << "\n";

  int ret=0;
  for (int i=0; i<N1; i++)
    for (int j=0; j<N2; j++)
      for (int k=0; k<N3; k++){
        outFile << (*this)(i,j,k) << " ";
        if (ret==2)
          outFile << "\n";
        ret=(ret+1)%3;
      }
  outFile << "\n";
  outFile << "object 'awesome' class field" << "\n";
  outFile.close();
}
//-----------------------------------------------------------------------------
void grid3D::writeToFileNc(const char* file){
  
  int ncid, x_dimid, y_dimid, z_dimid, varid;
  int dimids[3], retval, i, j, k;
  
  if (retval = nc_create(file, NC_CLOBBER, &ncid))
     ERR(retval);

  if (retval = nc_def_dim(ncid, "x", N1, &x_dimid))
    ERR(retval);
  if (retval = nc_def_dim(ncid, "y", N2, &y_dimid))
    ERR(retval);
  if (retval = nc_def_dim(ncid, "z", N3, &z_dimid))
    ERR(retval);
  
  dimids[0] = x_dimid;
  dimids[1] = y_dimid;
  dimids[2] = z_dimid;

  if (retval = nc_def_var(ncid, "density", NC_FLOAT, 3, dimids, &varid))
    ERR(retval);

  if (retval = nc_enddef(ncid))
    ERR(retval);
  
  //Because of boundaries on grid3D, write to the file as a series of vectors
  //along the z-dimension.
  for (size_t i=0;i<N1;++i){
    for (size_t j=0;j<N2;++j){
      size_t start[3]={i,j,0};
      size_t count[3]={1,1,(size_t)N3};
      if (retval = nc_put_vara_double(ncid, varid, start, count, &grid[i][j][0]))
        ERR(retval);
    }
  }
  
  if (retval = nc_close(ncid))
    ERR(retval);
}
//-----------------------------------------------------------------------------
void grid3D::operator=(const double value){
  #pragma omp parallel for
  fastBndryLoop{
      data[n]=value;
  }
}
//-----------------------------------------------------------------------------
void grid3D::operator=(grid3D& grid){
  #pragma omp parallel for
  fastBndryLoop{
      data[n]=grid.data[n];
  }
}
//-----------------------------------------------------------------------------
void grid3D::operator+=(grid3D& grid){
  #pragma omp parallel for
  fastBndryLoop{
      data[n]+=grid.data[n];
    }
}
//-----------------------------------------------------------------------------
void grid3D::operator-=(grid3D& grid){
  #pragma omp parallel for
  fastBndryLoop{
      data[n]-=grid.data[n];
  }
}
//-----------------------------------------------------------------------------
void grid3D::operator*=(const double value){
  #pragma omp parallel for
  fastBndryLoop{
      data[n]*=value;
  }
}
//-----------------------------------------------------------------------------
void grid3D::allocate(int N1, int N2, int N3, int boundary){
  int n1=N1+2*boundary;
  int n2=N2+2*boundary;
  int n3=N3+2*boundary;
  data = new double[n1*n2*n3];
  grid = new double**[n1];
  
  //Create pointers to the 3 dimensions of the grid, 
  //and set the origin at (boundary,boundary,boundary)
  for (int i=0; i<n1; ++i){
    grid[i] = new double*[n2];
    for (int j=0; j<n2; ++j){
      grid[i][j] = data+i*n2*n3+j*n3;
      grid[i][j]+=boundary;
    }
    grid[i]+=boundary;
  }
  grid+=boundary;
}
