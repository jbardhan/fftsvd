#include "FFTSVD.h"
unsigned int  num_GMRES_iter = 1000;
int main(int argc, char* argv[]) {

  Matrix A = Matrix_allocate(5,4);
  A[0][0]= 0.5;
  A[0][1]= 0.5;
  A[0][2]= 0.5;
  A[0][3]= 0.5;

  A[1][0]= 0.5;
  A[1][1]= -1.5;
  A[1][2]= 0.5;
  A[1][3]= -1.5;

  A[2][0]= 1.207107;
  A[2][1]= -0.5;
  A[2][2]= 0.207107;
  A[2][3]= -0.5;

  A[3][0]= 0.;
  A[3][1]= 2.414214;
  A[3][2]= 0.;
  A[3][3]= -0.414214;

  A[4][0]= 1.707107;
  A[4][1]= 0.707107;
  A[4][2]= 0.292893;
  A[4][3]= -0.707107;

  Matrix Q = Matrix_allocate(5,4);
  Matrix R = Matrix_allocate(5,4);
  Matrix_QR(A, Q, R, 5, 4);

  /* Matrix_writefile("R",R, 5, 4); */
  /* Matrix_writefile("Q",Q, 5, 4); */

  Matrix_free(A);
  Matrix_free(Q);
  Matrix_free(R);

  Matrix Asym = Matrix_allocate(4,4);
  Matrix V    = Matrix_allocate(4,4);
  Vector d    = Vector_allocate(4);
  Asym[0][0] = 2.8966;
  Asym[0][1] = 2.1881;
  Asym[0][2] = 1.1965;
  Asym[0][3] = 2.1551;
  Asym[1][0] = 2.1881;
  Asym[1][1] = 1.9966;
  Asym[1][2] = 0.6827;
  Asym[1][3] = 1.8861;
  Asym[2][0] = 1.1965;
  Asym[2][1] = 0.6827;
  Asym[2][2] = 0.7590;
  Asym[2][3] = 0.5348;
  Asym[3][0] = 2.1551;
  Asym[3][1] = 1.8861;
  Asym[3][2] = 0.5348;
  Asym[3][3] = 2.0955;

  Matrix_eigendecomposition(Asym, V, d, 4, 4);
  Matrix_writefile("V",V,4,4);
  Vector_writefile("dvec",d,4);

  Vector_free(d);
  Matrix_free(V);
  Matrix_free(Asym);

}
