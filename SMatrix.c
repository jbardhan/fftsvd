#include "FFTSVD.h"

/* Constructors and Destructors */

SMatrix SMatrix_allocate(unsigned int rows, unsigned int columns) {
   unsigned int i;
   SMatrix matrix = NULL;
   
   matrix = (SMatrix)calloc(rows, sizeof(SVector));
   matrix[0] = (SVector)calloc(rows * columns, sizeof(sreal));
   
#pragma ivdep
   for (i = 1; i < rows; i++)
      matrix[i] = matrix[0] + i * columns;
   
   return matrix;
}

void SMatrix_free(SMatrix matrix) {
   if (matrix != NULL)
      free(matrix[0]);
   free(matrix);
}

/* Initialization and Copying */

void SMatrix_copy(SMatrix matrixdest, SMatrix matrixsrc, unsigned int rows, unsigned int columns) {
   unsigned int i;   
   SVector pdest = matrixdest[0];
   SVector psrc = matrixsrc[0];
   
#pragma ivdep
   for (i = 0; i < rows * columns; i++)
      pdest[i] = psrc[i];
}

// this function is horribly inefficient but is only needed in a few
// minor cases, once or twice per run.
void SMatrix_copypiece(SMatrix matrixdest, unsigned int destStartRow, unsigned int destStartCol,
                      SMatrix matrixsrc, unsigned int srcStartRow, unsigned int srcStartCol,
                      unsigned int numrows, unsigned int numcols) {
   unsigned int i, j;
   for (i = 0; i < numrows; i++) {
      for (j = 0; j < numcols; j++) {
         matrixdest[destStartRow+i][destStartCol+j] = matrixsrc[srcStartRow+i][srcStartCol+j];
      }
   }
}

void SMatrix_zero(SMatrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i;   
   SVector pmatrix = matrix[0];
   
#pragma ivdep
   for (i = 0; i < rows * columns; i++)
      pmatrix[i] = (sreal)0;
}

void SMatrix_identity(SMatrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i;
   
   #pragma ivdep
   for (i = 0; i < uimin(rows, columns); i++)
      matrix[i][i] = (sreal)1;
}

void SMatrix_Matrix(SMatrix smatrix, Matrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         smatrix[i][j] = matrix[i][j];
}

/* Linear Algebra Operations */

void SMatrix_stats(SMatrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i;   
   SVector pmatrix = matrix[0];
   sreal min = 1e10, max = -1e10;
   
   #pragma ivdep
   for (i = 0; i < rows * columns; i++) {
      if (pmatrix[i] > max) max = pmatrix[i];
      if (pmatrix[i] < min) min = pmatrix[i];
   }

   printf("MIN: %e MAX: %e\n", min, max);
}

/* SMatrix Frobenius norm */

sreal SMatrix_norm(SMatrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i;
   sreal norm = (sreal)0;
   SMatrix temp = SMatrix_allocate(columns, columns);
   
   SMatrix_multiplymatrix_transpose(temp, matrix, matrix, rows, columns, columns);
   
#pragma ivdep
   for (i = 0; i < columns; i++)
      norm += temp[i][i];
   
   norm = sqrt(norm);
   
   SMatrix_free(temp);
   
   return norm;
}

/* b = A*x */
void SMatrix_multiplyvector(SVector b, SMatrix A, SVector x, unsigned int rows, unsigned int columns) {
   unsigned int i, j;

   for (i = 0; i < rows; i++) {
      SVector Ai = A[i];
#pragma ivdep
      for (j = 0; j < columns; j++)
         b[i] += Ai[j] * x[j];
   }
}

/* b = A'*x */
void SMatrix_multiplyvector_transpose(SVector b, SMatrix A, SVector x, unsigned int rows, unsigned int columns) {
   unsigned int i, j;

   for (i = 0; i < rows; i++) {
      SVector Ai = A[i];
      sreal xi = x[i];
#pragma ivdep
      for (j = 0; j < columns; j++)
         b[j] += Ai[j] * xi;
   }
}

/* B = A*X */
void SMatrix_multiplymatrix(SMatrix B, SMatrix A, SMatrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX) {
   unsigned int i, j, k;
   
   SMatrix_zero(B, rowsA, columnsX);
   
   for (i = 0; i < rowsA; i++) {
      SVector Bi = B[i];
      for (j = 0; j < columnsA; j++) {
         sreal Aij = A[i][j];
         SVector Xj = X[j];
#pragma ivdep
         for (k = 0; k < columnsX; k++)
            Bi[k] += Aij * Xj[k];
      }
   }
}

/* B = A'*X */
void SMatrix_multiplymatrix_transpose(SMatrix B, SMatrix A, SMatrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX) {
   unsigned int i, j, k;
   
   SMatrix_zero(B, columnsA, columnsX);
   
   for (i = 0; i < columnsA; i++) {
      SVector Bi = B[i];
      for (j = 0; j < rowsA; j++) {
         sreal Aji = A[j][i];
         SVector Xj = X[j];
#pragma ivdep
         for (k = 0; k < columnsX; k++)
            Bi[k] += Aji * Xj[k];
      }
   }
}

/* B = A*X' */
void SMatrix_multiplytranspose_matrix(SMatrix B, SMatrix A, SMatrix X, unsigned int rowsA, unsigned int columnsA, unsigned int rowsX) {
   unsigned int i, j, k;
   
   SMatrix_zero(B, rowsA, rowsX);
   
   for (i = 0; i < rowsA; i++) {
      SVector Bi = B[i];
      for (j = 0; j < columnsA; j++) {
         sreal Aij = A[i][j];
#pragma ivdep
         for (k = 0; k < rowsX; k++)
            Bi[k] += Aij * X[k][j];
      }
   }
}

void SMatrix_rowbasis_pivotedmgs(SMatrix* U, unsigned int* rowrank, SMatrix matrix, unsigned int rows, unsigned int columns, sreal epsilon) {
   SMatrix transpose = SMatrix_allocate(rows, columns);
   
   SMatrix_copy(transpose, matrix, rows, columns);
   
   SMatrix_transpose(&transpose, rows, columns);
   
   SMatrix_columnbasis_pivotedmgs(U, rowrank, transpose, columns, rows, epsilon);
   
   SMatrix_transpose(U, columns, *rowrank);
   
   SMatrix_free(transpose);
}

void SMatrix_columnbasis_pivotedmgs(SMatrix* V, unsigned int* columnrank, SMatrix matrix, unsigned int rows, unsigned int columns, sreal epsilon) {
   unsigned int i, j, k, l, reorth;
   SMatrix A = SMatrix_allocate(rows, columns);
   SMatrix Q = SMatrix_allocate(rows, columns);
   SMatrix R = SMatrix_allocate(rows, columns);
   sreal matrixnorm = SMatrix_norm(matrix, rows, columns);

   SMatrix_copy(A, matrix, rows, columns);
   
   for (k = 0; k < columns; k++) {
      sreal currentnorm, maxnorm = (sreal)0;
      unsigned int largestnormcolumn = 0;

      for (l = k; l < columns; l++) {
         currentnorm = (sreal)0;
         
         for (i = 0; i < rows; i++)
            currentnorm += A[i][l] * A[i][l];
         currentnorm = sqrt(currentnorm);
         
         if (currentnorm > maxnorm) {
            maxnorm = currentnorm;
            largestnormcolumn = l;
         }
      }
      
      if ((k >= 1) && (maxnorm < (epsilon * matrixnorm))) {
         *columnrank = k;
         break;
      }
      
      for (i = 0; i < rows; i++) {
         sreal temp = A[i][largestnormcolumn];
         A[i][largestnormcolumn] = A[i][k];
         A[i][k] = temp;
      }

      for (i = 0; i < rows; i++) {
         sreal temp = R[i][largestnormcolumn];
         R[i][largestnormcolumn] = R[i][k];
         R[i][k] = temp;
      }

      R[k][k] = maxnorm;
      
      for (i = 0; i < rows; i++)
         Q[i][k] = A[i][k] / maxnorm;
         
      if (k < (columns - 1)) {
         for (reorth = 0; reorth < 2; reorth++) {
            SVector Rk = R[k];
      
            for (j = k + 1; j < columns; j++)
               Rk[j] = (sreal)0;

            for (i = 0; i < rows; i++) {
               sreal Qik = Q[i][k];
               SVector Ai = A[i];
#pragma ivdep
               for (j = k + 1; j < columns; j++)
                  Rk[j] += Qik * Ai[j];
            }
         
            for (i = 0; i < rows; i++) {
               sreal Qik = Q[i][k];
               SVector Ai = A[i];
#pragma ivdep
               for (j = k + 1; j < columns; j++)
                  Ai[j] -= Qik * Rk[j];
            }
         }
      }
      
      if (k == (columns - 1))
         *columnrank = columns;
   }

   *V = SMatrix_allocate(rows, *columnrank);
   
   for (i = 0; i < rows; i++) {
      SVector Vi = (*V)[i];
      SVector Qi = Q[i];
      unsigned int r = *columnrank;
#pragma ivdep
      for (j = 0; j < r; j++)
         Vi[j] = Qi[j];
   }

   //SMatrix_columnbasis_check(*V, rows, *columnrank);
   
   SMatrix_free(A);
   SMatrix_free(Q);
   SMatrix_free(R);
}

void SMatrix_columnbasis_check(SMatrix V, unsigned int rows, unsigned int columns) {
   unsigned int c1, c2, r;
   sreal maxselferror = 0.0, maxcrosserror = 0.0, error;

   for (c1 = 0; c1 < columns; c1++)
      for (c2 = c1; c2 < columns; c2++) {
         sreal dot = 0.0;
         for (r = 0; r < rows; r++)
            dot += V[r][c1] * V[r][c2];
         if (c1 == c2) {
            error = fabs(1.0 - dot);
            if (error > maxselferror)
               maxselferror = error;
         }
         else {
            error = fabs(dot);
            if (error > maxcrosserror)
               maxcrosserror = error;
         }
      }

   printf("MAXSELFERROR: %g MAXCROSSERROR: %g\n", maxselferror, maxcrosserror);
}

void SMatrix_transpose(SMatrix* matrix, unsigned int rows, unsigned int columns) {
   unsigned int i, j;

   SMatrix transpose = SMatrix_allocate(columns, rows);
   
   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         transpose[j][i] = (*matrix)[i][j];
         
   SMatrix_free(*matrix);
   
   *matrix = transpose;
}

void SMatrix_diff(SMatrix X, SMatrix A, SMatrix B, unsigned int rows, unsigned int columns) {
   unsigned int i, j;
   
   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         X[i][j] = A[i][j] - B[i][j];
}

void SMatrix_pseudoinverse_droptol(SMatrix XI, SMatrix X, unsigned int rows, unsigned int columns, sreal droptol) {
   char JOBU = 'A', JOBVT = 'A';
   unsigned int m = rows, n = columns, LDA = rows, LDU = rows, LDVT = columns;
   unsigned int mindim = uimin(rows, columns);
   unsigned int LWORK = uimax(3*uimin(m,n)+uimax(m,n),5*uimin(m,n));
   int INFO;
   SMatrix Xcolumnmajor = SMatrix_allocate(columns, rows);
   SMatrix UT = SMatrix_allocate(rows, rows);
   SMatrix S = SMatrix_allocate(columns, rows);
   SMatrix V = SMatrix_allocate(columns, columns);
   SMatrix temp = SMatrix_allocate(columns, rows);
   sreal D[mindim], WORK[LWORK], tol;
   unsigned int i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         Xcolumnmajor[j][i] = X[i][j];

#ifdef SREAL_IS_DOUBLE
#ifdef OMP
#pragma omp critical (lapack)
#endif
#ifdef _AIX
   dgesvd(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
           UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
#else
   dgesvd_(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
           UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
#endif
#else
#ifdef SREAL_IS_FLOAT
#ifdef OMP
#pragma omp critical (lapack)
#endif
#ifdef _AIX
   sgesvd(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
           UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
#else
   sgesvd_(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
           UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
#endif
#endif
#endif

   if (INFO) {
      printf("s/dgesvd returned error code %d\n", INFO);
      exit(-2);
   }

   SMatrix_free(Xcolumnmajor);

   tol = D[0] * droptol;

   for (i = 0; i < mindim; i++)
      if (D[i] > tol)
         S[i][i] = (sreal)1 / D[i];

   SMatrix_multiplymatrix(temp, V, S, columns, columns, rows);

   SMatrix_free(V);
   SMatrix_free(S);

   SMatrix_multiplymatrix(XI, temp, UT, columns, rows, rows);

   SMatrix_free(temp);
   SMatrix_free(UT);
}

void SMatrix_pseudoinverse(SMatrix XI, SMatrix X, unsigned int rows, unsigned int columns) {
   char JOBU = 'A', JOBVT = 'A';
   unsigned int m = rows, n = columns, LDA = rows, LDU = rows, LDVT = columns;
   unsigned int mindim = uimin(rows, columns);
   unsigned int LWORK = uimax(3*uimin(m,n)+uimax(m,n),5*uimin(m,n));
   int INFO;
   SMatrix Xcolumnmajor = SMatrix_allocate(columns, rows);
   SMatrix UT = SMatrix_allocate(rows, rows);
   SMatrix S = SMatrix_allocate(columns, rows);
   SMatrix V = SMatrix_allocate(columns, columns);
   SMatrix temp = SMatrix_allocate(columns, rows);
   sreal D[mindim], WORK[LWORK], tol;
   unsigned int i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         Xcolumnmajor[j][i] = X[i][j];

#ifdef SREAL_IS_DOUBLE
#ifdef OMP
#pragma omp critical (lapack)
#endif
#ifdef _AIX
   dgesvd(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
           UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
#else
   dgesvd_(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
           UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
#endif
#else
#ifdef SREAL_IS_FLOAT
#ifdef OMP
#pragma omp critical (lapack)
#endif
#ifdef _AIX
   sgesvd(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
           UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
#else
   sgesvd_(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
           UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
#endif
#endif
#endif

   if (INFO) {
      printf("s/dgesvd returned error code %d\n", INFO);
      exit(-2);
   }

   SMatrix_free(Xcolumnmajor);

   /* These tolerances are the default for MATLAB's pinv routine */

#ifdef SREAL_IS_DOUBLE
   tol = uimax(rows, columns) * D[0] * DBL_EPSILON;
#else
#ifdef SREAL_IS_FLOAT
   tol = uimax(rows, columns) * D[0] * FLT_EPSILON;
#endif
#endif

   for (i = 0; i < mindim; i++)
      if (D[i] > tol)
         S[i][i] = (sreal)1 / D[i];

   SMatrix_multiplymatrix(temp, V, S, columns, columns, rows);

   SMatrix_free(V);
   SMatrix_free(S);

   SMatrix_multiplymatrix(XI, temp, UT, columns, rows, rows);

   SMatrix_free(temp);
   SMatrix_free(UT);
}

void SMatrix_lsq_solve(SVector x, SMatrix A, SVector b, unsigned int rows, unsigned int columns) {
   char TRANS = 'N';
   unsigned int M = rows, N = columns, NRHS = 1, LDA = rows, LDB = uimax(rows,columns);
   unsigned int LWORK = uimax(1, M*N + uimax(M*N, NRHS));
   int INFO;
   sreal WORK[LWORK];
   SMatrix Acolumnmajor = SMatrix_allocate(columns, rows);
   SVector bcopy = SVector_allocate(uimax(rows, columns));
   unsigned int i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         Acolumnmajor[j][i] = A[i][j];

   SVector_copy(bcopy, b, rows);

#ifdef SREAL_IS_DOUBLE
#ifdef _AIX
   dgels(&TRANS, &M, &N, &NRHS, Acolumnmajor[0], &LDA, bcopy, &LDB, WORK, &LWORK, &INFO);
#else
   dgels_(&TRANS, &M, &N, &NRHS, Acolumnmajor[0], &LDA, bcopy, &LDB, WORK, &LWORK, &INFO);
#endif
#else
#ifdef SREAL_IS_FLOAT
#endif
#ifdef _AIX
   sgels(&TRANS, &M, &N, &NRHS, Acolumnmajor[0], &LDA, bcopy, &LDB, WORK, &LWORK, &INFO);
#else
   sgels_(&TRANS, &M, &N, &NRHS, Acolumnmajor[0], &LDA, bcopy, &LDB, WORK, &LWORK, &INFO);
#endif
#endif

   if (INFO) {
      printf("s/dgels returned error code %d\n", INFO);
      exit(-2);
   }

   SVector_copy(x, bcopy, columns);

   SVector_free(bcopy);
   SMatrix_free(Acolumnmajor);
}

/* I/O */

void SMatrix_writefile(char* filename, SMatrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i, j;
   FILE* file = NULL;
   
   file = fopen(filename, "w");
   
   for (i = 0; i < rows; i++) {
      for (j = 0; j < columns; j++)
         fprintf(file, "%20.15e ", matrix[i][j]);
      fprintf(file, "\n");
   }
   
   fclose(file);
}

void SMatrix_print(SMatrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i, j;
   
   for (i = 0; i < rows; i++) {
      for (j = 0; j < columns; j++)
         printf("%20.15e ", matrix[i][j]);
      printf("\n");
   }
}

void SMatrix_writebinary(FILE* file, SMatrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int count = fwrite(matrix[0], sizeof(sreal), rows*columns, file);
   if (count != rows*columns)
      printf("Short Write!\n");
}

void SMatrix_readbinary(SMatrix matrix, FILE* file, unsigned int rows, unsigned int columns) {
/*
   char filename[1024], buf[64];
   sprintf(buf, "/proc/self/fd/%d", fileno(file));
   unsigned int len = readlink(buf, &filename[0], sizeof(filename));
   filename[len] = '\0';
   //printf("READING %u %u %u %s\n", ftell(file), rows, columns, filename);
*/
   if (feof(file))
      printf("End of File!\n");
   if (ferror(file))
      printf("File Error %d!\n", ferror(file));
   unsigned int count = fread(matrix[0], sizeof(sreal), rows*columns, file);
   if (count != rows*columns)
      printf("Short Read!\n");
}
