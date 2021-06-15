/* Intel Performance Primitives (IPP) prototypes */

/*
complex double* w7_ippsMalloc_64fc(int len);
complex float* w7_ippsMalloc_32fc(int len);
complex double* m7_ippsMalloc_64fc(int len);
complex float* m7_ippsMalloc_32fc(int len);

void w7_ippsFree(void* ptr);
void m7_ippsFree(void* ptr);

void w7_ippsAddProduct_64fc(complex double* v1, complex double* v2, complex double* vaddprod, unsigned int length);
void w7_ippsAddProduct_32fc(complex float* v1, complex float* v2, complex float* vaddprod, unsigned int length);
void m7_ippsAddProduct_64fc(complex double* v1, complex double* v2, complex double* vaddprod, unsigned int length);
void m7_ippsAddProduct_32fc(complex float* v1, complex float* v2, complex float* vaddprod, unsigned int length);
*/

/* LAPACK prototypes */

#ifdef _AIX
void dgesvd(char* JOBU, char* JOBVT, unsigned int* M, unsigned int* N, 
             double* A, unsigned int* LDA, double* S, double* U,
             unsigned int* LDU, double* VT, unsigned int* LDVT, double* WORK,
             unsigned int* LWORK, int* INFO);
void sgesvd(char* JOBU, char* JOBVT, unsigned int* M, unsigned int* N, 
             float* A, unsigned int* LDA, float* S, float* U,
             unsigned int* LDU, float* VT, unsigned int* LDVT, float* WORK,
             unsigned int* LWORK, int* INFO);
void dsyev(char* JOBZ, char* UPLO, unsigned int* N, double* A,
            unsigned int* LDA, double* W, double* WORK, unsigned int* LWORK,
            int* INFO);
void ssyev(char* JOBZ, char* UPLO, unsigned int* N, float* A,
            unsigned int* LDA, float* W, float* WORK, unsigned int* LWORK,
            int* INFO);
void dgels(char* TRANS, unsigned int* M, unsigned int* N, unsigned int* NRHS,
            double* A, unsigned int* LDA, double* B, unsigned int* LDB, 
            double* WORK, unsigned int* LWORK, int* INFO);
void sgels(char* TRANS, unsigned int* M, unsigned int* N, unsigned int* NRHS,
            float* A, unsigned int* LDA, float* B, unsigned int* LDB, 
            float* WORK, unsigned int* LWORK, int* INFO);
#else
void dgesvd(char* JOBU, char* JOBVT, unsigned int* M, unsigned int* N, 
             double* A, unsigned int* LDA, double* S, double* U,
             unsigned int* LDU, double* VT, unsigned int* LDVT, double* WORK,
             unsigned int* LWORK, int* INFO);
void sgesvd(char* JOBU, char* JOBVT, unsigned int* M, unsigned int* N, 
             float* A, unsigned int* LDA, float* S, float* U,
             unsigned int* LDU, float* VT, unsigned int* LDVT, float* WORK,
             unsigned int* LWORK, int* INFO);
void dsyev(char* JOBZ, char* UPLO, unsigned int* N, double* A,
            unsigned int* LDA, double* W, double* WORK, unsigned int* LWORK,
            int* INFO);
void ssyev(char* JOBZ, char* UPLO, unsigned int* N, float* A,
            unsigned int* LDA, float* W, float* WORK, unsigned int* LWORK,
            int* INFO);
void dgels(char* TRANS, unsigned int* M, unsigned int* N, unsigned int* NRHS,
            double* A, unsigned int* LDA, double* B, unsigned int* LDB, 
            double* WORK, unsigned int* LWORK, int* INFO);
void sgels(char* TRANS, unsigned int* M, unsigned int* N, unsigned int* NRHS,
            float* A, unsigned int* LDA, float* B, unsigned int* LDB, 
            float* WORK, unsigned int* LWORK, int* INFO);
#endif

/* Inline Functions */

static inline real intpow(real r, unsigned int p) {
   switch (p) {
      case 0: return 1.0;
      case 1: return r;
      case 2: return r*r;
      case 3: return r*r*r;
      case 4: { real s = r*r; return s*s; }
      case 5: { real s = r*r; return s*s*r; }
      case 6: { real s = r*r; return s*s*s; }
      default: {
         real ans = 1.0;
                                                                                
         while (p != 0) {
            if (p & 1) ans *= r;
            p >>= 1;
            r *= r;
         }
                                                                                
         return ans;
      }
   }
}

static inline unsigned int uimin(unsigned int a, unsigned int b) {
   if (a < b)
      return a;
   else
      return b;
}
                                                                                
static inline unsigned int uimax(unsigned int a, unsigned int b) {
   if (a > b)
      return a;
   else
      return b;
}
