#include "FFTSVD.h"

void givens(real a, real b, real* c, real* s) {
   if (b == 0.0) {
      *c = 1.0;
      *s = 0.0;
   }
   else {
      real r1 = fabs(a);
      real r2 = fabs(b);
      real length = sqrt(r1*r1 + r2*r2);
      *c = a / length;
      *s = -b / length;
   }
}

void givensrotate(real c, real s, real* a, real* b) {
   real a1 = *a;
   real b1 = *b;
   *a = a1 * c - b1 * s;
   *b = c * b1 + s * a1;
}

void GMRES_cap(Tree tree, Preconditioner preconditioner, Vector rhs, Vector sol, real tol) {
   unsigned int size = tree->numpanels;
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   real normr;
   unsigned int i;
   int j, k;
   real residual;

   r = Vector_allocate(size);

   Preconditioner_solve(r, preconditioner, rhs);

   normr = Vector_norm(r, size);

   x = Vector_allocate(size);

   c = Vector_allocate(MAXITERATIONS+1);
   s = Vector_allocate(MAXITERATIONS+1);
   g = Vector_allocate(MAXITERATIONS+1);
   y = Vector_allocate(MAXITERATIONS+1);
   H = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   
   P = Vector_allocate(size);

   g[0] = Vector_norm(r, size);
   bv[0] = Vector_allocate(size);
   Vector_copy(bv[0], r, size);
   Vector_scale(bv[0], 1.0 / g[0], size);

   for (i = 0; i < MAXITERATIONS; i++) {
      Tree_multiply(P, tree, bv[i], SINGLE_LAYER_INT);
      Preconditioner_solve(P, preconditioner, P);
      
      for (j = 0; j <= i; j++) {
         H[i][j] = Vector_dot(P, bv[j], size);
         Vector_subtractscaledvector(P, H[i][j], bv[j], size);
      }
      
      H[i][i+1] = Vector_norm(P, size);
      bv[i+1] = Vector_allocate(size);
      Vector_copy(bv[i+1], P, size);
      Vector_scale(bv[i+1], 1.0 / H[i][i+1], size);
      
      for (k = 0; k < i; k++)
         givensrotate(c[k], s[k], &H[i][k], &H[i][k+1]);
      givens(H[i][i], H[i][i+1], &c[i], &s[i]);
      givensrotate(c[i], s[i], &H[i][i], &H[i][i+1]);
      
      g[i+1] = 0.0;
      givensrotate(c[i], s[i], &g[i], &g[i+1]);
      
      residual = fabs(g[i+1]) / normr;

#ifdef PRINT_GMRES_RESIDUALS
      printf("Iteration: %u Residual: %g\n", i+1, residual);
#endif
      if (residual < tol)
         break;
   }

   for (k = 0; k <= i; k++)
      y[k] = g[k];
      
   for (k = i; k >= 0; k--) {
      y[k] /= H[k][k];
      for (j = k-1; j >= 0; j--)
         y[j] -= H[k][j] * y[k];
   }

   for (j = 0; j <= i; j++)
      Vector_addscaledvector(x, y[j], bv[j], size);

   Vector_copy(sol, x, size);

   Vector_free(x);
   Vector_free(r);
   Vector_free(P);
   Vector_free(c);
   Vector_free(s);
   Vector_free(g);
   Vector_free(y);
   Matrix_free(H);
   
   for (j = 0; j <= i+1; j++)
      Vector_free(bv[j]);
}

void Solv_multiply(Vector b, Tree tree, Vector x, real idiel, real odiel) {
   unsigned int realsize = tree->numpanels, i;
   
   Vector temp1 = Vector_allocate(realsize);
   Vector temp2 = Vector_allocate(realsize);
   
   Tree_multiply(temp1, tree, x+realsize, SINGLE_LAYER_INT);
   Tree_multiply(b, tree, x, DOUBLE_LAYER_INT);
   
#pragma ivdep
   for (i = 0; i < realsize; i++) {
      b[i+realsize] = b[i];
      temp2[i] = temp1[i];
   }
   
#pragma ivdep
   for (i = 0; i < realsize; i++) {
      b[i] -= temp1[i];
      b[i+realsize] = -b[i+realsize];
      temp2[i] *= (idiel / odiel);
      b[i+realsize] += temp2[i];
   }
   
   // DOUBLE LAYER CORRECTION

#pragma ivdep
   for (i = 0; i < realsize; i++) {
      b[i+realsize] += 4.0 * M_PI * x[i];// * tree->panels[i]->area;
   }
   
   Vector_free(temp1);
   Vector_free(temp2);
}

void Solv_ecf_multiply(Vector b, Tree tree, Vector x, real idiel, real odiel) {
   unsigned int size = tree->numpanels, i;

   Tree_multiply(b, tree, x, NORMDERIV_SINGLE_LAYER_INT);
   
#pragma ivdep
   for (i = 0; i < size; i++) {
      b[i] *= (odiel - idiel) / (4.0 * M_PI * idiel);
      b[i] += (1.0 + (odiel / idiel)) * 0.5 * x[i];
   }
}

void Solv_ecf_multiply_qual(Vector b, Tree tree, Vector x, real idiel, real odiel) {
   unsigned int size = tree->numpanels, i, j;

   Vector xa = Vector_allocate(size);
   Vector diag = Vector_allocate(size);

   for (i = 0; i < size; i++)
      xa[i] = tree->panels[i]->area * x[i];

   Tree_multiply_transpose(b, tree, xa, DOUBLE_LAYER_INT);

   Tree_extractdiagonal(diag, tree, DOUBLE_LAYER_INT);

#pragma ivdep
   for (i = 0; i < size; i++) {
      b[i] -= diag[i] * xa[i];  // Get rid of wrongo diagonal element
      b[i] *= 1.0 / (4.0 * M_PI * idiel);  // Scale off diagonal entries
      b[i] += (-odiel / ((odiel - idiel) * idiel) + (diag[i]) / (4.0*M_PI*idiel)) * xa[i]; // Correct diagonal
   }

   Vector_free(diag);
   Vector_free(xa);
}

void GMRES_solv(Tree tree, Preconditioner preconditioner, Vector rhs, Vector sol, real tol, real idiel, real odiel) {
   unsigned int size = tree->numpanels * 2;
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   real normr;
   unsigned int i;
   int j, k;
   real residual;

   r = Vector_allocate(size);

   Preconditioner_solve(r, preconditioner, rhs);

   normr = Vector_norm(r, size);

   x = Vector_allocate(size);

   c = Vector_allocate(MAXITERATIONS+1);
   s = Vector_allocate(MAXITERATIONS+1);
   g = Vector_allocate(MAXITERATIONS+1);
   y = Vector_allocate(MAXITERATIONS+1);
   H = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   
   P = Vector_allocate(size);

   g[0] = Vector_norm(r, size);
   bv[0] = Vector_allocate(size);
   Vector_copy(bv[0], r, size);
   Vector_scale(bv[0], 1.0 / g[0], size);

   for (i = 0; i < MAXITERATIONS; i++) {
      Solv_multiply(P, tree, bv[i], idiel, odiel);
      Preconditioner_solve(P, preconditioner, P);
      
      for (j = 0; j <= i; j++) {
         H[i][j] = Vector_dot(P, bv[j], size);
         Vector_subtractscaledvector(P, H[i][j], bv[j], size);
      }
      
      H[i][i+1] = Vector_norm(P, size);
      bv[i+1] = Vector_allocate(size);
      Vector_copy(bv[i+1], P, size);
      Vector_scale(bv[i+1], 1.0 / H[i][i+1], size);
      
      for (k = 0; k < i; k++)
         givensrotate(c[k], s[k], &H[i][k], &H[i][k+1]);
      givens(H[i][i], H[i][i+1], &c[i], &s[i]);
      givensrotate(c[i], s[i], &H[i][i], &H[i][i+1]);
      
      g[i+1] = 0.0;
      givensrotate(c[i], s[i], &g[i], &g[i+1]);
      
      residual = fabs(g[i+1]) / normr;

#ifdef PRINT_GMRES_RESIDUALS
      printf("Iteration: %u Residual: %g\n", i+1, residual);
#endif

      if (residual < tol)
         break;
   }

   for (k = 0; k <= i; k++)
      y[k] = g[k];
      
   for (k = i; k >= 0; k--) {
      y[k] /= H[k][k];
      for (j = k-1; j >= 0; j--)
         y[j] -= H[k][j] * y[k];
   }

   for (j = 0; j <= i; j++)
      Vector_addscaledvector(x, y[j], bv[j], size);

   Vector_copy(sol, x, size);

   Vector_free(x);
   Vector_free(r);
   Vector_free(P);
   Vector_free(c);
   Vector_free(s);
   Vector_free(g);
   Vector_free(y);
   Matrix_free(H);
   for (j = 0; j <= i+1; j++)
      Vector_free(bv[j]);
}

void GMRES_solv_ecf(Tree tree, Preconditioner preconditioner, Vector rhs, Vector sol, real tol, real idiel, real odiel) {
   unsigned int size = tree->numpanels;
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   real normr;
   unsigned int i;
   int j, k;
   real residual;

   r = Vector_allocate(size);

   Preconditioner_solve(r, preconditioner, rhs);

   normr = Vector_norm(r, size);

   x = Vector_allocate(size);

   c = Vector_allocate(MAXITERATIONS+1);
   s = Vector_allocate(MAXITERATIONS+1);
   g = Vector_allocate(MAXITERATIONS+1);
   y = Vector_allocate(MAXITERATIONS+1);
   H = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   
   P = Vector_allocate(size);

   g[0] = Vector_norm(r, size);
   bv[0] = Vector_allocate(size);
   Vector_copy(bv[0], r, size);
   Vector_scale(bv[0], 1.0 / g[0], size);

   for (i = 0; i < MAXITERATIONS; i++) {
	  Solv_ecf_multiply(P, tree, bv[i], idiel, odiel);
      Preconditioner_solve(P, preconditioner, P);
      
      for (j = 0; j <= i; j++) {
         H[i][j] = Vector_dot(P, bv[j], size);
         Vector_subtractscaledvector(P, H[i][j], bv[j], size);
      }
      
      H[i][i+1] = Vector_norm(P, size);
      bv[i+1] = Vector_allocate(size);
      Vector_copy(bv[i+1], P, size);
      Vector_scale(bv[i+1], 1.0 / H[i][i+1], size);
      
      for (k = 0; k < i; k++)
         givensrotate(c[k], s[k], &H[i][k], &H[i][k+1]);
      givens(H[i][i], H[i][i+1], &c[i], &s[i]);
      givensrotate(c[i], s[i], &H[i][i], &H[i][i+1]);
      
      g[i+1] = 0.0;
      givensrotate(c[i], s[i], &g[i], &g[i+1]);
      
      residual = fabs(g[i+1]) / normr;

#ifdef PRINT_GMRES_RESIDUALS
      printf("Iteration: %u Residual: %g\n", i+1, residual);
#endif

      if (residual < tol)
         break;
   }

   for (k = 0; k <= i; k++)
      y[k] = g[k];
      
   for (k = i; k >= 0; k--) {
      y[k] /= H[k][k];
      for (j = k-1; j >= 0; j--)
         y[j] -= H[k][j] * y[k];
   }

   for (j = 0; j <= i; j++)
      Vector_addscaledvector(x, y[j], bv[j], size);

   Vector_copy(sol, x, size);

   Vector_free(x);
   Vector_free(r);
   Vector_free(P);
   Vector_free(c);
   Vector_free(s);
   Vector_free(g);
   Vector_free(y);
   Matrix_free(H);
   for (j = 0; j <= i+1; j++)
      Vector_free(bv[j]);
}

void GMRES_solv_ecf_qual(Tree tree, Preconditioner preconditioner, Vector rhs, Vector sol, real tol, real idiel, real odiel) {
   unsigned int size = tree->numpanels;
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   real normr;
   unsigned int i;
   int j, k;
   real residual;

   r = Vector_allocate(size);

   Preconditioner_solve(r, preconditioner, rhs);

   normr = Vector_norm(r, size);

   x = Vector_allocate(size);

   c = Vector_allocate(MAXITERATIONS+1);
   s = Vector_allocate(MAXITERATIONS+1);
   g = Vector_allocate(MAXITERATIONS+1);
   y = Vector_allocate(MAXITERATIONS+1);
   H = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   
   P = Vector_allocate(size);

   g[0] = Vector_norm(r, size);
   bv[0] = Vector_allocate(size);
   Vector_copy(bv[0], r, size);
   Vector_scale(bv[0], 1.0 / g[0], size);

   for (i = 0; i < MAXITERATIONS; i++) {
      Solv_ecf_multiply_qual(P, tree, bv[i], idiel, odiel);
      Preconditioner_solve(P, preconditioner, P);
      
      for (j = 0; j <= i; j++) {
         H[i][j] = Vector_dot(P, bv[j], size);
         Vector_subtractscaledvector(P, H[i][j], bv[j], size);
      }
      
      H[i][i+1] = Vector_norm(P, size);
      bv[i+1] = Vector_allocate(size);
      Vector_copy(bv[i+1], P, size);
      Vector_scale(bv[i+1], 1.0 / H[i][i+1], size);
      
      for (k = 0; k < i; k++)
         givensrotate(c[k], s[k], &H[i][k], &H[i][k+1]);
      givens(H[i][i], H[i][i+1], &c[i], &s[i]);
      givensrotate(c[i], s[i], &H[i][i], &H[i][i+1]);
      
      g[i+1] = 0.0;
      givensrotate(c[i], s[i], &g[i], &g[i+1]);
      
      residual = fabs(g[i+1]) / normr;

#ifdef PRINT_GMRES_RESIDUALS
      printf("Iteration: %u Residual: %g\n", i+1, residual);
#endif

      if (residual < tol)
         break;
   }

   for (k = 0; k <= i; k++)
      y[k] = g[k];
      
   for (k = i; k >= 0; k--) {
      y[k] /= H[k][k];
      for (j = k-1; j >= 0; j--)
         y[j] -= H[k][j] * y[k];
   }

   for (j = 0; j <= i; j++)
      Vector_addscaledvector(x, y[j], bv[j], size);

   Vector_copy(sol, x, size);

   Vector_free(x);
   Vector_free(r);
   Vector_free(P);
   Vector_free(c);
   Vector_free(s);
   Vector_free(g);
   Vector_free(y);
   Matrix_free(H);
   for (j = 0; j <= i+1; j++)
      Vector_free(bv[j]);
}

void GMRES_SurfaceOperator(SurfaceOperator so, Preconditioner preconditioner, Vector rhs, Vector sol, unsigned int numpanels, real tol) {
   unsigned int size = numpanels;
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   real normr;
   unsigned int i;
   int j, k;
   real residual;

   r = Vector_allocate(size);

   Preconditioner_solve(r, preconditioner, rhs);

   normr = Vector_norm(r, size);

   if (normr < tol) { // hand back zero solution
      Vector_zero(sol, numpanels);
      Vector_free(r);
      return;
   }
   x = Vector_allocate(size);

   c = Vector_allocate(MAXITERATIONS+1);
   s = Vector_allocate(MAXITERATIONS+1);
   g = Vector_allocate(MAXITERATIONS+1);
   y = Vector_allocate(MAXITERATIONS+1);
   H = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   
   P = Vector_allocate(size);

   g[0] = Vector_norm(r, size);
   bv[0] = Vector_allocate(size);
   Vector_copy(bv[0], r, size);
   Vector_scale(bv[0], 1.0 / g[0], size);

   for (i = 0; i < MAXITERATIONS; i++) {
      SurfaceOperator_topmultiply(so, P, bv[i]);
      Preconditioner_solve(P, preconditioner, P);
      
      for (j = 0; j <= i; j++) {
         H[i][j] = Vector_dot(P, bv[j], size);
         Vector_subtractscaledvector(P, H[i][j], bv[j], size);
      }
      
      H[i][i+1] = Vector_norm(P, size);
      bv[i+1] = Vector_allocate(size);
      Vector_copy(bv[i+1], P, size);
      Vector_scale(bv[i+1], 1.0 / H[i][i+1], size);
      
      for (k = 0; k < i; k++)
         givensrotate(c[k], s[k], &H[i][k], &H[i][k+1]);
      givens(H[i][i], H[i][i+1], &c[i], &s[i]);
      givensrotate(c[i], s[i], &H[i][i], &H[i][i+1]);
      
      g[i+1] = 0.0;
      givensrotate(c[i], s[i], &g[i], &g[i+1]);
      
      residual = fabs(g[i+1]) / normr;

#ifdef PRINT_GMRES_RESIDUALS
      printf("Iteration: %u Residual: %g\n", i+1, residual);
#endif
      if (residual < tol)
         break;
   }

   for (k = 0; k <= i; k++)
      y[k] = g[k];

	num_GMRES_iter = i;

   for (k = i; k >= 0; k--) {
      y[k] /= H[k][k];
      for (j = k-1; j >= 0; j--)
         y[j] -= H[k][j] * y[k];
   }

   for (j = 0; j <= i; j++)
      Vector_addscaledvector(x, y[j], bv[j], size);

   Vector_copy(sol, x, size);

   Vector_free(x);
   Vector_free(r);
   Vector_free(P);
   Vector_free(c);
   Vector_free(s);
   Vector_free(g);
   Vector_free(y);
   Matrix_free(H);
   
   for (j = 0; j <= i+1; j++)
      Vector_free(bv[j]);
}

void GMRES_QualocationOperator(QualocationOperator qo, Preconditioner preconditioner, Vector rhs, Vector sol, unsigned int numpanels, real tol) {
   unsigned int size = numpanels;
	if (num_GMRES_iter > MAXITERATIONS ) {
	  num_GMRES_iter = MAXITERATIONS;
	}
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   real normr;
   unsigned int i;
   int j, k;
   real residual;
	
   r = Vector_allocate(size);

   Preconditioner_solve(r, preconditioner, rhs);

   normr = Vector_norm(r, size);

   if (normr < tol) {
      Vector_zero(sol, size);
      return;
   }

   x = Vector_allocate(size);

   c = Vector_allocate(MAXITERATIONS+1);
   s = Vector_allocate(MAXITERATIONS+1);
   g = Vector_allocate(MAXITERATIONS+1);
   y = Vector_allocate(MAXITERATIONS+1);
   H = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   
   P = Vector_allocate(size);

   g[0] = Vector_norm(r, size);
   bv[0] = Vector_allocate(size);
   Vector_copy(bv[0], r, size);
   Vector_scale(bv[0], 1.0 / g[0], size);
	printf("initial residual = %lf\n",g[0]);
	
	real saveHH; 
   for (i = 0; i < num_GMRES_iter; i++) {
      QualocationOperator_multiply(qo, P, bv[i]);
      Preconditioner_solve(P, preconditioner, P);
      
      for (j = 0; j <= i; j++) {
         H[i][j] = Vector_dot(P, bv[j], size);
         Vector_subtractscaledvector(P, H[i][j], bv[j], size);
      }
      
      H[i][i+1] = Vector_norm(P, size);
		saveHH = H[i][i+1];
      
      for (k = 0; k < i; k++)
         givensrotate(c[k], s[k], &H[i][k], &H[i][k+1]);
      givens(H[i][i], H[i][i+1], &c[i], &s[i]);
      givensrotate(c[i], s[i], &H[i][i], &H[i][i+1]);
      
      g[i+1] = 0.0;
      givensrotate(c[i], s[i], &g[i], &g[i+1]);
      
      residual = fabs(g[i+1]) / normr;

#ifdef PRINT_GMRES_RESIDUALS
      printf("Iteration: %u Residual: %g\n", i+1, residual);
#endif
      if (residual < tol)
         break;

      bv[i+1] = Vector_allocate(size);
      Vector_copy(bv[i+1], P, size);
      Vector_scale(bv[i+1], 1.0 / saveHH, size);

   }

	if (! (residual < tol) ) {
	  // we hit the limit on num iters
      for (j = 0; j <= i; j++) {
         H[i][j] = Vector_dot(P, bv[j], size);
         Vector_subtractscaledvector(P, H[i][j], bv[j], size);
      }
      
      H[i][i+1] = Vector_norm(P, size);
	}
   for (k = 0; k <= i; k++) {
      y[k] = g[k];
	}

	num_GMRES_iter = i;
	
   for (k = i; k >= 0; k--) {
      y[k] /= H[k][k];
      for (j = k-1; j >= 0; j--)
         y[j] -= H[k][j] * y[k];
   }

   for (j = 0; j <= i; j++)
      Vector_addscaledvector(x, y[j], bv[j], size);

   Vector_copy(sol, x, size);

   Vector_free(x);
   Vector_free(r);
   Vector_free(P);
   Vector_free(c);
   Vector_free(s);
   Vector_free(g);
   Vector_free(y);
   Matrix_free(H);
   
   for (j = 0; j < i+1; j++) // used to be j<=i+1 but that always fires terrible warnings
      Vector_free(bv[j]);

	// funny memory leak problem can happen... but doesn't always. never figured it out entirely.
	// see the allocation of bv[i+1] in the main Krylov loop.  compiler
	// optimization options screw everything up.  bv[i+1] never is an actual Vector in gdb.
}
