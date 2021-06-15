#include "FFTSVD.h"
 
/* Constructors and Destructors */

QualocationOperator QualocationOperator_allocate() {
   QualocationOperator qo = (QualocationOperator)calloc(1, sizeof(_QualocationOperator));

   return qo;
}

void QualocationOperator_free(QualocationOperator qo) {
   unsigned int i;

   if (qo->tree)
      Tree_free(qo->tree);
   if (qo->M1M3)
      Tree_free(qo->M1M3);

   if (qo->charges)
      Charge_free(qo->charges);

   free(qo);
}

void QualocationOperator_multiply(QualocationOperator qo, Vector dest, Vector src) {
   Solv_ecf_multiply_qual(dest, qo->tree, src, qo->innerdielectric, qo->outerdielectric);
}

void QualocationOperator_generatePreconditioner(QualocationOperator qo, Preconditioner P) {
   Preconditioner_fill_diagonal_solv_ecf_qual(P, qo->tree, qo->innerdielectric, qo->outerdielectric);
}

void QualocationOperator_makeRHS(QualocationOperator qo, Vector RHS) {
  //  if (0) {
  if (qo->M1M3) {
      unsigned int i;

      Tree_multiply_transpose(RHS, qo->M1M3, qo->charges->charges, DOUBLE_LAYER_INT);
      
      for (i = 0; i < qo->M1M3->numpanels; i++)
         RHS[i] *= -1.0 / (4.0 * M_PI * qo->innerdielectric);

/*       Vector RHS2 = Vector_allocate(qo->M1M3->numpanels); */

/*       Charge_makerhs_ecf_qual(RHS2, qo->charges, qo->tree->panels, qo->tree->numpanels, qo->innerdielectric, qo->outerdielectric); */

/*       for (i = 0; i < qo->M1M3->numpanels; i++) */
/*          printf("%f %f %g\n", RHS[i], RHS2[i], fabs(RHS[i]-RHS2[i])/fabs(RHS2[i])); */

/*       Vector_free(RHS2); */
   }
   else
      Charge_makerhs_ecf_qual(RHS, qo->charges, qo->tree->panels, qo->tree->numpanels, qo->innerdielectric, qo->outerdielectric);
}

void QualocationOperator_calculateReactionPotentials(QualocationOperator qo,
                                                 Vector surfacePotentials,
                                                 Vector pointPotentials) {
   unsigned int i, c;

   Vector_zero(pointPotentials, qo->charges->numcharges);

   if (qo->M1M3)
      Tree_multiply(pointPotentials, qo->M1M3, surfacePotentials, SINGLE_LAYER_INT);
   else {
      for (i = 0; i < qo->tree->numpanels; i++)
         for (c = 0; c < qo->charges->numcharges; c++) {
            real slp = Integration(qo->charges->points[c], qo->tree->panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
            pointPotentials[c] += surfacePotentials[i] * slp;
         }
   }
}

void QualocationOperator_writematlabfile(char* filename, QualocationOperator qo, unsigned int totalnumpanels) {
   unsigned int i, j;
   Vector x = Vector_allocate(totalnumpanels);
   Vector ans = Vector_allocate(totalnumpanels);
   FILE* file = NULL;

   file = fopen(filename, "w");

   //fprintf(file, "A = zeros(%u, %u);\n", 2 * totalnumpanels, 2 * totalnumpanels);

   //fprintf(file, "A = [\n");

   for (i = 0; i < totalnumpanels; i++) {
      Vector_zero(x, totalnumpanels);
      x[i] = 1.0;
      QualocationOperator_multiply(qo, ans, x);
      for (j = 0; j < totalnumpanels; j++)
         fprintf(file, "%f ", ans[j]);
      fprintf(file, "\n");
   }

   //fprintf(file, "]';\n");

   Vector_free(x);
   Vector_free(ans);
}

void QualocationOperator_makeRHS_fromVector(QualocationOperator qo, Vector RHS,
                                        Vector srcCharges, unsigned int totalnumcharges) {
   unsigned int i;

   for (i = 0; i < totalnumcharges; i++) {
	  qo->charges->charges[i] = srcCharges[i];
	}
   QualocationOperator_makeRHS(qo, RHS);
}

void QualocationOperator_solve(QualocationOperator qo, Vector qol, Vector RHS) {  
}

void QualocationOperator_collectPotentials_toVector(QualocationOperator qo, Vector phi_r, Vector qol) {
   QualocationOperator_calculateReactionPotentials(qo, qol, phi_r);
}
