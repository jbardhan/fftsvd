#include "FFTSVD.h"

/* Constructors and Destructors */

Preconditioner Preconditioner_allocate(unsigned int rows, unsigned int columns) {
   Preconditioner preconditioner = (Preconditioner)calloc(1, sizeof(_Preconditioner));
   
   preconditioner->rows = rows;
   preconditioner->columns = columns;
   preconditioner->P = (PreconditionerElement**)calloc(preconditioner->columns, sizeof(PreconditionerElement*));
   preconditioner->numelements = (unsigned int*)calloc(preconditioner->columns, sizeof(unsigned int));
   preconditioner->perm_r = (int*)calloc(preconditioner->rows, sizeof(int));
   preconditioner->perm_c = (int*)calloc(preconditioner->columns, sizeof(int));
   preconditioner->nnz = 0;
   
   StatInit(&preconditioner->stat);
   
   return preconditioner;
}

void Preconditioner_fill_identity(Preconditioner preconditioner) {
   unsigned int i;
   
   for (i = 0; i < preconditioner->columns; i++)
      Preconditioner_set(preconditioner, i, i, 1.0);
}

void Preconditioner_fill_diagonal_cap(Preconditioner preconditioner, Tree tree) {
   unsigned int i;

   Vector diag = Vector_allocate(tree->numpanels);
   Tree_extractdiagonal(diag, tree, SINGLE_LAYER_INT);
   
   for (i = 0; i < tree->numpanels; i++)
      Preconditioner_set(preconditioner, i, i, diag[i]);

   Vector_free(diag);
}

void Preconditioner_fill_diagonal_solv_juffersimple(Preconditioner preconditioner, Tree tree, real idiel, real odiel) {
  unsigned int i;

  for (i = 0; i < tree->numpanels; i++) {
	 Preconditioner_set(preconditioner, i, i, 0.5 * (1. + odiel / idiel));
  }
  
}

void Preconditioner_fill_diagonal_solv_ecf_qual(Preconditioner preconditioner, Tree tree, real idiel, real odiel) {
   unsigned int i;

   Vector diag = Vector_allocate(tree->numpanels);

   Tree_extractdiagonal(diag, tree, DOUBLE_LAYER_INT);

   for (i = 0; i < tree->numpanels; i++)
      Preconditioner_set(preconditioner, i, i, (-odiel / ((odiel - idiel) * idiel) + (diag[i]) / (4.0*M_PI*idiel)) * tree->panels[i]->area);

   Vector_free(diag);
}

void Preconditioner_fill_diagonal_solv_RUB(Preconditioner preconditioner, Tree tree, real idiel, real odiel) {
   unsigned int i;

   Vector diag = Vector_allocate(tree->numpanels);

   Tree_extractdiagonal(diag, tree, DOUBLE_LAYER_INT);
	real testval, val;

   for (i = 0; i < tree->numpanels; i++) {
	  testval = (-1.0 / ((odiel - idiel) * idiel)) * 1.0 + 1.0 / (2.0 * idiel);
	  val = -odiel / ((odiel - idiel) * idiel) + 2.0*M_PI / (4.0 * M_PI * idiel);
	  //	  printf("diag[%d] = %lf val = %lf  testval = %lf\n", i, diag[i], val, testval);
	  
	  Preconditioner_set(preconditioner, i, i, (-1./ ((odiel - idiel) * idiel)) * tree->panels[i]->area);
	}
	
   Vector_free(diag);
}
void Preconditioner_fill_diagonal_solv_born(Preconditioner preconditioner, Tree tree, real idiel, real odiel) {
   unsigned int i;

   Vector diag = Vector_allocate(tree->numpanels);

   Tree_extractdiagonal(diag, tree, DOUBLE_LAYER_INT);

   for (i = 0; i < tree->numpanels; i++)
	  Preconditioner_set(preconditioner, i, i, (-odiel / ((odiel - idiel) * idiel)/*goes here*/) * tree->panels[i]->area);
	// snipped out  "+ (diag[i]) / (4.0*M_PI*idiel)"
	// there's a double negative in this part of the ecf_qual code!! it matters for curved panels.
	
   Vector_free(diag);
}

void Preconditioner_fill_blockdiagonal_cap(Preconditioner preconditioner, Cube cube, BEMKernelType kerneltype, void* parameters) {
   unsigned int i, j, cx, cy, cz;

   if (cube->leaf) {
      for (i = 0; i < cube->numpointindices; i++)
         for (j = 0; j < cube->numpanelindices; j++)
            Preconditioner_set(preconditioner, cube->pointindices[i], cube->panelindices[j], cube->D_single[i][j]);
   }

   for (cx = 0; cx <= 1; cx++)
      for (cy = 0; cy <= 1; cy++)
         for (cz = 0; cz <= 1; cz++)
            if (cube->children[cx][cy][cz] != NULL)
               Preconditioner_fill_blockdiagonal_cap(preconditioner, cube->children[cx][cy][cz], kerneltype, parameters);
}

void Preconditioner_fill_blockdiagonal_solv(Preconditioner preconditioner, Tree tree, real idiel, real odiel) {
   unsigned int i;

   Vector diag = Vector_allocate(tree->numpanels);
   Tree_extractdiagonal(diag, tree, DOUBLE_LAYER_INT);
   
   for (i = 0; i < tree->numpanels; i++) {
      Preconditioner_set(preconditioner, i, i, diag[i]);
      Preconditioner_set(preconditioner, i + tree->numpanels, i, 4.0*M_PI-diag[i]);
   }

   Tree_extractdiagonal(diag, tree, SINGLE_LAYER_INT);

   for (i = tree->numpanels; i < tree->numpanels * 2; i++) {
      Preconditioner_set(preconditioner, i, i, (idiel / odiel) * diag[i-tree->numpanels]);
      Preconditioner_set(preconditioner, i - tree->numpanels, i, -diag[i-tree->numpanels]);
   }

   Vector_free(diag);
}

void Preconditioner_fill_blockdiagonal_solv_cavity(Preconditioner preconditioner,
                                                   Panel* panels1, unsigned int numpanels1,
                                                   Panel* panels2, unsigned int numpanels2,
                                                   BEMKernelType kerneltype, void* parameters) {
   unsigned int i, base;
   real slp, dlp;
   for (i = 0; i < numpanels1; i++) {
      //integrator(panels1[i]->centroid, panels1[i], parameters, &slp, &dlp);
      Preconditioner_set(preconditioner, i, i, 2.0 * M_PI);
      Preconditioner_set(preconditioner, i + numpanels1, i, 2.0 * M_PI);
   }
   for (i = numpanels1; i < numpanels1 * 2; i++) {
      slp = Integration(panels1[i-numpanels1]->centroid, panels1[i-numpanels1], kerneltype, parameters, SINGLE_LAYER_INT);
      Preconditioner_set(preconditioner, i, i, 20 * slp);
      Preconditioner_set(preconditioner, i - numpanels1, i, -slp);
   }   

   base = 2*numpanels1;
   for (i = 0; i < numpanels2; i++) {
      //integrator(panels2[i]->centroid, panels2[i], parameters, &slp, &dlp);
      Preconditioner_set(preconditioner, i+base, i+base, 2.0 * M_PI);
      Preconditioner_set(preconditioner, i+base + numpanels2, i+base, 2.0 * M_PI);
   }
   for (i = numpanels2; i < numpanels2 * 2; i++) {
      slp = Integration(panels2[i-numpanels2]->centroid, panels2[i-numpanels2], kerneltype, parameters, SINGLE_LAYER_INT);
      Preconditioner_set(preconditioner, i+base, i+base, 0.05 * slp);
      Preconditioner_set(preconditioner, i - numpanels2+base, i+base, -slp);
   }
}

void Preconditioner_free(Preconditioner preconditioner) {
   unsigned int i;
   
   for (i = 0; i < preconditioner->columns; i++)
      free(preconditioner->P[i]);
   free(preconditioner->P);
   free(preconditioner->numelements);
   Destroy_SuperNode_Matrix(&preconditioner->L);
   Destroy_CompCol_Matrix(&preconditioner->U);
   free(preconditioner->perm_r);
   free(preconditioner->perm_c);
   StatFree(&preconditioner->stat);
   
   free(preconditioner);
}

/* Operations */

void Preconditioner_set(Preconditioner preconditioner, unsigned int row, unsigned int column, real value) {
   /* If this column is empty, allocate the first one */
   if (preconditioner->numelements[column] == 0) {
      preconditioner->P[column] = (PreconditionerElement*)calloc(1, sizeof(PreconditionerElement));
   }
   /* Otherwise, increase the storage by one element, copy, remap, and set the
      value */
   else if ((preconditioner->numelements[column] > 0) && ((preconditioner->numelements[column] & (preconditioner->numelements[column] - 1)) == 0)) {
      PreconditionerElement* temp = (PreconditionerElement*)calloc(preconditioner->numelements[column] * 2, sizeof(PreconditionerElement));
      memcpy(temp, preconditioner->P[column], preconditioner->numelements[column] * sizeof(PreconditionerElement));
      free(preconditioner->P[column]);
      preconditioner->P[column] = temp;
   }

   preconditioner->P[column][preconditioner->numelements[column]].row = row;
   preconditioner->P[column][preconditioner->numelements[column]].value = value;
   preconditioner->numelements[column]++;
   
   preconditioner->nnz++;
}

void Preconditioner_factor(Preconditioner preconditioner) {
   SuperMatrix P, PC;
   int* etree, info;
   superlu_options_t options;   
   int* row_ind;
   int* col_ptr;
   real* val;
   unsigned int r, c, elementcount = 0;
   
   /* First, map preconditioner->P and preconditioner->numelements into
      column compressed format */
   
   col_ptr = (int*)calloc(preconditioner->columns + 1, sizeof(int));
   row_ind = (int*)calloc(preconditioner->nnz, sizeof(int));
   val = (real*)calloc(preconditioner->nnz, sizeof(real));
   
   for (c = 0; c < preconditioner->columns; c++) {
      col_ptr[c] = elementcount;
      for (r = 0; r < preconditioner->numelements[c]; r++) {
         row_ind[elementcount] = preconditioner->P[c][r].row;
         val[elementcount] = preconditioner->P[c][r].value;
         elementcount++;
      }
   }
   col_ptr[preconditioner->columns] = elementcount;
   
   /* Call the SuperLU routine to form the compressed column preconditioner
      matrix */
   
#ifdef REAL_IS_FLOAT
   sCreate_CompCol_Matrix(&P, preconditioner->rows, preconditioner->columns, preconditioner->nnz, val, row_ind, col_ptr, SLU_NC, SLU_S, SLU_GE);
#else
#ifdef REAL_IS_DOUBLE
   dCreate_CompCol_Matrix(&P, preconditioner->rows, preconditioner->columns, preconditioner->nnz, val, row_ind, col_ptr, SLU_NC, SLU_D, SLU_GE);
#endif
#endif

   /* Determine the permutation and the elimination tree */

   get_perm_c(2, &P, preconditioner->perm_c);
   
   etree = (int*)calloc(preconditioner->columns, sizeof(int));
   
   set_default_options(&options);

   sp_preorder(&options, &P, preconditioner->perm_c, etree, &PC);

   /* Perform the LU factorization */
   GlobalLU_t Glu;
   
#ifdef REAL_IS_FLOAT
   sgstrf(&options, &PC, sp_ienv(2), sp_ienv(1), etree, NULL, 0, preconditioner->perm_c, preconditioner->perm_r, &preconditioner->L, &preconditioner->U, &Glu, &preconditioner->stat, &info);
#else
#ifdef REAL_IS_DOUBLE
   dgstrf(&options, &PC, sp_ienv(2), sp_ienv(1), etree, NULL, 0, preconditioner->perm_c, preconditioner->perm_r, &preconditioner->L, &preconditioner->U, &Glu, &preconditioner->stat, &info);
#endif
#endif

   if (info) {
      printf("SuperLU factorization error: %d\n", info);
      exit(-999);
   }

   Destroy_SuperMatrix_Store(&P);
   Destroy_CompCol_Permuted(&PC);   
   free(etree);
   free(col_ptr);
   free(row_ind);
   free(val);
}

void Preconditioner_solve(Vector x, Preconditioner preconditioner, Vector b) {
   SuperMatrix X;
   int info;

   Vector bcopy = Vector_allocate(preconditioner->rows);
   
   Vector_copy(bcopy, b, preconditioner->rows);

#ifdef REAL_IS_FLOAT
   sCreate_Dense_Matrix(&X, preconditioner->rows, 1, bcopy, preconditioner->rows, SLU_DN, SLU_S, SLU_GE);
   sgstrs(NOTRANS, &preconditioner->L, &preconditioner->U, preconditioner->perm_c, preconditioner->perm_r, &X, &preconditioner->stat, &info);
#else
#ifdef REAL_IS_DOUBLE
   dCreate_Dense_Matrix(&X, preconditioner->rows, 1, bcopy, preconditioner->rows, SLU_DN, SLU_D, SLU_GE);
   dgstrs(NOTRANS, &preconditioner->L, &preconditioner->U, preconditioner->perm_c, preconditioner->perm_r, &X, &preconditioner->stat, &info);
#endif
#endif

   Vector_copy(x, (Vector)((DNformat*)X.Store)->nzval, preconditioner->columns);

   Vector_free(bcopy);
   Destroy_SuperMatrix_Store(&X);
}

unsigned int Preconditioner_memory(Preconditioner preconditioner) {
   mem_usage_t mem_usage;
   unsigned int Pmem;

#ifdef REAL_IS_FLOAT
   sQuerySpace(&preconditioner->L, &preconditioner->U, &mem_usage);
#else
#ifdef REAL_IS_DOUBLE
   dQuerySpace(&preconditioner->L, &preconditioner->U, &mem_usage);
#endif
#endif

   Pmem = sizeof(struct _Preconditioner) +
      preconditioner->columns * sizeof(PreconditionerElement*) + /* *P */
      preconditioner->nnz * sizeof(PreconditionerElement) + /* **P */
      preconditioner->columns * sizeof(unsigned int) + /* numelements */
      2 * preconditioner->columns * sizeof(int) + /* perm_r/c */
      mem_usage.for_lu;
                       
   return Pmem;
}

void Preconditioner_writematlabfile(char* filename, Preconditioner preconditioner) {
   unsigned int c, i;
   FILE* file = NULL;

   file = fopen(filename, "w");

   fprintf(file, "Pdata = [\n");

   for (c = 0; c < preconditioner->columns; c++)
      for (i = 0; i < preconditioner->numelements[c]; i++)
         fprintf(file, "%u %u %f\n", preconditioner->P[c][i].row, c, preconditioner->P[c][i].value);

   fprintf(file, "];\n");

   fprintf(file, "P = spalloc(%u, %u, %u);\n", preconditioner->rows, preconditioner->columns, preconditioner->nnz);
   fprintf(file, "for i = 1:%u\n", preconditioner->nnz);
   fprintf(file, "P(Pdata(i,1)+1,Pdata(i,2)+1) = Pdata(i,3);\n");
   fprintf(file, "end\n");

   fclose(file);
}

