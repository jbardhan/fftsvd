#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* quifile = NULL;
   QUI qui;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels;
   unsigned int numpanels, numconductors;
   Vector3D* centroids;
   unsigned int i, j, k, k2, index;
   Vector rhs;
   Vector sol;
   Tree tree;
   Preconditioner preconditioner;
   real *capacitance;

   if (argc != 6) {
      printf("Usage: %s [surface.qui] [svd tol] [gmres tol] [gridpoints] [maxpanels]\n", argv[0]);
      return -1;
   }

   qui = QUI_allocate();

   quifile = fopen(argv[1], "r");

   if (quifile == NULL) {
      perror("Error opening QUI file");
      return -2;
   }

   QUI_read(qui, quifile);

   fclose(quifile);

   svdtol = atof(argv[2]);
   gmrestol = atof(argv[3]);
   gridpoints = atoi(argv[4]);
   maxpanelsperfinestcube = atoi(argv[5]);

   numpanels = qui->numpanels;
   numconductors = qui->numconductors;
   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   QUI_getpanels(qui, panels);

   printf("done.\n");

   centroids = (Vector3D*)calloc(numpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels; i++) {
      centroids[i] = Vector3D_allocate();
      Vector3D_copy(centroids[i], panels[i]->centroid);
   }

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, centroids, numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_LAYER_INT, 0.0);

   printf("done.\n");

   printf("Determining local and interacting lists... ");
   fflush(stdout);

   Tree_lists(tree);

   printf("done.\n");

   printf("Filling tree structure... ");
   fflush(stdout);

   Tree_fill(tree);

   printf("done.\n");

   Tree_memory(tree);

   printf("Creating preconditioner... ");
   fflush(stdout);

   preconditioner = Preconditioner_allocate(numpanels, numpanels);
   //Preconditioner_fill_identity(preconditioner);
   Preconditioner_fill_diagonal_cap(preconditioner, tree);
   Preconditioner_factor(preconditioner);

   printf("done.\n");
 
   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));
   
#ifdef REAL_IS_DOUBLE
   capacitance = (double *)calloc(numconductors * numconductors, sizeof(double));
#else
#ifdef REAL_IS_FLOAT
   capacitance = (float *)calloc(numconductors * numconductors, sizeof(float));
#endif
#endif
   rhs = Vector_allocate(numpanels);
   sol = Vector_allocate(numpanels);
     
   for (k = 0; k < numconductors; k++) {
      Vector_zero(rhs, numpanels);

      for (i = 0; i < qui->conductorpanelcounts[k]; i++)
         rhs[qui->conductorpanelindices[k][i]] = 1.0;

      Vector_zero(sol, numpanels);
     
      GMRES_cap(tree, preconditioner, rhs, sol, gmrestol);

      for (k2 = 0; k2 < numconductors; k2++) {
         for (i = 0; i < qui->conductorpanelcounts[k2]; i++) {
            index = qui->conductorpanelindices[k2][i];
            capacitance[k2 * numconductors + k] += sol[index] * panels[index]->area;
         }
         capacitance[k2 * numconductors + k] *= 4.0 * M_PI * 8.854187818E-12 * 1e9;
      }
   }
   printf("Capacitance matrix (nF):\n");
   for (i = 0; i < numconductors; i++) {
      for (j = 0; j < numconductors; j++) {
#ifdef REAL_IS_DOUBLE
         printf("%.10lf   ",capacitance[i * numconductors + j]);
#else
#ifdef REAL_IS_FLOAT
         printf("%.10f   ",capacitance[i * numconductors + j]);
#endif
#endif
      }
      printf("\n");
   }
   free(capacitance);
   Preconditioner_free(preconditioner);
   Tree_free(tree);
   Vector_free(sol);
   Vector_free(rhs);
   for (i = 0; i < numpanels; i++)
      Vector3D_free(centroids[i]);
   free(centroids);
   for (i = 0; i < numpanels; i++)
      Panel_free(panels[i]);
   free(panels);
   QUI_free(qui);

   return 0;
}
