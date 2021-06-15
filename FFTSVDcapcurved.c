#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* gstfile = NULL;
   FILE* torfile = NULL;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   GST* gstpanels;
   TOR* torpanels;
   Panel* panels;
   unsigned int numgstpanels, numtorpanels, numpanels;
   Vector3D* centroids;
   unsigned int i, j, k, k2, index;
   Vector rhs;
   Vector sol;
   Tree tree;
   Preconditioner preconditioner;
   real capacitance = 0.0;

   if (argc != 7) {
      printf("Usage: %s [gst] [tor] [svd tol] [gmres tol] [gridpoints] [maxpanels]\n", argv[0]);
      return -1;
   }

   gstfile = fopen(argv[1], "r");

   if (gstfile == NULL) {
      perror("Error opening GST file");
      return -2;
   }

   GST_readfile(&numgstpanels, &gstpanels, gstfile, 0);

   fclose(gstfile);

   torfile = fopen(argv[2], "r");

   if (torfile == NULL) {
      perror("Error opening TOR file");
      return -2;
   }

   TOR_readfile(&numtorpanels, &torpanels, torfile, 0);

   fclose(torfile);

   svdtol = atof(argv[3]);
   gmrestol = atof(argv[4]);
   gridpoints = atoi(argv[5]);
   maxpanelsperfinestcube = atoi(argv[6]);

   numpanels = numgstpanels + numtorpanels;

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("NUMPANELS: %u\n", numpanels);

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   for (i = 0; i < numgstpanels; i++) {
      panels[i] = Panel_allocate();
      Panel_GST(panels[i], gstpanels[i], 0);
   }
   for (i = 0; i < numtorpanels; i++) {
      panels[i+numgstpanels] = Panel_allocate();
      Panel_TOR(panels[i+numgstpanels], torpanels[i], 0);
   }

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
   
   rhs = Vector_allocate(numpanels);
   sol = Vector_allocate(numpanels);
     
   for (i = 0; i < numpanels; i++)
      rhs[i] = 1.0;
     
   GMRES_cap(tree, preconditioner, rhs, sol, gmrestol);

   real totalarea = 0.0;

   for (i = 0; i < numpanels; i++) {
      capacitance  += sol[i] * panels[i]->area;
      totalarea += panels[i]->area;
   }

   printf("Surface Area: %f\n", totalarea);

   capacitance *= 4.0 * M_PI * 8.854187818E-12 * 1e9;

   printf("Capacitance: %f\n", capacitance);

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

   return 0;
}
