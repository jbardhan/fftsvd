#include "FFTSVD.h"

/*#define CONVERSION 332.06*/
#define CONVERSION 332.112

int main(int argc, char* argv[]) {
   FILE* gstfile = NULL;
   FILE* torfile = NULL;
   FILE* chargefile = NULL;
   Charge charge;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   GST* gstpanels;
   TOR* torpanels;
   Panel* panels;
   unsigned int numgstpanels = 0, numtorpanels = 0, numpanels = 0;
   Vector3D* centroids;
   unsigned int i, c;
   Vector rhs;
   Vector sol;
   Tree tree;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0, idiel, odiel;

   if (argc != 10) {
      printf("Usage: %s [surface.gst] [surface.tor] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels] [idiel] [odiel]\n", argv[0]);
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

   numpanels = numgstpanels + numtorpanels;
                                                                                
   printf("NUMPANELS: %u\n", numpanels);

   charge = Charge_allocate();
   
   chargefile = fopen(argv[3], "r");
   
   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2;
   }

   Charge_read(charge, chargefile);
   
   fclose(chargefile);

   svdtol = atof(argv[4]);
   gmrestol = atof(argv[5]);
   gridpoints = atoi(argv[6]);
   maxpanelsperfinestcube = atoi(argv[7]);
   idiel = atof(argv[8]);
   odiel = atof(argv[9]);

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

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

   rhs = Vector_allocate(numpanels);

   Charge_makerhs_ecf_qual(rhs, charge, panels, numpanels, idiel, odiel);

   sol = Vector_allocate(numpanels);

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, centroids, numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, DOUBLE_LAYER_INT, 0.0);

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

   printf("Constructing preconditioner... ");
   fflush(stdout);
   
   preconditioner = Preconditioner_allocate(numpanels, numpanels);

   Preconditioner_fill_diagonal_solv_ecf_qual(preconditioner, tree, idiel, odiel);

   //Preconditioner_fill_identity(preconditioner);
      
   Preconditioner_factor(preconditioner);

   printf("done.\n");

   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

   GMRES_solv_ecf_qual(tree, preconditioner, rhs, sol, gmrestol, idiel, odiel);
   
   phi = Vector_allocate(charge->numcharges);

   for (i = 0; i < numpanels; i++) {
      //printf("%f\n", rhs[i]);
      for (c = 0; c < charge->numcharges; c++) {
         real slp, dlp;
         
         slp = Integration(charge->points[c], panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
         
         phi[c] += CONVERSION * (sol[i] * slp) / idiel;
      }
   }

   for (c = 0; c < charge->numcharges; c++)
      energy += 0.5 * charge->charges[c] * phi[c];
   
   printf("Solvation Energy: %f kcal/mol\n", energy);


     {
     unsigned int i, j;
     Vector x = Vector_allocate(tree->numpanels);
     Vector ans = Vector_allocate(tree->numpoints);
     FILE* file = NULL;
                                                                                
     file = fopen("qual.m", "w");
                                                                                
     fprintf(file, "A = zeros(%u, %u);\n", tree->numpanels, tree->numpoints);
                                                                                
     fprintf(file, "A = [\n");

     for (i = 0; i < tree->numpanels; i++) {
     Vector_zero(x, tree->numpanels);
     x[i] = 1.0;
     Solv_ecf_multiply_qual(ans, tree, x, idiel, odiel);
     for (j = 0; j < tree->numpoints; j++)
     fprintf(file, "%f ", ans[j]);
     fprintf(file, "\n");
     }
                                                                                
     fprintf(file, "]';\n");
                                                                                
     Vector_free(x);
     Vector_free(ans);
     }   


   Preconditioner_free(preconditioner);
   Tree_free(tree);
   Vector_free(sol);
   Vector_free(rhs);
   Vector_free(phi);
   for (i = 0; i < numpanels; i++)
      Vector3D_free(centroids[i]);
   free(centroids);
   for (i = 0; i < numpanels; i++)
      Panel_free(panels[i]);
   free(panels);
   Charge_free(charge);

   return 0;
}
