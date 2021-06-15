#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* chargefile = NULL;
   VertFace vertface_d, vertface_c;
   Charge charge;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels_d;
   Panel* panels_c;
   Panel* totalpanels;
   unsigned int numpanels_d;
   unsigned int numpanels_c;
   unsigned int totalnumpanels;
   Vector3D* centroids_d;
   Vector3D* centroids_c;
   Vector3D* totalcentroids;
   unsigned int i, c;
   Vector rhs;
   Vector sol;
   SurfaceOperator sopo, sopi, sopc;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0;
   if (argc != 10) {
      printf("Usage: %s [surface.vert] [surface.face] [cavity.vert] [cavity.face] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels]\n", argv[0]);
      return -1;
   }

   vertface_d = VertFace_allocate();

   vertfile = fopen(argv[1], "r");

   if (vertfile == NULL) {
      perror("Error opening vertices file");
      return -2;
   }

   VertFace_readvert(vertface_d, vertfile);

   fclose(vertfile);

   facefile = fopen(argv[2], "r");
    
   if (facefile == NULL) {
      perror("Error opening face file");
      return -2;
   }
    
   VertFace_readface_flip(vertface_d, facefile);
    
   fclose(facefile);

   vertface_c = VertFace_allocate();

   vertfile = fopen(argv[3], "r");

   if (vertfile == NULL) {
      perror("Error opening vertices file");
      return -2;
   }

   VertFace_readvert(vertface_c, vertfile);

   fclose(vertfile);

   facefile = fopen(argv[4], "r");
    
   if (facefile == NULL) {
      perror("Error opening face file");
      return -2;
   }
    
   VertFace_readface_flip(vertface_c, facefile);
    
   fclose(facefile);

   charge = Charge_allocate();
   
   chargefile = fopen(argv[5], "r");
   
   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2;
   }

   Charge_read(charge, chargefile);
   charge->globalindexstart = 0;
   fclose(chargefile);

   svdtol = atof(argv[6]);
   gmrestol = atof(argv[7]);
   gridpoints = atoi(argv[8]);
   maxpanelsperfinestcube = atoi(argv[9]);

   numpanels_d = vertface_d->numfaces;
   numpanels_c = vertface_c->numfaces;

   totalnumpanels = numpanels_d + numpanels_c;

   panels_d = (Panel*)calloc(numpanels_d, sizeof(Panel));
   panels_c = (Panel*)calloc(numpanels_c, sizeof(Panel));
   totalpanels = (Panel*)calloc(totalnumpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface_d, panels_d);
   VertFace_getpanels(vertface_c, panels_c);

   for (i = 0; i < numpanels_d; i++)
      totalpanels[i] = panels_d[i];
   for (i = 0; i < numpanels_c; i++)
      totalpanels[i+numpanels_d] = panels_c[i];

   printf("done.\n");

   centroids_d = (Vector3D*)calloc(numpanels_d, sizeof(Vector3D));
   centroids_c = (Vector3D*)calloc(numpanels_c, sizeof(Vector3D));
   totalcentroids = (Vector3D*)calloc(totalnumpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels_d; i++) {
      centroids_d[i] = Vector3D_allocate();
      Vector3D_copy(centroids_d[i], panels_d[i]->centroid);
      totalcentroids[i] = centroids_d[i];
   }
   for (i = 0; i < numpanels_c; i++) {
      centroids_c[i] = Vector3D_allocate();
      Vector3D_copy(centroids_c[i], panels_c[i]->centroid);
      totalcentroids[i+numpanels_d] = centroids_c[i];
   }

   rhs = Vector_allocate(2 * totalnumpanels);
   sol = Vector_allocate(2 * totalnumpanels);

   printf("Constructing and allocating surface operator... ");
   fflush(stdout);

   // Surface Operator Poisson Outside

   sopo = SurfaceOperator_allocate();

   sopo->kernel = POISSON_KERNEL;
   sopo->epsilon = 80.0;
   sopo->kappa = 0.0;
   sopo->mynumpanels = 0;
   sopo->tree = Tree_allocate(panels_d, numpanels_d, centroids_d, numpanels_d, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   sopo->numchildren = 1;
   sopo->children = (SurfaceOperator*)calloc(1, sizeof(SurfaceOperator));
   sopo->resultInternal = Vector_allocate(numpanels_d);
   sopo->resultExternal = Vector_allocate(numpanels_d);
   sopo->parent = NULL;
   sopo->charges = NULL;

   // Surface Operator Poisson Inside

   sopi = SurfaceOperator_allocate();

   sopi->kernel = POISSON_KERNEL;
   sopi->epsilon = 4.0;
   sopi->kappa = 0.0;
   sopi->mynumpanels = numpanels_d;
   sopi->tree = Tree_allocate(totalpanels, totalnumpanels, totalcentroids, totalnumpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   sopi->numchildren = 1;
   sopi->children = (SurfaceOperator*)calloc(1, sizeof(SurfaceOperator));
   sopi->resultInternal = Vector_allocate(totalnumpanels);
   sopi->resultExternal = Vector_allocate(totalnumpanels);
   sopi->parent = sopo;
   sopi->charges = charge;
   sopo->children[0] = sopi;

   // Surface Operator Poisson Cavity

   sopc = SurfaceOperator_allocate();

   sopc->kernel = POISSON_KERNEL;
   sopc->epsilon = 80.0;
   sopc->kappa = 0.0;
   sopc->mynumpanels = numpanels_c;
   sopc->tree = Tree_allocate(panels_c, numpanels_c, centroids_c, numpanels_c, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   sopc->numchildren = 0;
   sopc->children = NULL;
   sopc->resultInternal = Vector_allocate(numpanels_c);
   sopc->resultExternal = Vector_allocate(numpanels_c);
   sopc->parent = sopi;
   sopc->charges = NULL;
   sopi->children[0] = sopc;

   unsigned int index = 0;

   SurfaceOperator_initStartIndices(sopo, &index);

   printf("done.\n");

   printf("Determining local and interacting lists... ");
   fflush(stdout);

   Tree_lists(sopo->tree);
   Tree_lists(sopi->tree);
   Tree_lists(sopc->tree);

   printf("done.\n");

   printf("Filling tree structure... ");
   fflush(stdout);

   Tree_fill(sopo->tree);
   Tree_fill(sopi->tree);
   Tree_fill(sopc->tree);

   printf("done.\n");

   Tree_memory(sopo->tree);
   Tree_memory(sopi->tree);
   Tree_memory(sopc->tree);

   printf("Constructing preconditioner... ");
   fflush(stdout);
   
   preconditioner = Preconditioner_allocate(2 * totalnumpanels, 2 * totalnumpanels);
   SurfaceOperator_generatePreconditioner(sopo, preconditioner);
   Preconditioner_factor(preconditioner);
   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

   printf("done.\n");

   SurfaceOperator_makeRHS(sopo, rhs);

   GMRES_SurfaceOperator(sopo, preconditioner, rhs, sol, 2*totalnumpanels, gmrestol);

#ifdef PRINT_OPERATOR
   printf("x = [");
   for (i = 0; i < 2 * totalnumpanels; i++)
      printf("%f\n", sol[i]);
   printf("];\n");

   printf("b = [");
   for (i = 0; i < 2 * totalnumpanels; i++)
      printf("%f\n", rhs[i]);
   printf("];\n");

   printf("A = [");
   for (i = 0; i < 2*totalnumpanels; i++) {
      Vector_zero(rhs, 2 * totalnumpanels);
      rhs[i] = 1.0;
      SurfaceOperator_topmultiply(sopo, sol, rhs);
      for (j = 0; j < 2 * totalnumpanels; j++)
         printf("%f  ", sol[j]);
      printf("\n");
   }
   printf("]';\n");
#endif
   phi = Vector_allocate(charge->numcharges);

   SurfaceOperator_calculateReactionPotentials(sopo, sol, phi);
   for (c = 0; c < charge->numcharges; c++)
      energy += 0.5 * charge->charges[c] * phi[c];

   energy *= 332.112 / 4.0;
   
   printf("Solvation Energy: %f kcal/mol\n", energy);

   SurfaceOperator_writematlabfile("solv.m", sopo, totalnumpanels);

   Preconditioner_free(preconditioner);
   SurfaceOperator_free(sopo);
   Vector_free(sol);
   Vector_free(rhs);
   Vector_free(phi);
   for (i = 0; i < numpanels_d; i++)
      Vector3D_free(centroids_d[i]);
   for (i = 0; i < numpanels_c; i++)
      Vector3D_free(centroids_c[i]);
   free(centroids_d);
   free(centroids_c);
   free(totalcentroids);
   for (i = 0; i < numpanels_d; i++)
      Panel_free(panels_d[i]);
   for (i = 0; i < numpanels_c; i++)
      Panel_free(panels_c[i]);
   free(panels_d);
   free(panels_c);
   free(totalpanels);
   Charge_free(charge);
   VertFace_free(vertface_d);
   VertFace_free(vertface_c);

   return 0;
}
