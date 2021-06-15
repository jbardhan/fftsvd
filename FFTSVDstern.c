#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* chargefile = NULL;
   VertFace vertface_d, vertface_s;
   Charge charge;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels_d;
   Panel* panels_s;
   Panel* totalpanels;
   unsigned int numpanels_d;
   unsigned int numpanels_s;
   unsigned int totalnumpanels;
   Vector3D* centroids_d;
   Vector3D* centroids_s;
   Vector3D* totalcentroids;
   unsigned int i, c;
   Vector rhs;
   Vector sol;
   SurfaceOperator sopi, sopo, soho;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0;
   if (argc != 10) {
      printf("Usage: %s [surface.vert] [surface.face] [stern.vert] [stern.face] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels]\n", argv[0]);
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

   vertface_s = VertFace_allocate();

   vertfile = fopen(argv[3], "r");

   if (vertfile == NULL) {
      perror("Error opening vertices file");
      return -2;
   }

   VertFace_readvert(vertface_s, vertfile);

   fclose(vertfile);

   facefile = fopen(argv[4], "r");
    
   if (facefile == NULL) {
      perror("Error opening face file");
      return -2;
   }
    
   VertFace_readface_flip(vertface_s, facefile);
    
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

   VertFace_fix(vertface_d, 0);
   VertFace_fix(vertface_s, 1);

   numpanels_d = vertface_d->numfaces;
   numpanels_s = vertface_s->numfaces;

   totalnumpanels = numpanels_d + numpanels_s;

   panels_d = (Panel*)calloc(numpanels_d, sizeof(Panel));
   panels_s = (Panel*)calloc(numpanels_s, sizeof(Panel));
   totalpanels = (Panel*)calloc(totalnumpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface_d, panels_d);
   VertFace_getpanels(vertface_s, panels_s);

   for (i = 0; i < numpanels_s; i++)
      totalpanels[i] = panels_s[i];
   for (i = 0; i < numpanels_d; i++)
      totalpanels[i+numpanels_s] = panels_d[i];

   printf("done.\n");

   centroids_d = (Vector3D*)calloc(numpanels_d, sizeof(Vector3D));
   centroids_s = (Vector3D*)calloc(numpanels_s, sizeof(Vector3D));
   totalcentroids = (Vector3D*)calloc(totalnumpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels_s; i++) {
      centroids_s[i] = Vector3D_allocate();
      Vector3D_copy(centroids_s[i], panels_s[i]->centroid);
      totalcentroids[i] = centroids_s[i];
   }
   for (i = 0; i < numpanels_d; i++) {
      centroids_d[i] = Vector3D_allocate();
      Vector3D_copy(centroids_d[i], panels_d[i]->centroid);
      totalcentroids[i+numpanels_s] = centroids_d[i];
   }

   rhs = Vector_allocate(2 * totalnumpanels);
   sol = Vector_allocate(2 * totalnumpanels);

   printf("Constructing and allocating surface operator... ");
   fflush(stdout);

   real kappa = sqrt(0.145) / 3.047;

   soho = SurfaceOperator_allocate();

   soho->kernel = HELMHOLTZ_KERNEL;
   soho->epsilon = 80.0;
   soho->kappa = kappa;
   soho->mynumpanels = 0;
   soho->tree = Tree_allocate(panels_s, numpanels_s, centroids_s, numpanels_s, maxpanelsperfinestcube, HELMHOLTZ_KERNEL, &kappa, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   soho->numchildren = 1;
   soho->children = (SurfaceOperator*)calloc(1, sizeof(SurfaceOperator));
   soho->resultInternal = Vector_allocate(numpanels_s);
   soho->resultExternal = Vector_allocate(numpanels_s);
   soho->parent = NULL;
   soho->charges = NULL;

   sopo = SurfaceOperator_allocate();

   sopo->kernel = POISSON_KERNEL;
   sopo->epsilon = 80.0;
   sopo->kappa = 0.0;
   sopo->mynumpanels = numpanels_s;
   sopo->tree = Tree_allocate(totalpanels, totalnumpanels, totalcentroids, totalnumpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   sopo->numchildren = 1;
   sopo->children = (SurfaceOperator*)calloc(1, sizeof(SurfaceOperator));
   sopo->resultInternal = Vector_allocate(totalnumpanels);
   sopo->resultExternal = Vector_allocate(totalnumpanels);
   sopo->parent = soho;
   sopo->charges = NULL;
   soho->children[0] = sopo;

   sopi = SurfaceOperator_allocate();

   sopi->kernel = POISSON_KERNEL;
   sopi->epsilon = 4.0;
   sopi->kappa = 0.0;
   sopi->mynumpanels = numpanels_d;
   sopi->tree = Tree_allocate(panels_d, numpanels_d, centroids_d, numpanels_d, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   sopi->numchildren = 0;
   sopi->children = NULL;
   sopi->resultInternal = Vector_allocate(numpanels_d);
   sopi->resultExternal = Vector_allocate(numpanels_d);
   sopi->parent = sopo;
   sopi->charges = charge;
   sopo->children[0] = sopi;

   unsigned int index = 0;

   SurfaceOperator_initStartIndices(soho, &index);

   printf("done.\n");

   printf("Determining local and interacting lists... ");
   fflush(stdout);

   Tree_lists(soho->tree);
   Tree_lists(sopo->tree);
   Tree_lists(sopi->tree);

   printf("done.\n");

   printf("Filling tree structure... ");
   fflush(stdout);

   Tree_fill(soho->tree);
   Tree_fill(sopo->tree);
   Tree_fill(sopi->tree);

   printf("done.\n");

   Tree_memory(soho->tree);
   Tree_memory(sopo->tree);
   Tree_memory(sopi->tree);

   printf("Constructing preconditioner... ");
   fflush(stdout);
   
   preconditioner = Preconditioner_allocate(2 * totalnumpanels, 2 * totalnumpanels);
   SurfaceOperator_generatePreconditioner(soho, preconditioner);
   //SurfaceOperator_generatePreconditioner_block(soho, preconditioner);
   Preconditioner_factor(preconditioner);
   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

   printf("done.\n");

   SurfaceOperator_makeRHS(soho, rhs);

   GMRES_SurfaceOperator(soho, preconditioner, rhs, sol, 2*totalnumpanels, gmrestol);

   phi = Vector_allocate(charge->numcharges);

   SurfaceOperator_calculateReactionPotentials(soho, sol, phi);
   for (c = 0; c < charge->numcharges; c++)
      energy += 0.5 * charge->charges[c] * phi[c];

   energy *= 332.112 / 4.0;
   
   printf("Solvation Energy: %f kcal/mol\n", energy);

   Preconditioner_free(preconditioner);
   SurfaceOperator_free(soho);
   Vector_free(sol);
   Vector_free(rhs);
   Vector_free(phi);
   for (i = 0; i < numpanels_d; i++)
      Vector3D_free(centroids_d[i]);
   for (i = 0; i < numpanels_s; i++)
      Vector3D_free(centroids_s[i]);
   free(centroids_d);
   free(centroids_s);
   free(totalcentroids);
   for (i = 0; i < numpanels_d; i++)
      Panel_free(panels_d[i]);
   for (i = 0; i < numpanels_s; i++)
      Panel_free(panels_s[i]);
   free(panels_d);
   free(panels_s);
   free(totalpanels);
   Charge_free(charge);
   VertFace_free(vertface_d);
   VertFace_free(vertface_s);

   return 0;
}
