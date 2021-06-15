#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* chargefile = NULL;
   VertFace vertface_d;
   Charge charge;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels_d;
   unsigned int numpanels_d;
   Vector3D* centroids_d;
   unsigned int i, c;
   Vector rhs;
   Vector sol;
   SurfaceOperator sopi, soho;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0;
   if (argc != 8) {
      printf("Usage: %s [surface.vert] [surface.face] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels]\n", argv[0]);
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

   svdtol = atof(argv[4]);
   gmrestol = atof(argv[5]);
   gridpoints = atoi(argv[6]);
   maxpanelsperfinestcube = atoi(argv[7]);

   numpanels_d = vertface_d->numfaces;

   panels_d = (Panel*)calloc(numpanels_d, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface_d, panels_d);

   printf("done.\n");

   centroids_d = (Vector3D*)calloc(numpanels_d, sizeof(Vector3D));
   
   for (i = 0; i < numpanels_d; i++) {
      centroids_d[i] = Vector3D_allocate();
      Vector3D_copy(centroids_d[i], panels_d[i]->centroid);
   }

   charge = Charge_allocate();
   
   chargefile = fopen(argv[3], "r");
   
   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2;
   }

   Charge_read_centroids(charge, chargefile, centroids_d, numpanels_d);
   charge->globalindexstart = 0;
   fclose(chargefile);

   rhs = Vector_allocate(2 * numpanels_d);
   sol = Vector_allocate(2 * numpanels_d);

   printf("Constructing and allocating surface operator... ");
   fflush(stdout);

   real kappa = sqrt(0.145) / 3.047;

   soho = SurfaceOperator_allocate();

   soho->kernel = HELMHOLTZ_KERNEL;
   soho->epsilon = 80.0;
   soho->kappa = kappa;
   soho->mynumpanels = 0;
   soho->tree = Tree_allocate(panels_d, numpanels_d, centroids_d, numpanels_d, maxpanelsperfinestcube, HELMHOLTZ_KERNEL, &kappa, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   soho->numchildren = 1;
   soho->children = (SurfaceOperator*)calloc(1, sizeof(SurfaceOperator));
   soho->resultInternal = Vector_allocate(numpanels_d);
   soho->resultExternal = Vector_allocate(numpanels_d);
   soho->parent = NULL;
   soho->charges = NULL;

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
   sopi->parent = soho;
   sopi->charges = charge;
   soho->children[0] = sopi;

   unsigned int index = 0;

   SurfaceOperator_initStartIndices(soho, &index);

   printf("done.\n");

   printf("Determining local and interacting lists... ");
   fflush(stdout);

   Tree_lists(soho->tree);
   Tree_lists(sopi->tree);

   printf("done.\n");

   printf("Filling tree structure... ");
   fflush(stdout);

   Tree_fill(soho->tree);
   Tree_fill(sopi->tree);

   printf("done.\n");

   Tree_memory(soho->tree);
   Tree_memory(sopi->tree);

   printf("Constructing preconditioner... ");
   fflush(stdout);
   
   preconditioner = Preconditioner_allocate(2 * numpanels_d, 2 * numpanels_d);
   SurfaceOperator_generatePreconditioner(soho, preconditioner);
   //SurfaceOperator_generatePreconditioner_block(soho, preconditioner);
   Preconditioner_factor(preconditioner);
   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

   printf("done.\n");

   SurfaceOperator_makeRHS(soho, rhs);

   GMRES_SurfaceOperator(soho, preconditioner, rhs, sol, 2*numpanels_d, gmrestol);

   phi = Vector_allocate(charge->numcharges);

   SurfaceOperator_calculateReactionPotentials(soho, sol, phi);
   for (c = 0; c < charge->numcharges; c++)
      energy += 0.5 * charge->charges[c] * phi[c];

   energy *= 332.112 / 4.0;
   
   printf("Solvation Energy: %f kcal/mol\n", energy);

   /*
     printf("Surface Potentials:\n");

     for (i = 0; i < numpanels_d; i++)
     printf("%f %f\n", 332.112 * rhs[i] / 4.0, 332.112 * phi[charge->numcharges-numpanels_d+i] / 4.0);
   */

   Preconditioner_free(preconditioner);
   SurfaceOperator_free(soho);
   Vector_free(sol);
   Vector_free(rhs);
   Vector_free(phi);
   for (i = 0; i < numpanels_d; i++)
      Vector3D_free(centroids_d[i]);
   free(centroids_d);
   for (i = 0; i < numpanels_d; i++)
      Panel_free(panels_d[i]);
   free(panels_d);
   Charge_free(charge);
   VertFace_free(vertface_d);

   return 0;
}
