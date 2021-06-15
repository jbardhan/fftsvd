#include "FFTSVD.h"

/*#define CONVERSION 332.06*/
#define CONVERSION 332.112

int main(int argc, char* argv[]) {
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* chargefile = NULL;
   VertFace vertface;
   Charge charge;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels;
   unsigned int numpanels;
   Vector3D* centroids;
   unsigned int i, c;
   Vector rhs;
   Vector sol;
   Tree tree;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0, idiel, odiel;

   if (argc != 10) {
      printf("Usage: %s [surface.vert] [surface.face] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels] [idiel] [odiel]\n", argv[0]);
      return -1;
   }

   vertface = VertFace_allocate();

   vertfile = fopen(argv[1], "r");

   if (vertfile == NULL) {
      perror("Error opening vertices file");
      return -2;
   }

   VertFace_readvert(vertface, vertfile);

   fclose(vertfile);

   facefile = fopen(argv[2], "r");
    
   if (facefile == NULL) {
      perror("Error opening face file");
      return -2;
   }
    
   VertFace_readface_flip(vertface, facefile);
   //VertFace_readface(vertface, facefile);
    
   fclose(facefile);

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
   idiel = atoi(argv[8]);
   odiel = atoi(argv[9]);

   if (odiel < 1.1)
      VertFace_fix(vertface, 1);
   else
      VertFace_fix(vertface, 0);

   numpanels = vertface->numfaces;

   printf("NUMPANELS: %u\n", numpanels);

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface, panels);

   printf("done.\n");

   real area = 0.0;
   for (i = 0; i < numpanels; i++)
      area += panels[i]->area;

   printf("Surface Area: %f\n", area);

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
   /*
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
     fprintf(file, "%e ", ans[j]);
     fprintf(file, "\n");
     }

     fprintf(file, "]';\n");

     fprintf(file, "rhs = [\n");
     for (i = 0; i < tree->numpanels; i++)
     fprintf(file, "%f\n", rhs[i]);
     fprintf(file, "];\n");
                                                                                
     Vector_free(x);
     Vector_free(ans);
     }
   */
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
   VertFace_free(vertface);

   return 0;
}
