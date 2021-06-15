#include "FFTSVD.h"

void Solv_ecf_multiply(Vector b, Tree tree, Vector x, real idiel, real odiel);
unsigned int num_GMRES_iter;
/*#define CONVERSION 332.06*/
#define CONVERSION 332.112
void saveECFColloc(Tree tree, char *filename) {
  unsigned int i,j;
  unsigned int numpanels = tree->numpanels;
  Matrix M = Matrix_allocate(numpanels, numpanels);
  Vector x = Vector_allocate(numpanels);
  Vector A2_x = Vector_allocate(numpanels);

  for (i = 0; i < numpanels; i++) {
	 Vector_zero(x, numpanels);
	 Vector_zero(A2_x, numpanels);
	 x[i] = 1.0;

	 Solv_ecf_multiply(A2_x, tree, x, 4.0, 80.0);
	 
	 for (j = 0; j < numpanels; j++) {
		M[j][i] = A2_x[j];
	 }
  }

  Matrix_writefile(filename, M, numpanels, numpanels);

  Vector_free(A2_x);
  Vector_free(x);
  Matrix_free(M);
}

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
   real energy = 0.0;

   if (argc != 9) {
      printf("Usage: %s [surface.vert] [surface.face] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels] A2_colloc.m\n", argv[0]);
      return -1;
   }

   printf("WARNING: This program does not work right now.\n");
   printf("         It needs a special modification to the integrator which is\n");
   printf("         not general.\n");

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
    
   // NOT USING _flip HERE!

   VertFace_readface(vertface, facefile);
    
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

   numpanels = vertface->numfaces;

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface, panels);

   printf("done.\n");

   centroids = (Vector3D*)calloc(2 * numpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels; i++) {
      centroids[i] = Vector3D_allocate();
      Vector3D_copy(centroids[i], panels[i]->centroid);
		centroids[i+numpanels] = Vector3D_allocate();
		Vector3D_copy(centroids[i+numpanels], panels[i]->normal);
   }

   rhs = Vector_allocate(numpanels);

   Charge_makerhs_ecf(rhs, charge, panels, numpanels);

   sol = Vector_allocate(numpanels);

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, centroids, numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, NORMDERIV_SINGLE_LAYER_INT, 0.0);

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

   Preconditioner_fill_identity(preconditioner);
      
   Preconditioner_factor(preconditioner);

   printf("done.\n");

   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

	saveECFColloc(tree, argv[8]);
	
	exit(-1);

   GMRES_solv_ecf(tree, preconditioner, rhs, sol, gmrestol, 4.0, 80.0);
   
   phi = Vector_allocate(charge->numcharges);

   for (i = 0; i < numpanels; i++)
      for (c = 0; c < charge->numcharges; c++) {
         real slp, dlp;
         
         slp = Integration(charge->points[c], panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
         
         phi[c] += CONVERSION * (sol[i] * slp) / 4.0;
      }
   
   for (c = 0; c < charge->numcharges; c++)
      energy += 0.5 * charge->charges[c] * phi[c];
   
   printf("Solvation Energy: %f kcal/mol\n", energy);

   Preconditioner_free(preconditioner);
   Tree_free(tree);
   Vector_free(sol);
   Vector_free(rhs);
   Vector_free(phi);
   for (i = 0; i <2* numpanels; i++)
      Vector3D_free(centroids[i]);
   free(centroids);
   for (i = 0; i < numpanels; i++)
      Panel_free(panels[i]);
   free(panels);
   Charge_free(charge);
   VertFace_free(vertface);

   return 0;
}
