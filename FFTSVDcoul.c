#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* chargefile = NULL;
   Charge charge;
   real svdtol, energy;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels;
   unsigned int numpanels;
   unsigned int i;
   Tree tree;
   Vector phi;

   if (argc != 5) {
      printf("Usage: %s [molecule.charge] [svd tol] [gridpoints] [maxpanels]\n", argv[0]);
      return -1;
   }

   charge = Charge_allocate();

   chargefile = fopen(argv[1], "r");

   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2;
   }

   Charge_read(charge, chargefile);

   fclose(chargefile);

   svdtol = atof(argv[2]);
   gridpoints = atoi(argv[3]);
   maxpanelsperfinestcube = atoi(argv[4]);

   numpanels = charge->numcharges;
   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   for (i = 0; i < numpanels; i++) {
      panels[i] = Panel_allocate();
      Panel_Vector3D(panels[i], charge->points[i]);
   }

   printf("done.\n");

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, charge->points, numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_LAYER_INT, 0.0);

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

   phi = Vector_allocate(numpanels);

   Tree_multiply(phi, tree, charge->charges, SINGLE_LAYER_INT);

   energy = Vector_dot(phi, charge->charges, numpanels);

   energy *= 332.112 * 0.5 / 4.0;

   printf("Coulombic Energy: %f\n", energy);

   Tree_free(tree);
   Vector_free(phi);
   for (i = 0; i < numpanels; i++)
      Panel_free(panels[i]);
   free(panels);
   Charge_free(charge);

   return 0;
}
