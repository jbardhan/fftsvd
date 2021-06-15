#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* ljfile = NULL;
   VertFace vertface;
   LJparameters ljparameters;
   real svdtol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels;
   unsigned int numpanels;
   unsigned int i, j;
   Vector rhs;
   Vector sol12, sol6;
   Tree tree12, tree6;
   real totalIntegral = 0.0;
   real parameter = 1.0;

   if (argc != 7) {
      printf("Usage: %s [surface.vert] [surface.face] [lj] [svd tol] [gridpoints] [maxpanels]\n", argv[0]);
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

   fclose(facefile);

   ljparameters = LJparameters_allocate();

   ljfile = fopen(argv[3], "r");

   if (ljfile == NULL) {
      perror("Error opening LJ file");
      return -2;
   }

   LJparameters_read(ljparameters, ljfile);

   fclose(ljfile);

   svdtol = atof(argv[4]);
   gridpoints = atoi(argv[5]);
   maxpanelsperfinestcube = atoi(argv[6]);

   VertFace_fix(vertface, 0);

   numpanels = vertface->numfaces;

   printf("NUMPANELS: %u\n", numpanels);

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface, panels);
   
   printf("done.\n");

   /* test LJ integral */

   real sa = 0.0;

   for (i = 0; i < numpanels; i++)
      sa += panels[i]->area;

   printf("SURFACE AREA: %f\n", sa);

   totalIntegral = 0.0;
   double six = 0.0, twelve = 0.0;

   for (i = 0; i < numpanels; i++) {
      real slp, dlp;      

      for (j = 0; j < ljparameters->numatoms; j++) {      
         dlp = Integration(ljparameters->coordinates[j], panels[i], LJ12_KERNEL, &parameter, DOUBLE_LAYER_INT);
         totalIntegral += ljparameters->R12terms[j] * dlp;
         twelve += ljparameters->R12terms[j] * dlp;

         dlp = Integration(ljparameters->coordinates[j], panels[i], LJ6_KERNEL, &parameter, DOUBLE_LAYER_INT);
         totalIntegral += -ljparameters->R6terms[j] * dlp;
         six += ljparameters->R6terms[j] * dlp;
      }

      if (i % 1000 == 0) printf("%u\n", i);
   }

   printf("direct integral is %20.10f\n", -totalIntegral);
   printf("%f %f\n", six, twelve);
   exit(999);

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree12 = Tree_allocate(panels, numpanels, ljparameters->coordinates, 
                          ljparameters->numatoms, maxpanelsperfinestcube,
                          LJ12_KERNEL, &parameter, 
                          gridpoints, svdtol, DOUBLE_LAYER_INT, 2.0);
   tree6 = Tree_allocate(panels, numpanels, ljparameters->coordinates, 
                         ljparameters->numatoms, maxpanelsperfinestcube,
                         LJ6_KERNEL, &parameter, 
                         gridpoints, svdtol, DOUBLE_LAYER_INT, 2.0);

   printf("done.\n");

   printf("Determining local and interacting lists... ");
   fflush(stdout);

   Tree_lists(tree12);
   Tree_lists(tree6);

   printf("done.\n");

   printf("Filling tree structure... ");
   fflush(stdout);

   Tree_fill(tree12);
   Tree_fill(tree6);

   printf("done.\n");

   Tree_memory(tree12);
   Tree_memory(tree6);

   rhs = Vector_allocate(numpanels);
   LJparameters_makerhs(rhs, ljparameters, panels, numpanels);

   sol12 = Vector_allocate(ljparameters->numatoms);
   sol6 = Vector_allocate(ljparameters->numatoms);

   Tree_multiply(sol12, tree12, rhs, DOUBLE_LAYER_INT);
   Tree_multiply(sol6, tree6, rhs, DOUBLE_LAYER_INT);

   totalIntegral = 0.0;

   six = twelve = 0.0;

   for (i = 0; i < ljparameters->numatoms; i++) {
      totalIntegral += ljparameters->R12terms[i] * sol12[i] -
                       ljparameters->R6terms[i] * sol6[i];
      six -= ljparameters->R6terms[i] * sol6[i];
      twelve += ljparameters->R12terms[i] * sol12[i];
   }

   printf("fast method: %f\n", totalIntegral);
   printf("%f %f\n", six, twelve);

   //Tree_writematlabfile("six.m", tree6, DOUBLE_LAYER_INT);
   //Tree_writematlabfile("twelve.m", tree12, DOUBLE_LAYER_INT);

   Tree_free(tree12);
   Tree_free(tree6);
   
   for (i = 0; i < numpanels; i++)
      Panel_free(panels[i]);
   free(panels);
   VertFace_free(vertface);

   return 0;
}
