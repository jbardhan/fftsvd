#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* gstfile = NULL;
   FILE* ljfile = NULL;
   LJparameters ljparameters;
   real svdtol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   GST* gstpanels;
   unsigned int numgstpanels;
   Panel* panels;
   unsigned int numpanels;
   unsigned int i, j;
   Vector rhs;
   Vector sol12, sol6;
   Tree tree12, tree6;
   real totalIntegral = 0.0;
   real parameter = 1.0;

   if (argc != 6) {
      printf("Usage: %s [surface.gst] [lj] [svd tol] [gridpoints] [maxpanels]\n", argv[0]);
      return -1;
   }

   gstfile = fopen(argv[1], "r");

   if (gstfile == NULL) {
      perror("Error opening GST file");
      return -2;
   }

   GST_readfile(&numgstpanels, &gstpanels, gstfile, 0);

   fclose(gstfile);

   ljparameters = LJparameters_allocate();

   ljfile = fopen(argv[2], "r");

   if (ljfile == NULL) {
      perror("Error opening LJ file");
      return -2;
   }

   LJparameters_read(ljparameters, ljfile);

   fclose(ljfile);

   svdtol = atof(argv[3]);
   gridpoints = atoi(argv[4]);
   maxpanelsperfinestcube = atoi(argv[5]);

   numpanels = numgstpanels;

   printf("NUMPANELS: %u\n", numpanels);

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   for (i = 0; i < numgstpanels; i++) {
      panels[i] = Panel_allocate();
      Panel_GST(panels[i], gstpanels[i], 0);
   }
   
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

      if (i % 100 == 0) printf("%u\n", i);
   }

   printf("direct integral is %20.10f\n", totalIntegral);
   printf("%f %f\n", six, twelve);
   exit(999);

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree12 = Tree_allocate(panels, numpanels, ljparameters->coordinates, 
                        ljparameters->numatoms, maxpanelsperfinestcube,
	                LJ12_KERNEL, &parameter, 
                        gridpoints, svdtol, DOUBLE_LAYER_INT, 4.0);
   tree6 = Tree_allocate(panels, numpanels, ljparameters->coordinates, 
                        ljparameters->numatoms, maxpanelsperfinestcube,
	                LJ6_KERNEL, &parameter, 
                        gridpoints, svdtol, DOUBLE_LAYER_INT, 4.0);

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

   for (i = 0; i < ljparameters->numatoms; i++)
      totalIntegral += ljparameters->R12terms[i] * sol12[i] -
         ljparameters->R6terms[i] * sol6[i];

   printf("fast method: %f\n", totalIntegral);

   //Tree_writematlabfile("six.m", tree6, DOUBLE_LAYER_INT);
   //Tree_writematlabfile("twelve.m", tree12, DOUBLE_LAYER_INT);

   Tree_free(tree12);
   Tree_free(tree6);
   
   for (i = 0; i < numpanels; i++)
      Panel_free(panels[i]);
   free(panels);

   return 0;
}
