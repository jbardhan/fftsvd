#include "FFTSVD.h"

int main(int argc, char* argv[]) {
   FILE* gstfile = NULL;
   FILE* torfile = NULL;
   FILE* chargefile = NULL;
   FILE* facefile = NULL;
   GST* gstpanels = NULL;
   TOR* torpanels = NULL;
   Charge charge = NULL;
   Panel* panels = NULL;
   unsigned int numgstpanels = 0, numtorpanels = 0, numpanels = 0;
   unsigned int i, j;

   if (argc != 4) {
      printf("Usage: %s [gstfile] [torfile] [chargefile]\n", argv[0]);
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

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   for (i = 0; i < numgstpanels; i++) {
      panels[i] = Panel_allocate();
      Panel_GST(panels[i], gstpanels[i], 0);
   }

   for (i = 0; i < numtorpanels; i++) {
      panels[i+numgstpanels] = Panel_allocate();
      Panel_TOR(panels[i+numgstpanels], torpanels[i], 0);
   }

   chargefile = fopen(argv[3], "r");
    
   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2;
   }

   charge = Charge_allocate();

   Charge_read(charge, chargefile);

   for (i = 0; i < numpanels; i++) {
      for (j = 0; j < charge->numcharges; j++) {
         real ans = Integration(charge->points[j], panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
         printf("%e\n", ans);
      }
   }

   return 0;
}
