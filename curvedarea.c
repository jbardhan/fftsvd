#include "FFTSVD.h"

/*#define CONVERSION 332.06*/
#define CONVERSION 332.112

int main(int argc, char* argv[]) {
   FILE* gstfile = NULL;
   FILE* torfile = NULL;
   GST* gstpanels;
   TOR* torpanels;
   Panel* panels;
   unsigned int numgstpanels = 0, numtorpanels = 0, numpanels = 0, i;

   if (argc != 3) {
      printf("Usage: %s [gstfile] [torfile]\n", argv[0]);
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

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   for (i = 0; i < numgstpanels; i++) {
      panels[i] = Panel_allocate();
      Panel_GST(panels[i], gstpanels[i], 0);
   }

   for (i = 0; i < numtorpanels; i++) {
      panels[i+numgstpanels] = Panel_allocate();
      Panel_TOR(panels[i+numgstpanels], torpanels[i], 0);
   }

   real caparea = 0.0;
   real pitarea = 0.0;
   real torusarea = 0.0;
   real dqcaparea = 0.0;
   real dqpitarea = 0.0;
   real dqtorusarea = 0.0;

   Vector3D point = Vector3D_allocate();

   point->x = 1000000.0;

   for (i = 0; i < numgstpanels; i++) {
      real slp;
      Integration_general_GST_single(point, gstpanels[i], CONSTANT_KERNEL, NULL, &slp);
      if (gstpanels[i]->pit) {
         pitarea += slp;
         dqpitarea += Integration(point, panels[i], CONSTANT_KERNEL, NULL, SINGLE_LAYER_INT);
      }
      else {
         caparea += slp;
         real value = Integration(point, panels[i], CONSTANT_KERNEL, NULL, SINGLE_LAYER_INT);
         printf("%g\n", value);
         dqcaparea += value;
      }
   }

   for (i = 0; i < numtorpanels; i++) {
      real slp;
      Integration_general_TOR_single(point, torpanels[i], CONSTANT_KERNEL, NULL, &slp);
      torusarea += slp;
      dqtorusarea += Integration(point, panels[numgstpanels+i], CONSTANT_KERNEL, NULL, SINGLE_LAYER_INT);
   }

   printf("\n");

   printf("Cap Area: %f\n", caparea);
   printf("Pit Area: %f\n", pitarea);
   printf("Torus Area: %f\n", torusarea);
   printf("Surface Area: %f\n", caparea+pitarea+torusarea);

   printf("\n");

   printf("DQ Cap Area: %f\n", dqcaparea);
   printf("DQ Pit Area: %f\n", dqpitarea);
   printf("DQ Torus Area: %f\n", dqtorusarea);
   printf("DQ Surface Area: %e\n", dqcaparea+dqpitarea+dqtorusarea);

   return 0;
}
