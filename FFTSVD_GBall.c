#include "FFTSVD.h"

/*#define CONVERSION 332.06*/
#define CONVERSION 332.112

extern unsigned int allcount;
extern unsigned int quadcount;

int main(int argc, char* argv[]) {
   FILE* gstfile = NULL;
   FILE* torfile = NULL;
   FILE* chargefile = NULL;
   Charge charge;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   real *parameters;
   real *lesyngParam;
   GST* gstpanels;
   TOR* torpanels;
   Panel* panels;
   unsigned int numgstpanels = 0, numtorpanels = 0, numpanels = 0;
   Vector3D* centroids;
   unsigned int i, c, j;
   real *selfenergiesGhosh,*selfenergiesGrycuk, *selfenergiesLesyng;
   real *Born_radiiGhosh, *Born_radiiGrycuk, *Born_radiiLesyng;
   Vector rhs;
   Vector sol;
   Tree tree;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0;
   real idiel, odiel;

   if (argc != 6) {
      printf("Usage: %s [gstfile] [torfile] [chargefile] [idiel] [odiel]\n", argv[0]);
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

   charge = Charge_allocate();
   
   chargefile = fopen(argv[3], "r");
   
   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2;
   }

   Charge_read(charge, chargefile);
   
   fclose(chargefile);

   idiel = atof(argv[4]);
   odiel = atof(argv[5]);

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   printf("done.\n");

   centroids = (Vector3D*)calloc(numpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels; i++) {
      centroids[i] = Vector3D_allocate();
      Vector3D_copy(centroids[i], panels[i]->centroid);
   }

   rhs = Vector_allocate(2 * numpanels);

   sol = Vector_allocate(2 * numpanels);

   real caparea = 0.0;
   real pitarea = 0.0;
   real torarea = 0.0;
   for (i = 0; i < numgstpanels; i++) {
      if (gstpanels[i]->pit)
         pitarea += panels[i]->area;
      else
         caparea += panels[i]->area;
   } 
   for (i = numgstpanels; i < numpanels; i++)
      torarea += panels[i]->area;
   printf("Cap Area: %f\n", caparea);
   printf("Pit Area: %f\n", pitarea);
   printf("Tor Area: %f\n", torarea);
   printf("Surface Area: %f\n", caparea+pitarea+torarea);

   fflush(stdout);
   selfenergiesGhosh = (real *)calloc(charge->numcharges, sizeof(real));
   Born_radiiGhosh = (real *)calloc(charge->numcharges, sizeof(real));
   selfenergiesGrycuk = (real *)calloc(charge->numcharges, sizeof(real));
   Born_radiiGrycuk = (real *)calloc(charge->numcharges, sizeof(real));
   selfenergiesLesyng = (real *)calloc(charge->numcharges, sizeof(real));
   Born_radiiLesyng = (real *)calloc(charge->numcharges, sizeof(real));
   parameters = (real *)calloc(1, sizeof(real));
   parameters[0] = 1.0;
   lesyngParam = (real *)calloc(1, sizeof(real));
   lesyngParam[0] = 4.32 / pow((idiel / odiel) + 0.33, 0.3);  // set to 6 to check lesyng integration! 
   printf("E_Ghosh  R_Ghosh  E_Grycuk  R_Grycuk  E_Lesyng  R_Lesyng\n");
   for (j = 0; j < charge->numcharges; j++) {
      for (i = 0; i < numpanels; i++) {
         selfenergiesGhosh[j] += CONVERSION/(8 * M_PI)
            * (1/odiel - 1/idiel) * Integration(charge->points[j], panels[i], GHOSH_KERNEL, parameters, DOUBLE_LAYER_INT);
         selfenergiesGrycuk[j] -= Integration(charge->points[j], panels[i], GRYCUK_KERNEL, parameters, DOUBLE_LAYER_INT);
         selfenergiesLesyng[j] -= Integration(charge->points[j], panels[i], LESYNG_KERNEL, lesyngParam, DOUBLE_LAYER_INT);
      }
      Born_radiiGhosh[j] = 1 / (selfenergiesGhosh[j] * 2 / (CONVERSION * (1/odiel - 1/idiel)));

      Born_radiiGrycuk[j] =  1. / pow( (3. / (4. * M_PI)) * selfenergiesGrycuk[j], 1./3.);
      selfenergiesGrycuk[j] = CONVERSION / 2 * (1/odiel - 1/idiel) / Born_radiiGrycuk[j];

      Born_radiiLesyng[j] = 1. /
         pow((lesyngParam[0] - 3)/(4 *M_PI) * selfenergiesLesyng[j]
             , 1./(lesyngParam[0] - 3));
      selfenergiesLesyng[j] = CONVERSION / 2 * (1/odiel - 1/idiel) / Born_radiiLesyng[j];
      printf("atom  #%d  %lf  %lf  %lf  %lf %lf  %lf\n", j,
             selfenergiesGhosh[j], Born_radiiGhosh[j],
             selfenergiesGrycuk[j], Born_radiiGrycuk[j],
             selfenergiesLesyng[j], Born_radiiLesyng[j]);
   }
   
   free(Born_radiiGhosh);
   free(Born_radiiGrycuk);
   free(Born_radiiLesyng);
   free(selfenergiesGhosh);
   free(selfenergiesGrycuk);
   free(selfenergiesLesyng);
   Vector_free(sol);
   Vector_free(rhs);
   //   Vector_free(phi);
   for (i = 0; i < numpanels; i++)
      Vector3D_free(centroids[i]);
   free(centroids);
   for (i = 0; i < numpanels; i++)
      Panel_free(panels[i]);
   free(panels);
   free(gstpanels);
   free(torpanels);
   Charge_free(charge);

   return 0;
}
