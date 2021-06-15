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
   GST* gstpanels;
   TOR* torpanels;
   Panel* panels;
   unsigned int numgstpanels = 0, numtorpanels = 0, numpanels = 0;
   Vector3D* centroids;
   unsigned int i, c;
   Vector rhs;
   Vector sol;
   Tree tree;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0;
   real idiel, odiel;

   if (argc != 10) {
      printf("Usage: %s [gstfile] [torfile] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels] [idiel] [odiel]\n", argv[0]);
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

   svdtol = atof(argv[4]);
   gmrestol = atof(argv[5]);
   gridpoints = atoi(argv[6]);
   maxpanelsperfinestcube = atoi(argv[7]);
   idiel = atof(argv[8]);
   odiel = atof(argv[9]);

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   printf("done.\n");

   centroids = (Vector3D*)calloc(numpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels; i++) {
      centroids[i] = Vector3D_allocate();
      Vector3D_copy(centroids[i], panels[i]->centroid);
   }

   rhs = Vector_allocate(2 * numpanels);

   Charge_makerhs(rhs, charge, panels, numpanels);

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

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, centroids, numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);

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

/*
     {
     FILE* file = fopen("O", "w");

     Vector x = Vector_allocate(2*numpanels);
     Vector ans = Vector_allocate(2*numpanels);

     for (i = 0; i < 2*numpanels; i++) {
     Vector_zero(x, 2*numpanels);
     x[i] = 1.0;
     Solv_multiply(ans, tree, x, idiel, odiel);
     unsigned int j;
     for (j = 0; j < 2*numpanels; j++)
     fprintf(file, "%f ", ans[j]);
     fprintf(file, "\n");
     }

     fclose(file);
     }

     //exit(-999);
*/
   printf("Constructing preconditioner... ");
   fflush(stdout);
   
   preconditioner = Preconditioner_allocate(2 * numpanels, 2 * numpanels);
   //Preconditioner_fill_identity(preconditioner);
      
   Preconditioner_fill_blockdiagonal_solv(preconditioner, tree, idiel, odiel);
   Preconditioner_factor(preconditioner);

   printf("done.\n");

   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

   GMRES_solv(tree, preconditioner, rhs, sol, gmrestol, idiel, odiel);
   
   phi = Vector_allocate(charge->numcharges);

   /* If this is set to 1, the M3 operator will be "accelerated".  The
      resulting Tree operator is actually very inaccurate.  This is because
      almost nothing is done densely, since the charges are usually quite
      far from the boundary */
      
   if (0) {
      Tree M3;
      Vector dlpphi = Vector_allocate(charge->numcharges);
   
      M3 = Tree_allocate(panels, numpanels, charge->points, charge->numcharges, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
      Tree_lists(M3);
      Tree_fill(M3);
      Tree_memory(M3);
   
      Tree_multiply(phi, M3, sol+numpanels, SINGLE_LAYER_INT);
      Tree_multiply(dlpphi, M3, sol, DOUBLE_LAYER_INT);
   
      Vector_subtractvector(phi, dlpphi, charge->numcharges);
      Vector_scale(phi, CONVERSION / idiel, charge->numcharges);

      Vector_free(dlpphi);
      Tree_free(M3);
   } else {
      /* Otherwise, we just compute the elements of M3 and multiply with them
         as they are needed */
   
      for (i = 0; i < numpanels; i++)
         for (c = 0; c < charge->numcharges; c++) {
            real slp, dlp;
         
            slp = Integration(charge->points[c], panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
            dlp = Integration(charge->points[c], panels[i], POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);
         
            phi[c] += CONVERSION * (sol[i+numpanels] * slp - sol[i] * dlp) / idiel;
         }
   }
   
   for (c = 0; c < charge->numcharges; c++)
      energy += 0.5 * charge->charges[c] * phi[c];
   
   printf("Solvation Energy: %f kcal/mol\n", energy);

   printf("ALLCOUNT: %u\n", allcount);
   printf("QUADCOUNT: %u\n", quadcount);

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
   free(gstpanels);
   free(torpanels);
   Charge_free(charge);

   return 0;
}
