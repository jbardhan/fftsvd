#include "FFTSVD.h"

/*#define CONVERSION 332.06*/
#define CONVERSION 332.112

int main(int argc, char* argv[]) {
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* chargefile = NULL;
   VertFace vertface;
   Charge charge_src, charge_dest;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels;
   unsigned int numpanels;
   Vector3D* centroids;
   unsigned int i, c, colcount;
   Vector rhs;
   Vector sol;
   Vector curchargevec;
   Tree tree;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0;

   if (argc != 9) {
      printf("Usage: %s [surface.vert] [surface.face] [srcchargefile] [destchargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels]\n", argv[0]);
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

   charge_src = Charge_allocate();
   
   chargefile = fopen(argv[3], "r");
   
   if (chargefile == NULL) {
      perror("Error opening src charge file");
      return -2;
   }

   Charge_read(charge_src, chargefile);
   
   fclose(chargefile);

   charge_dest = Charge_allocate();
   chargefile = fopen(argv[4], "w");
   if (chargefile == NULL) {
      perror("Error opening dest charge file");
      return -2;
   }
   Charge_read(charge_dest, chargefile);
   fclose(chargefile);
   
   svdtol = atof(argv[5]);
   gmrestol = atof(argv[6]);
   gridpoints = atoi(argv[7]);
   maxpanelsperfinestcube = atoi(argv[8]);

   numpanels = vertface->numfaces;

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface, panels);

   printf("done.\n");

   centroids = (Vector3D*)calloc(numpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels; i++) {
      centroids[i] = Vector3D_allocate();
      Vector3D_copy(centroids[i], panels[i]->centroid);
   }


   sol = Vector_allocate(2 * numpanels);

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, centroids, numpanels, maxpanelsperfinestcube, GreensFunction_oneoverr, Integration_oneoverr, NULL, gridpoints, svdtol, SINGLE_LAYER | DOUBLE_LAYER);

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
   
   preconditioner = Preconditioner_allocate(2 * numpanels, 2 * numpanels);
   //Preconditioner_fill_identity(preconditioner);
      
   Preconditioner_fill_blockdiagonal_solv(preconditioner, panels, numpanels, Integration_oneoverr, NULL);
   Preconditioner_factor(preconditioner);

   printf("done.\n");

   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

   rhs = Vector_allocate(2 * numpanels);

   Charge_makerhs(rhs, charge_src, panels, numpanels);
   
   GMRES_solv(tree, preconditioner, rhs, sol, gmrestol);
   
   phi = Vector_allocate(charge_dest->numcharges);
   
   /* If this is set to 1, the M3 operator will be "accelerated".  The
      resulting Tree operator is actually very inaccurate.  This is because
      almost nothing is done densely, since the charges are usually quite
      far from the boundary */
   
   if (0) {
      Tree M3;
      Vector dlpphi = Vector_allocate(charge_dest->numcharges);
     
      M3 = Tree_allocate(panels, numpanels, charge_dest->points, charge_dest->numcharges, maxpanelsperfinestcube, GreensFunction_oneoverr, Integration_oneoverr, NULL, gridpoints, svdtol, SINGLE_LAYER | DOUBLE_LAYER);
      Tree_lists(M3);
      Tree_fill(M3);
      Tree_memory(M3);
     
      Tree_multiply(phi, M3, sol+numpanels, SINGLE_LAYER);
      Tree_multiply(dlpphi, M3, sol, DOUBLE_LAYER);
     
      Vector_subtractvector(phi, dlpphi, charge_dest->numcharges);
      Vector_scale(phi, CONVERSION / 4.0, charge_dest->numcharges);
     
      Vector_free(dlpphi);
      Tree_free(M3);
   } else {
      /* Otherwise, we just compute the elements of M3 and multiply with them
         as they are needed */

      for (i = 0; i < numpanels; i++)
         for (c = 0; c < charge_dest->numcharges; c++) {
            real slp, dlp;
         
            Integration_oneoverr(charge_dest->points[c], panels[i], NULL, &slp, &dlp);
         
            phi[c] += CONVERSION * (sol[i+numpanels] * slp - sol[i] * dlp) / 4.0;
         }
   }
   printf("pot = [");
   for (i = 0; i < charge_dest->numcharges; i++)
      printf("%f\n", phi[i]);
   printf("];\n");
   //   printf("% Solvation Energy: %f kcal/mol\n", energy);

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
   Charge_free(charge_src);
   Charge_free(charge_dest);
   
   VertFace_free(vertface);

   return 0;
}
