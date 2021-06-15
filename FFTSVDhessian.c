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
   unsigned int i, c, colcount;
   Vector rhs;
   Vector sol;
   Vector curchargevec;
   Tree tree;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0;

   if (argc != 8) {
      printf("Usage: %s [surface.vert] [surface.face] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels]\n", argv[0]);
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

   centroids = (Vector3D*)calloc(numpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels; i++) {
      centroids[i] = Vector3D_allocate();
      Vector3D_copy(centroids[i], panels[i]->centroid);
   }


   sol = Vector_allocate(2 * numpanels);

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, centroids, numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, NONE, 0.0);

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
   curchargevec = Vector_allocate(charge->numcharges);
   printf("L = [");
   for (colcount = 0; colcount < charge->numcharges; colcount++) {
      Vector_zero(curchargevec, charge->numcharges);
      curchargevec[colcount] = 1.0;
        
      Charge_makerhs_fromVector(rhs, charge, curchargevec, panels, numpanels);
      //    Charge_makerhs(rhs, charge, panels, numpanels);
     
      GMRES_solv(tree, preconditioner, rhs, sol, gmrestol);
     
      phi = Vector_allocate(charge->numcharges);
     
      /* If this is set to 1, the M3 operator will be "accelerated".  The
         resulting Tree operator is actually very inaccurate.  This is because
         almost nothing is done densely, since the charges are usually quite
         far from the boundary */
     
      if (0) {
         Tree M3;
         Vector dlpphi = Vector_allocate(charge->numcharges);
       
         M3 = Tree_allocate(panels, numpanels, charge->points, charge->numcharges, maxpanelsperfinestcube, GreensFunction_oneoverr, Integration_oneoverr, NULL, gridpoints, svdtol, SINGLE_LAYER | DOUBLE_LAYER, NONE);
         Tree_lists(M3);
         Tree_fill(M3);
         Tree_memory(M3);
       
         Tree_multiply(phi, M3, sol+numpanels, SINGLE_LAYER);
         Tree_multiply(dlpphi, M3, sol, DOUBLE_LAYER);
       
         Vector_subtractvector(phi, dlpphi, charge->numcharges);
         Vector_scale(phi, CONVERSION / 4.0, charge->numcharges);
       
         Vector_free(dlpphi);
         Tree_free(M3);
      } else {
         /* Otherwise, we just compute the elements of M3 and multiply with them
            as they are needed */
       
         for (i = 0; i < numpanels; i++)
            for (c = 0; c < charge->numcharges; c++) {
               real slp, dlp;
           
               Integration_oneoverr(charge->points[c], panels[i], NULL, &slp, &dlp);
           
               phi[c] += CONVERSION * (sol[i+numpanels] * slp - sol[i] * dlp) / 4.0;
            }
      }
     
      for (c = 0; c < charge->numcharges; c++)
         printf("\t%f", phi[c]);
      printf(";\n");
   }
   printf("]';\n");
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
   Charge_free(charge);
   VertFace_free(vertface);

   return 0;
}
