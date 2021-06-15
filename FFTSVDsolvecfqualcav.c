#include "FFTSVD.h"
#include <unistd.h>

/*#define CONVERSION 332.06*/
#define CONVERSION 332.112
#undef M1M3

int main(int argc, char* argv[]) {
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* chargefile = NULL;
   VertFace vertface;
   VertFace* cavvertface;
   Charge charge;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels;
   unsigned int numdielpanels;
   unsigned int numcavities = 0;
   unsigned int numpanels;
   Vector3D* centroids;
   unsigned int i, c;
   Vector rhs, rhs2;
   Vector sol;
   Tree tree;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0, idiel, odiel;
#ifdef M1M3
   Tree m1m3;
#endif

   if (argc != 11) {
      printf("Usage: %s [surface.vert] [surface.face] [cavitybase] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels] [idiel] [odiel]\n", argv[0]);
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

   for (i = 1; ; i++) {
      char filename[256];
      sprintf(filename, "%s_%u.face", argv[3], i);
      if (!access(filename, R_OK))
         numcavities++;
      else
         break;
   }

   printf("NUMCAVITIES: %u\n", numcavities);

   cavvertface = (VertFace*)calloc(numcavities, sizeof(VertFace));

   for (i = 0; i < numcavities; i++) {
      char filename[256];
      sprintf(filename, "%s_%u.vert", argv[3], i+1);
      cavvertface[i] = VertFace_allocate();
      FILE* cavvertfile = fopen(filename, "r");
      VertFace_readvert(cavvertface[i], cavvertfile);
      fclose(cavvertfile);
      sprintf(filename, "%s_%u.face", argv[3], i+1);
      FILE* cavfacefile = fopen(filename, "r");
      VertFace_readface_flip(cavvertface[i], cavfacefile);
      fclose(cavfacefile);
   }

   charge = Charge_allocate();
   
   chargefile = fopen(argv[4], "r");
   
   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2;
   }

   Charge_read(charge, chargefile);
   
   fclose(chargefile);

   svdtol = atof(argv[5]);
   gmrestol = atof(argv[6]);
   gridpoints = atoi(argv[7]);
   maxpanelsperfinestcube = atoi(argv[8]);
   idiel = atoi(argv[9]);
   odiel = atoi(argv[10]);

   if (odiel < 1.1) {
      VertFace_fix(vertface, 1);
      for (i = 0; i < numcavities; i++)
         VertFace_fix(cavvertface[i], 1);
   }
   else {
      VertFace_fix(vertface, 0);
      for (i = 0; i < numcavities; i++)
         VertFace_fix(cavvertface[i], 0);
   }

   numdielpanels = vertface->numfaces;
   numpanels = numdielpanels;
   for (i = 0; i < numcavities; i++)
      numpanels += cavvertface[i]->numfaces;

   printf("NUMPANELS: %u\n", numpanels);

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface, panels);
   c = 0;
   for (i = 0; i < numcavities; i++) {
      VertFace_getpanels(cavvertface[i], panels + numdielpanels + c);
      c += cavvertface[i]->numfaces;
   }

   printf("done.\n");

   real area = 0.0;
   for (i = 0; i < numpanels; i++)
      area += panels[i]->area;

   printf("Surface Area: %f\n", area);

   centroids = (Vector3D*)calloc(numpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels; i++) {
      centroids[i] = Vector3D_allocate();
      Vector3D_copy(centroids[i], panels[i]->centroid);
   }

   rhs = Vector_allocate(numpanels);
   rhs2 = Vector_allocate(numpanels);

#ifdef M1M3
   m1m3 = Tree_allocate(panels, numpanels, charge->points, charge->numcharges, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   Tree_lists(m1m3);
   Tree_fill(m1m3);
   Tree_memory(m1m3);
   Vector q = Vector_allocate(charge->numcharges);
   for (i = 0; i < charge->numcharges; i++)
      q[i] = -charge->charges[i] / (4.0 * M_PI * idiel);
   Tree_multiply_transpose(rhs, m1m3, q, DOUBLE_LAYER_INT);
   Vector_free(q);
   Charge_makerhs_ecf_qual(rhs2, charge, panels, numpanels, idiel, odiel);
   for (i = 0; i < numpanels; i++)
      printf("%e %e %e\n", rhs2[i], rhs[i], fabs(rhs2[i]-rhs[i])/fabs(rhs2[i]));
   return 999;
#else
   Charge_makerhs_ecf_qual(rhs, charge, panels, numpanels, idiel, odiel);
#endif

   sol = Vector_allocate(numpanels);

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, centroids, numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, DOUBLE_LAYER_INT, 0.0);

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

   Preconditioner_fill_diagonal_solv_ecf_qual(preconditioner, tree, idiel, odiel);

   //Preconditioner_fill_identity(preconditioner);
      
   Preconditioner_factor(preconditioner);

   printf("done.\n");

   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

   GMRES_solv_ecf_qual(tree, preconditioner, rhs, sol, gmrestol, idiel, odiel);
   
   phi = Vector_allocate(charge->numcharges);

#ifdef M1M3
   Tree_multiply(phi, m1m3, sol, SINGLE_LAYER_INT);
   Vector_scale(phi, CONVERSION / idiel, charge->numcharges);
#else
   for (i = 0; i < numpanels; i++) {
      //printf("%f\n", rhs[i]);
      for (c = 0; c < charge->numcharges; c++) {
         real slp, dlp;
         
         slp = Integration(charge->points[c], panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
         
         phi[c] += CONVERSION * (sol[i] * slp) / idiel;
      }
   }
#endif

   for (c = 0; c < charge->numcharges; c++)
      energy += 0.5 * charge->charges[c] * phi[c];
   
   printf("Solvation Energy: %f kcal/mol\n", energy);
   /*
     {
     unsigned int i, j;
     Vector x = Vector_allocate(tree->numpanels);
     Vector ans = Vector_allocate(tree->numpoints);
     FILE* file = NULL;
     FILE* areafile = NULL;
     FILE* rhsfile = NULL;
     FILE* solfile = NULL;
                                                                                
     file = fopen("qual", "w");
     areafile = fopen("area", "w");
     rhsfile = fopen("rhs", "w");
     solfile = fopen("sol", "w");
                                                                                
     for (i = 0; i < tree->numpanels; i++) {
     Vector_zero(x, tree->numpanels);
     x[i] = 1.0;
     Solv_ecf_multiply_qual(ans, tree, x, idiel, odiel);
     for (j = 0; j < tree->numpoints; j++)
     fprintf(file, "%e ", ans[j]);
     fprintf(file, "\n");
     fprintf(areafile, "%f\n", panels[i]->area);
     fprintf(rhsfile, "%f\n", rhs[i]);
     fprintf(solfile, "%f\n", sol[i]);
     }

     Vector_free(x);
     Vector_free(ans);

     fclose(file);
     fclose(areafile);
     fclose(rhsfile);
     fclose(solfile);
     }
   */
   Preconditioner_free(preconditioner);
   Tree_free(tree);
#ifdef M1M3
   Tree_free(m1m3);
#endif
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
   for (i = 0; i < numcavities; i++)
      VertFace_free(cavvertface[i]);
   free(cavvertface);

   return 0;
}
