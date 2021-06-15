#include "FFTSVD.h"
#ifdef MVTIME
#include <sys/time.h>
#include <sys/resource.h>
#endif
#ifdef OMP
#include <omp.h>
#endif

/* Constructors and Destructors */

Tree Tree_allocate(Panel* panels, unsigned int numpanels, Vector3D* points, unsigned int numpoints, unsigned int maxpanelsperfinestcube, BEMKernelType kerneltype, void* parameters, unsigned int gridpoints, real epsilon, BEMLayerType layertype, real minimumcubesize) {
   Tree tree = (Tree)calloc(1, sizeof(_Tree));
   unsigned int* panelindices = (unsigned int*)calloc(numpanels, sizeof(unsigned int));
   unsigned int* pointindices = (unsigned int*)calloc(numpoints, sizeof(unsigned int));
   unsigned int indices[3] = {0, 0, 0};
   Vector3D bounds[2];
   unsigned int p;
   int i;
   real minx = 1000000.0, miny = 1000000.0, minz = 1000000.0;
   real maxx = -1000000.0, maxy = -1000000.0, maxz = -1000000.0;
   real centerx, centery, centerz, maximumextent = 0.0;
   
   tree->panels = panels;
   tree->numpanels = numpanels;
   tree->points = points;
   tree->numpoints = numpoints;
   tree->maxpanelsperfinestcube = maxpanelsperfinestcube;
   tree->kerneltype = kerneltype;
   tree->parameters = parameters;
   tree->gridpoints = gridpoints;
   tree->epsilon = epsilon;
   tree->layertype = layertype;
   tree->minimumcubesize = minimumcubesize;

	if (tree->layertype == NORMDERIV_SINGLE_LAYER_INT)
	  	tree->normals = &(points[numpoints]);  // only used if we have NORMDERIV_SINGLE_LAYER_INT used. jpb 9/21/08

   /* Initialize quadrature points */

   switch (tree->gridpoints) {
   case 3: tree->numquadraturepoints = 25;
      tree->quadraturepoints = quadrature25;
      break;
   case 4: tree->numquadraturepoints = 49;
      tree->quadraturepoints = quadrature49;
      break;
   case 5: tree->numquadraturepoints = 81;
      tree->quadraturepoints = quadrature81;
      break;
   default: printf("No equivalent density quadrature rule associated with %u gridpoints\n", tree->gridpoints);
      exit(-999);
   }

   /* Initialize cube hierarchy */

   for (p = 0; p < numpanels; p++)
      panelindices[p] = p;

   for (p = 0; p < numpoints; p++)
      pointindices[p] = p;

   bounds[0] = Vector3D_allocate();
   bounds[1] = Vector3D_allocate();

   for (p = 0; p < numpanels; p++) {
      if (panels[p]->centroid->x < minx) minx = panels[p]->centroid->x;
      if (panels[p]->centroid->y < miny) miny = panels[p]->centroid->y;
      if (panels[p]->centroid->z < minz) minz = panels[p]->centroid->z;
      if (panels[p]->centroid->x > maxx) maxx = panels[p]->centroid->x;
      if (panels[p]->centroid->y > maxy) maxy = panels[p]->centroid->y;
      if (panels[p]->centroid->z > maxz) maxz = panels[p]->centroid->z;
   }

   for (p = 0; p < numpoints; p++) {
      if (points[p]->x < minx) minx = points[p]->x;
      if (points[p]->y < miny) miny = points[p]->y;
      if (points[p]->z < minz) minz = points[p]->z;
      if (points[p]->x > maxx) maxx = points[p]->x;
      if (points[p]->y > maxy) maxy = points[p]->y;
      if (points[p]->z > maxz) maxz = points[p]->z;
   }

   /* This gives rectangular cubes */
   /*
     bounds[0]->x = minx;
     bounds[0]->y = miny;
     bounds[0]->z = minz;
     bounds[1]->x = maxx;
     bounds[1]->y = maxy;
     bounds[1]->z = maxz;
   */
    
   /* This gives cubic cubes */

   centerx = (minx + maxx) * 0.5;
   centery = (miny + maxy) * 0.5;
   centerz = (minz + maxz) * 0.5;
   
   if ((maxx - minx) > maximumextent)
      maximumextent = maxx - minx;
   if ((maxy - miny) > maximumextent)
      maximumextent = maxy - miny;
   if ((maxz - minz) > maximumextent)
      maximumextent = maxz - minz;

   bounds[0]->x = centerx - maximumextent * 0.5;
   bounds[0]->y = centery - maximumextent * 0.5;
   bounds[0]->z = centerz - maximumextent * 0.5;
   bounds[1]->x = centerx + maximumextent * 0.5;
   bounds[1]->y = centery + maximumextent * 0.5;
   bounds[1]->z = centerz + maximumextent * 0.5;
#ifndef SCATTER
   /* Determine optimal number of partitioning levels */
   for (tree->partitioningdepth = 2; ; tree->partitioningdepth++) {
	  tree->root = Cube_allocate(panels, panelindices, numpanels, points, pointindices, numpoints, 0, indices, bounds, NULL, tree, tree->partitioningdepth);

	  if ((numpanels == 0) || (numpoints == 0))
		 break;
	  
	  if (((bounds[1]->x - bounds[0]->x) / pow(2.0, tree->partitioningdepth)) < tree->minimumcubesize) {
         printf("Minimum cube size reached\n");
         break;
	  }
	  
      unsigned int numleaves = 0, maxpanels = 0, numleaveswithinpanellimit = 0;

      Cube_leafpanelstats(tree->root, &numleaves, &maxpanels, &numleaveswithinpanellimit, tree->maxpanelsperfinestcube);

      /* FastCap criterion, 90% of the cubes must have fewer panels than 
         maxpanelsperfinestcube */
      if (((real)numleaveswithinpanellimit / (real)numleaves) > 0.90)
         break;

      Cube_free(tree->root);
   }
#else
	Cube_free(tree->root);
	tree->partitioningdepth = resolution;
	tree->root = Cube_allocate(panels, panelindices, numpanels, points, pointindices, numpoints, 0, indices, bounds, NULL, tree, tree->partitioningdepth);
#endif
   printf("Using %u partitioning levels\n", tree->partitioningdepth);

   /* Initialize gridpointsperlevel */

   tree->gridpointsperlevel = (unsigned int*)calloc(tree->partitioningdepth, sizeof(unsigned int));

   /* Right now, set every level to the user specified value */
/*
     for (i = 0; i < tree->partitioningdepth; i++)
        tree->gridpointsperlevel[i] = tree->gridpoints;
*/
   /* Or, use what maltman found is reasonable.  Use the user-specified 
      grid size for the finest two levels, and increase one grid size for
      each level coarser than that.  This should cause no noticable 
      performance or memory difference, but increase accuracy significantly. */

   tree->gridpointsperlevel[tree->partitioningdepth-1] = tree->gridpoints;

   if (tree->partitioningdepth >= 2)
      tree->gridpointsperlevel[tree->partitioningdepth-2] = tree->gridpoints;

   if (tree->partitioningdepth >= 3) {
      p = 1;

      for (i = tree->partitioningdepth-3; i >= 0; i--) {
         tree->gridpointsperlevel[i] = tree->gridpoints + p;
         p++;
      }
   }

   free(panelindices);
   free(pointindices);
   Vector3D_free(bounds[0]);
   Vector3D_free(bounds[1]);
   
   return tree;
}

void Tree_free(Tree tree) {
   unsigned int dimension = 3 + 4 * LOCAL;
   unsigned int i, j, k, d;

   Cube_free(tree->root);

#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
#ifdef POLYNOMIAL
   for (d = 0; d < tree->partitioningdepth; d++)
      Matrix_free(tree->pinv_Fgridpoints[d]);
   free(tree->pinv_Fgridpoints);
#else
   for (d = 0; d < tree->partitioningdepth; d++)
      Matrix_free(tree->pinv_Q2Sgridpoints[d]);
   free(tree->pinv_Q2Sgridpoints);
#endif
#endif
   free(tree->gridpointsperlevel);

#ifndef _DO_NOT_GENERATE_SPARSIFIED_OPERATOR_
   for (d = 0; d < tree->partitioningdepth; d++) {
      for (i = 0; i < dimension; i++) {
         for (j = 0; j < dimension; j++) {
            for (k = 0; k < dimension; k++)
               ComplexSVector_free(tree->Tprecomputed[d][i][j][k]);
            free(tree->Tprecomputed[d][i][j]);
         }
         free(tree->Tprecomputed[d][i]);
      }
      free(tree->Tprecomputed[d]);
   }
#endif
   free(tree->Tprecomputed);

#ifdef SERIALIZE
   char filename[1024];

   for (i = 0; i < 64; i++) {
      sprintf(filename, "/tmp/fftsvd/%lu_D_single_%u", tree, i);
      fclose(tree->Dfiles_single[i]);
      unlink(filename);
      sprintf(filename, "/tmp/fftsvd/%lu_D_double_%u", tree, i);
      fclose(tree->Dfiles_double[i]);
      unlink(filename);

      sprintf(filename, "/tmp/fftsvd/%lu_VT_%u", tree, i);
      fclose(tree->VTfiles[i]);
      unlink(filename);
      sprintf(filename, "/tmp/fftsvd/%lu_U_%u", tree, i);
      fclose(tree->Ufiles[i]);
      unlink(filename);

      sprintf(filename, "/tmp/fftsvd/%lu_PV_single_%u", tree, i);
      fclose(tree->PVfiles_single[i]);
      unlink(filename);
      sprintf(filename, "/tmp/fftsvd/%lu_PV_double_%u", tree, i);
      fclose(tree->PVfiles_double[i]);
      unlink(filename);

      sprintf(filename, "/tmp/fftsvd/%lu_UTI_%u", tree, i);
      fclose(tree->UTIfiles[i]);
      unlink(filename);
   }
#endif

   free(tree);
}

/* Operations */

void Tree_lists(Tree tree) {
   Cube_lists(tree->root);
}

void Tree_fill_D(Tree tree) {
   int i;

   printf("D... ");
   fflush(stdout);

   if (tree->partitioningdepth == 2) {
#ifdef OMP
#pragma omp parallel private(i)
#pragma omp for schedule(dynamic, 1)
#endif
      for (i = 0; i < 8; i++) {
#ifdef OMP
         printf("Thread %d does %u\n", omp_get_thread_num(), i);
#endif
         if (tree->root->children[0][0][i])
            Cube_fill_D(tree->root->children[0][0][i]);
#ifdef OMP
         printf("Thread %d finished %u\n", omp_get_thread_num(), i);
#endif
      }
   }
   else {
#ifdef OMP
#pragma omp parallel private(i)
#pragma omp for schedule(dynamic, 1)
#endif
      for (i = 0; i < 64; i++) {
#ifdef OMP
         printf("Thread %d does %u\n", omp_get_thread_num(), i);
#endif
         if (tree->root->children[0][0][i/8])
            if (tree->root->children[0][0][i/8]->children[0][0][i%8])
               Cube_fill_D(tree->root->children[0][0][i/8]->children[0][0][i%8]);
#ifdef OMP
         printf("Thread %d finished %u\n", omp_get_thread_num(), i);
#endif
      }
   }
}

void Tree_fill_pinv_Fitting(Tree tree) {
   unsigned int i, j, k, d;
   real rootboxsize = tree->root->bounds[1]->x - tree->root->bounds[0]->x;

   /* First, generate pinv_???gridpoints for every level */

#ifdef POLYNOMIAL
   tree->pinv_Fgridpoints = (Matrix*)calloc(tree->partitioningdepth, sizeof(Matrix));
#else
   tree->pinv_Q2Sgridpoints = (Matrix*)calloc(tree->partitioningdepth, sizeof(Matrix));
#endif

   for (d = 0; d < tree->partitioningdepth; d++) {
      real boxsize = rootboxsize / pow(2.0, d);
      real gridspacing = boxsize / (tree->gridpointsperlevel[d] - 1);
      unsigned int gp3 = tree->gridpointsperlevel[d]*tree->gridpointsperlevel[d]*tree->gridpointsperlevel[d];
      unsigned int padgridsize = (2*tree->gridpointsperlevel[d]-1)*(2*tree->gridpointsperlevel[d]-1)*((2*tree->gridpointsperlevel[d]-1)/2+1);

      /* Create a Vector3D array of the gridpoints for this depth*/

      Vector3D* g_points = (Vector3D*)calloc(gp3, sizeof(Vector3D));

      unsigned int count = 0;

      for (i = 0; i < tree->gridpointsperlevel[d]; i++)
         for (j = 0; j < tree->gridpointsperlevel[d]; j++)
            for (k = 0; k < tree->gridpointsperlevel[d]; k++) {
               g_points[count] = Vector3D_allocate();
               g_points[count]->x = i * gridspacing;
               g_points[count]->y = j * gridspacing;
               g_points[count]->z = k * gridspacing;
               count++;
            }

      Vector3D center = Vector3D_allocate();

      center->x = 0.5 * boxsize;
      center->y = 0.5 * boxsize;
      center->z = 0.5 * boxsize;

#ifdef POLYNOMIAL
      unsigned int numcoeffs = tree->gridpointsperlevel[d]*(tree->gridpointsperlevel[d]+1)*(tree->gridpointsperlevel[d]+2)/6;

      Matrix Fgridpoints = Polynomial_F_points(center, boxsize,
                                               g_points, gp3,
                                               tree->gridpointsperlevel[d]);

      tree->pinv_Fgridpoints[d] = Matrix_allocate(gp3, numcoeffs);

      Matrix_pseudoinverse(tree->pinv_Fgridpoints[d], Fgridpoints, numcoeffs, gp3);

      Matrix_free(Fgridpoints);
#else
      Matrix Q2Sgridpoints = EquivDensity_Q2S_points(center, boxsize,
                                                     g_points, gp3,
                                                     tree->quadraturepoints,
                                                     tree->numquadraturepoints,
                                                     tree->kerneltype,
                                                     tree->parameters);

      tree->pinv_Q2Sgridpoints[d] = Matrix_allocate(gp3, tree->numquadraturepoints);

      Matrix_pseudoinverse(tree->pinv_Q2Sgridpoints[d], Q2Sgridpoints, tree->numquadraturepoints, gp3);

      Matrix_free(Q2Sgridpoints);
#endif

      Vector3D_free(center);

      for (i = 0; i < gp3; i++)
         Vector3D_free(g_points[i]);

      free(g_points);
   }
}

void Tree_fill_P_I(Tree tree) {
   int i;

   printf("P... I... ");
   fflush(stdout);

#ifdef OMP
#pragma omp parallel private(i)
#pragma omp for schedule(dynamic, 1)
#endif
   for (i = 0; i < 64; i++)
      if (tree->root->children[0][0][i/8])
         if (tree->root->children[0][0][i/8]->children[0][0][i%8])
            Cube_fill_P_I(tree->root->children[0][0][i/8]->children[0][0][i%8], 1);
}

void Tree_fill_T(Tree tree) {
   unsigned int dimension = 3 + 4 * LOCAL;
   int halfdimension = 1 + 2 * LOCAL;
   unsigned int i, j, k, d;
   int x, y, z;
   real rootboxsize = tree->root->bounds[1]->x - tree->root->bounds[0]->x;

   /* Precompute the translation Fourier shift */

   printf("T... ");
   fflush(stdout);
   tree->Tprecomputed = (ComplexSVector****)calloc(tree->partitioningdepth, sizeof(ComplexSVector***));

   for (d = 0; d < tree->partitioningdepth; d++) {
      tree->Tprecomputed[d] = (ComplexSVector***)calloc(dimension, sizeof(ComplexSVector**));
      for (i = 0; i < dimension; i++) {
         tree->Tprecomputed[d][i] = (ComplexSVector**)calloc(dimension, sizeof(ComplexSVector*));
         for (j = 0; j < dimension; j++)
            tree->Tprecomputed[d][i][j] = (ComplexSVector*)calloc(dimension, sizeof(ComplexSVector));
      }
   }

   for (d = 0; d < tree->partitioningdepth; d++) {
      real gridspacing = rootboxsize / pow(2.0, d) / (tree->gridpointsperlevel[d] - 1);
      unsigned int padgridsize = (2*tree->gridpointsperlevel[d]-1)*(2*tree->gridpointsperlevel[d]-1)*((2*tree->gridpointsperlevel[d]-1)/2+1);
      for (x = -halfdimension; x <= halfdimension; x++)
         for (y = -halfdimension; y <= halfdimension; y++)
            for (z = -halfdimension; z <= halfdimension; z++) {
               if ((x < -LOCAL) || (x > LOCAL) ||
                   (y < -LOCAL) || (y > LOCAL) ||
                   (z < -LOCAL) || (z > LOCAL)) {
                  Vector3D shift = Vector3D_allocate();                 

                  shift->x = (real)x * (tree->gridpointsperlevel[d] - 1);
                  shift->y = (real)y * (tree->gridpointsperlevel[d] - 1);
                  shift->z = (real)z * (tree->gridpointsperlevel[d] - 1);

                  tree->Tprecomputed[d]
                     [x+halfdimension]
                     [y+halfdimension]
                     [z+halfdimension] = 
                     ComplexSVector_allocate(padgridsize);

                  FFT_calcDiagonalTranslationOperator(d,
                                                      tree->gridpointsperlevel[d],
                                                      shift,
                                                      gridspacing,
                                                      tree->kerneltype,
                                                      tree->parameters,
                                                      tree->Tprecomputed[d]
                                                      [x+halfdimension]
                                                      [y+halfdimension]
                                                      [z+halfdimension]);

                  Vector3D_free(shift);
               }
            }
   }
}

void Tree_fill_U_VT(Tree tree) {
   int i;

   printf("U... VT... ");
   fflush(stdout);

#ifdef OMP
#pragma omp parallel private(i)
#pragma omp for schedule(dynamic, 1)
#endif
   for (i = 0; i < 64; i++) {
#ifdef OMP
      printf("Thread %d does %u\n", omp_get_thread_num(), i);
#endif
      if (tree->root->children[0][0][i/8])
         if (tree->root->children[0][0][i/8]->children[0][0][i%8])
            Cube_fill_U_VT(tree->root->children[0][0][i/8]->children[0][0][i%8]);
#ifdef OMP
      printf("Thread %d finished %u\n", omp_get_thread_num(), i);
#endif
   }
}

void Tree_fill_PV_UTI(Tree tree) {
   int i;

   printf("PV... UTI... ");
   fflush(stdout);

#ifdef OMP
#pragma omp parallel private(i)
#pragma omp for schedule(dynamic, 1)
#endif
   for (i = 0; i < 64; i++) {
#ifdef OMP
      printf("Thread %d does %u\n", omp_get_thread_num(), i);
#endif
      if (tree->root->children[0][0][i/8])
         if (tree->root->children[0][0][i/8]->children[0][0][i%8])
            Cube_fill_PV_UTI(tree->root->children[0][0][i/8]->children[0][0][i%8]);
#ifdef OMP
      printf("Thread %d finished %u\n", omp_get_thread_num(), i);
#endif
   }
}

#ifdef ADAPTIVE

void Tree_fill_K(Tree tree) {
   int i;

   printf("K... ");
   fflush(stdout);

#ifdef OMP
#pragma omp parallel private(i)
#pragma omp for schedule(dynamic, 1)
#endif
   for (i = 0; i < 64; i++)
      if (tree->root->children[0][0][i/8])
         if (tree->root->children[0][0][i/8]->children[0][0][i%8])
            Cube_fill_K(tree->root->children[0][0][i/8]->children[0][0][i%8]);
}

#endif

void Tree_fill(Tree tree) {
#ifdef SERIALIZE
   unsigned int i;
   char filename[1024];

   mkdir("/tmp/fftsvd", 0755);

   for (i = 0; i < 64; i++) {
      sprintf(filename, "/tmp/fftsvd/%lu_D_single_%u", tree, i);
      tree->Dfiles_single[i] = fopen(filename, "w+");
      sprintf(filename, "/tmp/fftsvd/%lu_D_double_%u", tree, i);
      tree->Dfiles_double[i] = fopen(filename, "w+");

      sprintf(filename, "/tmp/fftsvd/%lu_VT_%u", tree, i);
      tree->VTfiles[i] = fopen(filename, "w+");
      sprintf(filename, "/tmp/fftsvd/%lu_U_%u", tree, i);
      tree->Ufiles[i] = fopen(filename, "w+");

      sprintf(filename, "/tmp/fftsvd/%lu_PV_single_%u", tree, i);
      tree->PVfiles_single[i] = fopen(filename, "w+");
      sprintf(filename, "/tmp/fftsvd/%lu_PV_double_%u", tree, i);
      tree->PVfiles_double[i] = fopen(filename, "w+");

      sprintf(filename, "/tmp/fftsvd/%lu_UTI_%u", tree, i);
      tree->UTIfiles[i] = fopen(filename, "w+");
   }
#endif

   Tree_fill_D(tree);
   Tree_fill_pinv_Fitting(tree);
   Tree_fill_T(tree);

#ifdef ACCELERATED_SAMPLING
   Tree_fill_P_I(tree);
#endif

   Tree_fill_U_VT(tree);

#ifdef SERIALIZE
   for (i = 0; i < 64; i++) {
      rewind(tree->VTfiles[i]);
      rewind(tree->Ufiles[i]);
   }
#endif

   Tree_fill_PV_UTI(tree);

#ifdef ADAPTIVE
   Tree_fill_K(tree);
#endif

   printf("\n");
}

void Tree_memory(Tree tree) {
   unsigned int padgridsize = (2*tree->gridpoints-1)*(2*tree->gridpoints-1)*((2*tree->gridpoints-1)/2+1);
   unsigned int Tmem = tree->partitioningdepth *
      ((3 + 4 * LOCAL) * (3 + 4 * LOCAL) * (3 + 4 * LOCAL) -
       (1 + 2 * LOCAL) * (1 + 2 * LOCAL) * (1 + 2 * LOCAL)) *
      padgridsize * 
      sizeof(complexreal);

   unsigned int Cubemem = 0, Umem = 0, VTmem = 0, UTImem = 0, PVmem = 0, Kmem = 0, Dmem = 0;
   unsigned int Panelmem = 0, i;

   for (i = 0; i < tree->numpanels; i++)
      Panelmem += Panel_memory(tree->panels[i], tree->layertype);

   Cube_memory(tree->root, &Cubemem, &Umem, &VTmem, &UTImem, &PVmem, &Kmem, &Dmem);

   printf("Memory Use For Panels: %u\n", Panelmem);
   printf("Memory Use For Cubes: %u\n", Cubemem);
   printf("Memory Use For U: %u\n", Umem);
   printf("Memory Use For VT: %u\n", VTmem);
   printf("Memory Use For UTI: %u\n", UTImem);
   printf("Memory Use For PV: %u\n", PVmem);
   printf("Memory Use For T: %u\n", Tmem);
#ifdef ADAPTIVE
   printf("Memory Use For K: %u\n", Kmem);
#endif
   printf("Memory Use For D: %u\n", Dmem);
   printf("Memory Use Total: %u\n", Panelmem+Cubemem+Umem+VTmem+UTImem+PVmem+Tmem+Kmem+Dmem);
}

void Tree_multiplyboth(Vector b, Tree tree, Vector xs, Vector xd) {
#ifdef MVTIME   
   struct rusage ruse;
   struct timeval tval;
   real starttime, endtime, startwalltime, endwalltime;
#endif
#ifdef OMP
   int i;
#endif

#ifdef MVTIME
   getrusage(RUSAGE_SELF, &ruse);
   starttime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   startwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;
#endif

#ifdef SERIALIZE
   unsigned int f;
   for (f = 0; f < 64; f++) {
      rewind(tree->Dfiles_single[f]);
      rewind(tree->Dfiles_double[f]);
      rewind(tree->VTfiles[f]);
      rewind(tree->Ufiles[f]);
      rewind(tree->PVfiles_single[f]);
      rewind(tree->PVfiles_double[f]);
      rewind(tree->UTIfiles[f]);
   }
#endif

   Vector_zero(b, tree->numpoints);

#ifdef OMP
   #pragma omp parallel private(i)
   {
   #pragma omp for schedule(dynamic, 1)
   for (i = 0; i < 8; i++)
      if (tree->root->children[0][0][i]) {
         Cube cube = tree->root->children[0][0][i];
         Cube_clear_tempvectors(cube);
         Cube_multiplyboth_FFT_PV_VT(cube, xs, xd);
      }

   #pragma omp barrier

   #pragma omp for schedule(dynamic, 1)
   for (i = 0; i < 8; i++)
      if (tree->root->children[0][0][i]) {
         Cube cube = tree->root->children[0][0][i];
         Cube_multiply_T(cube);
#ifdef ADAPTIVE
         Cube_multiply_K(cube, SINGLE_LAYER_INT);
         Cube_multiply_K(cube, DOUBLE_LAYER_INT);
#endif
         Cube_multiply_U_UTI_IFFT(b, cube);
         Cube_multiply_D(b, cube, xs, SINGLE_LAYER_INT);
         Cube_multiply_D(b, cube, xd, DOUBLE_LAYER_INT);
      }
   }
#else
   Cube_clear_tempvectors(tree->root);
   Cube_multiplyboth_FFT_PV_VT(tree->root, xs, xd);
   Cube_multiply_T(tree->root);
#ifdef ADAPTIVE
   Cube_multiply_K(tree->root, SINGLE_LAYER_INT);
   Cube_multiply_K(tree->root, DOUBLE_LAYER_INT);
#endif
   Cube_multiply_U_UTI_IFFT(b, tree->root);
   Cube_multiply_D(b, tree->root, xs, SINGLE_LAYER_INT);
   Cube_multiply_D(b, tree->root, xd, DOUBLE_LAYER_INT);
#endif

#ifdef MVTIME
   getrusage(RUSAGE_SELF, &ruse);
   endtime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   endwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;

#ifdef OMP
   printf("MV TIME: %.2f s (%.2f)\n", (endtime - starttime) / omp_get_max_threads(), endwalltime - startwalltime);
#else
   printf("MV TIME: %.2f s (%.2f)\n", endtime - starttime, endwalltime - startwalltime);
#endif
#endif
}

void Tree_multiply(Vector b, Tree tree, Vector x, BEMLayerType layertype) {
#ifdef MVTIME   
   struct rusage ruse;
   struct timeval tval;
   real starttime, endtime, startwalltime, endwalltime;
#endif
#ifdef OMP
   int i;
#endif

#ifdef MVTIME
   getrusage(RUSAGE_SELF, &ruse);
   starttime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   startwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;
#endif

#ifdef SERIALIZE
   unsigned int f;
   for (f = 0; f < 64; f++) {
      rewind(tree->Dfiles_single[f]);
      rewind(tree->Dfiles_double[f]);
      rewind(tree->VTfiles[f]);
      rewind(tree->Ufiles[f]);
      rewind(tree->PVfiles_single[f]);
      rewind(tree->PVfiles_double[f]);
      rewind(tree->UTIfiles[f]);
   }
#endif

   Vector_zero(b, tree->numpoints);

#ifdef OMP
   #pragma omp parallel private(i)
   {
   #pragma omp for schedule(dynamic, 1)
   for (i = 0; i < 8; i++)
      if (tree->root->children[0][0][i]) {
         Cube cube = tree->root->children[0][0][i];
         Cube_clear_tempvectors(cube);
         Cube_multiply_FFT_PV_VT(cube, x, layertype);
      }

   #pragma omp barrier

   #pragma omp for schedule(dynamic, 1)
   for (i = 0; i < 8; i++)
      if (tree->root->children[0][0][i]) {
         Cube cube = tree->root->children[0][0][i];
         Cube_multiply_T(cube);
#ifdef ADAPTIVE
         Cube_multiply_K(cube, layertype);
#endif
         Cube_multiply_U_UTI_IFFT(b, cube);
         Cube_multiply_D(b, cube, x, layertype);
      }
   }
#else
   Cube_clear_tempvectors(tree->root);
   Cube_multiply_FFT_PV_VT(tree->root, x, layertype);
   Cube_multiply_T(tree->root);
#ifdef ADAPTIVE
   Cube_multiply_K(tree->root, layertype);
#endif
   Cube_multiply_U_UTI_IFFT(b, tree->root);
   Cube_multiply_D(b, tree->root, x, layertype);
#endif

#ifdef MVTIME
   getrusage(RUSAGE_SELF, &ruse);
   endtime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   endwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;

#ifdef OMP
   printf("MV TIME: %.2f s (%.2f)\n", (endtime - starttime) / omp_get_max_threads(), endwalltime - startwalltime);
#else
   printf("MV TIME: %.2f s (%.2f)\n", endtime - starttime, endwalltime - startwalltime);
#endif
#endif
}

void Tree_multiply_transpose(Vector b, Tree tree, Vector x, BEMLayerType layertype) {
#ifdef MVTIME   
   struct rusage ruse;
   struct timeval tval;
   real starttime, endtime, startwalltime, endwalltime;
#endif
#ifdef OMP
   int i;
#endif

#ifdef MVTIME
   getrusage(RUSAGE_SELF, &ruse);
   starttime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   startwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;
#endif

   Vector_zero(b, tree->numpanels);

#ifdef OMP
   #pragma omp parallel private(i)
   {
   #pragma omp for schedule(dynamic, 1)
   for (i = 0; i < 8; i++)
      if (tree->root->children[0][0][i]) {
         Cube cube = tree->root->children[0][0][i];
         Cube_clear_tempvectors(cube);
         Cube_multiply_FFT_UTIT_UT(cube, x);
      }

   #pragma omp barrier

   #pragma omp for schedule(dynamic, 1)
   for (i = 0; i < 8; i++)
      if (tree->root->children[0][0][i]) {
         Cube cube = tree->root->children[0][0][i];
         Cube_multiply_T_transpose(cube);
         Cube_multiply_V_PVT_IFFT(b, cube, layertype);
         Cube_multiply_DT(b, cube, x, layertype);
      }
   }
#else
   Cube_clear_tempvectors(tree->root);
   Cube_multiply_FFT_UTIT_UT(tree->root, x);
   Cube_multiply_T_transpose(tree->root);
   Cube_multiply_V_PVT_IFFT(b, tree->root, layertype);
   Cube_multiply_DT(b, tree->root, x, layertype);
#endif

#ifdef MVTIME
   getrusage(RUSAGE_SELF, &ruse);
   endtime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   endwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;

#ifdef OMP
   printf("MV TIME: %.2f s (%.2f)\n", (endtime - starttime) / omp_get_max_threads(), endwalltime - startwalltime);
#else
   printf("MV TIME: %.2f s (%.2f)\n", endtime - starttime, endwalltime - startwalltime);
#endif
#endif
}

void Tree_writematlabfile(char* filename, Tree tree, BEMLayerType layertype) {
   unsigned int i, j;
   Vector x = Vector_allocate(tree->numpanels);
   Vector ans = Vector_allocate(tree->numpoints);
   FILE* file = NULL;

   file = fopen(filename, "w");

   fprintf(file, "A = zeros(%u, %u);\n", tree->numpanels, tree->numpoints);

   fprintf(file, "A = [\n");

   for (i = 0; i < tree->numpanels; i++) {
      Vector_zero(x, tree->numpanels);
      x[i] = 1.0;
      Tree_multiply(ans, tree, x, layertype);
      for (j = 0; j < tree->numpoints; j++)
         fprintf(file, "%e ", ans[j]);
      fprintf(file, "\n");
   }

   fprintf(file, "]';\n");

   Vector_free(x);
   Vector_free(ans);
}

void Tree_extractdiagonal(Vector d, Tree tree, BEMLayerType layertype) {
#ifdef SERIALIZE
   unsigned int f;
   if (layertype == SINGLE_LAYER_INT) {
      for (f = 0; f < 64; f++)
         rewind(tree->Dfiles_single[f]);
   }
   else if (layertype == DOUBLE_LAYER_INT) {
      for (f = 0; f < 64; f++)
         rewind(tree->Dfiles_double[f]);
   }
#endif

   Cube_extractdiagonal(d, tree->root, layertype);
}
