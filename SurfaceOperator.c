#include "FFTSVD.h"
 
/* Constructors and Destructors */

SurfaceOperator SurfaceOperator_allocate() {
   SurfaceOperator so = (SurfaceOperator)calloc(1, sizeof(_SurfaceOperator));

   return so;
}

void SurfaceOperator_free(SurfaceOperator so) {
   unsigned int i;

   if (so->tree)
      Tree_free(so->tree);
   if (so->M3)
      Tree_free(so->M3);

   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_free(so->children[i]);

   free(so->children);

   if (so->charges)
      Charge_free(so->charges);

   free(so);
}

// we still need: a SurfaceOperator_initialize() or _load() function!

// initstartindices uses a preorder tree walk (CLR)
// malty, is there a better ordering--could we see some sort of thrashing
// given how i've implemented split and join these vectors?? i don't think
// we do split/join more than 2 * num surfs times
//  ie, one alternative would be to split this into two loops
//  one for loop does the index setting on children
//  the second for loop calls initStartIndices() on children
void SurfaceOperator_initStartIndices(SurfaceOperator so, unsigned int *index) {
   unsigned int i;
   for (i = 0; i < so->numchildren; i++) {
      so->children[i]->startphi=*index;
      so->children[i]->startdphi= *index + so->children[i]->mynumpanels;
      *index = *index + 2 * so->children[i]->mynumpanels;
      SurfaceOperator_initStartIndices(so->children[i], index);
   }
}

// joinpotentials assembles the appropriate phi, dphi vectors for the
// surfaceoperator object's Tree_multiply function
void SurfaceOperator_joinpotentials(SurfaceOperator so, Vector src, Vector phi, Vector dphi, unsigned int scaleDielectric) {
   unsigned int i;
   unsigned int curIndex = 0;
   Vector_copypiece(phi, curIndex, src, so->startphi, so->mynumpanels);
   Vector_copyscaledpiece(dphi, curIndex, src, so->startdphi, so->mynumpanels, -1.0);
   curIndex += so->mynumpanels;

   for (i = 0; i < so->numchildren; i++) {
      Vector_copyscaledpiece(phi, curIndex, src, so->children[i]->startphi,
                             so->children[i]->mynumpanels, -1.0);
      Vector_copyscaledpiece(dphi, curIndex, src, so->children[i]->startdphi,
                             so->children[i]->mynumpanels, so->children[i]->epsilon / so->epsilon);
      curIndex += so->children[i]->mynumpanels;
   }
}

// splitpotentials takes the result of so->tree's multiply and puts it into
// the appropriate result.  i don't know why i arranged it this way, it would
// be way easier to copy it into the global destination vector once, sparing
// ourselves the double copying...
// notice the use of two different sets of result vectors: resultInternal and
// resultExternal.  when we write green's thm for a region, we have to take the
// limit as our point goes to the outer surface, and to any inner bounding surfaces.
// in terms of our matrix equation, the rows corresponding to the outer surface
// limit are numbered the same as the _potential_ vector on the outer surface;
// the rows corresponding to inner surfaces' limits are numbered the same as
// the _potential_norm_derivative_ vectors on the appropriate inner surfaces.
void SurfaceOperator_splitpotentials(SurfaceOperator so, Vector v) {
   unsigned int i, curIndex = 0;
   Vector_copypiece(so->resultInternal, 0, v, curIndex, so->mynumpanels);
   curIndex += so->mynumpanels;
   for (i = 0; i < so->numchildren; i++) {
      Vector_copypiece(so->children[i]->resultExternal, 0, v, curIndex,
                       so->children[i]->mynumpanels);
      curIndex += so->children[i]->mynumpanels;
   }
}

// this is the top level function!!! it is NOT recursive.  call this for GMRES
void SurfaceOperator_topmultiply(SurfaceOperator so,
                                 Vector dest,
                                 Vector src)
{
   SurfaceOperator_resultclear(so);
   SurfaceOperator_multiply(so, src);
   SurfaceOperator_joinglobal(so, dest);
}

// assemble entire output vector from top level multiply
void SurfaceOperator_joinglobal(SurfaceOperator so, Vector v)
{
   unsigned int i;
   Vector_copypiece(v, so->startphi, so->resultInternal, 0, so->mynumpanels);
   Vector_copypiece(v, so->startdphi, so->resultExternal, 0, so->mynumpanels);
   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_joinglobal(so->children[i], v);
}

// the _recursive_ surface multiply.  remember that our top level surface is
// the universe, whose bounding surface has zero panels, and the surrounded
// surfaces are really the Stern layers of all the molecules in our system (if
// we're using Stern layers) or the molecular surfaces for all molecules in our
// system.
void SurfaceOperator_multiply(SurfaceOperator so, Vector src)
{
   // first initialize everything.
   unsigned int i;
   real factor;
   Vector joinedSurfacesSrcPhi, joinedSurfacesSrcDPhi, joinedSurfacesDest;
   joinedSurfacesSrcPhi = Vector_allocate(so->tree->numpanels);
   joinedSurfacesSrcDPhi = Vector_allocate(so->tree->numpanels);
   joinedSurfacesDest = Vector_allocate(so->tree->numpanels);
   Vector_zero(joinedSurfacesDest, so->tree->numpanels);

   // join different surface potentials to produce correctly sized vector for tree!
   SurfaceOperator_joinpotentials(so, src,
                                  joinedSurfacesSrcPhi,
                                  joinedSurfacesSrcDPhi,
                                  1);
  
   // this is why we zero'd joinedSurfacesDest: tree_multiply doesn't clear
   // notice i'm relying on a combined double layer/single layer integration
   Tree_multiplyboth(joinedSurfacesDest, so->tree,
                     joinedSurfacesSrcDPhi, // associated with single layer kernel
                     joinedSurfacesSrcPhi); // associated with double layer kernel
  
   // do double layer correction;  ONLY required as we go to a surface from
   // the OUTSIDE.  this means, skip correction on SurfOper so itself, but
   // do children's
   //  printf("doing double layer correction for surface with %d panels and %d children", so->mynumpanels, so->numchildren);
   if (so->kernel == HELMHOLTZ_KERNEL)
      factor = -4.0;
   else
      factor = -4.0;

   for (i = so->mynumpanels; i < so->tree->numpanels; i++) 
      joinedSurfacesDest[i] += factor * M_PI * joinedSurfacesSrcPhi[i];
  
   // now split to the operators' result vectors
   SurfaceOperator_splitpotentials(so, joinedSurfacesDest);
  
   Vector_free(joinedSurfacesSrcPhi);
   Vector_free(joinedSurfacesSrcDPhi);
   Vector_free(joinedSurfacesDest);

   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_multiply(so->children[i], src);
}

void SurfaceOperator_resultclear(SurfaceOperator so)
{
   unsigned int i;
   Vector_zero(so->resultInternal, so->mynumpanels);
   Vector_zero(so->resultExternal, so->mynumpanels);
   for (i = 0; i < so->numchildren; i++) {
      SurfaceOperator_resultclear(so->children[i]);
   }
}

void SurfaceOperator_generatePreconditioner(SurfaceOperator so, Preconditioner P)
{
   unsigned int i;
   real relepsilon = 1.0;
   real parameters = so->kappa;

   if (so->parent != NULL) {
      relepsilon = so->epsilon / so->parent->epsilon;
      parameters = so->parent->kappa;
   }

   Vector diag = Vector_allocate(so->tree->numpanels);
  
   Tree_extractdiagonal(diag, so->tree, DOUBLE_LAYER_INT);
   for (i = 0; i < so->mynumpanels; i++) {
      Preconditioner_set(P, so->startphi + i, so->startphi + i, diag[i]);
      Preconditioner_set(P, so->startdphi + i, so->startphi + i, diag[i]);
   }

   Tree_extractdiagonal(diag, so->tree, SINGLE_LAYER_INT);
   for (i = 0; i < so->mynumpanels; i++) {
      Preconditioner_set(P, so->startphi + i, so->startdphi + i, -diag[i]);
      Preconditioner_set(P, so->startdphi + i, so->startdphi + i, relepsilon * diag[i]);
   }
  
   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_generatePreconditioner(so->children[i], P);
}

void SurfaceOperator_generatePreconditioner_diagonal(SurfaceOperator so, Preconditioner P)
{
   unsigned int i;
   real relepsilon = 1.0;
   real parameters = so->kappa;

   if (so->parent != NULL) {
      relepsilon = so->epsilon / so->parent->epsilon;
      parameters = so->parent->kappa;
   }

   Vector diag = Vector_allocate(so->tree->numpanels);
  
   Tree_extractdiagonal(diag, so->tree, DOUBLE_LAYER_INT);
   for (i = 0; i < so->mynumpanels; i++) {
      Preconditioner_set(P, so->startphi + i, so->startphi + i, diag[i]);
   }

   Tree_extractdiagonal(diag, so->tree, SINGLE_LAYER_INT);
   for (i = 0; i < so->mynumpanels; i++) {
      Preconditioner_set(P, so->startdphi + i, so->startdphi + i, relepsilon * diag[i]);
   }
  
   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_generatePreconditioner_diagonal(so->children[i], P);
}

void SurfaceOperator_generatePreconditioner_identity(SurfaceOperator so, Preconditioner P)
{
   unsigned int i;
   real relepsilon = 1.0;
   real parameters = so->kappa;

   if (so->parent != NULL) {
      relepsilon = so->epsilon / so->parent->epsilon;
      parameters = so->parent->kappa;
   }

   Vector diag = Vector_allocate(so->tree->numpanels);
  
   for (i = 0; i < so->mynumpanels; i++) {
      Preconditioner_set(P, so->startphi + i, so->startphi + i, 1.0);
   }

   for (i = 0; i < so->mynumpanels; i++) {
      Preconditioner_set(P, so->startdphi + i, so->startdphi + i, 1.0);
   }
  
   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_generatePreconditioner_identity(so->children[i], P);
}

void SurfaceOperator_generatePreconditioner_block_recursive(Preconditioner P, unsigned int startphi, unsigned int startdphi, real relepsilon, Cube cube) {
   unsigned int cx, cy, cz, i, j;

   if (cube->leaf) {
      for (i = 0; i < cube->numpointindices; i++)
         for (j = 0; j < cube->numpanelindices; j++) {
            if (i != j)
               continue;

            Preconditioner_set(P, startphi+cube->pointindices[i], startphi+cube->panelindices[j], cube->D_double[i][j]);
            Preconditioner_set(P, startdphi+cube->pointindices[i], startphi+cube->panelindices[j], -cube->D_double[i][j]);
            Preconditioner_set(P, startphi+cube->pointindices[i], startdphi+cube->panelindices[j], -cube->D_single[i][j]);
            Preconditioner_set(P, startdphi+cube->pointindices[i], startdphi+cube->panelindices[j], relepsilon * cube->D_single[i][j]);
         }
   }

   /* Recurse over children */

   for (cx = 0; cx <= 1; cx++)
      for (cy = 0; cy <= 1; cy++)
         for (cz = 0; cz <= 1; cz++)
            if (cube->children[cx][cy][cz] != NULL)
               SurfaceOperator_generatePreconditioner_block_recursive(P, startphi, startdphi, relepsilon, cube->children[cx][cy][cz]);
}

void SurfaceOperator_generatePreconditioner_block(SurfaceOperator so, Preconditioner P)
{
   unsigned int i;
   real relepsilon = 1.0;

   if (so->parent != NULL)
      relepsilon = so->epsilon / so->parent->epsilon;

   if (so->startphi != so->startdphi)
      SurfaceOperator_generatePreconditioner_block_recursive(P, so->startphi, so->startdphi, relepsilon, so->tree->root);
  
   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_generatePreconditioner_block(so->children[i], P);
}

void SurfaceOperator_makeRHS(SurfaceOperator so, Vector RHS) 
{
   unsigned int i, j, k;
   // charges in my region affect not only My Exerior Bounding Surface,
   // but also all my children's exterior bounding surfaces--that is,
   // my bounding interior surfaces
   if (so->charges != NULL) {
  
      for (i = 0; i < so->charges->numcharges; i++) {
         if (so->charges->charges[i] == 0.0)
            continue;

         for (j = 0; j < so->mynumpanels; j++) {
            RHS[so->startphi + j] += so->charges->charges[i] /
               Vector3D_distance(so->charges->points[i], so->tree->panels[j]->centroid);
         }
         for (k = 0; k < so->numchildren; k++) {
            for (j = 0; j < so->children[k]->mynumpanels; j++) {
               RHS[so->children[k]->startdphi + j] += so->charges->charges[i] /
                  Vector3D_distance(so->charges->points[i], so->children[k]->tree->panels[j]->centroid);
            }
         }
      }

   }

   // loop over all children
   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_makeRHS(so->children[i], RHS);
  
}

void SurfaceOperator_calculateReactionPotentials(SurfaceOperator so,
                                                 Vector surfacePotentials,
                                                 Vector pointPotentials)
{
   unsigned int i;

   // this is the fast solver approach: one way (stupid)--add reactTree
   // tree to SurfaceOperator --alternative is to add them on to the
   // end of the standard tree, add "totalnumelements" and use that
   // here in place of numpanels for joinedSurfacesReact--then just
   // grab from end of the resulting vector....  the other thing to add
   // is a "globalindexstart" to Charge
   /*   Vector joinedSurfacesPhi, joinedSurfacesDPhi, joinedSurfacesReact; */
   /*   joinedSurfacesPhi = Vector_allocate(so->tree->numpanels); */
   /*   joinedSurfacesDPhi = Vector_allocate(so->tree->numpanels); */
   /*   joinedSurfacesReact = Vector_allocate(so->charges->numcharges); // would be replace with totalnumelements */
   /*   SurfaceOperator_joinpotentials(so, src, joinedSurfacesSrcPhi, joinedSurfacesSrcDPhi, 1); */
   /*   Tree_multiplyboth(joinedSurfacesReact, so->reactTree, joinedSurfacesDPhi,//would replace with so->tree */
   /*                     joinedSurfacesPhi); */
   /*   for (i = 0; i < so->charges->numcharges; i++) { */
   /*     pointPotentials[so->charges->globalindexstart + i] = joinedSurfacesReact[i]; //would replace with numpanels+i */
   /*   } */
   /*   Vector_free(joinedSurfacesPhi); */
   /*   Vector_free(joinedSurfacesDPhi); */
   /*   Vector_free(joinedSurfacesReact); */

   if (so->M3) {
      if (so->charges != NULL)  {
         Vector joinedSurfacesPhi, joinedSurfacesDPhi, joinedSurfacesReact;
         joinedSurfacesPhi = Vector_allocate(so->tree->numpanels);
         joinedSurfacesDPhi = Vector_allocate(so->tree->numpanels);
         joinedSurfacesReact = Vector_allocate(so->charges->numcharges);
         SurfaceOperator_joinpotentials(so, surfacePotentials, joinedSurfacesPhi, joinedSurfacesDPhi, 1);
         Tree_multiplyboth(joinedSurfacesReact, so->M3, joinedSurfacesDPhi, joinedSurfacesPhi);
         for (i = 0; i < so->charges->numcharges; i++)
            pointPotentials[so->charges->globalindexstart + i] = -joinedSurfacesReact[i];
         Vector_free(joinedSurfacesPhi);
         Vector_free(joinedSurfacesDPhi);
         Vector_free(joinedSurfacesReact);
      }
   }
   else {
      // this is the direct approach since trees are inaccurate
      int j,k;
      real slp, dlp;
      if (so->charges != NULL)  {
#ifdef OMP
#pragma omp parallel private(i, j, slp, dlp)
#pragma omp for schedule(dynamic, 1)
#endif
         for (j = 0; j < so->charges->numcharges; j++) { // reverse indices for parallel
            for (i = 0; i < so->mynumpanels; i++) {
               slp = Integration(so->charges->points[j], so->tree->panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
               dlp = Integration(so->charges->points[j], so->tree->panels[i], POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);
               pointPotentials[so->charges->globalindexstart+j] +=
                  (surfacePotentials[so->startdphi+i] * slp - surfacePotentials[so->startphi+i] * dlp);
               if (Vector3D_equal(so->charges->points[j], so->tree->panels[i]->centroid))
                  pointPotentials[so->charges->globalindexstart+j] -= 4.0 * M_PI * surfacePotentials[so->startphi+i];
            }
         }
         // may have inadvertently inverted so->epsilon/so->child->epsilon, don't think so though
         for (k = 0; k < so->numchildren; k++) {
            for (i = 0; i < so->children[k]->mynumpanels; i++) {
               for (j = 0; j < so->charges->numcharges; j++) {
                  slp = Integration(so->charges->points[j], so->children[k]->tree->panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
                  dlp = Integration(so->charges->points[j], so->children[k]->tree->panels[i], POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);
                  pointPotentials[so->charges->globalindexstart+j] +=
                     (-so->children[k]->epsilon/so->epsilon * surfacePotentials[so->children[k]->startdphi+i] * slp
                      + surfacePotentials[so->children[k]->startphi+i] * dlp);
               }
            }
         }
      }
   }

   // recurse over children
   for (i = 0; i < so->numchildren; i++)
      SurfaceOperator_calculateReactionPotentials(so->children[i], surfacePotentials, pointPotentials);
}

void SurfaceOperator_writematlabfile(char* filename, SurfaceOperator so, unsigned int totalnumpanels) {
   unsigned int i, j;
   Vector x = Vector_allocate(2 * totalnumpanels);
   Vector ans = Vector_allocate(2 * totalnumpanels);
   FILE* file = NULL;

   file = fopen(filename, "w");

   //fprintf(file, "A = zeros(%u, %u);\n", 2 * totalnumpanels, 2 * totalnumpanels);

   //fprintf(file, "A = [\n");

   for (i = 0; i < 2 * totalnumpanels; i++) {
      Vector_zero(x, 2 * totalnumpanels);
      x[i] = 1.0;
      SurfaceOperator_topmultiply(so, ans, x);
      for (j = 0; j < 2 * totalnumpanels; j++)
         fprintf(file, "%f ", ans[j]);
      fprintf(file, "\n");
   }

	fclose(file);
   //fprintf(file, "]';\n");

   Vector_free(x);
   Vector_free(ans);
}

// CVS keeps screwing me!! argh
void SurfaceOperator_makeRHS_fromVector(SurfaceOperator so, Vector RHS,
                                        Vector srcCharges, unsigned int totalnumcharges) {
   unsigned int i, j, curLowDielectricBody, offset;
   offset = 0;
   real curchargevalue;
   SurfaceOperator curDielectricSurface;
   SurfaceOperator relevantOuterSurface;
   Vector3D curchargeposition;
   curLowDielectricBody = 0;
	// for right now this is BROKEN -- will ONLY WORK when there exists ONE and ONLY ONE
/* 	// salt-exclusion region!! */
/* 	if (so->kernel == HELMHOLTZ_KERNEL) { */
/* 	  curDielectricSurface = so->children[0]->children[curLowDielectricBody]; */
/* 	} else { */
/* 	  curDielectricSurface = so->children[curLowDielectricBody]; */
/* 	} */
/*    Vector3D curchargeposition; */

   if (so->children[curLowDielectricBody]->charges == NULL)
      relevantOuterSurface = so->children[0];
   else
      relevantOuterSurface = so;

   curDielectricSurface = relevantOuterSurface->children[curLowDielectricBody];
   offset = 0;

   for (i = 0; i < totalnumcharges; i++) {
      curchargevalue = srcCharges[i];
      if (i - offset > curDielectricSurface->charges->numcharges) {
         curLowDielectricBody++;
         if (curLowDielectricBody >= relevantOuterSurface->numchildren) { 
            printf("could not place charges!! dying.\n");
            exit(-1);
         }
         curDielectricSurface = relevantOuterSurface->children[curLowDielectricBody];
         offset = i;
      }
      curchargeposition = curDielectricSurface->charges->points[i-offset];
      for (j = 0; j < curDielectricSurface->mynumpanels; j++) {
         RHS[curDielectricSurface->startphi + j] += curchargevalue /
            Vector3D_distance(curchargeposition, curDielectricSurface->tree->panels[j]->centroid);
      }
   }
}

void SurfaceOperator_solve(SurfaceOperator so, Vector sol, Vector RHS) {
  
}

// it is expected that phi_r is of the right length!!  total number of reaction field points
void SurfaceOperator_collectPotentials_toVector(SurfaceOperator so, Vector phi_r, Vector sol) {
   SurfaceOperator_calculateReactionPotentials(so, sol, phi_r);
}
