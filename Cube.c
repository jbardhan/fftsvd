#include "FFTSVD.h"

/* Constructors and Destructors */
#ifdef SCATTER

int edgeIntersectsTriangle(Vector3D v1, Vector3D v2, Vector3D normal, Vector3D *triangle)
{
  static int initialized = 0;
  static Vector3D u, v, w0, w, dir, testIntersection;
  if (! initialized) {
	 u = Vector3D_allocate();
	 v = Vector3D_allocate();
	 w0 = Vector3D_allocate();
	 w = Vector3D_allocate();
	 dir = Vector3D_allocate();
	 testIntersection = Vector3D_allocate();
	 initialized = 1;
  }

  real r, a, b;
  Vector3D_sub(u, triangle[1], triangle[0]);
  Vector3D_sub(v, triangle[2], triangle[0]);
  Vector3D_sub(w0, v1, triangle[0]);
  Vector3D_sub(dir, v2, v1);
  a = -Vector3D_dot(normal, w0);
  b = Vector3D_dot(normal, dir);

  r = a/b;

  if ((r < 0.0) || (r > 1.0)) {
/* 	 Vector3D_free(u); */
/* 	 Vector3D_free(v); */
/* 	 Vector3D_free(w0); */
/* 	 Vector3D_free(w); */
/* 	 Vector3D_free(dir); */
	 return 0;
  }
  Vector3D_addscaled(testIntersection, v1, r, dir);

/*   if (trackCube) { */
/* 	 printf("r=%f    testInt = %f %f %f\n",r, testIntersection->x, testIntersection->y, testIntersection->z); */
/*   } */
  
  real uu, uv, vv, wu, wv, D;
  uu = Vector3D_dot(u,u);
  uv = Vector3D_dot(u,v);
  vv = Vector3D_dot(v,v);
  Vector3D_sub(w, testIntersection, triangle[0]);
  wu = Vector3D_dot(w,u);
  wv = Vector3D_dot(w,v);
  D = uv * uv - uu * vv;

  real s, t;
  s = (uv * wv - vv * wu) / D;
  if ((s < 0.0) || (s > 1.0)) {
/* 	 Vector3D_free(u); */
/* 	 Vector3D_free(v); */
/* 	 Vector3D_free(w0); */
/* 	 Vector3D_free(w); */
/* 	 Vector3D_free(testIntersection); */
/* 	 Vector3D_free(dir); */
	 return 0;
  }
  t = (uv * wu - uu * wv) / D;
  if ((t < 0.0) || ((s + t) > 1.0)) {
/* 	 Vector3D_free(u); */
/* 	 Vector3D_free(v); */
/* 	 Vector3D_free(w0); */
/* 	 Vector3D_free(w); */
/* 	 Vector3D_free(testIntersection); */
/* 	 Vector3D_free(dir); */
	 return 0;
  }
/*   if (trackCube) { */
/* 	 printf("v1 = %f, %f, %f\n", v1->x, v1->y, v1->z); */
/* 	 printf("v2 = %f, %f, %f\n", v2->x, v2->y, v2->z); */
/*   } */

/*   Vector3D_free(w); */
/*   Vector3D_free(w0); */
/*   Vector3D_free(v); */
/*   Vector3D_free(u); */
/*   Vector3D_free(testIntersection); */
/*   Vector3D_free(dir); */
  return 1;

}

int edgeIntersectsSquare(Vector3D v2, Vector3D v1, 
								 Vector3D s1, Vector3D s2, Vector3D s3, Vector3D s4)
{

  static int initialized = 0;
  static Vector3D triangle1[3];
  static Vector3D triangle2[3];
  static Vector3D normal, d1, d2;
  unsigned int i;
  int doesIntersect = 0;

  if (! initialized) {
	 for (i = 0; i < 3; i++) {
		triangle1[i] = Vector3D_allocate();
		triangle2[i] = Vector3D_allocate();
	 }
	 normal = Vector3D_allocate();
	 d1 = Vector3D_allocate();
	 d2 = Vector3D_allocate();
	 initialized = 1;
  }
  
  Vector3D_copy(triangle1[0], s1);
  Vector3D_copy(triangle1[1], s2);
  Vector3D_copy(triangle1[2], s3);

  Vector3D_copy(triangle2[0], s1);
  Vector3D_copy(triangle2[1], s3);
  Vector3D_copy(triangle2[2], s4);

  Vector3D_sub(d1, s2, s1);
  Vector3D_sub(d2, s3, s1);

  Vector3D_cross(normal, d1, d2); Vector3D_normalize(normal);
  
  // assuming clockwise but it don't matter to jesus
  if (edgeIntersectsTriangle(v1, v2, normal, triangle1)
		|| edgeIntersectsTriangle(v1, v2, normal, triangle2)) {
	 doesIntersect = 1;
  } else {
	 doesIntersect = 0;
  }

/*   for (i = 0; i < 3; i++) { */
/* 	 Vector3D_free(triangle1[i]); */
/* 	 Vector3D_free(triangle2[i]); */
/*   } */
/*   Vector3D_free(normal); */
/*   Vector3D_free(d1); */
/*   Vector3D_free(d2); */

  return doesIntersect;
/*   Vector3D normal = Vector3D_allocate(); */
/*   Vector3D_copy(normal, inplane); */
/*   real scaleFactor = 1.0; */
/*   if (inplane == min) */
/* 	 scaleFactor = -1.0; */
/*   Vector3D_scale(normal, scaleFactor); */
/*   real junk1 = -(Vector3D_dot(normal, inplane) - Vector3D_dot(normal, v2)); */
/*   real junk2 = Vector3D_dot(normal, edgeVector); */
/*   Vector3D_free(normal); */
/*   if (fabs((float)junk2) < 1e-8) {  */
/* 	 return 0; */
/*   } */
/*   if (! ((junk1/junk2 >= 0) & (junk1/junk2 <= 1))) { */
/* 	 return 0; */
	 
/*   } */
/*   real s = junk1/junk2; */
/*   Vector3D_addscaled(edgeVector, v1, s, edgeVector); */
/*   if ((face->x > 0)  */
/* 		&& ((edgeVector->y > min->y) && (edgeVector->y < max->y)) */
/* 		&& ((edgeVector->z > min->z) && (edgeVector->z < max->z)) ){ */
/* 	 doesIntersect = 1; */
/*   } else if ((face->y > 0)  */
/* 				 && ((edgeVector->x > min->x) && (edgeVector->x < max->x))  */
/* 				 && ((edgeVector->z > min->z) && (edgeVector->z < max->z)) ){ */
/* 	 doesIntersect = 1; */
/*   } else if ((face->z > 0) */
/* 				 && ((edgeVector->y > min->y) && (edgeVector->y < max->y)) */
/* 				 && ((edgeVector->x > min->x) && (edgeVector->x < max->x)) ){ */
/* 	 doesIntersect = 1; */
/*   } */

/*   return doesIntersect; */

  
}

int panelIntersectsChildCube(Panel panel, Cube parent, unsigned int cx, unsigned int cy, unsigned int cz) {

  static int initialized = 0;
  unsigned int i;
  FlatPanel fp = (FlatPanel)panel->realpanel;
  int doesIntersect = 0;
  static Vector3D delta, childcenter, diff, childmin, childmax, X, Y, Z, edge, v[8];
  if (! initialized) {
	 delta = Vector3D_allocate();
	 childcenter = Vector3D_allocate();
	 diff = Vector3D_allocate();
	 childmin = Vector3D_allocate();
	 childmax = Vector3D_allocate();
	 edge = Vector3D_allocate();
	 X = Vector3D_allocate(); X->x = 1.0; X->y = 0.0; X->z = 0.0;
	 Y = Vector3D_allocate(); Y->x = 0.0; Y->y = 1.0; Y->z = 0.0;
	 Z = Vector3D_allocate(); Z->x = 0.0; Z->y = 0.0; Z->z = 1.0;
	 for (i = 0; i < 8; i++) {
		v[i] = Vector3D_allocate();
	 }
	 initialized = 1;
  }
  Vector3D_sub(delta, parent->bounds[1],parent->bounds[0]);
  real parentEdgeLength = delta->x;
  real childEdgeLength = 0.5 * parentEdgeLength;

  Vector3D_copy(childcenter, parent->center);
  childcenter->x += ((cx > 0)?+1:-1) * childEdgeLength / 2;
  childcenter->y += ((cy > 0)?+1:-1) * childEdgeLength / 2;
  childcenter->z += ((cz > 0)?+1:-1) * childEdgeLength / 2;

  // first check whether encompassing sphere can't intersect plane of triangle
  Vector3D_sub(diff, childcenter, fp->vertex[0]);
  real distanceToPlane = Vector3D_dot(diff, panel->normal);
  if (distanceToPlane*distanceToPlane > 3 * childEdgeLength * childEdgeLength) {
/* 	 Vector3D_free(diff); */
/* 	 Vector3D_free(childcenter); */
	 doesIntersect = 0;
	 return doesIntersect;
  }
  
  // then check whether any vertex is inside.  if so, boom!
  childmin->x = childcenter->x - childEdgeLength/2;
  childmin->y = childcenter->y - childEdgeLength/2;
  childmin->z = childcenter->z - childEdgeLength/2;
  childmax->x = childcenter->x + childEdgeLength/2;
  childmax->y = childcenter->y + childEdgeLength/2;
  childmax->z = childcenter->z + childEdgeLength/2;
/*   if (trackCube) { */
/*   printf("childcube (%d,%d,%d):\n", cx, cy, cz); */
/*   printf("min: %f %f %f\n", childmin->x,childmin->y, childmin->z); */
/*   printf("max: %f %f %f\n", childmax->x,childmax->y, childmax->z); */
/*   } */
  for (i = 0; i < 3; i++) {
	 if ((fp->vertex[i]->x > childmin->x) &&
		  (fp->vertex[i]->x < childmax->x) &&
		  (fp->vertex[i]->y > childmin->y) &&
		  (fp->vertex[i]->y < childmax->y) &&
		  (fp->vertex[i]->z > childmin->z) &&
		  (fp->vertex[i]->z < childmax->z)) {
		doesIntersect = 1;
/* 		Vector3D_free(diff); */
/* 		Vector3D_free(childmin); */
/* 		Vector3D_free(childmax); */
/* 		Vector3D_free(childcenter); */

/* 		if (trackCube) { */
/* 		  printf("panel has a vertex inside the cube centered at %f %f %f!\n",childcenter->x,childcenter->y,childcenter->z); */
/* 		  printf("vertex = %f, %f, %f\n", fp->vertex[i]->x, fp->vertex[i]->y,fp->vertex[i]->z); */
/* 		  printf("min = %f, %f, %f\n", childmin->x, childmin->y,childmin->z); */
/* 		  printf("max = %f, %f, %f\n", childmax->x, childmax->y,childmax->z); */
/* 		} */
		return doesIntersect;
	 }
  }

  // check whether any triangle edge is intersecting one of the sides
  unsigned int next;
/*   for (i = 0; i < 3; i++) { */
/* 	 next = i + 1; */
/* 	 if (next > 2) { */
/* 		next = 0; */
/* 	 } */
/* 	 Vector3D_sub(edge, fp->vertex[next], fp->vertex[i]); */
/* 	 if (edgeIntersectsSquare(childmin, childmax, X, childmin, edge, fp->vertex[next], fp->vertex[i]) */
/* 		  || edgeIntersectsSquare(childmin, childmax, X, childmax, edge, fp->vertex[next], fp->vertex[i]) */
/* 		  || edgeIntersectsSquare(childmin, childmax, Y, childmin, edge, fp->vertex[next], fp->vertex[i]) */
/* 		  || edgeIntersectsSquare(childmin, childmax, Y, childmax, edge, fp->vertex[next], fp->vertex[i]) */
/* 		  || edgeIntersectsSquare(childmin, childmax, Z, childmin, edge, fp->vertex[next], fp->vertex[i]) */
/* 		  || edgeIntersectsSquare(childmin, childmax, Z, childmax, edge, fp->vertex[next], fp->vertex[i])) { */
/* 		doesIntersect += 4; */
/* 		printf("panel has an edge inside the cube!\n"); */
/* 		i = 3; // break */
/* 	 } */
/*   } */


  unsigned int index, ix, iy, iz;
  for (ix = 0; ix < 2; ix++) {
	 edge->x = ix * childEdgeLength;
	 for (iy = 0; iy < 2; iy++) {
		edge->y = iy * childEdgeLength;
		for (iz = 0; iz < 2; iz++) {
		  edge->z = iz * childEdgeLength;
		  index = ix + iy*2  + iz*4;
		  Vector3D_copy(v[index], childmin);
		  Vector3D_add(v[index], edge, v[index]);
/* 		  if (trackCube) { */
/* 			 printf("v[%d] = %f %f %f\n", index, v[index]->x, v[index]->y, v[index]->z); */
/* 		  } */
		}
	 }
  }

  for (i = 0; i < 3; i++) {
	 next = (i+1) % 3;
	 if (   edgeIntersectsSquare(fp->vertex[next],fp->vertex[i], v[0], v[1], v[3], v[2])
		  || edgeIntersectsSquare(fp->vertex[next],fp->vertex[i], v[4], v[5], v[7], v[6])
		  || edgeIntersectsSquare(fp->vertex[next],fp->vertex[i], v[0], v[1], v[5], v[4])
		  || edgeIntersectsSquare(fp->vertex[next],fp->vertex[i], v[1], v[3], v[7], v[5])
		  || edgeIntersectsSquare(fp->vertex[next],fp->vertex[i], v[3], v[2], v[6], v[7])
		  || edgeIntersectsSquare(fp->vertex[next],fp->vertex[i], v[2], v[0], v[4], v[6])) {
/* 		if (trackCube) {  */
/* 		  printf("triangle edge %d intersects the panel face!\n", i); */
/* 		  printf("edge v1 = %f, %f, %f\nedge v2 = %f, %f, %f\n", fp->vertex[i]->x, fp->vertex[i]->y, fp->vertex[i]->z, fp->vertex[next]->x, fp->vertex[next]->y, fp->vertex[next]->z); */
/* 		} */
		doesIntersect += 4;
	 }
  }

/*   if (trackCube) { */
/* 	 printf("vertex0 = %f, %f, %f\n", fp->vertex[0]->x, fp->vertex[0]->y, fp->vertex[0]->z); */
/* 	 printf("vertex1 = %f, %f, %f\n", fp->vertex[1]->x, fp->vertex[1]->y, fp->vertex[1]->z); */
/* 	 printf("vertex2 = %f, %f, %f\n", fp->vertex[2]->x, fp->vertex[2]->y, fp->vertex[2]->z); */
/* 	 printf("normal = %f, %f, %f\n", panel->normal->x, panel->normal->y, panel->normal->z); */
/* 	 printf("panelaxis = %f, %f, %f\n", fp->panelaxis[2]->x, fp->panelaxis[2]->y, fp->panelaxis[2]->z); */
/*   } */
  // check whether any childcube edge intersects the triangle
  if (   edgeIntersectsTriangle(v[0], v[1], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[1], v[3], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[3], v[2], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[2], v[0], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[0], v[4], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[1], v[5], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[2], v[6], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[3], v[7], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[4], v[5], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[5], v[7], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[7], v[6], panel->normal, fp->vertex)
		|| edgeIntersectsTriangle(v[6], v[4], panel->normal, fp->vertex) ) {
	 doesIntersect += 2;
  }
/*   for (i = 0; i < 8; i++) { */
/* 	 Vector3D_free(v[i]); */
/*   } */

/*   if (trackCube) { */
/*   if (doesIntersect == 2) { */
/* 	 printf("cube just intersects the panel in its middle and nowhere else!\n"); */
/*   } */
/*   } */
/*   free(v); */
/*   Vector3D_free(X); */
/*   Vector3D_free(Y); */
/*   Vector3D_free(Z); */
/*   Vector3D_free(edge); */
/*   Vector3D_free(diff); */
/*   Vector3D_free(childmin); */
/*   Vector3D_free(childmax); */
/*   Vector3D_free(childcenter); */
  return (doesIntersect > 0);
}
#endif

Cube Cube_allocate(Panel* panels, unsigned int* panelindices, unsigned int numpanelindices, Vector3D* points, unsigned int* pointindices, unsigned int numpointindices, unsigned int level, unsigned int indices[3], Vector3D bounds[2], Cube parent, Tree tree, unsigned int partitioningdepth) {
   Cube cube = (Cube)calloc(1, sizeof(_Cube));
   unsigned int i;
   unsigned int cx, cy, cz;

   cube->level = level;
   cube->indices[0] = indices[0];
   cube->indices[1] = indices[1];
   cube->indices[2] = indices[2];
   cube->bounds[0] = Vector3D_allocate();
   cube->bounds[1] = Vector3D_allocate();
   Vector3D_copy(cube->bounds[0], bounds[0]);
   Vector3D_copy(cube->bounds[1], bounds[1]);
   cube->center = Vector3D_allocate();
   Vector3D_add(cube->center, cube->bounds[0], cube->bounds[1]);
   Vector3D_scale(cube->center, 0.5);

#ifdef SCATTER
	cube->isInsideProtein = 0;
	cube->densityInsideCube = 0.;
#endif
	
   cube->panelindices = (unsigned int*)calloc(numpanelindices, sizeof(unsigned int));
   for (i = 0; i < numpanelindices; i++)
      cube->panelindices[i] = panelindices[i];
   cube->numpanelindices = numpanelindices;

   cube->pointindices = (unsigned int*)calloc(numpointindices, sizeof(unsigned int));
   for (i = 0; i < numpointindices; i++)
      cube->pointindices[i] = pointindices[i];
   cube->numpointindices = numpointindices;

   cube->parent = parent;
   cube->tree = tree;

#ifdef SERIALIZE
   if (cube->level == 0)
      cube->serialid = 0;
   else if (cube->level == 1)
      cube->serialid = 4*cube->indices[0] + 2*cube->indices[1] + cube->indices[2];
   else if (cube->level == 2)
      cube->serialid = 16*cube->indices[0] + 4*cube->indices[1] + cube->indices[2];
   else
      cube->serialid = parent->serialid;
#endif

   if (cube->level == (partitioningdepth-1)) {
      cube->leaf = 1;
      return cube;
   }

   /* Allocate the children */

   for (cx = 0; cx <= 1; cx++)
      for (cy = 0; cy <= 1; cy++)
         for (cz = 0; cz <= 1; cz++) {
            unsigned int* childpanelindices;
            unsigned int numchildpanelindices = 0;
            unsigned int* childpointindices;
            unsigned int numchildpointindices = 0;
            unsigned int childindices[3];
            Vector3D childbounds[2];
            unsigned int count = 0;
            
            for (i = 0; i < numpanelindices; i++) {
#ifndef SCATTER				  
               if ((((cx == 0) && (panels[panelindices[i]]->centroid->x < cube->center->x)) ||
                    ((cx == 1) && (panels[panelindices[i]]->centroid->x >= cube->center->x))) &&
                   (((cy == 0) && (panels[panelindices[i]]->centroid->y < cube->center->y)) ||
                    ((cy == 1) && (panels[panelindices[i]]->centroid->y >= cube->center->y))) &&
                   (((cz == 0) && (panels[panelindices[i]]->centroid->z < cube->center->z)) ||
                    ((cz == 1) && (panels[panelindices[i]]->centroid->z >= cube->center->z))))
					  numchildpanelindices++;
#else
					if (panelIntersectsChildCube(panels[panelindices[i]], cube, cx, cy, cz)) 
					  numchildpanelindices++;
#endif
            }

            for (i = 0; i < numpointindices; i++) {
               if ((((cx == 0) && (points[pointindices[i]]->x < cube->center->x)) ||
                    ((cx == 1) && (points[pointindices[i]]->x >= cube->center->x))) &&
                   (((cy == 0) && (points[pointindices[i]]->y < cube->center->y)) ||
                    ((cy == 1) && (points[pointindices[i]]->y >= cube->center->y))) &&
                   (((cz == 0) && (points[pointindices[i]]->z < cube->center->z)) ||
                    ((cz == 1) && (points[pointindices[i]]->z >= cube->center->z))))
                  numchildpointindices++;
            }

#ifndef SCATTER
            if ((numchildpanelindices == 0) && (numchildpointindices == 0)) {
               cube->children[cx][cy][cz] = NULL;
               continue;
            }
#else
				int shouldBeLeaf = 0;
				if (numchildpanelindices == 0) {
				  //				  printf("level = %d cube should be leaf\n", cube->level+1);
				  shouldBeLeaf = 1;
				}
#endif
            childpanelindices = (unsigned int*)calloc(numchildpanelindices, sizeof(unsigned int));

            count = 0;

            for (i = 0; i < numpanelindices; i++) {
#ifndef SCATTER
               if ((((cx == 0) && (panels[panelindices[i]]->centroid->x < cube->center->x)) ||
                    ((cx == 1) && (panels[panelindices[i]]->centroid->x >= cube->center->x))) &&
                   (((cy == 0) && (panels[panelindices[i]]->centroid->y < cube->center->y)) ||
                    ((cy == 1) && (panels[panelindices[i]]->centroid->y >= cube->center->y))) &&
                   (((cz == 0) && (panels[panelindices[i]]->centroid->z < cube->center->z)) ||
                    ((cz == 1) && (panels[panelindices[i]]->centroid->z >= cube->center->z)))) {
                  childpanelindices[count] = panelindices[i];
                  count++;
               }
#else
					if (panelIntersectsChildCube(panels[panelindices[i]], cube, cx, cy, cz)) {
					  childpanelindices[count] = panelindices[i];
					  count++;
					}
#endif
            }

            childpointindices = (unsigned int*)calloc(numchildpointindices, sizeof(unsigned int));

            count = 0;

            for (i = 0; i < numpointindices; i++) {
               if ((((cx == 0) && (points[pointindices[i]]->x < cube->center->x)) ||
                    ((cx == 1) && (points[pointindices[i]]->x >= cube->center->x))) &&
                   (((cy == 0) && (points[pointindices[i]]->y < cube->center->y)) ||
                    ((cy == 1) && (points[pointindices[i]]->y >= cube->center->y))) &&
                   (((cz == 0) && (points[pointindices[i]]->z < cube->center->z)) ||
                    ((cz == 1) && (points[pointindices[i]]->z >= cube->center->z)))) {
                  childpointindices[count] = pointindices[i];
                  count++;
               }
            }
            
            childbounds[0] = Vector3D_allocate();
            childbounds[1] = Vector3D_allocate();
            
            Vector3D_copy(childbounds[0], cube->bounds[0]);
            Vector3D_copy(childbounds[1], cube->center);
            
            if (cx == 1) {
               real diff = cube->center->x - childbounds[0]->x;
               childbounds[0]->x += diff;
               childbounds[1]->x += diff;
            }
            if (cy == 1) {
               real diff = cube->center->y - childbounds[0]->y;
               childbounds[0]->y += diff;
               childbounds[1]->y += diff;
            }
            if (cz == 1) {
               real diff = cube->center->z - childbounds[0]->z;
               childbounds[0]->z += diff;
               childbounds[1]->z += diff;
            }
            
            childindices[0] = 2 * indices[0]; 
            childindices[1] = 2 * indices[1];
            childindices[2] = 2 * indices[2];
            
            if (cx == 1)
               childindices[0]++;
            if (cy == 1)
               childindices[1]++;
            if (cz == 1)
               childindices[2]++;

#ifndef SCATTER
            cube->children[cx][cy][cz] = Cube_allocate(panels, childpanelindices, numchildpanelindices, points, childpointindices, numchildpointindices, level+1, childindices, childbounds, cube, tree, partitioningdepth);
#else
				if (shouldBeLeaf) {
				  cube->children[cx][cy][cz] = Cube_allocate(panels, childpanelindices, numchildpanelindices, points, childpointindices, numchildpointindices, partitioningdepth-1, childindices, childbounds, cube, tree, partitioningdepth);
				} else {
            cube->children[cx][cy][cz] = Cube_allocate(panels, childpanelindices, numchildpanelindices, points, childpointindices, numchildpointindices, level+1, childindices, childbounds, cube, tree, partitioningdepth);
				}
#endif            

            free(childpanelindices);
            free(childpointindices);
            Vector3D_free(childbounds[0]);
            Vector3D_free(childbounds[1]);
         }

   return cube;
}

void Cube_free(Cube cube) {
   unsigned int i;

   if (cube != NULL) {
      Vector3D_free(cube->bounds[0]);
      Vector3D_free(cube->bounds[1]);
      Vector3D_free(cube->center);
      free(cube->panelindices);
      free(cube->pointindices);
      free(cube->localcubes);
      free(cube->interactingcubes);
      SMatrix_free(cube->VTsrc);
      SMatrix_free(cube->Udest);
#ifdef ADAPTIVE
      SVector_free(cube->VT_q);
      SVector_free(cube->sum_K_VT_q);
#endif
      ComplexSVector_free(cube->FFT_PV_VT_q);
      ComplexSVector_free(cube->sum_T_FFT_PV_VT_q);
      //NOTRANS
		ComplexSVector_free(cube->FFT_UTIT_UT_q);
      ComplexSVector_free(cube->sum_T_FFT_UTIT_UT_q);
      SMatrix_free(cube->UTI);
      SMatrix_free(cube->PV_single);
      SMatrix_free(cube->PV_double);
#ifdef ADAPTIVE
      if (cube->K_single)
         for (i = 0; i < cube->numinteractingcubes; i++)
            if (cube->K_single[i])
               SMatrix_free(cube->K_single[i]);
      if (cube->K_double)
         for (i = 0; i < cube->numinteractingcubes; i++)
            if (cube->K_double[i])
               SMatrix_free(cube->K_double[i]);
#endif

      if (cube->leaf) {
         SMatrix_free(cube->D_single);
         SMatrix_free(cube->D_double);
      } else {
         Cube_free(cube->children[0][0][0]);
         Cube_free(cube->children[0][0][1]);
         Cube_free(cube->children[0][1][0]);
         Cube_free(cube->children[0][1][1]);
         Cube_free(cube->children[1][0][0]);
         Cube_free(cube->children[1][0][1]);
         Cube_free(cube->children[1][1][0]);
         Cube_free(cube->children[1][1][1]);
      }
   }

   free(cube);
}

/* Operations */

void Cube_lists(Cube cube) {
   unsigned int cx, cy, cz, i, count;

   /* Initialize */
   
   cube->numlocalcubes = 0;
   cube->localcubes = NULL;
   cube->numinteractingcubes = 0;
   cube->interactingcubes = NULL;

   /* Siblings are guaranteed to be local cubes */
   
   if (cube->parent != NULL)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->parent->children[cx][cy][cz] != NULL)
                  if (cube->parent->children[cx][cy][cz] != cube)
                     cube->numlocalcubes++;
   
   /* Other local cubes are children of the parent's local cubes that are
      within the LOCAL cutoff */
   
   if (cube->parent != NULL)
      for (i = 0; i < cube->parent->numlocalcubes; i++)
         for (cx = 0; cx <= 1; cx++)
            for (cy = 0; cy <= 1; cy++)
               for (cz = 0; cz <= 1; cz++)
                  if (cube->parent->localcubes[i]->children[cx][cy][cz] != NULL)
                     if (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[0] - cube->indices[0]) <= LOCAL)
                        if (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[1] - cube->indices[1]) <= LOCAL)
                           if (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[2] - cube->indices[2]) <= LOCAL)
                              cube->numlocalcubes++;

   /* Allocate */

   if (cube->numlocalcubes > 0)
      cube->localcubes = (Cube*)calloc(cube->numlocalcubes, sizeof(Cube*));

   count = 0;

   /* Siblings are guaranteed to be local cubes */
   
   if (cube->parent != NULL)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->parent->children[cx][cy][cz] != NULL)
                  if (cube->parent->children[cx][cy][cz] != cube) {
                     cube->localcubes[count] = cube->parent->children[cx][cy][cz];
                     count++;
                  }
   
   /* Other local cubes are children of the parent's local cubes that are
      within the LOCAL cutoff */
   
   if (cube->parent != NULL)
      for (i = 0; i < cube->parent->numlocalcubes; i++)
         for (cx = 0; cx <= 1; cx++)
            for (cy = 0; cy <= 1; cy++)
               for (cz = 0; cz <= 1; cz++)
                  if (cube->parent->localcubes[i]->children[cx][cy][cz] != NULL)
                     if (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[0] - cube->indices[0]) <= LOCAL)
                        if (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[1] - cube->indices[1]) <= LOCAL)
                           if (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[2] - cube->indices[2]) <= LOCAL) {
                              cube->localcubes[count] = cube->parent->localcubes[i]->children[cx][cy][cz];
                              count++;
                           }

   /* Interacting cubes are children of the parent's local cubes that are
      outside the LOCAL cutoff and have at least one point in them */
   
   if (cube->parent != NULL)
      for (i = 0; i < cube->parent->numlocalcubes; i++)
         for (cx = 0; cx <= 1; cx++)
            for (cy = 0; cy <= 1; cy++)
               for (cz = 0; cz <= 1; cz++)
                  if (cube->parent->localcubes[i]->children[cx][cy][cz] != NULL)
                     if ((abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[0] - cube->indices[0]) > LOCAL) ||
                         (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[1] - cube->indices[1]) > LOCAL) ||
                         (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[2] - cube->indices[2]) > LOCAL))
                        cube->numinteractingcubes++;

   /* Allocate */

   if (cube->numinteractingcubes > 0)
      cube->interactingcubes = (Cube*)calloc(cube->numinteractingcubes, sizeof(Cube*));

   count = 0;

   /* Interacting cubes are children of the parent's local cubes that are
      outside the LOCAL cutoff */
   
   if (cube->parent != NULL)
      for (i = 0; i < cube->parent->numlocalcubes; i++)
         for (cx = 0; cx <= 1; cx++)
            for (cy = 0; cy <= 1; cy++)
               for (cz = 0; cz <= 1; cz++)
                  if (cube->parent->localcubes[i]->children[cx][cy][cz] != NULL)
                     if ((abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[0] - cube->indices[0]) > LOCAL) ||
                         (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[1] - cube->indices[1]) > LOCAL) ||
                         (abs(cube->parent->localcubes[i]->children[cx][cy][cz]->indices[2] - cube->indices[2]) > LOCAL)) {
                        cube->interactingcubes[count] = cube->parent->localcubes[i]->children[cx][cy][cz];
                        count++;
                     }

   /* Recurse over children */

   for (cx = 0; cx <= 1; cx++)
      for (cy = 0; cy <= 1; cy++)
         for (cz = 0; cz <= 1; cz++)
            if (cube->children[cx][cy][cz] != NULL)
               Cube_lists(cube->children[cx][cy][cz]);
}

void Cube_fill_D(Cube cube) {
   unsigned int i, j, r, c;
   unsigned int cx, cy, cz;
   real slp, dlp;   
   BEMLayerType layertype = cube->tree->layertype;
   BEMKernelType kerneltype = cube->tree->kerneltype;
   void* parameters = cube->tree->parameters;
   Vector3D* points = cube->tree->points;
	Vector3D* normals = cube->tree->normals;
   Panel* panels = cube->tree->panels;

   /* If this cube is a leaf, we do self+local directly */

   if (cube->leaf) {
      cube->dcolumns = cube->numpanelindices;
      cube->drows = cube->numpointindices;
   
      for (i = 0; i < cube->numlocalcubes; i++)
         cube->dcolumns += cube->localcubes[i]->numpanelindices;

      if ((cube->drows > 0) && (cube->dcolumns > 0)) {
         if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
            cube->D_single = SMatrix_allocate(cube->drows, cube->dcolumns);
         if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
            cube->D_double = SMatrix_allocate(cube->drows, cube->dcolumns);
			if (layertype == NORMDERIV_SINGLE_LAYER_INT)
			  cube->D_single =  SMatrix_allocate(cube->drows, cube->dcolumns);

         c = 0;
      
         for (i = 0; i < cube->numpanelindices; i++) {
            for (r = 0; r < cube->numpointindices; r++) {
               if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
                  cube->D_single[r][c] = Integration(points[cube->pointindices[r]], panels[cube->panelindices[c]], kerneltype, parameters, SINGLE_LAYER_INT);
               if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
                  cube->D_double[r][c] = Integration(points[cube->pointindices[r]], panels[cube->panelindices[c]], kerneltype, parameters, DOUBLE_LAYER_INT);
               if (layertype == NORMDERIV_SINGLE_LAYER_INT)
                  cube->D_single[r][c] = Integration(points[cube->pointindices[r]], panels[cube->panelindices[c]], kerneltype, (void *)(normals[cube->pointindices[r]]), NORMDERIV_SINGLE_LAYER_INT);

            }
            c++;
         }
         for (i = 0; i < cube->numlocalcubes; i++)
            for (j = 0; j < cube->localcubes[i]->numpanelindices; j++) {
               for (r = 0; r < cube->numpointindices; r++) {
                  if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
                     cube->D_single[r][c] = Integration(points[cube->pointindices[r]], panels[cube->localcubes[i]->panelindices[j]], kerneltype, parameters, SINGLE_LAYER_INT);
                  if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
                     cube->D_double[r][c] = Integration(points[cube->pointindices[r]], panels[cube->localcubes[i]->panelindices[j]], kerneltype, parameters, DOUBLE_LAYER_INT);
						if (layertype == NORMDERIV_SINGLE_LAYER_INT)
                     cube->D_single[r][c] = Integration(points[cube->pointindices[r]], panels[cube->localcubes[i]->panelindices[j]], kerneltype, (void *)(normals[cube->pointindices[r]]), NORMDERIV_SINGLE_LAYER_INT);
               }
               c++;
            }
#ifdef SERIALIZE
         if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
            Matrix_writebinary(cube->tree->Dfiles_single[cube->serialid], cube->D_single, cube->drows, cube->dcolumns);
            Matrix_free(cube->D_single);
            cube->D_single = NULL;
         }
         if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
            Matrix_writebinary(cube->tree->Dfiles_double[cube->serialid], cube->D_double, cube->drows, cube->dcolumns);
            Matrix_free(cube->D_double);
            cube->D_double = NULL;
         }
			if (layertype == NORMDERIV_SINGLE_LAYER_INT) {
			  printf("SERIALIZATION IS NOT SET UP FOR NORMDERIV_SINGLE_LAYER_INT!\nDying.\n");
			  exit(-1);
			}
#endif
      }
   }
   else /* Recurse over childern */
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_fill_D(cube->children[cx][cy][cz]);
}

void Cube_fill_P_I(Cube cube, unsigned int recurse) {
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];
   unsigned int count = 0;
   unsigned int i, j, k;
   real boxlength = cube->bounds[1]->x - cube->bounds[0]->x;
   BEMLayerType layertype = cube->tree->layertype;
   real (*quadraturepoints)[4] = cube->tree->quadraturepoints;
   unsigned int numquadraturepoints = cube->tree->numquadraturepoints;
   unsigned int numcoeffs = gridpoints[cube->level]*(gridpoints[cube->level]+1)*(gridpoints[cube->level]+2)/6;
   BEMKernelType kerneltype = cube->tree->kerneltype;
   void* parameters = cube->tree->parameters;
   Panel* panels = cube->tree->panels;
   Vector3D* points = cube->tree->points;
#ifdef POLYNOMIAL
   Matrix* pinv_Fgridpoints = cube->tree->pinv_Fgridpoints;
#else
   Matrix* pinv_Q2Sgridpoints = cube->tree->pinv_Q2Sgridpoints;
#endif

   /* If there is nothing on the interation list, simply recurse over
      the children (or quit if a leaf) */

   if (cube->numinteractingcubes == 0) {
      if (recurse)
         if (!cube->leaf)
            for (cx = 0; cx <= 1; cx++)
               for (cy = 0; cy <= 1; cy++)
                  for (cz = 0; cz <= 1; cz++)
                     if (cube->children[cx][cy][cz] != NULL)
                        Cube_fill_P_I(cube->children[cx][cy][cz], recurse);

      return;
   }

   /* P = Projection (gridpoints^3 * numpanels) */
   /* P only exists if there are panels in this cube */
   if (cube->numpanelindices > 0) {
      if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT) || (layertype == NORMDERIV_SINGLE_LAYER_INT)) {
         cube->P_single = SMatrix_allocate(gp3, cube->numpanelindices);

#ifdef POLYNOMIAL
         Matrix Fpanels = Polynomial_F_panels(cube, boxlength, panels,
                                              gridpoints[cube->level],
                                              SINGLE_LAYER_INT);

         Matrix P_single = Matrix_allocate(gp3, cube->numpanelindices);

         Matrix_multiplymatrix(P_single, pinv_Fgridpoints[cube->level], Fpanels,
                               gp3, numcoeffs, cube->numpanelindices);

         SMatrix_Matrix(cube->P_single, P_single, gp3, cube->numpanelindices);
         
         Matrix_free(P_single);

         Matrix_free(Fpanels);
#else
         Matrix Q2Spanels = EquivDensity_Q2S_panels(cube, boxlength, panels,
                                                    quadraturepoints, 
                                                    numquadraturepoints, 
                                                    kerneltype, parameters,
                                                    SINGLE_LAYER_INT);

         Matrix P_single = Matrix_allocate(gp3, cube->numpanelindices);

         Matrix_multiplymatrix(P_single, pinv_Q2Sgridpoints[cube->level], Q2Spanels, 
                               gp3, numquadraturepoints, cube->numpanelindices);

         SMatrix_Matrix(cube->P_single, P_single, gp3, cube->numpanelindices);
         
         Matrix_free(P_single);
         Matrix_free(Q2Spanels);
#endif
      }

      if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
         cube->P_double = SMatrix_allocate(gp3, cube->numpanelindices);

#ifdef POLYNOMIAL
         Matrix Fpanels = Polynomial_F_panels(cube, boxlength, panels,
                                              gridpoints[cube->level],
                                              DOUBLE_LAYER_INT);

         Matrix P_double = Matrix_allocate(gp3, cube->numpanelindices);

         Matrix_multiplymatrix(P_double, pinv_Fgridpoints[cube->level], Fpanels,
                               gp3, numcoeffs, cube->numpanelindices);

         SMatrix_Matrix(cube->P_double, P_double, gp3, cube->numpanelindices);

         Matrix_fre(P_double);
         Matrix_free(Fpanels);
#else
         Matrix Q2Spanels = EquivDensity_Q2S_panels(cube, boxlength, panels,
                                                    quadraturepoints, 
                                                    numquadraturepoints, 
                                                    kerneltype, parameters,
                                                    DOUBLE_LAYER_INT);

         Matrix P_double = Matrix_allocate(gp3, cube->numpanelindices);

         Matrix_multiplymatrix(P_double, pinv_Q2Sgridpoints[cube->level], Q2Spanels, gp3, numquadraturepoints, cube->numpanelindices);

         SMatrix_Matrix(cube->P_double, P_double, gp3, cube->numpanelindices);

         Matrix_free(P_double);
         Matrix_free(Q2Spanels);
#endif
      }
   }

   /* I = Interpolation (numpoints * gridpoints^3) */
   /* I only exists if there are points in this cube */

   if (cube->numpointindices > 0) {
      Vector3D* evalpoints = (Vector3D*)calloc(cube->numpointindices, sizeof(Vector3D));

      for (i = 0; i < cube->numpointindices; i++) {
         evalpoints[i] = Vector3D_allocate();
         Vector3D_copy(evalpoints[i], points[cube->pointindices[i]]);
      }

#ifdef POLYNOMIAL
      Matrix Fevalpoints = Polynomial_F_points(cube->center, boxlength,
                                        evalpoints, cube->numpointindices,
                                        gridpoints[cube->level]);
#else
      Matrix Q2Sevalpoints = EquivDensity_Q2S_points(cube->center, boxlength,
                                              evalpoints, cube->numpointindices, 
                                              quadraturepoints,
                                              numquadraturepoints, 
                                              kerneltype,
                                              parameters);
#endif

 
      cube->I_ = SMatrix_allocate(gp3, cube->numpointindices);

#ifdef POLYNOMIAL
      Matrix I_ = Matrix_allocate(gp3, cube->numpointindices);

      Matrix_multiplymatrix(I_, pinv_Fgridpoints[cube->level], Fevalpoints, gp3, numcoeffs, cube->numpointindices);
     
      SMatrix_Matrix(cube->I_, I_, gp3, cube->numpointindices);

      Matrix_free(I_);
      Matrix_free(Fevalpoints);
#else     
      Matrix I_ = Matrix_allocate(gp3, cube->numpointindices);

      Matrix_multiplymatrix(I_, pinv_Q2Sgridpoints[cube->level], Q2Sevalpoints, gp3, numquadraturepoints, cube->numpointindices);
     
      SMatrix_Matrix(cube->I_, I_, gp3, cube->numpointindices);

      Matrix_free(I_);
      Matrix_free(Q2Sevalpoints);
#endif
     
      SMatrix_transpose(&cube->I_, gp3, cube->numpointindices);

		if (cube->tree->layertype == NORMDERIV_SINGLE_LAYER_INT) {
#ifdef POLYNOMIAL
		  printf("NORMDERIV_SINGLE_LAYER_INT is not set up for POLYNOMIAL!\n");
#endif
		  Vector3D* evalnormals = (Vector3D*)calloc(cube->numpointindices, sizeof(Vector3D));
		  Vector3D* normals = cube->tree->normals;
		  for (i = 0; i < cube->numpointindices; i++) {
			 evalnormals[i] = Vector3D_allocate();
			 Vector3D_copy(evalnormals[i], normals[cube->pointindices[i]]);
		  }
		  
		  cube->I_double = SMatrix_allocate(gp3, cube->numpointindices);
		  Matrix I_double = Matrix_allocate(gp3, cube->numpointindices);
		  Matrix Q2Sevalpoints_double = EquivDensity_Q2S_points_deriv(cube->center, boxlength,
																						  evalpoints, evalnormals, cube->numpointindices, 
																						  quadraturepoints,
																						  numquadraturepoints, 
																						  kerneltype,
																						  parameters);
		  
		  Matrix_multiplymatrix(I_double, pinv_Q2Sgridpoints[cube->level], Q2Sevalpoints_double, gp3, numquadraturepoints, cube->numpointindices);
		  
		  SMatrix_Matrix(cube->I_double, I_double, gp3, cube->numpointindices);
		  Matrix_free(I_double);
		  Matrix_free(Q2Sevalpoints_double);
		  SMatrix_transpose(&cube->I_double, gp3, cube->numpointindices);
		  for (i = 0; i < cube->numpointindices; i++)
			 Vector3D_free(evalnormals[i]);
		  free(evalnormals);
		}
		
		for (i = 0; i < cube->numpointindices; i++)
		  Vector3D_free(evalpoints[i]);
		free(evalpoints);
		
   }

   /* Recurse over children */

   if (recurse)
      if (!cube->leaf)
         for (cx = 0; cx <= 1; cx++)
            for (cy = 0; cy <= 1; cy++)
               for (cz = 0; cz <= 1; cz++)
                  if (cube->children[cx][cy][cz] != NULL)
                     Cube_fill_P_I(cube->children[cx][cy][cz], recurse);
}

void Sampling_multiply_PT_IFFT_T_FFT_IT(SVector ans, Cube pointcube, Cube cube, SVector x, BEMLayerType layertype) {
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];
   unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
   SVector tempgp3 = SVector_allocate(gp3);
   ComplexSVector temppadgridsize = ComplexSVector_allocate(padgridsize);
   ComplexSVector temppadgridsize2 = ComplexSVector_allocate(padgridsize);
   unsigned int i;
   ComplexSVector**** Tprecomputed = cube->tree->Tprecomputed;

   // IT

   SMatrix_multiplyvector_transpose(tempgp3, pointcube->I_, x, pointcube->numpointindices, gp3);

   // FFT

   FFT_forwardGridTransform(cube->level, gridpoints[cube->level], tempgp3, temppadgridsize);

   // T

   unsigned int halfdimension = 1 + 2 * LOCAL;
   int dx = cube->indices[0] - pointcube->indices[0];
   int dy = cube->indices[1] - pointcube->indices[1];
   int dz = cube->indices[2] - pointcube->indices[2];

   ComplexSVector T = Tprecomputed[cube->level]
      [dx+halfdimension]
      [dy+halfdimension]
      [dz+halfdimension];

   ComplexSVector_addelementmultiplyvector(temppadgridsize2, T, temppadgridsize, padgridsize);

   // IFFT

   SVector_zero(tempgp3, gp3);

   FFT_backwardGridTransform(cube->level, gridpoints[cube->level], temppadgridsize2, tempgp3);

   // PT

   if (layertype == SINGLE_LAYER_INT)
      SMatrix_multiplyvector_transpose(ans, cube->P_single, tempgp3, gp3, cube->numpanelindices);
   else if (layertype == DOUBLE_LAYER_INT)
      SMatrix_multiplyvector_transpose(ans, cube->P_double, tempgp3, gp3, cube->numpanelindices);

   ComplexSVector_free(temppadgridsize2);
   ComplexSVector_free(temppadgridsize);
   SVector_free(tempgp3);
}

void Sampling_multiply_I_IFFT_T_FFT_P(SVector ans, Cube cube, Cube panelcube, SVector x, BEMLayerType layertype) {
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];
   unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
   SVector tempgp3 = SVector_allocate(gp3);
   ComplexSVector temppadgridsize = ComplexSVector_allocate(padgridsize);
   ComplexSVector temppadgridsize2 = ComplexSVector_allocate(padgridsize);
   unsigned int i;
   ComplexSVector**** Tprecomputed = cube->tree->Tprecomputed;

   // P

   if (layertype == SINGLE_LAYER_INT)
      SMatrix_multiplyvector(tempgp3, panelcube->P_single, x, gp3, panelcube->numpanelindices);
   else if (layertype == DOUBLE_LAYER_INT)
      SMatrix_multiplyvector(tempgp3, panelcube->P_double, x, gp3, panelcube->numpanelindices);

   // FFT

   FFT_forwardGridTransform(cube->level, gridpoints[cube->level], tempgp3, temppadgridsize);

   // T

   unsigned int halfdimension = 1 + 2 * LOCAL;
   int dx = cube->indices[0] - panelcube->indices[0];
   int dy = cube->indices[1] - panelcube->indices[1];
   int dz = cube->indices[2] - panelcube->indices[2];

   ComplexSVector T = Tprecomputed[cube->level]
      [dx+halfdimension]
      [dy+halfdimension]
      [dz+halfdimension];

   ComplexSVector_addelementmultiplyvector(temppadgridsize2, T, temppadgridsize, padgridsize);

   // IFFT

   SVector_zero(tempgp3, gp3);

   FFT_backwardGridTransform(cube->level, gridpoints[cube->level], temppadgridsize2, tempgp3);

   // I

   SMatrix_multiplyvector(ans, cube->I_, tempgp3, cube->numpointindices, gp3);

   ComplexSVector_free(temppadgridsize2);
   ComplexSVector_free(temppadgridsize);
   SVector_free(tempgp3);
}

void Cube_fill_U_VT(Cube cube) {
   unsigned int i, p, pass, count;
   unsigned int cx, cy, cz;
   unsigned int r, c;
   real slp, dlp;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
   unsigned int* interactingcubeswithpanels;
   unsigned int numinteractingcubeswithpanels = 0;
   unsigned int* interactingcubeswithpoints;
   unsigned int numinteractingcubeswithpoints = 0;
   BEMLayerType layertype = cube->tree->layertype;
   BEMKernelType kerneltype = cube->tree->kerneltype;
   void* parameters = cube->tree->parameters;
   Panel* panels = cube->tree->panels;
   Vector3D* points = cube->tree->points;

	Vector3D* normals = NULL;
	if (cube->tree->layertype == NORMDERIV_SINGLE_LAYER_INT) 
	  normals = cube->tree->normals;

   real epsilon = cube->tree->epsilon;

   /* Since a cube can contain points, panels or both, in order to do sampling
      properly, we must determine which interacting cubes have what in them */

   for (i = 0; i < cube->numinteractingcubes; i++)
      if (cube->interactingcubes[i]->numpanelindices > 0)
         numinteractingcubeswithpanels++;

   interactingcubeswithpanels = (unsigned int*)calloc(numinteractingcubeswithpanels, sizeof(unsigned int));
   
   count = 0;
   
   for (i = 0; i < cube->numinteractingcubes; i++)
      if (cube->interactingcubes[i]->numpanelindices > 0) {
         interactingcubeswithpanels[count] = i;
         count++;
      }

   for (i = 0; i < cube->numinteractingcubes; i++)
      if (cube->interactingcubes[i]->numpointindices > 0)
         numinteractingcubeswithpoints++;

   interactingcubeswithpoints = (unsigned int*)calloc(numinteractingcubeswithpoints, sizeof(unsigned int));
   
   count = 0;
   
   for (i = 0; i < cube->numinteractingcubes; i++)
      if (cube->interactingcubes[i]->numpointindices > 0) {
         interactingcubeswithpoints[count] = i;
         count++;
      }

   /* Sampling version
   
   Here is the sampling strategy:

   For VTsrc, we will generate interactingcubeswithpoints full rows by
   evaluating every source panel at collocation points in each
   interacting cube with at least one point.

   For Udest, we will generate interactingcubeswithpanels full columns by
   evaluating every source point at panels in each interacting cube
   with at least one panel.
   */

   if ((cube->numpanelindices > 0) && (numinteractingcubeswithpoints > 0)) {
      for (pass = 1; ; pass++) {
         Matrix R;
      
         if (layertype == SINGLE_AND_DOUBLE_LAYER_INT)
            R = Matrix_allocate(2 * numinteractingcubeswithpoints * pass, cube->numpanelindices);
         else
            R = Matrix_allocate(numinteractingcubeswithpoints * pass, cube->numpanelindices);

#ifdef ACCELERATED_SAMPLING
         /* Generate the rows via FFT acceleration */
         /* NOTE: NOT HANDLING OUTER NORMAL DERIVATIVE! */

         for (p = 0; p < pass; p++)
            for (r = 0; r < numinteractingcubeswithpoints; r++) {
               Cube pointcube = cube->interactingcubes[interactingcubeswithpoints[r]];
               unsigned int pointindex = p * pointcube->numpointindices / pass;

               Vector pointvector = Vector_allocate(pointcube->numpointindices);
               Vector samplerow = Vector_allocate(cube->numpanelindices);
               pointvector[pointindex] = 1.0;

               // Accelerated multiply to get row

               if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
                  Sampling_multiply_PT_IFFT_T_FFT_IT(samplerow, pointcube, cube, pointvector, SINGLE_LAYER_INT);
                  for (c = 0; c < cube->numpanelindices; c++)
                     R[p*numinteractingcubeswithpoints+r][c] = samplerow[c];
                  if (layertype == SINGLE_AND_DOUBLE_LAYER_INT) {
                     Vector_zero(samplerow, cube->numpanelindices);
                     Sampling_multiply_PT_IFFT_T_FFT_IT(samplerow, pointcube, cube, pointvector, DOUBLE_LAYER_INT);
                     for (c = 0; c < cube->numpanelindices; c++)
                        R[p*numinteractingcubeswithpoints+r+numinteractingcubeswithpoints*pass][c] = samplerow[c];
                  }
               }
               else if (layertype == DOUBLE_LAYER_INT) {
                  Sampling_multiply_PT_IFFT_T_FFT_IT(samplerow, pointcube, cube, pointvector, DOUBLE_LAYER_INT);
                  for (c = 0; c < cube->numpanelindices; c++)
                     R[p*numinteractingcubeswithpoints+r][c] = samplerow[c];
               }

               Vector_free(pointvector);
               Vector_free(samplerow);
            }
#else
         /* Generate the rows via integration */

         for (p = 0; p < pass; p++)
            for (r = 0; r < numinteractingcubeswithpoints; r++) {
               Cube pointcube = cube->interactingcubes[interactingcubeswithpoints[r]];
               unsigned int pointindex = p * pointcube->numpointindices / pass;

               if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
                  for (c = 0; c < cube->numpanelindices; c++)
                     R[p*numinteractingcubeswithpoints+r][c] = Integration(points[pointcube->pointindices[pointindex]], panels[cube->panelindices[c]], kerneltype, parameters, SINGLE_LAYER_INT);
                  if (layertype == SINGLE_AND_DOUBLE_LAYER_INT) {
                     for (c = 0; c < cube->numpanelindices; c++)
                        R[p*numinteractingcubeswithpoints+r+numinteractingcubeswithpoints*pass][c] = Integration(points[pointcube->pointindices[pointindex]], panels[cube->panelindices[c]], kerneltype, parameters, DOUBLE_LAYER_INT);
                  }
               }
               else if (layertype == DOUBLE_LAYER_INT) {
                  for (c = 0; c < cube->numpanelindices; c++)
                     R[p*numinteractingcubeswithpoints+r][c] = Integration(points[pointcube->pointindices[pointindex]], panels[cube->panelindices[c]], kerneltype, parameters, DOUBLE_LAYER_INT);
               }
					else if (layertype == NORMDERIV_SINGLE_LAYER_INT) {
					  for (c = 0; c < cube->numpanelindices; c++)
						 R[p*numinteractingcubeswithpoints+r][c] = Integration(points[pointcube->pointindices[pointindex]], panels[cube->panelindices[c]], kerneltype, (void *)(normals[pointcube->pointindices[pointindex]]), NORMDERIV_SINGLE_LAYER_INT);
					}
					
            }

#endif

         /* Once we have R, we construct a basis for these to make V */

         Matrix VTsrc;

         if (layertype == SINGLE_AND_DOUBLE_LAYER_INT)
            Matrix_rowbasis_pivotedmgs(&VTsrc, &cube->rowrank, R, 2 * numinteractingcubeswithpoints * pass, cube->numpanelindices, epsilon);
         else
            Matrix_rowbasis_pivotedmgs(&VTsrc, &cube->rowrank, R, numinteractingcubeswithpoints * pass, cube->numpanelindices, epsilon);
      
         Matrix_free(R);

         /* If the sampling was good, exit the loop, otherwise delete and resample */

         if (layertype == SINGLE_AND_DOUBLE_LAYER_INT) {
            if (cube->rowrank < (2*pass*numinteractingcubeswithpoints)/2) {
               cube->VTsrc = SMatrix_allocate(cube->rowrank, cube->numpanelindices);
               SMatrix_Matrix(cube->VTsrc, VTsrc, cube->rowrank, cube->numpanelindices);
               Matrix_free(VTsrc);
               break;
            }
         } else  {
            if (cube->rowrank < (pass*numinteractingcubeswithpoints)/2) {
               cube->VTsrc = SMatrix_allocate(cube->rowrank, cube->numpanelindices);
               SMatrix_Matrix(cube->VTsrc, VTsrc, cube->rowrank, cube->numpanelindices);
               Matrix_free(VTsrc);
               break;
            }
         }

         Matrix_free(VTsrc);
      }

#ifdef SERIALIZE
      SMatrix_writebinary(cube->tree->VTfiles[cube->serialid], cube->VTsrc, cube->rowrank, cube->numpanelindices);
      SMatrix_free(cube->VTsrc); 
      cube->VTsrc = NULL;
#endif
   }
   
   if ((cube->numpointindices > 0) && (numinteractingcubeswithpanels > 0)) {
      for (pass = 1; ; pass++) {
         Matrix C;
      
         if (layertype == SINGLE_AND_DOUBLE_LAYER_INT)
            C = Matrix_allocate(cube->numpointindices, 2 * numinteractingcubeswithpanels * pass);
         else
            C = Matrix_allocate(cube->numpointindices, numinteractingcubeswithpanels * pass);

#ifdef ACCELERATED_SAMPLING
         /* Generate the columns via FFT acceleration */

         for (p = 0; p < pass; p++)
            for (c = 0; c < numinteractingcubeswithpanels; c++) {
               Cube panelcube = cube->interactingcubes[interactingcubeswithpanels[c]];
               unsigned int panelindex = p * panelcube->numpanelindices / pass;

               Vector panelvector = Vector_allocate(panelcube->numpanelindices);
               Vector samplecolumn = Vector_allocate(cube->numpointindices);
               panelvector[panelindex] = 1.0;

               // Accelerated multiply to get column

               if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
                  Sampling_multiply_I_IFFT_T_FFT_P(samplecolumn, cube, panelcube, panelvector, SINGLE_LAYER_INT);
                  for (r = 0; r < cube->numpointindices; r++)
                     C[r][p*numinteractingcubeswithpanels+c] = samplecolumn[r];
                  if (layertype == SINGLE_AND_DOUBLE_LAYER_INT) {
                     Vector_zero(samplecolumn, cube->numpointindices);
                     Sampling_multiply_I_IFFT_T_FFT_P(samplecolumn, cube, panelcube, panelvector, DOUBLE_LAYER_INT);
                     for (r = 0; r < cube->numpointindices; r++)
                        C[r][p*numinteractingcubeswithpanels+c+numinteractingcubeswithpanels*pass] = samplecolumn[r];
                  }
               }
               else if (layertype == DOUBLE_LAYER_INT) {
                  Sampling_multiply_I_IFFT_T_FFT_P(samplecolumn, cube, panelcube, panelvector, DOUBLE_LAYER_INT);
                  for (r = 0; r < cube->numpointindices; r++)
                     C[r][p*numinteractingcubeswithpanels+c] = samplecolumn[r];
               }

               Vector_free(panelvector);
               Vector_free(samplecolumn);
            }
#else
         /* Generate the columns via integration */

         for (p = 0; p < pass; p++)
            for (c = 0; c < numinteractingcubeswithpanels; c++) {
               Cube panelcube = cube->interactingcubes[interactingcubeswithpanels[c]];
               unsigned int panelindex = p * panelcube->numpanelindices / pass;

               if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
                  for (r = 0; r < cube->numpointindices; r++)
                     C[r][p*numinteractingcubeswithpanels+c] = Integration(points[cube->pointindices[r]], panels[panelcube->panelindices[panelindex]], kerneltype, parameters, SINGLE_LAYER_INT);
                  if (layertype == SINGLE_AND_DOUBLE_LAYER_INT) {
                     for (r = 0; r < cube->numpointindices; r++)
                        C[r][p*numinteractingcubeswithpanels+c+numinteractingcubeswithpanels*pass] = Integration(points[cube->pointindices[r]], panels[panelcube->panelindices[panelindex]], kerneltype, parameters, DOUBLE_LAYER_INT);
                  }
               }
               else if (layertype == DOUBLE_LAYER_INT) {
					  for (r = 0; r < cube->numpointindices; r++)
                     C[r][p*numinteractingcubeswithpanels+c] = Integration(points[cube->pointindices[r]], panels[panelcube->panelindices[panelindex]], kerneltype, parameters, DOUBLE_LAYER_INT);
               }
					else if (layertype == NORMDERIV_SINGLE_LAYER_INT) {
					  C[r][p*numinteractingcubeswithpanels+c] = Integration(points[cube->pointindices[r]], panels[panelcube->panelindices[panelindex]], kerneltype, (void *)(normals[cube->pointindices[r]]), NORMDERIV_SINGLE_LAYER_INT);
					}
            }
#endif

         /* Once we have C, we construct a basis for these to make U */

         Matrix Udest;

         if (layertype == SINGLE_AND_DOUBLE_LAYER_INT)
            Matrix_columnbasis_pivotedmgs(&Udest, &cube->columnrank, C, cube->numpointindices, 2 * numinteractingcubeswithpanels * pass, epsilon);
         else
            Matrix_columnbasis_pivotedmgs(&Udest, &cube->columnrank, C, cube->numpointindices, numinteractingcubeswithpanels * pass, epsilon);
      
         Matrix_free(C);

         /* If the sampling was good, exit the loop, otherwise delete and resample */

         if (layertype == SINGLE_AND_DOUBLE_LAYER_INT) {
            if (cube->columnrank < (2*pass*numinteractingcubeswithpanels)/2) {
               cube->Udest = SMatrix_allocate(cube->numpointindices, cube->columnrank);
               SMatrix_Matrix(cube->Udest, Udest, cube->numpointindices, cube->columnrank);
               Matrix_free(Udest);
               break;
            }
         } else {
            if (cube->columnrank < (pass*numinteractingcubeswithpanels)/2) {
               cube->Udest = SMatrix_allocate(cube->numpointindices, cube->columnrank);
               SMatrix_Matrix(cube->Udest, Udest, cube->numpointindices, cube->columnrank);
               Matrix_free(Udest);
               break;
            }
         }
      
         Matrix_free(Udest);
      }

#ifdef SERIALIZE
      Matrix_writebinary(cube->tree->Ufiles[cube->serialid], cube->Udest, cube->numpointindices, cube->columnrank);
      Matrix_free(cube->Udest); 
      cube->Udest = NULL;
#endif
   }
   
   free(interactingcubeswithpanels);
   free(interactingcubeswithpoints);

#ifdef ADAPTIVE
   cube->VT_q = SVector_allocate(cube->rowrank);
   cube->sum_K_VT_q = SVector_allocate(cube->columnrank);
#endif
   cube->FFT_PV_VT_q = ComplexSVector_allocate(padgridsize);
   cube->sum_T_FFT_PV_VT_q = ComplexSVector_allocate(padgridsize);
   //NOTRANS
	cube->FFT_UTIT_UT_q = ComplexSVector_allocate(padgridsize);
   cube->sum_T_FFT_UTIT_UT_q = ComplexSVector_allocate(padgridsize);

   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_fill_U_VT(cube->children[cx][cy][cz]);
}

void Cube_fill_PV_UTI(Cube cube) {
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];
   unsigned int count = 0;
   unsigned int i, j, k;
   real boxlength = cube->bounds[1]->x - cube->bounds[0]->x;
   BEMLayerType layertype = cube->tree->layertype;

   /* If there is nothing on the interation list, simply recurse over
      the children (or quit if a leaf) */

   if (cube->numinteractingcubes == 0) {
      if (!cube->leaf)
         for (cx = 0; cx <= 1; cx++)
            for (cy = 0; cy <= 1; cy++)
               for (cz = 0; cz <= 1; cz++)
                  if (cube->children[cx][cy][cz] != NULL)
                     Cube_fill_PV_UTI(cube->children[cx][cy][cz]);

      return;
   }

#ifndef ACCELERATED_SAMPLING
   Cube_fill_P_I(cube, 0);
#endif

   /* PV = Projection * V (gridpoints^3 * rowrank) */
   /* PV only exists if there are panels in this cube */
   if ((cube->numpanelindices > 0) && (cube->rowrank > 0)) {
#ifdef SERIALIZE
      cube->VTsrc = SMatrix_allocate(cube->rowrank, cube->numpanelindices);
      SMatrix_readbinary(cube->VTsrc, cube->tree->VTfiles[cube->serialid], cube->rowrank, cube->numpanelindices);
#endif

      if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
         cube->PV_single = SMatrix_allocate(gp3, cube->rowrank);

         SMatrix_multiplytranspose_matrix(cube->PV_single, cube->P_single, cube->VTsrc, gp3, cube->numpanelindices, cube->rowrank);

         SMatrix_free(cube->P_single);
         cube->P_single = NULL;
      }

      if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
         cube->PV_double = SMatrix_allocate(gp3, cube->rowrank);

         SMatrix_multiplytranspose_matrix(cube->PV_double, cube->P_double, cube->VTsrc, gp3, cube->numpanelindices, cube->rowrank);

         SMatrix_free(cube->P_double);
         cube->P_double = NULL;
      }

#ifdef SERIALIZE
      SMatrix_free(cube->VTsrc);
      cube->VTsrc = NULL;
#endif

#ifdef SERIALIZE
      if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
         SMatrix_writebinary(cube->tree->PVfiles_single[cube->serialid], cube->PV_single, gp3, cube->rowrank);
         SMatrix_free(cube->PV_single);
         cube->PV_single = NULL;
      }
      if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
         SMatrix_writebinary(cube->tree->PVfiles_double[cube->serialid], cube->PV_double, gp3, cube->rowrank);
         SMatrix_free(cube->PV_double);
         cube->PV_double = NULL;
      }
#endif
   }

   /* UTI = UT * Interpolation (columnrank * gridpoints^3) */
   /* UTI only exists if there are points in this cube */

   if ((cube->numpointindices > 0) && (cube->columnrank > 0)) {
#ifdef SERIALIZE
      cube->Udest = SMatrix_allocate(cube->numpointindices, cube->columnrank);
      SMatrix_readbinary(cube->Udest, cube->tree->Ufiles[cube->serialid], cube->numpointindices, cube->columnrank);
#endif

      cube->UTI = SMatrix_allocate(cube->columnrank, gp3);

      SMatrix_multiplymatrix_transpose(cube->UTI, cube->Udest, cube->I_, cube->numpointindices, cube->columnrank, gp3);

      SMatrix_free(cube->I_);
      cube->I_ = NULL;

#ifdef SERIALIZE
      SMatrix_free(cube->Udest);
      cube->Udest = NULL;
#endif

#ifdef SERIALIZE
      SMatrix_writebinary(cube->tree->UTIfiles[cube->serialid], cube->UTI, cube->columnrank, gp3);
      SMatrix_free(cube->UTI);
      cube->UTI = NULL;
#endif
   }

   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_fill_PV_UTI(cube->children[cx][cy][cz]);
}

#ifdef ADAPTIVE

void Cube_fill_K(Cube cube) {
   unsigned int cx, cy, cz;
   unsigned int i, r, c, allK = 1;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];
   unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
   BEMLayerType layertype = cube->tree->layertype;
   BEMKernelType kerneltype = cube->tree->kerneltype;
   void* parameters = cube->tree->parameters;
   Vector3D* points = cube->tree->points;
   Panel* panels = cube->tree->panels;
   ComplexVector**** Tprecomputed = cube->tree->Tprecomputed;

   /* If there is nothing on the interation list, simply recurse over
      the children (or quit if a leaf) */

   if (cube->numinteractingcubes == 0) {
      if (!cube->leaf)
         for (cx = 0; cx <= 1; cx++)
            for (cy = 0; cy <= 1; cy++)
               for (cz = 0; cz <= 1; cz++)
                  if (cube->children[cx][cy][cz] != NULL)
                     Cube_fill_K(cube->children[cx][cy][cz]);

      return;
   }

   if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
      cube->K_single = (Matrix*)calloc(cube->numinteractingcubes, sizeof(Matrix));

   if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
      cube->K_double = (Matrix*)calloc(cube->numinteractingcubes, sizeof(Matrix));

   for (i = 0; i < cube->numinteractingcubes; i++) {
      /* Is it worth making this interaction adaptive?  To find out, we 
         count multiplies */

      unsigned int acceleratedops = gp3 * cube->interactingcubes[i]->rowrank +
         4 * padgridsize +
         cube->columnrank * gp3;
      acceleratedops = 4 * padgridsize;  /* is this right? */
      unsigned int Kops = cube->columnrank * cube->interactingcubes[i]->rowrank;

      if (Kops < acceleratedops) {
         /* Compute K */

         unsigned int halfdimension = 1 + 2 * LOCAL;
         int dx = cube->indices[0] - cube->interactingcubes[i]->indices[0];
         int dy = cube->indices[1] - cube->interactingcubes[i]->indices[1];
         int dz = cube->indices[2] - cube->interactingcubes[i]->indices[2];

         ComplexVector T = Tprecomputed[cube->level]
            [dx+halfdimension]
            [dy+halfdimension]
            [dz+halfdimension];

         if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
            cube->K_single[i] = Matrix_allocate(cube->columnrank,  cube->interactingcubes[i]->rowrank);
         if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
            cube->K_double[i] = Matrix_allocate(cube->columnrank,  cube->interactingcubes[i]->rowrank);

         /* Compute by FFTing, translating, and inverse FFTing every row 
            of PV, followed by multiplying by UTI */

         if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
            Matrix IFFT_T_FFT_PV = Matrix_allocate(gp3, cube->interactingcubes[i]->rowrank);

            for (c = 0; c < cube->interactingcubes[i]->rowrank; c++) {
               ComplexVector FFT_PV = ComplexVector_allocate(padgridsize);
               ComplexVector T_FFT_PV = ComplexVector_allocate(padgridsize);
               Vector tempgp3 = Vector_allocate(gp3);

               for (r = 0; r < gp3; r++)
                  tempgp3[r] = cube->interactingcubes[i]->PV_single[r][c];

               FFT_forwardGridTransform(cube->level, gridpoints[cube->level], tempgp3, FFT_PV);

               ComplexVector_addelementmultiplyvector(T_FFT_PV, T, FFT_PV, padgridsize);

               FFT_backwardGridTransform(cube->level, gridpoints[cube->level], T_FFT_PV, tempgp3);

               for (r = 0; r < gp3; r++)
                  IFFT_T_FFT_PV[r][c] = tempgp3[r];

               Vector_free(tempgp3);
               ComplexVector_free(T_FFT_PV);
               ComplexVector_free(FFT_PV);
            }

            Matrix_multiplymatrix(cube->K_single[i], cube->UTI, IFFT_T_FFT_PV, cube->columnrank, gp3, cube->interactingcubes[i]->rowrank);

            Matrix_free(IFFT_T_FFT_PV);
         }

         if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT)) {
            Matrix IFFT_T_FFT_PV = Matrix_allocate(gp3, cube->interactingcubes[i]->rowrank);

            for (c = 0; c < cube->interactingcubes[i]->rowrank; c++) {
               ComplexVector FFT_PV = ComplexVector_allocate(padgridsize);
               ComplexVector T_FFT_PV = ComplexVector_allocate(padgridsize);
               Vector tempgp3 = Vector_allocate(gp3);

               for (r = 0; r < gp3; r++)
                  tempgp3[r] = cube->interactingcubes[i]->PV_double[r][c];

               FFT_forwardGridTransform(cube->level, gridpoints[cube->level], tempgp3, FFT_PV);

               ComplexVector_addelementmultiplyvector(T_FFT_PV, T, FFT_PV, padgridsize);

               FFT_backwardGridTransform(cube->level, gridpoints[cube->level], T_FFT_PV, tempgp3);

               for (r = 0; r < gp3; r++)
                  IFFT_T_FFT_PV[r][c] = tempgp3[r];

               Vector_free(tempgp3);
               ComplexVector_free(T_FFT_PV);
               ComplexVector_free(FFT_PV);
            }

            Matrix_multiplymatrix(cube->K_double[i], cube->UTI, IFFT_T_FFT_PV, cube->columnrank, gp3, cube->interactingcubes[i]->rowrank);

            Matrix_free(IFFT_T_FFT_PV);
         }
      }
      else
         allK = 0;
   }

   if (allK) {
      /* If this cube is only using K matrices, we can dispose of all
         matrices associated with the MV product */
   }

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_fill_K(cube->children[cx][cy][cz]);
}

#endif

void Cube_memory(Cube cube, unsigned int* Cubemem, unsigned int* Umem, unsigned int* VTmem, unsigned int* UTImem, unsigned int* PVmem, unsigned int* Kmem, unsigned int* Dmem) {
   unsigned int cx, cy, cz;
   unsigned int i;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];
   unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
   BEMLayerType layertype = cube->tree->layertype;

   *Cubemem += sizeof(struct _Cube);
   *Cubemem += 3 * sizeof(struct _Vector3D);  // bounds, center
   *Cubemem += cube->numpanelindices * sizeof(unsigned int);  // panelindices
   *Cubemem += cube->numpointindices * sizeof(unsigned int);  // pointindices
   *Cubemem += cube->numlocalcubes * sizeof(Cube*);  // localcubes
   *Cubemem += cube->numinteractingcubes * sizeof(Cube*);  // interactingcubes
   *Cubemem += 2 * padgridsize * sizeof(complexsreal);  /* FFT_PV_VT_q/sum_T_FFT_PV_VT_q */
   //NOTRANS
   *Cubemem += 2 * padgridsize * sizeof(complexsreal);  /* FFT_UTIT_UT_q/sum_T_FFT_UTIT_UT_q */

   if (cube->leaf) {
      if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
         if (cube->D_single)
            *Dmem += cube->drows * cube->dcolumns * sizeof(sreal);
      if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
         if (cube->D_double)
            *Dmem += cube->drows * cube->dcolumns * sizeof(sreal);
   }

   if (cube->VTsrc)
      *VTmem += cube->rowrank * cube->numpanelindices * sizeof(sreal);  /* VTsrc */
   if (cube->Udest)
      *Umem += cube->numpointindices * cube->columnrank * sizeof(sreal);  /* Udest */
   if (cube->UTI)
      *UTImem += cube->columnrank * gp3 * sizeof(sreal);
   if (cube->PV_single)
      if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
         *PVmem += gp3 * cube->rowrank * sizeof(sreal);
   if (cube->PV_double)
      if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
         *PVmem += gp3 * cube->rowrank * sizeof(sreal);

#ifdef ADAPTIVE
   if (cube->K_single)
      if ((layertype == SINGLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
         for (i = 0; i < cube->numinteractingcubes; i++)
            if (cube->K_single[i])
               *Kmem += cube->interactingcubes[i]->rowrank * cube->columnrank * sizeof(sreal);
   if (cube->K_double)
      if ((layertype == DOUBLE_LAYER_INT) || (layertype == SINGLE_AND_DOUBLE_LAYER_INT))
         for (i = 0; i < cube->numinteractingcubes; i++)
            if (cube->K_double[i])
               *Kmem += cube->interactingcubes[i]->rowrank * cube->columnrank * sizeof(sreal);
#endif

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_memory(cube->children[cx][cy][cz], Cubemem, Umem, VTmem, UTImem, PVmem, Kmem, Dmem);
}

void Cube_multiply_D(Vector b, Cube cube, Vector x, BEMLayerType layertype) {
   unsigned int i, j, count;
   unsigned int cx, cy, cz;

   /* If the cube is a leaf, we need to multiply with the dense matrix */

   if ((cube->leaf) && (cube->drows > 0) && (cube->dcolumns > 0)) {
      SVector subx = SVector_allocate(cube->dcolumns);
      SVector subb = SVector_allocate(cube->drows);

      /* Copy into subx the entries of x we need to multiply by */
   
      count = 0;
   
      for (i = 0; i < cube->numpanelindices; i++) {
         subx[count] = x[cube->panelindices[i]];
         count++;
      }
      for (i = 0; i < cube->numlocalcubes; i++)
         for (j = 0; j < cube->localcubes[i]->numpanelindices; j++) {
            subx[count] = x[cube->localcubes[i]->panelindices[j]];
            count++;
         }
      
      /* Perform the dense matrix multiply
         subb = D * subx */

#ifdef SERIALIZE
      if (layertype == SINGLE_LAYER_INT) {
         cube->D_single = Matrix_allocate(cube->drows, cube->dcolumns);
         Matrix_readbinary(cube->D_single, cube->tree->Dfiles_single[cube->serialid], cube->drows, cube->dcolumns);
      }
      else if (layertype == DOUBLE_LAYER_INT) {
         cube->D_double = Matrix_allocate(cube->drows, cube->dcolumns);
         Matrix_readbinary(cube->D_double, cube->tree->Dfiles_double[cube->serialid], cube->drows, cube->dcolumns);
      }
#endif
   
      if (layertype == SINGLE_LAYER_INT)
         SMatrix_multiplyvector(subb, cube->D_single, subx, cube->drows, cube->dcolumns);
      else if (layertype == DOUBLE_LAYER_INT)
         SMatrix_multiplyvector(subb, cube->D_double, subx, cube->drows, cube->dcolumns);

#ifdef SERIALIZE
      if (layertype == SINGLE_LAYER_INT) {
         Matrix_free(cube->D_single);
         cube->D_single = NULL;
      }
      else if (layertype == DOUBLE_LAYER_INT) {
         Matrix_free(cube->D_double);
         cube->D_double = NULL;
      }
#endif

      /* Add subb to the appropriate entries of b */
   
      for (i = 0; i < cube->numpointindices; i++)
         b[cube->pointindices[i]] += subb[i];

      SVector_free(subx);
      SVector_free(subb);
   }
   else /* Recurse over children */
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_D(b, cube->children[cx][cy][cz], x, layertype);
}

void Cube_multiply_DT(Vector b, Cube cube, Vector x, BEMLayerType layertype) {
   unsigned int i, j, count;
   unsigned int cx, cy, cz;

   /* If the cube is a leaf, we need to multiply with the dense matrix */

   if ((cube->leaf) && (cube->drows > 0) && (cube->dcolumns > 0)) {
      SVector subx = SVector_allocate(cube->drows);
      SVector subb = SVector_allocate(cube->dcolumns);

      /* Copy into subx the entries of x we need to multiply by */
   
      for (i = 0; i < cube->numpointindices; i++)
         subx[i] = x[cube->pointindices[i]];
   
      /* Perform the dense matrix multiply
         subb = DT * subx */
   
      if (layertype == SINGLE_LAYER_INT)
         SMatrix_multiplyvector_transpose(subb, cube->D_single, subx, cube->drows, cube->dcolumns);
      if (layertype == DOUBLE_LAYER_INT)
         SMatrix_multiplyvector_transpose(subb, cube->D_double, subx, cube->drows, cube->dcolumns);

      /* Add subb to the appropriate entries of b */

      count = 0;
   
      for (i = 0; i < cube->numpanelindices; i++) {
         b[cube->panelindices[i]] += subb[count];
         count++;
      }
      for (i = 0; i < cube->numlocalcubes; i++)
         for (j = 0; j < cube->localcubes[i]->numpanelindices; j++) {
            b[cube->localcubes[i]->panelindices[j]] += subb[count];
            count++;
         }

      SVector_free(subx);
      SVector_free(subb);
   }
   else /* Recurse over children */
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_DT(b, cube->children[cx][cy][cz], x, layertype);
}

void Cube_multiply_FFT_PV_VT(Cube cube, Vector x, BEMLayerType layertype) {
   unsigned int i;
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];

   if ((cube->numinteractingcubes > 0) && (cube->numpanelindices > 0) && (cube->rowrank > 0)) {
      SVector subx = SVector_allocate(cube->numpanelindices);
#ifndef ADAPTIVE
      SVector VT_q = SVector_allocate(cube->rowrank);
#endif
      SVector PV_VT_q = SVector_allocate(gp3);

      /* Copy into subx the entries of x we need to multiply by */

      for (i = 0; i < cube->numpanelindices; i++)
         subx[i] = x[cube->panelindices[i]];

      /* Project subx with VTsrc */

#ifdef SERIALIZE
      cube->VTsrc = Matrix_allocate(cube->rowrank, cube->numpanelindices);
      SMatrix_readbinary(cube->VTsrc, cube->tree->VTfiles[cube->serialid], cube->rowrank, cube->numpanelindices);
#endif

#ifdef ADAPTIVE
      SMatrix_multiplyvector(cube->VT_q, cube->VTsrc, subx, cube->rowrank, cube->numpanelindices);
#else
      SMatrix_multiplyvector(VT_q, cube->VTsrc, subx, cube->rowrank, cube->numpanelindices);
#endif

      SVector_free(subx);

#ifdef SERIALIZE
      SMatrix_free(cube->VTsrc);
      cube->VTsrc = NULL;
#endif

      /* Multiply VT_q by PV */

#ifdef SERIALIZE
      if (layertype == SINGLE_LAYER_INT) {
         cube->PV_single = SMatrix_allocate(gp3, cube->rowrank);
         SMatrix_readbinary(cube->PV_single, cube->tree->PVfiles_single[cube->serialid], gp3, cube->rowrank);
      }
      else if (layertype == DOUBLE_LAYER_INT) {
         cube->PV_double = SMatrix_allocate(gp3, cube->rowrank);
         SMatrix_readbinary(cube->PV_double, cube->tree->PVfiles_double[cube->serialid], gp3, cube->rowrank);
      }
#endif

#ifdef ADAPTIVE
      if (layertype == SINGLE_LAYER_INT)
         SMatrix_multiplyvector(PV_VT_q, cube->PV_single, cube->VT_q, gp3, cube->rowrank);
      else if (layertype == DOUBLE_LAYER_INT)
         SMatrix_multiplyvector(PV_VT_q, cube->PV_double, cube->VT_q, gp3, cube->rowrank);
#else
      if (layertype == SINGLE_LAYER_INT)
         SMatrix_multiplyvector(PV_VT_q, cube->PV_single, VT_q, gp3, cube->rowrank);
      else if (layertype == DOUBLE_LAYER_INT)
         SMatrix_multiplyvector(PV_VT_q, cube->PV_double, VT_q, gp3, cube->rowrank);
#endif

#ifdef SERIALIZE
      if (layertype == SINGLE_LAYER_INT) {
         SMatrix_free(cube->PV_single);
         cube->PV_single = NULL;
      }
      else if (layertype == DOUBLE_LAYER_INT) {
         SMatrix_free(cube->PV_double);
         cube->PV_double = NULL;
      }
#endif

      /* Apply the forward FFT */

      FFT_forwardGridTransform(cube->level, gridpoints[cube->level], PV_VT_q, cube->FFT_PV_VT_q);

#ifndef ADAPTIVE
      SVector_free(VT_q);
#endif
      SVector_free(PV_VT_q);
   }
         
   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_FFT_PV_VT(cube->children[cx][cy][cz], x, layertype);
}

void Cube_multiply_FFT_UTIT_UT(Cube cube, Vector x) {
   unsigned int i;
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];

   if ((cube->numinteractingcubes > 0) && (cube->numpointindices > 0) && (cube->columnrank > 0)) {
      SVector subx = SVector_allocate(cube->numpointindices);
      SVector UT_q = SVector_allocate(cube->columnrank);
      SVector UTIT_UT_q = SVector_allocate(gp3);

      /* Copy into subx the entries of x we need to multiply by */

      for (i = 0; i < cube->numpointindices; i++)
         subx[i] = x[cube->pointindices[i]];

      /* Project subx with UTdest */

      SMatrix_multiplyvector_transpose(UT_q, cube->Udest, subx, cube->numpointindices, cube->columnrank);

      SVector_free(subx);

      /* Multiply UT_q by UTIT */

      SMatrix_multiplyvector_transpose(UTIT_UT_q, cube->UTI, UT_q, cube->columnrank, gp3);

      /* Apply the forward FFT */

      FFT_forwardGridTransform(cube->level, gridpoints[cube->level], UTIT_UT_q, cube->FFT_UTIT_UT_q);

      SVector_free(UT_q);
      SVector_free(UTIT_UT_q);
   }
         
   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_FFT_UTIT_UT(cube->children[cx][cy][cz], x);
}

void Cube_multiplyboth_FFT_PV_VT(Cube cube, Vector xs, Vector xd) {
   unsigned int i;
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];

   if ((cube->numinteractingcubes > 0) && (cube->numpanelindices > 0) && (cube->rowrank > 0)) {
      SVector subxs = SVector_allocate(cube->numpanelindices);
      SVector subxd = SVector_allocate(cube->numpanelindices);
      SVector sVT_q = SVector_allocate(cube->rowrank);
      SVector dVT_q = SVector_allocate(cube->rowrank);
      SVector PV_VT_q = SVector_allocate(gp3);

      /* Copy into subxs/d the entries of xs/d we need to multiply by */

      for (i = 0; i < cube->numpanelindices; i++) {
         subxs[i] = xs[cube->panelindices[i]];
         subxd[i] = xd[cube->panelindices[i]];
      }

      /* Project subxs/d with VTsrc */

#ifdef SERIALIZE
      cube->VTsrc = SMatrix_allocate(cube->rowrank, cube->numpanelindices);
      SMatrix_readbinary(cube->VTsrc, cube->tree->VTfiles[cube->serialid], cube->rowrank, cube->numpanelindices);
#endif

      SMatrix_multiplyvector(sVT_q, cube->VTsrc, subxs, cube->rowrank, cube->numpanelindices);
      SMatrix_multiplyvector(dVT_q, cube->VTsrc, subxd, cube->rowrank, cube->numpanelindices);

      SVector_free(subxs);
      SVector_free(subxd);

#ifdef SERIALIZE
      SMatrix_free(cube->VTsrc);
      cube->VTsrc = NULL;
#endif

      /* Multiply s/dVT_q by PV */

#ifdef SERIALIZE
      cube->PV_single = SMatrix_allocate(gp3, cube->rowrank);
      SMatrix_readbinary(cube->PV_single, cube->tree->PVfiles_single[cube->serialid], gp3, cube->rowrank);
      cube->PV_double = SMatrix_allocate(gp3, cube->rowrank);
      SMatrix_readbinary(cube->PV_double, cube->tree->PVfiles_double[cube->serialid], gp3, cube->rowrank);
#endif

      SMatrix_multiplyvector(PV_VT_q, cube->PV_single, sVT_q, gp3, cube->rowrank);
      SMatrix_multiplyvector(PV_VT_q, cube->PV_double, dVT_q, gp3, cube->rowrank);

#ifdef SERIALIZE
      SMatrix_free(cube->PV_single);
      cube->PV_single = NULL;
      SMatrix_free(cube->PV_double);
      cube->PV_double = NULL;
#endif

      /* Apply the forward FFT */

      FFT_forwardGridTransform(cube->level, gridpoints[cube->level], PV_VT_q, cube->FFT_PV_VT_q);

      SVector_free(sVT_q);
      SVector_free(dVT_q);
      SVector_free(PV_VT_q);
   }
         
   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiplyboth_FFT_PV_VT(cube->children[cx][cy][cz], xs, xd);
}

void Cube_multiply_T(Cube cube) {
   unsigned int c;
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
   ComplexSVector**** Tprecomputed = cube->tree->Tprecomputed;

   if ((cube->numinteractingcubes > 0) && (cube->numpointindices > 0)) {
      for (c = 0; c < cube->numinteractingcubes; c++) {
         /* If the interacting cube has no panels, it can't contribute to me */

         if (cube->interactingcubes[c]->numpanelindices == 0)
            continue;

#ifdef ADAPTIVE
         /* If this interaction is handled by a K-matrix, dont translate */

         if ((cube->K_single && cube->K_single[c]) ||
             (cube->K_double && cube->K_double[c]))
            continue;
#endif

         /* Determine the correct FFT translation to use from the
            precomputed set */
      
         unsigned int halfdimension = 1 + 2 * LOCAL;
         int dx = cube->indices[0] - cube->interactingcubes[c]->indices[0];
         int dy = cube->indices[1] - cube->interactingcubes[c]->indices[1];
         int dz = cube->indices[2] - cube->interactingcubes[c]->indices[2];
      
         ComplexSVector T = Tprecomputed[cube->level]
            [dx+halfdimension]
            [dy+halfdimension]
            [dz+halfdimension];

         /* Translate cubes with the appropriate T vector */
            
         ComplexSVector_addelementmultiplyvector(cube->sum_T_FFT_PV_VT_q, T, cube->interactingcubes[c]->FFT_PV_VT_q, padgridsize);
      }
   }
   
   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_T(cube->children[cx][cy][cz]);
}

void Cube_multiply_T_transpose(Cube cube) {
   unsigned int c;
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
   ComplexSVector**** Tprecomputed = cube->tree->Tprecomputed;

   if ((cube->numinteractingcubes > 0) && (cube->numpanelindices > 0)) {
      for (c = 0; c < cube->numinteractingcubes; c++) {
         /* If the interacting cube has no points, it can't contribute to me */

         if (cube->interactingcubes[c]->numpointindices == 0)
            continue;

         /* Determine the correct FFT translation to use from the
            precomputed set */
      
         unsigned int halfdimension = 1 + 2 * LOCAL;
         int dx = cube->indices[0] - cube->interactingcubes[c]->indices[0];
         int dy = cube->indices[1] - cube->interactingcubes[c]->indices[1];
         int dz = cube->indices[2] - cube->interactingcubes[c]->indices[2];

         ComplexSVector T = Tprecomputed[cube->level]
            [dx+halfdimension]
            [dy+halfdimension]
            [dz+halfdimension];

         /* Translate cubes with the appropriate T vector */
            
         ComplexSVector_addelementmultiplyvector(cube->sum_T_FFT_UTIT_UT_q, T, cube->interactingcubes[c]->FFT_UTIT_UT_q, padgridsize);
      }
   }
   
   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_T_transpose(cube->children[cx][cy][cz]);
}

void Cube_multiply_U_UTI_IFFT(Vector b, Cube cube) {
   unsigned int i;
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];

   if ((cube->numinteractingcubes > 0) && (cube->numpointindices > 0) && (cube->columnrank > 0)) {
      SVector subb = SVector_allocate(cube->numpointindices);
      SVector IFFT_sum_T_FFT_PV_VT_q = SVector_allocate(gp3);
      SVector UTI_IFFT_sum_T_FFT_PV_VT_q = SVector_allocate(cube->columnrank);

      /* Apply the backward IFFT */

      FFT_backwardGridTransform(cube->level, gridpoints[cube->level], cube->sum_T_FFT_PV_VT_q, IFFT_sum_T_FFT_PV_VT_q);

      /* Multiply by UTI */

#ifdef SERIALIZE
      cube->UTI = SMatrix_allocate(cube->columnrank, gp3);
      SMatrix_readbinary(cube->UTI, cube->tree->UTIfiles[cube->serialid], cube->columnrank, gp3);
#endif

      SMatrix_multiplyvector(UTI_IFFT_sum_T_FFT_PV_VT_q, cube->UTI, IFFT_sum_T_FFT_PV_VT_q, cube->columnrank, gp3);

      SVector_free(IFFT_sum_T_FFT_PV_VT_q);

#ifdef SERIALIZE
      SMatrix_free(cube->UTI);
      cube->UTI = NULL;
#endif

#ifdef ADAPTIVE
      /* Add in contribution from K matrices */

      Vector_addvector(UTI_IFFT_sum_T_FFT_PV_VT_q, cube->sum_K_VT_q, cube->columnrank);
#endif

      /* Project again with the interacting Udest */

#ifdef SERIALIZE
      cube->Udest = SMatrix_allocate(cube->numpointindices, cube->columnrank);
      SMatrix_readbinary(cube->Udest, cube->tree->Ufiles[cube->serialid], cube->numpointindices, cube->columnrank);
#endif

      SMatrix_multiplyvector(subb, cube->Udest, UTI_IFFT_sum_T_FFT_PV_VT_q, cube->numpointindices, cube->columnrank);

      SVector_free(UTI_IFFT_sum_T_FFT_PV_VT_q);

#ifdef SERIALIZE
      Matrix_free(cube->Udest);
      cube->Udest = NULL;
#endif

      /* Add subb to the appropriate entries of b */
      
      for (i = 0; i < cube->numpointindices; i++)
         b[cube->pointindices[i]] += subb[i];
         
      SVector_free(subb);
   }

   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_U_UTI_IFFT(b, cube->children[cx][cy][cz]);
}

void Cube_multiply_V_PVT_IFFT(Vector b, Cube cube, BEMLayerType layertype) {
   unsigned int i;
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];

   if ((cube->numinteractingcubes > 0) && (cube->numpanelindices > 0) && (cube->rowrank > 0)) {
      SVector subb = SVector_allocate(cube->numpanelindices);
      SVector IFFT_sum_T_FFT_UTIT_UT_q = SVector_allocate(gp3);
      SVector PVT_IFFT_sum_T_FFT_UTIT_UT_q = SVector_allocate(cube->rowrank);

      /* Apply the backward IFFT */

      FFT_backwardGridTransform(cube->level, gridpoints[cube->level], cube->sum_T_FFT_UTIT_UT_q, IFFT_sum_T_FFT_UTIT_UT_q);

      /* Multiply by PVT */

      if (layertype == SINGLE_LAYER_INT)
         SMatrix_multiplyvector_transpose(PVT_IFFT_sum_T_FFT_UTIT_UT_q, cube->PV_single, IFFT_sum_T_FFT_UTIT_UT_q, gp3, cube->rowrank);
      if (layertype == DOUBLE_LAYER_INT)
         SMatrix_multiplyvector_transpose(PVT_IFFT_sum_T_FFT_UTIT_UT_q, cube->PV_double, IFFT_sum_T_FFT_UTIT_UT_q, gp3, cube->rowrank);

      SVector_free(IFFT_sum_T_FFT_UTIT_UT_q);

      /* Project again with the interacting VTsrc */

      SMatrix_multiplyvector_transpose(subb, cube->VTsrc, PVT_IFFT_sum_T_FFT_UTIT_UT_q, cube->rowrank, cube->numpanelindices);

      SVector_free(PVT_IFFT_sum_T_FFT_UTIT_UT_q);

      /* Add subb to the appropriate entries of b */
      
      for (i = 0; i < cube->numpanelindices; i++)
         b[cube->panelindices[i]] += subb[i];
         
      SVector_free(subb);
   }

   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_V_PVT_IFFT(b, cube->children[cx][cy][cz], layertype);
}

#ifdef ADAPTIVE

void Cube_multiply_K(Cube cube, BEMLayerType layertype) {
   unsigned int cx, cy, cz;
   unsigned int i;

   if (layertype == SINGLE_LAYER_INT) {
      for (i = 0; i < cube->numinteractingcubes; i++)
         if (cube->K_single && cube->K_single[i])
            Matrix_multiplyvector(cube->sum_K_VT_q, cube->K_single[i], cube->interactingcubes[i]->VT_q, cube->columnrank, cube->interactingcubes[i]->rowrank);
   }
   else if (layertype == DOUBLE_LAYER_INT) {
      for (i = 0; i < cube->numinteractingcubes; i++)
         if (cube->K_double && cube->K_double[i])
            Matrix_multiplyvector(cube->sum_K_VT_q, cube->K_double[i], cube->interactingcubes[i]->VT_q, cube->columnrank, cube->interactingcubes[i]->rowrank);
   }

   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_multiply_K(cube->children[cx][cy][cz], layertype);
}

#endif

void Cube_clear_tempvectors(Cube cube) {
   unsigned int cx, cy, cz;
   unsigned int* gridpoints = cube->tree->gridpointsperlevel;
   unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);

   if (cube->numinteractingcubes > 0) {
#ifdef ADAPTIVE
      Vector_zero(cube->VT_q, cube->rowrank);
      Vector_zero(cube->sum_K_VT_q, cube->columnrank);
#endif
      ComplexSVector_zero(cube->FFT_PV_VT_q, padgridsize);
      ComplexSVector_zero(cube->sum_T_FFT_PV_VT_q, padgridsize);
      //NOTRANS
      ComplexSVector_zero(cube->FFT_UTIT_UT_q, padgridsize);
      ComplexSVector_zero(cube->sum_T_FFT_UTIT_UT_q, padgridsize);
   }
   
   /* Recurse over children */

   if (!cube->leaf)
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_clear_tempvectors(cube->children[cx][cy][cz]);
}

void Cube_leafpanelstats(Cube cube, unsigned int* numleaves, unsigned int* maxpanels, unsigned int* numleaveswithinlimitpanels, unsigned int panellimit) {
   unsigned int cx, cy, cz;

   if (cube->leaf) {
      (*numleaves)++;
      if (cube->numpanelindices > *maxpanels)
         *maxpanels = cube->numpanelindices;
      if (cube->numpanelindices <= panellimit)
         (*numleaveswithinlimitpanels)++;
   }
   else
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_leafpanelstats(cube->children[cx][cy][cz], numleaves, maxpanels, numleaveswithinlimitpanels, panellimit);
}

void Cube_extractdiagonal(Vector diag, Cube cube, BEMLayerType layertype) {
   unsigned int cx, cy, cz;
   unsigned int i;

   if (cube->leaf) {
#ifdef SERIALIZE
      if (layertype == SINGLE_LAYER_INT) {
         cube->D_single = Matrix_allocate(cube->drows, cube->dcolumns);
         Matrix_readbinary(cube->D_single, cube->tree->Dfiles_single[cube->serialid], cube->drows, cube->dcolumns);
      }
      else if (layertype == DOUBLE_LAYER_INT) {
         cube->D_double = Matrix_allocate(cube->drows, cube->dcolumns);
         Matrix_readbinary(cube->D_double, cube->tree->Dfiles_double[cube->serialid], cube->drows, cube->dcolumns);
      }
#endif

      for (i = 0; i < cube->numpanelindices; i++) {
         if (layertype == SINGLE_LAYER_INT)
            diag[cube->panelindices[i]] = cube->D_single[i][i];
         else if (layertype == DOUBLE_LAYER_INT)
            diag[cube->panelindices[i]] = cube->D_double[i][i];
      }

#ifdef SERIALIZE
      if (layertype == SINGLE_LAYER_INT) {
         Matrix_free(cube->D_single);
         cube->D_single = NULL;
      }
      else if (layertype == DOUBLE_LAYER_INT) {
         Matrix_free(cube->D_double);
         cube->D_double = NULL;
      }
#endif
   }
   else {
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL)
                  Cube_extractdiagonal(diag, cube->children[cx][cy][cz], layertype);
   }
}
