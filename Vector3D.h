/* Typedef */

typedef struct _Vector3D {
   real x, y, z;
} _Vector3D;

typedef _Vector3D* Vector3D;

/* Constructors and Destructors */

Vector3D Vector3D_allocate();
void Vector3D_free(Vector3D v);

/* Initialization and Copying */
 
void Vector3D_copy(Vector3D vectordest, Vector3D vectorsrc);
  
/* Arithmetic Operations */

void Vector3D_add(Vector3D sum, Vector3D v1, Vector3D v2);
void Vector3D_sub(Vector3D diff, Vector3D v1, Vector3D v2);
void Vector3D_scale(Vector3D v, real scale);
void Vector3D_cross(Vector3D cross, Vector3D v1, Vector3D v2);
real Vector3D_dot(Vector3D v1, Vector3D v2);
real Vector3D_length(Vector3D v);
void Vector3D_normalize(Vector3D v);
real Vector3D_distance(Vector3D v1, Vector3D v2);
unsigned int Vector3D_equal(Vector3D v1, Vector3D v2);
void Vector3D_transform(Vector3D dest, real *A, Vector3D src);
void Vector3D_addscaled(Vector3D sum, Vector3D v1, real scale, Vector3D v2);
void Vector3D_transformVecs(Vector3D Ax, Vector3D A1, Vector3D A2, Vector3D A3,Vector3D x);
void Vector3D_transformVecs_inverse(Vector3D Ax, Vector3D A1, Vector3D A2, Vector3D A3, Vector3D x);

void Vector3D_print(Vector3D v);
