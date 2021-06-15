#ifndef FFTSVDpbeAPI_H
#define FFTSVDpbeAPI_H

#include "FFTSVD.h"
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>

#define VERSION "0.01"

#define KT_CONVERSION 561.0  /* For compatibility with DelPhi */
#define DEBYE_CONSTANT 3.047  /* For compatibility with DelPhi */

#define GRID_SIZE 4
#define SVD_ERROR 1e-5

#define MAX_PANELS_PER_FINEST_CUBE 32

/* Important Structures */

typedef struct {
   char record[7];
   unsigned int atomnumber;
   char atomname[5];
   char alternatelocation;
   char residuename[4];
   char chain;
   unsigned int residuenumber;
   char residueinsertion;
   real x;
   real y;
   real z;
   real occupancy;
   real temperature;
   unsigned int footnotenumber;

   real radius;
   real charge;
   real potential;
} PDBentry;

typedef struct {
   char atomlabel[7];
   char residuelabel[4];
   real radius;
} SIZentry;

typedef struct {
   char atomlabel[7];
   char residuelabel[4];
   unsigned int residuenumber;
   char chain;
   real charge;
} CRGentry;

/* Important Filenames */

extern const char* delphiparameterfilename;
extern const char* delphisizfilename;
extern const char* delphicrgfilename;
extern const char* delphipdbinputfilename;
extern const char* delphiphioutputfilename;
extern const char* delphipdbcollectfilename;
extern const char* delphifrcfilename;
extern const char* delphiphiinputfilename;
extern const char* delphipdboutputfilename;
extern const char* delphisrffilename;

/* Parameters */

extern real dielectricvertexdensity;
extern real saltvertexdensity;
extern real innerdielectric;
extern real outerdielectric;
extern real lambda;  // nonlocal length scale parameter
extern real inftydielectric; // nonlocal short-range limit of solvent dielectric
extern real ionicstrength;
extern real ionexclusionradius;
extern real proberadius;
extern real errortolerance;
extern char energycalculation;
extern unsigned int outputmodpdb;
extern unsigned int outputfrc;
extern unsigned int addcoulombic;
extern unsigned int accelerateM3;
extern unsigned int usequalocation;
extern unsigned int usecurved;

/* Utility */

void error(char* message, ...);
void warning(char* message, ...);
void removeWhitespace(char* s);
int yesno(char* s);

/* Input */

void readParams(const char* filename);
void checkParams();
void readPDB(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries);
void readCRD(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries);
void readXYZR(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries);
void readSIZ(const char* filename, unsigned int* numSIZentries, SIZentry** SIZentries);
void readCRG(const char* filename, unsigned int* numCRGentries, CRGentry** CRGentries);

/* Output */

void writeMODPDB(const char* filename, PDBentry* PDBentries, unsigned int numPDBentries);
void writeFRC(const char* filename, PDBentry* PDBentries, unsigned int numPDBentries);

/* Molecule Handling */

void assignRadiiCharges(PDBentry* PDBentries, unsigned int numPDBentries, SIZentry* SIZentries, unsigned int numSIZentries, CRGentry* CRGentries, unsigned int numCRGentries);
void generateFlatSRF(const char* pdbfilename, const char* sizfilename, const char* srffilename);
void generateCurvedSRF(const char* pdbfilename, const char* sizfilename, const char* srffilename);
void readSRF(const char* filename, Panel*** saltpanels, unsigned int** numsaltpanels, unsigned int* numsalts,
Panel*** dielectricpanels, unsigned int** numdielectricpanels, unsigned int* numdielectrics, unsigned int** dielectricparent,
Panel*** dielectriccavitypanels, unsigned int** numdielectriccavitypanels, unsigned int* numdielectriccavities, unsigned int** dielectriccavityparent,
Panel*** saltcavitypanels, unsigned int** numsaltcavitypanels, unsigned int* numsaltcavities, unsigned int** saltcavityparent, unsigned int* numtotalpanels);

/* Surface Operator */

void generateSurfaceOperator(SurfaceOperator* pbesurfaceoperator,
PDBentry* PDBentries, unsigned int numPDBentries,
Panel** saltpanels, unsigned int* numsaltpanels, unsigned int numsalts, 
Panel** dielectricpanels, unsigned int* numdielectricpanels, unsigned int numdielectrics, unsigned int* dielectricparent,
Panel** dielectriccavitypanels, unsigned int* numdielectriccavitypanels, unsigned int numdielectriccavities, unsigned int* dielectriccavityparent, 
Panel** saltcavitypanels, unsigned int* numsaltcavitypanels, unsigned int numsaltcavities, unsigned int* saltcavityparent,
unsigned int numtotalpanels);

void generateQualocationOperator(QualocationOperator* qualocationoperator,
PDBentry* PDBentries, unsigned int numPDBentries,
Panel** dielectricpanels, unsigned int* numdielectricpanels, unsigned int numdielectrics,
Panel** dielectriccavitypanels, unsigned int* numdielectriccavitypanels, unsigned int numdielectriccavities, unsigned int* dielectriccavityparent,
unsigned int numtotalpanels);

void generateQualocationOperatorMinimal(QualocationOperator* qualocationoperator,
													 PDBentry* PDBentries, unsigned int numPDBentries,
													 Panel** dielectricpanels, unsigned int* numdielectricpanels, unsigned int numdielectrics,
													 Panel** dielectriccavitypanels, unsigned int* numdielectriccavitypanels,
													 unsigned int numdielectriccavities, unsigned int* dielectriccavityparent,
													 unsigned int numtotalpanels);
void generateQualocationOperatorA1A3(QualocationOperator qualocationoperator);

void generateSurfaceOperatorPreconditioner(Preconditioner* preconditioner, SurfaceOperator pbesurfaceoperator, unsigned int numtotalpanels);
void generateQualocationOperatorPreconditioner(Preconditioner* preconditioner, QualocationOperator qualocationoperator, unsigned int numtotalpanels);

void generateRHS(Vector* rhs, SurfaceOperator pbesurfaceoperator, unsigned int numtotalpanels);
void generateQualocationRHS(Vector* rhs, QualocationOperator qualocationoperator, unsigned int numtotalpanels);

void solveSurfaceOperator(Vector* solution, SurfaceOperator pbesurfaceoperator, Preconditioner preconditioner, Vector rhs, unsigned int numtotalpanels);
void solveQualocationOperator(Vector* solution, QualocationOperator qualocationoperator, Preconditioner preconditioner, Vector rhs, unsigned int numtotalpanels);

void collectPotentials(PDBentry* PDBcollectentries, unsigned int numPDBcollectentries, PDBentry* PDBinputentries, unsigned int numPDBinputentries, SurfaceOperator pbesurfaceoperator, Vector solution);
void collectQualocationPotentials(PDBentry* PDBcollectentries, unsigned int numPDBcollectentries, PDBentry* PDBinputentries, unsigned int numPDBinputentries, QualocationOperator qualocationoperator, Vector solution);

/* for GAMESS_FFTSVD */

void solveBEMproblem(unsigned int *numpanels, real *paneldata,
							real *innerdielectric, real *outerdielectric,
							real *rhs, real *output);
#endif
