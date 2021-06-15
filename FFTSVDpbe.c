#include "FFTSVDpbeAPI.h"

/* Global Variables */

/* Surface Related */

unsigned int numsalts = 0;
Panel** saltpanels = NULL;
unsigned int* numsaltpanels = NULL;

unsigned int numdielectrics = 0;
Panel** dielectricpanels = NULL;
unsigned int* numdielectricpanels = NULL;
unsigned int* dielectricparent = NULL;

unsigned int numdielectriccavities = 0;
Panel** dielectriccavitypanels = NULL;
unsigned int* numdielectriccavitypanels = NULL;
unsigned int* dielectriccavityparent = NULL;

unsigned int numsaltcavities = 0;
Panel** saltcavitypanels = NULL;
unsigned int* numsaltcavitypanels = NULL;
unsigned int* saltcavityparent = NULL;

unsigned int numtotalpanels = 0;

/* Surface Operator Related */

SurfaceOperator pbesurfaceoperator;  /* The big momma */
QualocationOperator pbequalocationoperator;  /* The big qual momma */
Preconditioner preconditioner;
real* rhs;  /* The RHS for the whole problem, length 2*numtotalpanels */
real* solution;  /* The solution for the whole problem, length 2*numtotalpanels */

/* PBE Related */

real energy = 0.0f;
real kappa = 0.0f;

/* Molecule Related */

PDBentry* PDBinputentries;
unsigned int numPDBinputentries;
PDBentry* PDBcollectentries;
unsigned int numPDBcollectentries;
SIZentry* SIZentries;
unsigned int numSIZentries;
CRGentry* CRGentries;
unsigned int numCRGentries;

void hessian() {
   FILE* h = fopen("hessian", "w");
   Charge charges = pbesurfaceoperator->children[0]->charges;

   if (charges == NULL)
      charges = pbesurfaceoperator->children[0]->children[0]->charges;

   if (charges == NULL) {
      printf("FUTZ\n");
      exit(-1);
   }

   unsigned int i, j;

   for (i = 0; i < charges->numcharges; i++) {
      if (PDBinputentries[i].chain != 'D')
         continue;

      printf("Computing row %u\n", i+1);

      for (j = 0; j < charges->numcharges; j++)
         charges->charges[j] = 0.0;

      charges->charges[i] = 1.0;

      generateRHS(&rhs, pbesurfaceoperator, numtotalpanels);
      solveSurfaceOperator(&solution, pbesurfaceoperator, preconditioner, rhs, numtotalpanels);
      collectPotentials(PDBinputentries, numPDBinputentries, PDBinputentries, numPDBinputentries, pbesurfaceoperator, solution);

      for (j = 0; j < numPDBinputentries; j++) {
         if (PDBinputentries[j].chain != 'D')
            continue;

         fprintf(h, "%f ", PDBinputentries[j].potential * 0.592 * 0.5);
      }
      fprintf(h, "\n");

      free(rhs);
      free(solution);
   }

   fclose(h);

   FILE* frc = fopen(delphifrcfilename, "w");
   fclose(frc);

   exit(0);
}

int main() {
   setlinebuf(stdout);
   accelerateM3 = 0;

   printf("**************************************************\n");
   printf("*FFTSVDpbe %4s                                  *\n", VERSION);
   printf("*FFTSVD Boundary Element Poisson-Boltzmann Solver*\n");
   printf("**************************************************\n");

   printf("\n");
   printf("Parameter Setup\n");
   printf("---------------\n");
   printf("\n");

   printf("Reading parameters from %s\n", delphiparameterfilename);
   readParams(delphiparameterfilename);

   printf("Checking for valid parameters\n");
   checkParams();

   printf("Reading input PDB from %s\n", delphipdbinputfilename);
   readPDB(delphipdbinputfilename, &numPDBinputentries, &PDBinputentries);

   printf("Reading radii from %s\n", delphisizfilename);
   readSIZ(delphisizfilename, &numSIZentries, &SIZentries);

   printf("Reading charges from %s\n", delphicrgfilename);
   readCRG(delphicrgfilename, &numCRGentries, &CRGentries);

   printf("\n");
   printf("Molecule Setup\n");
   printf("--------------\n");
   printf("\n");

   printf("Assigning radii and charges to atoms\n");
   assignRadiiCharges(PDBinputentries, numPDBinputentries, SIZentries, numSIZentries, CRGentries, numCRGentries);

   if (outputmodpdb) {
      printf("Writing modified PDB to %s\n", delphipdboutputfilename);
      writeMODPDB(delphipdboutputfilename, PDBinputentries, numPDBinputentries);
   }

   printf("\n");
   printf("Surface Generation\n");
   printf("------------------\n");
   printf("\n");

   printf("Generating salt and dielectric surfaces\n\n");

   if (usecurved)
      generateCurvedSRF(delphipdbinputfilename, delphisizfilename, "operator.srf");
   else
      generateFlatSRF(delphipdbinputfilename, delphisizfilename, "operator.srf");

   printf("\n");

   printf("Processing SRF file\n\n");
   readSRF("operator.srf", &saltpanels, &numsaltpanels, &numsalts,
           &dielectricpanels, &numdielectricpanels, &numdielectrics, &dielectricparent,
           &dielectriccavitypanels, &numdielectriccavitypanels, &numdielectriccavities, &dielectriccavityparent,
           &saltcavitypanels, &numsaltcavitypanels, &numsaltcavities, &saltcavityparent,
           &numtotalpanels);

   printf("\n");
   printf("Surface Operator Setup\n");
   printf("----------------------\n");
   printf("\n");

   if (usequalocation) {
      printf("Generating qualocation operator\n");
      generateQualocationOperator(&pbequalocationoperator,
                                  PDBinputentries, numPDBinputentries,
                                  dielectricpanels, numdielectricpanels, numdielectrics,
                                  dielectriccavitypanels, numdielectriccavitypanels, numdielectriccavities, dielectriccavityparent,
                                  numtotalpanels);
   }
   else {
      printf("Generating surface operator\n");
      generateSurfaceOperator(&pbesurfaceoperator,
                       PDBinputentries, numPDBinputentries,
                       saltpanels, numsaltpanels, numsalts,
                       dielectricpanels, numdielectricpanels, numdielectrics, dielectricparent,
                       dielectriccavitypanels, numdielectriccavitypanels, numdielectriccavities, dielectriccavityparent,
                       saltcavitypanels, numsaltcavitypanels, numsaltcavities, saltcavityparent,
                       numtotalpanels);
   }

   if (usequalocation) {
      printf("\nGenerating qualocation operator preconditioner\n");
      generateQualocationOperatorPreconditioner(&preconditioner, pbequalocationoperator, numtotalpanels);
   }
   else {
      printf("\nGenerating surface operator preconditioner\n");
      generateSurfaceOperatorPreconditioner(&preconditioner, pbesurfaceoperator, numtotalpanels);
   }

   // remove for normal action
   //hessian();

   if (usequalocation) {
      printf("Generating qualocation operator RHS\n");
      generateQualocationRHS(&rhs, pbequalocationoperator, numtotalpanels);
   }
   else {
      printf("Generating surface operator RHS\n");
      generateRHS(&rhs, pbesurfaceoperator, numtotalpanels);
   }

   printf("\n");
   printf("FFTSVD BEM GMRES Solve\n");
   printf("----------------------\n");
   printf("\n");

   if (usequalocation)
      solveQualocationOperator(&solution, pbequalocationoperator, preconditioner, rhs, numtotalpanels);
   else
      solveSurfaceOperator(&solution, pbesurfaceoperator, preconditioner, rhs, numtotalpanels);

   printf("\n");
   printf("Potentials\n");
   printf("----------\n");
   printf("\n");

   printf("Reading collect PDB from %s\n", delphipdbcollectfilename);
   readPDB(delphipdbcollectfilename, &numPDBcollectentries, &PDBcollectentries);

   printf("Assigning radii and charges to atoms\n");
   assignRadiiCharges(PDBcollectentries, numPDBcollectentries, SIZentries, numSIZentries, CRGentries, numCRGentries);

   printf("Collecting potentials at atom centers\n");
   if (usequalocation)
      collectQualocationPotentials(PDBcollectentries, numPDBcollectentries, PDBinputentries, numPDBinputentries, pbequalocationoperator, solution);
   else
      collectPotentials(PDBcollectentries, numPDBcollectentries, PDBinputentries, numPDBinputentries, pbesurfaceoperator, solution);

   if (outputfrc) {
      printf("Writing potentials at atom centers to %s\n", delphifrcfilename);
      writeFRC(delphifrcfilename, PDBcollectentries, numPDBcollectentries);
   }

   //SurfaceOperator_writematlabfile("operator.txt", pbesurfaceoperator, numtotalpanels);
   //QualocationOperator_writematlabfile("operator.txt", pbequalocationoperator, numtotalpanels);

   //SurfaceOperator_free(pbesurfaceoperator);

   printf("\n");
   printf("*********************\n");
   printf("*FFTSVDpbe Complete!*\n");
   printf("*********************\n");

   // MEMORY SNIPPET
   FILE* status = fopen("/proc/self/status", "r");
   unsigned int i, kb;
   char line[256], garbage1[256], garbage2[256];

   for (i = 0; i < 12; i++) {
      fgets(line, 256, status);
   }
   fgets(line, 256, status);
   sscanf(line, "%s %u %s", garbage1, &kb, garbage2);

   printf("%s\n", garbage1);
   printf("MEMORY USE: %u\n", kb);

   fclose(status);

   return 0;
}
