/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(ordern,FixOrderN)

#else

#ifndef LMP_FIX_ORDERN_H
#define LMP_FIX_ORDERN_H

#include <cstdio>
#include "fix.h"

namespace LAMMPS_NS {

class FixOrderN : public Fix {
 public:
  FixOrderN(class LAMMPS *, int, char **);  // constructor
  ~FixOrderN(); // destructor
  int setmask();  // the code is needed at the end of step
  void init();  // initialization
  void setup(int);
  void end_of_step();

  #define MAXGROUPS 32
  #define SQR(x) ((x)*(x))

  //void write_restart(FILE *);
  //void restart(char *);
  //double memory_usage();

 private:

  int me;   // the ID of the local processor
  int mode; // which type of transport properties (defined as enum)
  int nvalues;  // the number of input values (it must be 1)
  int nfreq;   // the frequency of writing files to disk
  int icompute; // the ID of compute for the tranpsort property
  int nrows;  // number of data passed to fix ordern via compute
  int startstep;  // the first timestep that it starts sampling
  int tnb;   // total number of blocks (different time scales)
  int tnbe;  // total number of block elements (elemets of similar timescales)
  int vecsize;  // number of vector elements constructed from nrows
  int sampsize; // number of vector elements for constructing sample arrays

  bigint nvalid;        // the next(current) timestep to call end_of_step()
  bigint nvalid_last;   // the previous timestep that end_of_step() was called

  double deltat;    // timeinterval (nevery*dt)
  double time;    // correlation time
  double boltz;   // Boltzmann constant
  double nktv2p;  // conversion factors
  double volume;  // volume of the simulation box
	double coef; // conversion coefficient


  long filepos1; // the location of header lines for file 1
  long filepos2; // the location of header lines for file 2

  double *recdata; // data passed from compute to fix;

  char *filename1;  // filename for self-diffusion, shear viscosity, and conductivity
  char *filename2;  // filename for onsager coefficients and bulk viscosity

  char *format; // the default format of writing data to files (%g)
  char *format_user;  // user-defined format of writing data to files
  char *title; // optional text for the first line of the text
  char *idcompute;   // the ID of compute 

  FILE *fp1;  // output file 1 (self-diffusion, shear viscosity, thermal conductivity)
  FILE *fp2;  // output file 2 (onsager coefficients, bulk viscosity)


  // General variables
  int count;			// how many cycles have been passed
  int icount;     // how many nevery, used for integration
  // ORDER-N ALGORITHM variables
  double ***samp;		// samples of int(P^2) & int(P): 7+1
  double **nsamp;	// total number of samples
  double ***oldint;	// The lowest bound of the integral
  int *nbe;	// (BlockLength) nuber of active elements of blocks
  int cnb;  // (NumberOfBlock) current active number of blocks
  int cnbe;	//  (CurrentBlockLength) current active number of elemets of the block 


  // DIFFUSIVITY vairables
  //int numgroup; // number of groups
  // double ****BlockDATA;
  //double **TmpPos;
  double ***PosC_ii;
  double ****PosC_ij;
  double ****PosCorrSum;
  //int **Groups; 		// The group id of each atom
  // Finding corresponding groups to each atom
  int **groupinfo;
  int **atomingroup;
  int atomgroup;  // the ID of each group
  int ngroup ;  // total # of groups
  int tnatom;  // total # of atoms in system
  int natom;    // total # of atoms in groups
  int sortID;   // A virtual ID of the atom in a for loop (to tnatom)
  int ID;   // the real ID of an atom on all cores (to tnatom)
  int atomID; // the ID of an atom inside a group (to natom)
  int atommask; // the group mask of an atom 

  // VISCOSITY vairables
  double sumP, numP, avgP;   // (avgpressure, sumpressure, numpressure) hydrostatic pressure
  // order of data: dP_xx, dP_yy, dP_zz, P_xy, P_xz, P_yz, Phydro
  double *data;  // accumulated integrals
  double *rint;   // running integral
  double *simpf0, *simpf1; // Simpson's rule of integration
  double dist;			// (integral) The integral we need to sample, i.e., "accint-oldint"


  // THERMCOND variables
  
  
  // METHODS
  bigint nextvalid();
  void integrate();
  void write_diffusivity();
  void write_viscosity();
  void write_thermcond();




  // MAYBE LATER

  

  //int noff,overwrite;
  //double *vector;
  //double *vector_total;
  //double **vector_list;
  //double *column;
  //double **array;
  //double **array_total;
  //double ***array_list;

  //int column_length(int);
  //void invoke_scalar(bigint);
  //void invoke_vector(bigint);
  //void options(int, int, char **);
  //void allocate_arrays();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: No values in fix ave/time command

Self-explanatory.

E: Invalid fix ave/time off column

Self-explanatory.

E: Compute ID for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time compute does not calculate a scalar

Self-explanatory.

E: Fix ave/time compute does not calculate a vector

Self-explanatory.

E: Fix ave/time compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ave/time compute does not calculate an array

Self-explanatory.

E: Fix ave/time compute array is accessed out-of-range

An index for the array is out of bounds.

E: Fix ID for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time fix does not calculate a scalar

Self-explanatory.

E: Fix ave/time fix does not calculate a vector

Self-explanatory.

E: Fix ave/time fix vector cannot be variable length

Self-explanatory.

E: Fix ave/time fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix for fix ave/time not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/time
is requesting a value on a non-allowed timestep.

E: Fix ave/time fix does not calculate an array

Self-explanatory.

E: Fix ave/time fix array cannot be variable length

Self-explanatory.

E: Fix ave/time fix array is accessed out-of-range

An index for the array is out of bounds.

E: Variable name for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time variable is not equal-style variable

Self-explanatory.

E: Fix ave/time variable is not vector-style variable

UNDOCUMENTED

E: Fix ave/time mode vector variable cannot be indexed

UNDOCUMENTED

E: Error writing file header

Something in the output to the file triggered an error.

E: Fix ave/time cannot set output array intensive/extensive from these inputs

One of more of the vector inputs has individual elements which are
flagged as intensive or extensive.  Such an input cannot be flagged as
all intensive/extensive when turned into an array by fix ave/time.

E: Invalid timestep reset for fix ave/time

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing out time averaged data

Something in the output to the file triggered an error.

E: Fix ave/time vector-style variable changed length

UNDOCUMENTED

E: Fix ave/time columns are inconsistent lengths

Self-explanatory.

E: Cannot open fix ave/time file %s

The specified file cannot be opened.  Check that the path and name are
correct.

U: Fix ave/time cannot use variable with vector mode

Variables produce scalar values.

*/
