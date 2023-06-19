/* ----------------------------------------------------------------------
   fix ordern is a child class of "fix", developed based on two classes
   of "fix ave/time" and "fix ave/correlate/long", provided by LAMMPS.
   This command is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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

#include "fix.h"

#include <cstdio>

namespace LAMMPS_NS {

class FixOrderN : public Fix {
 public:
  FixOrderN(class LAMMPS *, int, char **);  // constructor
  ~FixOrderN(); // destructor
  int setmask();  // the code is needed at the end of step
  void init();  // initialization
  void setup(int);
  void end_of_step();
  void invoke_scalar(bigint);

  void write_restart(FILE *);
  void restart(char *);

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
  int vecsize;  // number of vector elements constructed from nrows (exp. diff)
  int sampsize; // number of vector elements for constructing sample arrays
  int flag_Dxyz; // if components of diffusivities in x, y, & z should be written
  int flag_TCconv; // if convective components of thermal conductivity should be written
  int restart_continue; // checks to not sample twice in a timestep due to restarting

  bigint nvalid;        // the next(current) timestep to call end_of_step()
  bigint nvalid_last;   // the previous timestep that end_of_step() was called

  double deltat;    // timeinterval (nevery*dt)
  double time;    // correlation time (only for writting data)
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
  double **nsamp;	// total number of samples
  double ***oldint;	// The lowest bound of the integral
  double *rint;   // running integral
  int *nbe;	// (BlockLength) number of active elements of blocks
  int cnb;  // (NumberOfBlock) current active number of blocks
  int cnbe;	//  (CurrentBlockLength) current active number of elemets of the block 


  // DIFFUSIVITY vairables
  double ***PosC_ii;
  double ***PosC_iix;
  double ***PosC_iiy;
  double ***PosC_iiz;
  double ****PosC_ij;
  double ****PosC_ijx;
  double ****PosC_ijy;
  double ****PosC_ijz;
  double ****PosCorrSum;
  int **groupinfo;
  int **atomingroup;
  int atomgroup;  // the ID of each group
  int tngroup;  // total # of defined groups (tngroup >= ngroup)
  int ngroup;  // total # of groups
  int tnatom;  // total # of atoms in system
  int natom;    // total # of atoms in groups
  int sortID;   // A virtual ID of the atom in a for loop (to tnatom)
  int ID;   // the real ID of an atom on all cores (to tnatom)
  int atomID; // the ID of an atom inside a group (to natom)
  int atommask; // the group mask of an atom 

  // VISCOSITY/THERMCOND vairables
  double ***samp;		// samples of MSDs
  double *data;  // accumulated integrals
  double sumP, numP, avgP;   // hydrostatic pressure (sum/num = avg)
  double *simpf0, *simpf1; // Simpson's rule of integration
  double dist;			// (integral) The integral we need to sample, i.e., "accint-oldint"  
  
  // PRIVATE METHODS
  bigint nextvalid();
  void integrate();
  void write_diffusivity();
  void write_viscosity();
  void write_thermcond();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect number of inputs for fix ordern command

There should be only one compute ID.

E: All components of the vector are required for fix ordern command

No use of brackets []

E: Illegal fix ordern command: illegal numbers for fix ordern command

Self-explanatory.

E: Illegal fix ordern command: nevery is not a factor of nfreq

In case of viscosity and thermal conductivity, 2*nevery should be a factor of nfreq

E: Illegal fix ordern command: nevery is not a factor of start

Self-explanatory.

E: No compute ID for fix ordern command

Self-explanatory.

E: No global compute vector is computed for fix ordern command

Self-explanatory.

E: Input vector for fix ordern command has variable size

Self-explanatory.

E: Cannot open fix ordern files

Self-explanatory.

E: Error in writing file header for fix ordern command

Self-explanatory.

E: Invalid timestep reset for fix order

Self-explanatory.

E: restart is not possible with different parameters

Self-explanatory.

E: restart error: change in the number of group

Self-explanatory.

E: restart error: change in the number of atoms

Self-explanatory.

*/
