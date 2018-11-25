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
  ~FixOrderN();                             // destructor
  int setmask();                            // the code is needed at the end of step
  void init();
  void setup(int);
  void end_of_step();

  void write_restart(FILE *);
  void restart(char *);
  double compute_scalar();
  double compute_vector(int);
  double compute_array(int,int);

 private:

  int me;   // the ID of the local processor
  int mode; // which type of transport properties (defined as enum)
  int nvalues;  // the number of input values (it must be 1)
  int nfreq;   // the frequency of writing files to disk
  int icompute; // the ID of compute for the tranpsort property
  int nrows;  // number of data passed to fix ordern via compute

  bigint nvalid;        // WHAT IS IT???
  bigint nvalid_last;   // WHAT IS IT???


  long filepos1; // the location of header lines for file 1
  long filepos2; // the location of header lines for file 2

  char *filename1;  // filename for self-diffusion, shear viscosity, and conductivity
  char *filename2;  // filename for onsager coefficients and bulk viscosity

  char *format; // the default format of writing data to files (%g)
  char *format_user;  // user-defined format of writing data to files
  char *title; // optional text for the first line of the text
  char *idcompute;   // the ID of compute 

  FILE *fp1;  // output file 1 (self-diffusion, shear viscosity, thermal conductivity)
  FILE *fp2;  // output file 2 (onsager coefficients, bulk viscosity)


  // DIFFUSION inputs
  int numgroup; // number of groups



  
  //int me,nvalues;
  //int nrepeat,nfreq,irepeat;
  //bigint nvalid,nvalid_last;
  //int *which,*argindex,*value2index,*offcol;
  //int *varlen;               // 1 if value is from variable-length compute (NOT NEEDED)
  //char **ids;
  //FILE *fp;
  //int any_variable_length;
  //int all_variable_length;
  //int lockforever;

  //int ave,nwindow,mode;
  //int noff,overwrite;
  //int *offlist;
  //char *format,*format_user;
  //char *title1,*title2,*title3;

  //int norm,iwindow,window_limit;
  //double *vector;
  //double *vector_total;
  //double **vector_list;
  //double *column;
  //double **array;
  //double **array_total;
  //double ***array_list;

  int column_length(int);
  void invoke_scalar(bigint);
  void invoke_vector(bigint);
  void options(int, int, char **);
  void allocate_arrays();
  bigint nextvalid();
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
