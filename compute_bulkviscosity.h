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

#ifdef COMPUTE_CLASS

ComputeStyle(bulkviscosity,ComputeBulkViscosity)

#else

#ifndef LMP_COMPUTE_BULK_VISCOSITY_H
#define LMP_COMPUTE_BULK_VISCOSITY_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeBulkViscosity : public Compute {
 public:
  ComputeBulkViscosity(class LAMMPS *, int, char **);
  virtual ~ComputeBulkViscosity();
  void init();
  double compute_scalar();
  void reset_extra_compute_fix(const char *);
  // SJ: Start editing
  #define NUMBER_OF_BLOCKS		10
  #define NUMBER_OF_BLOCKELEMENTS	10
  // SJ: End editing

 protected:
  double boltz,nktv2p,inv_volume;
  int nvirial,dimension;
  double **vptr;
  double *kspace_virial;
  Compute *temperature;
  char *id_temp;
  double virial[6];
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int fixflag,kspaceflag;

  void virial_compute(int, int);
  // SJ: Start editing
  int me;				// The rank of each processor
  int nprocs; 				// The total number of processors
  int samplerate;			// This is dt. Every 2dt an integration is carried out
  double timeinterval;
  int count, countint;			// how many cycles have been passed (count = countint*2)
  double Pxx, Pyy, Pzz, Pxy, Pxz, Pyz;	//All stress tensor components 
  double Pxx_s, Pyy_s, Pzz_s;		// traceless components
  double pressure;			// The hydrostatic pressure: (Pxx+Pyy+Pzz)/3
  double avgpressure;				// The average hydrostatic pressure
  double sumpressure;				// The sum of hydrostatic pressure
  double numpressure;				// The number of hydrostatic pressure
  double tmpa[7], tmpb[7], tmpc[7];	// Temporary arrays for simpson's rule integration
  double lastint[7];			// The last integral between t-2dt and t
  double accint[7];			// The accumlative integral upto the newest timestep
  double oldint[NUMBER_OF_BLOCKS][NUMBER_OF_BLOCKELEMENTS][7];	// The lowest bound of the integral
  double integral;			// The integral we need to sample, i.e., "accint-oldint"
  //double integrals[NUMBER_OF_BLOCKS][NUMBER_OF_BLOCKELEMENTS][6];	// Integrals to t
  double samples[NUMBER_OF_BLOCKS][NUMBER_OF_BLOCKELEMENTS][8];		// samples of int(P^2) & int(P): 7+1
  double samplecount[NUMBER_OF_BLOCKS][NUMBER_OF_BLOCKELEMENTS][8];	// total number of samples: 7+1 
  int BlockLength[NUMBER_OF_BLOCKS];	// Length of each active block
  int NumberOfBlock;			// Up to how which block, integration has been done
  int CurrentBlockLength;		// 
  int WriteFileEvery;			// Writing the viscosity file every
  FILE *FilePtr;			// Pointer to the file
  // SJ: end editing
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute pressure must use group all

Virial contributions computed by potentials (pair, bond, etc) are
computed on all atoms.

E: Could not find compute pressure temperature ID

The compute ID for calculating temperature does not exist.

E: Compute pressure temperature ID does not compute temperature

The compute ID assigned to a pressure computation must compute
temperature.

E: Compute pressure requires temperature ID to include kinetic energy

The keflag cannot be used unless a temperature compute is provided.

E: Virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

E: Must use 'kspace_modify pressure/scalar no' for tensor components with kspace_style msm

Otherwise MSM will compute only a scalar pressure.  See the kspace_modify
command for details on this setting.

*/
