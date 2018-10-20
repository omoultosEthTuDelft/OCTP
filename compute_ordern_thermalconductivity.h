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

/*
    compute_ordern_thermalconductivity is based on the compute_heatflux command.
*/

#ifdef COMPUTE_CLASS

ComputeStyle(thermalconductivity/ordern,ComputeOrderNThermalConductivity)

#else

#ifndef LMP_COMPUTE_ORDERN_THERMALCONDUCTIVITY_H
#define LMP_COMPUTE_ORDERN_THERMALCONDUCTIVITY_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeOrderNThermalConductivity : public Compute {
 public:
  ComputeOrderNThermalConductivity(class LAMMPS *, int, char **);
  ~ComputeOrderNThermalConductivity();
  void init();
  void compute_vector();
  #define NUMBER_OF_BLOCKS		10
  #define NUMBER_OF_BLOCKELEMENTS	10

 protected:
  char *id_ke,*id_pe,*id_stress;
  class Compute *c_ke,*c_pe,*c_stress;
  
  double boltz,inv_volume;
  int dimension;
  int me;				// The rank of each processor
  int nprocs; 				// The total number of processors
  int samplerate;			// This is dt. Every 2dt an integration is carried out
  double timeinterval;
  int count, countint;			// how many cycles have been passed (count = countint*2)
  double Jtx, Jty, Jtz, Jcx, Jcy, Jcz;	//The components of total (t) and convective (c) flux
  double pressure;			// The hydrostatic pressure: (Pxx+Pyy+Pzz)/3
  double tmpa[6], tmpb[6], tmpc[6];	// Temporary arrays for simpson's rule integration
  double lastint[6];			// The last integral between t-2dt and t
  double oldint[NUMBER_OF_BLOCKS][NUMBER_OF_BLOCKELEMENTS][6];	// The lowest bound of the integral
  double accint[6];			// The accumlative integral upto the newest timestep
  double integral;			// The integral we need to sample, i.e., "accint-oldint"
  //double integrals[NUMBER_OF_BLOCKS][NUMBER_OF_BLOCKELEMENTS][6];	// Integrals to t
  double samples[NUMBER_OF_BLOCKS][NUMBER_OF_BLOCKELEMENTS][6];		// samples of int(P^2)
  double samplecount[NUMBER_OF_BLOCKS][NUMBER_OF_BLOCKELEMENTS][6];	// total number of samples 
  int BlockLength[NUMBER_OF_BLOCKS];	// Length of each active block
  int NumberOfBlock;			// Up to how which block, integration has been done
  int CurrentBlockLength;		// 
  int WriteFileEvery;			// Writing the viscosity file every
  FILE *FilePtr;			// Pointer to the file
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute heat/flux compute ID

Self-explanatory.

E: Compute heat/flux compute ID does not compute ke/atom

Self-explanatory.

E: Compute heat/flux compute ID does not compute pe/atom

Self-explanatory.

E: Compute heat/flux compute ID does not compute stress/atom

Self-explanatory.

*/
