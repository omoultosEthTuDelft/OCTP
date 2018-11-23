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
    compute_ordern_diffusion is based on the compute_msd command.
*/
#ifdef COMPUTE_CLASS

ComputeStyle(diffusion/ordern,ComputeOrdernDiffusion)

#else

#ifndef LMP_COMPUTE_ORDERN_DIFFUSION_H
#define LMP_COMPUTE_ORDERN_DIFFUSION_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeOrdernDiffusion : public Compute {
 public:
  ComputeOrdernDiffusion(class LAMMPS *, int, char **);
  virtual ~ComputeOrdernDiffusion();
  void init();
  virtual void compute_vector();
  void set_arrays(int);
  #define MAX_NUMBER_OF_BLOCKS            10
  #define MAX_NUMBER_OF_BLOCKELEMENTS     10
  #define MAX_NUMBER_OF_ATOMS         100000
  #define MAX_NUMBER_OF_TYPES             32
  #define SQR(x) ((x)*(x))
 
 protected:
  int comflag;   // comflag = 1 if reference moves with center of mass
  int avflag;    // avflag = 1 if using average position as reference
  int naverage;  // number of samples for average position
  bigint nmsd;
  double masstotal;
  char *id_fix;
  class FixStore *fix;
  // The parameters related to MPI and gathering the positions
  int me; // This is for the definition of rank of processor in LAMMPS
  int nprocs; // This is for the number of processors
  int *tmprecvcnts; // Number of atoms per each processor [array]
  int *recvcnts; // Number of atoms per each processor [array]
  int *displs; // Displacement array for receiving the correct data set from each processor [array]
  double sendbuff [MAX_NUMBER_OF_ATOMS*5]; // The sending array from all processors
  double recvbuff [MAX_NUMBER_OF_ATOMS*5]; // The receiving array from all processors
  int groupatom1, groupatom2;
  double timeinterval;
  int samplerate;
  int count , NumberOfBlocks ;
  int realatom;	// The tag of each atom in the global simulation system

  // Storing the main data for the order-n algorithm
  double ****BlockDATA;  
      // The "sum" arrays  have the total sum from all cores
  double PosC_ii[MAX_NUMBER_OF_BLOCKS][MAX_NUMBER_OF_BLOCKELEMENTS][MAX_NUMBER_OF_TYPES];
  double PosC_ij[MAX_NUMBER_OF_BLOCKS][MAX_NUMBER_OF_BLOCKELEMENTS][MAX_NUMBER_OF_TYPES][MAX_NUMBER_OF_TYPES];
  double PosCorrSum[MAX_NUMBER_OF_BLOCKS][MAX_NUMBER_OF_BLOCKELEMENTS][MAX_NUMBER_OF_TYPES][3];
  double TmpPos [MAX_NUMBER_OF_ATOMS][4];
  int BlockLength[MAX_NUMBER_OF_BLOCKS];
  double count_samples[MAX_NUMBER_OF_BLOCKS][MAX_NUMBER_OF_BLOCKELEMENTS];
  int Groups[MAX_NUMBER_OF_ATOMS][3]; 		// The group id of each atom
  double  sum_dself ,  sum_cij , msdcounting;   // Sum of msd
  // Finding corresponding groups to each atom
  int tmpgroup[32][2];
  int tmpnumgroup ;
  int tmpfoundgroup ;
  //int WriteFileEvery ;
  int WriteFileEvery ;
  FILE *FilePtrii;	
  FILE *FilePtrij;	

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute msd fix ID

Self-explanatory.

*/
