/* ----------------------------------------------------------------------
   compute rdf/ext is a child class of "compute", developed based on two 
   classes of "compute msd" and "compute rdf", provided by LAMMPS.
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

#ifdef COMPUTE_CLASS

ComputeStyle(rdf/ext,ComputeRDFExtend)

#else

#ifndef LMP_COMPUTE_RDF__Extend_H
#define LMP_COMPUTE_RDF__Extend_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRDFExtend : public Compute {
 public:
  ComputeRDFExtend(class LAMMPS *, int, char **);
  virtual ~ComputeRDFExtend();
  void init();
  virtual double compute_scalar();  // remember what this returns

 protected:
  #define MAXGROUPS 32
  #define SQR(x) ((x)*(x))
  #define CUBE(x) ((x)*(x)*(x))
  // The parameters related to MPI and gathering the positions
  int me; // This is for the definition of rank of processor in LAMMPS
  int nprocs; // This is for the number of processors
  int *tmprecvcnts; // Number of atoms per each processor [array]
  int *recvcnts; // Number of atoms per each processor [array]
  int *displs; // Displacement array for receiving the correct data set from each processor [array]
  double *sendbuff;
  double *recvbuff;

  // New static arrays or dynamic arrays
  double **TmpPos;
  int **Groups;

  // RDF
  int binsize;
  double r, r2, sf, pairs;
  double rin , rout, vin, vout;
  int whichbin;
  double NumberOfParticles;
  double ***Gg;
  double ***Gglocal;
  double boxsize;
  double maxdist;
  double hBox;
  double Delta;
  double Delta_1;
  double Cube_deltaboxsize;
  double Sqr_deltaboxsize;
  double Ggt;
  double dr[3]; 
  // RDF Correction
  double ***PartSum;
  double SphereVolFrac;
  double gr_correction;
  // The things in common with the other one
  double timeinterval;
  int count ;
  int WriteFileEvery;
  FILE *FilePtr;		// The file handle for writing the calculated RDF 
  int samplerate;
  int realatom;
  char *filename;  // filename
  
  // Finding corresponding groups to each atom
  int tmpgroup[MAXGROUPS][3];  // group_id ; groupbits ; number_of_particles
  int tmpnumgroup ;
  int tmpfoundgroup ;
  int groupatom1, groupatom2; 
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  

E: Not a cubic simulation box for rdf/ext

Self-explanatory.

*/
