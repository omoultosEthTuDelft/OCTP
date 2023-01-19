/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(rdf/force,ComputeRDFForce);
// clang-format on
#else

#ifndef LMP_COMPUTE_RDFFORCE_H
#define LMP_COMPUTE_RDFFORCE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRDFForce : public Compute {
 public:
  ComputeRDFForce(class LAMMPS *, int, char **);
  ~ComputeRDFForce() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;

 private:
  int nbin;                // # of rdf bins
  int cutflag;             // user cutoff flag
  int npairs;              // # of rdf pairs
  int icompute_temp;        // compute index for temperature
  double delr, delrinv;    // bin width and its inverse
  double cutoff_user;      // user-specified cutoff
  double mycutneigh;       // user-specified cutoff + neighbor skin
  int ***rdfpair;          // map 2 type pair to rdf pair for each histo
  int **nrdfpair;          // # of histograms for each type pair
  int *ilo, *ihi, *jlo, *jhi;
  double **hist;       // histogram bins
  double **histall;    // summed histogram bins across all procs

  int *typecount;
  int *icount, *jcount;
  int *duplicates;

  class NeighList *list;    // half neighbor list
  void init_norm();
  bigint natoms_old;
};

}    // namespace LAMMPS_NS

#endif
#endif
