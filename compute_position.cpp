/* ----------------------------------------------------------------------
   compute position is a child class of "compute", developed based on the
   "compute msd", provided by LAMMPS.
   This command is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_position.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePosition::ComputePosition(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute position command");

  vector_flag = 1;
  size_vector = atom->natoms*5;
  extvector = 0;
  create_attribute = 0;

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  MPI_Comm_rank(world, &me); // Assigning the rank of a molecule for each core
  MPI_Comm_size(world, &nprocs); // The number of processors 
  tmprecvcnts = new int[nprocs]; 
  recvcnts = new int[nprocs];
  displs = new int[nprocs];
  sendbuff = new double[size_vector];
  vector = new double[size_vector]; // The output global vector
}

/* ---------------------------------------------------------------------- */

ComputePosition::~ComputePosition()
{
  delete [] vector;
  delete [] tmprecvcnts;
  delete [] recvcnts;
  delete [] displs;
  delete [] sendbuff;
}

/* ---------------------------------------------------------------------- */

void ComputePosition::init()
{
  for (int ii = 0; ii < size_vector; ii++)
  {
      vector[ii]=-1.0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePosition::compute_vector()
{
  invoked_vector = update->ntimestep;

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int xbox,ybox,zbox;
  double xtmp, ytmp, ztmp;

  if (domain->triclinic == 0) 
  {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) 
      {  
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
	      xtmp = x[i][0] + xbox*xprd;
	      ytmp = x[i][1] + ybox*yprd;
	      ztmp = x[i][2] + zbox*zprd;
	      sendbuff[5*i] = xtmp;
	      sendbuff[5*i+1] = ytmp;
	      sendbuff[5*i+2] = ztmp;
	      sendbuff[5*i+3] = (double) ((atom->tag[i]) - 1);
	      sendbuff[5*i+4] = (double) mask[i];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) 
      {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        xtmp = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
        ytmp = x[i][1] + h[1]*ybox + h[3]*zbox;
        ztmp = x[i][2] + h[2]*zbox;
	      sendbuff[5*i] = xtmp;
	      sendbuff[5*i+1] = ytmp;
	      sendbuff[5*i+2] = ztmp;
	      sendbuff[5*i+3] = (double) ((atom->tag[i]) - 1);
	      sendbuff[5*i+4] = (double) mask[i];
      }
  }
  // Sending the position of all atoms from all cores (0, 1, ..., n) to all cores
  for (int ii = 0; ii < nprocs; ii++) 
    tmprecvcnts[ii] = (ii == me) ?  5*nlocal : 0 ;
  MPI_Allreduce(tmprecvcnts, recvcnts, nprocs, MPI_INT , MPI_SUM, world);

  for (int ii = 0; ii < nprocs; ii++)
  {
    displs[ii] = 0;
  }
  for (int ii = 1; ii < nprocs; ii++)
    for (int jj = ii; jj < nprocs ; jj++)
    {
      displs[jj] += recvcnts[ii-1];
    }
  MPI_Allgatherv(sendbuff,5*nlocal, MPI_DOUBLE, vector, recvcnts, displs, MPI_DOUBLE, world);
}
