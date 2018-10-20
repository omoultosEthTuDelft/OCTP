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

/*
    compute_ordern_viscosity is based on the compute_pressure command.
*/


#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "compute_ordern_viscosity.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeOrdernViscosity::ComputeOrdernViscosity(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute pressure command");
  if (igroup) error->all(FLERR,"Compute pressure must use group all");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;
  pressflag = 1;
  timeflag = 1;

  // store temperature ID used by pressure computation
  // insure it is valid for temperature computation

  if (strcmp(arg[3],"NULL") == 0) id_temp = NULL;
  else {
    int n = strlen(arg[3]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[3]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute pressure temperature ID");
    if (modify->compute[icompute]->tempflag == 0)
      error->all(FLERR,"Compute pressure temperature ID does not "
                 "compute temperature");
  }

  // process optional args

  if (narg == 4) {
    keflag = 1;
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = fixflag = 1;
  } else {
    keflag = 0;
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"ke") == 0) keflag = 1;
      else if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedralflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) improperflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg],"fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg],"virial") == 0) {
        pairflag = 1;
        bondflag = angleflag = dihedralflag = improperflag = 1;
        kspaceflag = fixflag = 1;
      } else error->all(FLERR,"Illegal compute pressure command");
      iarg++;
    }
  }

  // error check

  if (keflag && id_temp == NULL)
    error->all(FLERR,"Compute pressure requires temperature ID "
	       "to include kinetic energy");

  vector = new double[6];
  nvirial = 0;
  vptr = NULL;
  // SJ: Start editing
  MPI_Comm_rank(world, &me);    
  MPI_Comm_size(world, &nprocs);
  // SJ: end editing
}

/* ---------------------------------------------------------------------- */

ComputeOrdernViscosity::~ComputeOrdernViscosity()
{
  delete [] id_temp;
  delete [] vector;
  delete [] vptr;
}

/* ---------------------------------------------------------------------- */

void ComputeOrdernViscosity::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;

  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  if (keflag) {
    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute pressure temperature ID");
    temperature = modify->compute[icompute];
  }

  // detect contributions to virial
  // vptr points to all virial[6] contributions

  delete [] vptr;
  nvirial = 0;
  vptr = NULL;

  if (pairflag && force->pair) nvirial++;
  if (bondflag && atom->molecular && force->bond) nvirial++;
  if (angleflag && atom->molecular && force->angle) nvirial++;
  if (dihedralflag && atom->molecular && force->dihedral) nvirial++;
  if (improperflag && atom->molecular && force->improper) nvirial++;
  if (fixflag)
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->virial_flag) nvirial++;

  if (nvirial) {
    vptr = new double*[nvirial];
    nvirial = 0;
    if (pairflag && force->pair) vptr[nvirial++] = force->pair->virial;
    if (bondflag && force->bond) vptr[nvirial++] = force->bond->virial;
    if (angleflag && force->angle) vptr[nvirial++] = force->angle->virial;
    if (dihedralflag && force->dihedral)
      vptr[nvirial++] = force->dihedral->virial;
    if (improperflag && force->improper)
      vptr[nvirial++] = force->improper->virial;
    if (fixflag)
      for (int i = 0; i < modify->nfix; i++)
        if (modify->fix[i]->virial_flag)
          vptr[nvirial++] = modify->fix[i]->virial;
  }

  // flag Kspace contribution separately, since not summed across procs

  if (kspaceflag && force->kspace) kspace_virial = force->kspace->virial;
  else kspace_virial = NULL;

  // SJ: Start editing
  WriteFileEvery = 10000; 
  count = 0;
  numpressure = 0;
  sumpressure = 0;
  avgpressure = 0;
  countint = 1;
  NumberOfBlock = 1;
  for (int i = 0; i < NUMBER_OF_BLOCKS; i++)
  {
    BlockLength[i]=1;
  }
  // SJ: end editing
}

/* ----------------------------------------------------------------------
   compute pressure tensor
   assume KE tensor has already been computed
------------------------------------------------------------------------- */

double ComputeOrdernViscosity::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  if (force->kspace && kspace_virial && force->kspace->scalar_pressure_flag)
    error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' for "
	       "tensor components with kspace_style msm");

  // invoke temperature if it hasn't been already

  double t;
  double *ke_tensor;
  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
    {
      temperature->compute_vector();
      t = temperature->compute_scalar();
    }
    else
    {
      t = temperature->scalar;
    }
    ke_tensor = temperature->vector;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6,3);
    if (keflag) {
      for (int i = 0; i < 6; i++)
        vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
    } else
      for (int i = 0; i < 6; i++)
        vector[i] = virial[i] * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4,2);
    if (keflag) {
      vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
      vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
      vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    } else {
      vector[0] = virial[0] * inv_volume * nktv2p;
      vector[1] = virial[1] * inv_volume * nktv2p;
      vector[3] = virial[3] * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    }
  }
  // Start: SJ editing
  pressure = (vector[0]+vector[1]+vector[2])/3.0;
  Pxx = vector[0];
  Pyy = vector[1];
  Pzz = vector[2];
  Pxy = vector[3];
  Pxz = vector[4];
  Pyz = vector[5];
  Pxx_s = Pxx - pressure;
  Pyy_s = Pyy - pressure;
  Pzz_s = Pzz - pressure; 
  if (me == 0)		// Do everything on the master processor
  {
    sumpressure += pressure;
    numpressure += 1.0;
    if (count == 0)
    {
      tmpa[0] = Pxx_s;
      tmpa[1] = Pyy_s;
      tmpa[2] = Pzz_s;
      tmpa[3] = Pxy;
      tmpa[4] = Pxz;
      tmpa[5] = Pyz;
      tmpa[6] = pressure;
      accint[0] = accint[1] = accint[2] = accint[3] = accint[4] = accint[5] = accint[6] = 0.0;
      for (int currentblock = 0; currentblock<NUMBER_OF_BLOCKS; currentblock++)
      {
        for (int k = 0; k < NUMBER_OF_BLOCKELEMENTS; k++)
	{
	  for (int i = 0; i < 7; i++)
	  {
	    oldint[currentblock][k][i] = 0.0;
	    samples[currentblock][k][i] = 0.0;
	    samplecount[currentblock][k][i] = 0.0;
	  }
	  samples[currentblock][k][7] = 0.0;
	  samplecount[currentblock][k][7] = 0.0;
	}
      }
    } 
    else if ( count % 2 == 1)
    {
      if (count == 1)
      {
        timeinterval = (double) (update->ntimestep - 0)*(update->dt);	// CHECK IT LATER FOR BETTER SAMPLING #############
	samplerate = (update->ntimestep);
      }
      tmpc[0] = Pxx_s;
      tmpc[1] = Pyy_s;
      tmpc[2] = Pzz_s;
      tmpc[3] = Pxy;
      tmpc[4] = Pxz;
      tmpc[5] = Pyz;
      tmpc[6] = pressure;
    }
    else
    {
      tmpb[0] = Pxx_s;
      tmpb[1] = Pyy_s;
      tmpb[2] = Pzz_s;
      tmpb[3] = Pxy;
      tmpb[4] = Pxz;
      tmpb[5] = Pyz;
      tmpb[6] = pressure;
      // Simpson's rule to integrate int(P.dt) from "t-2dt" to "t"
      for (int i = 0; i < 7; i++)
      {
        lastint[i] = (2.0*timeinterval)*(tmpa[i]+4*tmpc[i]+tmpb[i])/6.0;
      }
      // Updating the values of tmpa, tmpb, and tmpc
      for (int i = 0; i < 7; i++)
      {
        tmpa[i] = tmpb[i];
        tmpb[i] = 0.0;
        tmpc[i] = 0.0;
      }
      // Add previous integral to the whole from "0" to "t"
      for (int i = 0; i < 7; i++)
      {
        accint[i] += lastint[i];
      }
      // loop over all blocks
      int i = countint/(pow(NUMBER_OF_BLOCKELEMENTS,NumberOfBlock));
      while (i != 0)
      {
        NumberOfBlock++;
        i /= NUMBER_OF_BLOCKELEMENTS;
      }
      for(int currentblock = 0; currentblock < NumberOfBlock; currentblock++)
      {
        if ((countint)%((int)pow(NUMBER_OF_BLOCKELEMENTS,currentblock))==0)
	{
	  CurrentBlockLength = MIN(BlockLength[currentblock],NUMBER_OF_BLOCKELEMENTS);
	  for (int k=NUMBER_OF_BLOCKELEMENTS-1; k>=NUMBER_OF_BLOCKELEMENTS-CurrentBlockLength; k--)
	  {
            //if ((countint)%((k+1)*((int)pow(NUMBER_OF_BLOCKELEMENTS,currentblock)))==0)
	    {
	      for (int i = 0; i < 7; i++)
              {
	        // correction for dimensions in total, and adding kb and volume (WITHOUT TEMPERATURE)
		integral = accint[i]-oldint[currentblock][k][i];
	        //samples[currentblock][k][i] += (integral*integral)*(1.0/inv_volume/2.0/boltz)*(1.0/nktv2p);
	        samples[currentblock][k][i] += (integral*integral);	// ADD THE COEFFICIENT LATER
	        samplecount[currentblock][k][i] += 1.0;
		if (i == 6)
		{
		  samples[currentblock][k][7] += (integral);
	          samplecount[currentblock][k][7] += 1.0;
		}
	      }
	    }
	  }
	  BlockLength[currentblock]++;
	  for (int i=0; i < 7; i++)
	  {
	    for (int k=1; k < NUMBER_OF_BLOCKELEMENTS; k++)
	      oldint[currentblock][k-1][i] = oldint[currentblock][k][i] ;
	    oldint[currentblock][NUMBER_OF_BLOCKELEMENTS-1][i] = accint[i]; 
	  }
	}
      }
      // Writing in a file
      if ( ((countint)%WriteFileEvery==0) || (((update->endstep)-(update->ntimestep))<samplerate) )
      {
        FilePtr = fopen("bulkviscosity.dat","w");
	fprintf(FilePtr,"#  NOTE: the coefficients are without temperature in the denominator\n");
	fprintf(FilePtr,"#  Time \tmu_xx\tmu_yy\tmu_zz\tmu_xy\tmu_xz\tmu_yz\tmu_shear\tmu_total\tsample_num\tbulk_visc\n");
	double coefficient = (1.0/inv_volume/2.0/boltz)*(1.0/nktv2p);
	for (int currentblock=0; currentblock<MIN(NUMBER_OF_BLOCKS,NumberOfBlock); currentblock++)
	{
	  if (currentblock == MIN(NUMBER_OF_BLOCKS,NumberOfBlock)-1)
	    CurrentBlockLength = MIN(BlockLength[currentblock],NUMBER_OF_BLOCKELEMENTS)-1;
	  else
	    CurrentBlockLength = MIN(BlockLength[currentblock],NUMBER_OF_BLOCKELEMENTS);
	  for (int k = 1; k <= CurrentBlockLength; k++)	// Just neglect the first data on the right (k=2)
	  {
	    double Deltat = ((1.0*k)*(2.0*timeinterval)*pow(NUMBER_OF_BLOCKELEMENTS,currentblock));
            fprintf(FilePtr,"%-9g\t",Deltat);
	    // SHEAR VISCOSITY
	    double totalstress = 0.0;
	    double totalshear = 0.0;
	    for (int i = 0; i< 6; i++)
	    {
	      double stresscomp = coefficient*samples[currentblock][NUMBER_OF_BLOCKELEMENTS-k][i]/samplecount[currentblock][NUMBER_OF_BLOCKELEMENTS-k][i];
	      fprintf(FilePtr,"%-10g\t",stresscomp);
	      if (i<3)
	      {
	        totalstress += 1.0*1.0*stresscomp/10.0;		// should it be 1.0 or 4.0/3.0 ?? ########
		//// Tenney2010: 4/3 is the IMPLICIT contribution factor #################################
	      }
	      else
	      {
	        totalstress += 2.0*1.0*stresscomp/10.0;
	        totalshear += stresscomp/3.0;
	      }
	    }
	    fprintf(FilePtr,"%-10g\t%-10g\t",totalshear,totalstress);
	    // BULK VISCOSITY
	    double intp2 = samples[currentblock][NUMBER_OF_BLOCKELEMENTS-k][6]/samplecount[currentblock][NUMBER_OF_BLOCKELEMENTS-k][6];
	    double intp = samples[currentblock][NUMBER_OF_BLOCKELEMENTS-k][7]/samplecount[currentblock][NUMBER_OF_BLOCKELEMENTS-k][7];
	    avgpressure = sumpressure / numpressure;
	    double bulkvis = coefficient * ( intp2 - 2.0*intp*(avgpressure*Deltat) + (avgpressure*Deltat)*(avgpressure*Deltat) );
	    fprintf(FilePtr,"%-10g\t",bulkvis);
	    fprintf(FilePtr,"%-10g\n",samplecount[currentblock][NUMBER_OF_BLOCKELEMENTS-k][0]);
	   

	  }
	}
	fclose(FilePtr);
      }
      countint++;
    }
  count++;
  }
  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    if (keflag)
      scalar = (temperature->dof * boltz * t + virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
  } else {
    if (keflag)
      scalar = (temperature->dof * boltz * t + virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
  }

  return scalar;
  
  // End: SJ editing
}

/* ---------------------------------------------------------------------- */

void ComputeOrdernViscosity::virial_compute(int n, int ndiag)
{
  int i,j;
  double v[6],*vcomponent;

  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from forces and fixes

  for (j = 0; j < nvirial; j++) {
    vcomponent = vptr[j];
    for (i = 0; i < n; i++) v[i] += vcomponent[i];
  }

  // sum virial across procs

  MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,world);

  // KSpace virial contribution is already summed across procs

  if (kspace_virial)
    for (i = 0; i < n; i++) virial[i] += kspace_virial[i];

  // LJ long-range tail correction, only if pair contributions are included

  if (force->pair && pairflag && force->pair->tail_flag)
    for (i = 0; i < ndiag; i++) virial[i] += force->pair->ptail * inv_volume;
}

/* ---------------------------------------------------------------------- */

void ComputeOrdernViscosity::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_temp;
  int n = strlen(id_new) + 1;
  id_temp = new char[n];
  strcpy(id_temp,id_new);
}
