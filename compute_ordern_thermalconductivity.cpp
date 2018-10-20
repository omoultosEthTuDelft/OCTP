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

/* ----------------------------------------------------------------------
   Contributing authors: German Samolyuk (ORNL) and
                         Mario Pinto (Computational Research Lab, Pune, India)
------------------------------------------------------------------------- */

/*
    compute_ordern_thermalconductivity is based on the compute_heatflux command.
*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "compute_ordern_thermalconductivity.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

ComputeOrderNThermalConductivity::ComputeOrderNThermalConductivity(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal compute thermalconductivity command");

  vector_flag = 1;
  size_vector = 6;
  extvector = 1;

  // store ke/atom, pe/atom, stress/atom IDs used by heat flux computation
  // insure they are valid for these computations

  int n = strlen(arg[3]) + 1;
  id_ke = new char[n];
  strcpy(id_ke,arg[3]);

  n = strlen(arg[4]) + 1;
  id_pe = new char[n];
  strcpy(id_pe,arg[4]);

  n = strlen(arg[5]) + 1;
  id_stress = new char[n];
  strcpy(id_stress,arg[5]);

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  int istress = modify->find_compute(id_stress);
  if (ike < 0 || ipe < 0 || istress < 0)
    error->all(FLERR,"Could not find compute thermalconductivity compute ID");
  if (strcmp(modify->compute[ike]->style,"ke/atom") != 0)
    error->all(FLERR,"Compute thermalconductivity compute ID does not compute ke/atom");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Compute thermalconductivity compute ID does not compute pe/atom");
  if (modify->compute[istress]->pressatomflag == 0)
    error->all(FLERR,"Compute thermalconductivity compute ID does not compute stress/atom");

  vector = new double[6];
  
  // SJ: Start editing
  MPI_Comm_rank(world, &me);    
  MPI_Comm_size(world, &nprocs);
  // SJ: end editing
}

/* ---------------------------------------------------------------------- */

ComputeOrderNThermalConductivity::~ComputeOrderNThermalConductivity()
{
  delete [] id_ke;
  delete [] id_pe;
  delete [] id_stress;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeOrderNThermalConductivity::init()
{
  // error checks
  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  int istress = modify->find_compute(id_stress);
  if (ike < 0 || ipe < 0 || istress < 0)
    error->all(FLERR,"Could not find compute heat/flux compute ID");

  c_ke = modify->compute[ike];
  c_pe = modify->compute[ipe];
  c_stress = modify->compute[istress];
  
  // SJ: Start editing
  boltz = force->boltz;
  dimension = domain->dimension;
  
  WriteFileEvery = 10000; 
  count = 0;
  countint = 1;
  NumberOfBlock = 1;
  for (int i = 0; i < NUMBER_OF_BLOCKS; i++)
  {
    BlockLength[i]=1;
  }
  // SJ: end editing  
  
  
}

/* ---------------------------------------------------------------------- */

void ComputeOrderNThermalConductivity::compute_vector()
{
  invoked_vector = update->ntimestep;

  // invoke 3 computes if they haven't been already

  if (!(c_ke->invoked_flag & INVOKED_PERATOM)) {
    c_ke->compute_peratom();
    c_ke->invoked_flag |= INVOKED_PERATOM;
  }
  if (!(c_pe->invoked_flag & INVOKED_PERATOM)) {
    c_pe->compute_peratom();
    c_pe->invoked_flag |= INVOKED_PERATOM;
  }
  if (!(c_stress->invoked_flag & INVOKED_PERATOM)) {
    c_stress->compute_peratom();
    c_stress->invoked_flag |= INVOKED_PERATOM;
  }

  // heat flux vector = jc[3] + jv[3]
  // jc[3] = convective portion of heat flux = sum_i (ke_i + pe_i) v_i[3]
  // jv[3] = virial portion of heat flux = sum_i (stress_tensor_i . v_i[3])
  // normalization by volume is not included

  double *ke = c_ke->vector_atom;
  double *pe = c_pe->vector_atom;
  double **stress = c_stress->array_atom;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double jc[3] = {0.0,0.0,0.0};
  double jv[3] = {0.0,0.0,0.0};
  double eng;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      eng = pe[i] + ke[i];
      jc[0] += eng*v[i][0];
      jc[1] += eng*v[i][1];
      jc[2] += eng*v[i][2];
      jv[0] -= stress[i][0]*v[i][0] + stress[i][3]*v[i][1] + stress[i][4]*v[i][2];
      jv[1] -= stress[i][3]*v[i][0] + stress[i][1]*v[i][1] + stress[i][5]*v[i][2];
      jv[2] -= stress[i][4]*v[i][0] + stress[i][5]*v[i][1] + stress[i][2]*v[i][2];
    }
  }

  // convert jv from stress*volume to energy units via nktv2p factor

  double nktv2p = force->nktv2p;
  jv[0] /= nktv2p;
  jv[1] /= nktv2p;
  jv[2] /= nktv2p;

  // sum across all procs
  // 1st 3 terms are total heat flux
  // 2nd 3 terms are just conductive portion

  double data[6] = {jc[0]+jv[0],jc[1]+jv[1],jc[2]+jv[2],jc[0],jc[1],jc[2]};
  MPI_Allreduce(data,vector,6,MPI_DOUBLE,MPI_SUM,world);
  
  
  
  
  // Start: SJ editing
  if (dimension == 3)
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  else
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
  Jtx = vector[0];
  Jty = vector[1];
  Jtz = vector[2];
  Jcx = vector[3];
  Jcy = vector[4];
  Jcz = vector[5];
  if (me == 0)        // Do everything on the master processor
  {
    if (count == 0)
    {
      tmpa[0] = Jtx;
      tmpa[1] = Jty;
      tmpa[2] = Jtz;
      tmpa[3] = Jcx;
      tmpa[4] = Jcy;
      tmpa[5] = Jcz;
      accint[0] = accint[1] = accint[2] = accint[3] = accint[4] = accint[5] = 0.0;
      for (int currentblock = 0; currentblock<NUMBER_OF_BLOCKS; currentblock++)
      {
        for (int k = 0; k < NUMBER_OF_BLOCKELEMENTS; k++)
        {
          for (int i = 0; i < 6; i++)
          {
            oldint[currentblock][k][i] = 0.0;
            samples[currentblock][k][i] = 0.0;
            samplecount[currentblock][k][i] = 0.0;
          }
        }
      }
    } 
    else if ( count % 2 == 1)
    {
      if (count == 1)
      {
        timeinterval = (double) (update->ntimestep - 0)*(update->dt);
        samplerate = (update->ntimestep);
      }
      tmpc[0] = Jtx;
      tmpc[1] = Jty;
      tmpc[2] = Jtz;
      tmpc[3] = Jcx;
      tmpc[4] = Jcy;
      tmpc[5] = Jcz;
    }
    else
    {
      tmpb[0] = Jtx;
      tmpb[1] = Jty;
      tmpb[2] = Jtz;
      tmpb[3] = Jcx;
      tmpb[4] = Jcy;
      tmpb[5] = Jcz;
      // Simpson's rule to integrate int(P.dt) from "t-2dt" to "t"
      for (int i = 0; i < 6; i++)
      {
        lastint[i] = (2.0*timeinterval)*(tmpa[i]+4*tmpc[i]+tmpb[i])/6.0;
      }
      // Updating the values of tmpa, tmpb, and tmpc
      for (int i = 0; i < 6; i++)
      {
        tmpa[i] = tmpb[i];
        tmpb[i] = 0.0;
        tmpc[i] = 0.0;
      }
      // Add previous integral to the whole from "0" to "t"
      for (int i = 0; i < 6; i++)
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
            {
              for (int i = 0; i < 6; i++)
              {
                integral = accint[i]-oldint[currentblock][k][i];
                samples[currentblock][k][i] += (integral*integral);	// ADD THE COEFFICIENT LATER
                samplecount[currentblock][k][i] += 1.0;
              }
            }
          }
          BlockLength[currentblock]++;
          for (int i=0; i < 6; i++)
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
        // coefficient = 1/(2*V*kB). NOTE: add 1/(T^2) during the post processing
        double coefficient = ((1.0*inv_volume)/2.0/boltz);
        FilePtr = fopen("thermalconductivity.dat","w");
        fprintf(FilePtr,"# NOTE: the coefficients are without temperature^2 in the denominator\n");
        fprintf(FilePtr,"# Time \tJtot_x\tJtot_y\tJtot_z\tJconv_x\tJconv_y\tJconv_z\tJtotal\tJconvective\tsample_num\n");
        for (int currentblock=0; currentblock<MIN(NUMBER_OF_BLOCKS,NumberOfBlock); currentblock++)
        {
          if (currentblock == MIN(NUMBER_OF_BLOCKS,NumberOfBlock)-1)
            CurrentBlockLength = MIN(BlockLength[currentblock],NUMBER_OF_BLOCKELEMENTS)-1;
          else
            CurrentBlockLength = MIN(BlockLength[currentblock],NUMBER_OF_BLOCKELEMENTS);
          for (int k = 1; k <= CurrentBlockLength; k++)    // Just neglect the first data on the right (k=2)
          {
            double Deltat = ((1.0*k)*(2.0*timeinterval)*pow(NUMBER_OF_BLOCKELEMENTS,currentblock));
            fprintf(FilePtr,"%-9g\t",Deltat);
            double Jtotal = 0.0;
            double Jconvective = 0.0;
            for (int i = 0; i< 6; i++)
            {
              double fluxcomp = coefficient*samples[currentblock][NUMBER_OF_BLOCKELEMENTS-k][i]/samplecount[currentblock][NUMBER_OF_BLOCKELEMENTS-k][i];
              fprintf(FilePtr,"%-10g\t",fluxcomp);
              if (i<3)
              {
                Jtotal += 1.0*fluxcomp/3.0; 
              }
              else
              {
                Jconvective += 1.0*fluxcomp/3.0;
              }
            }
            fprintf(FilePtr,"%-10g\t%-10g\t%-10g\n",Jtotal,Jconvective,samplecount[currentblock][NUMBER_OF_BLOCKELEMENTS-k][0]);
          }
        }
        fclose(FilePtr);
      }
      countint++;
    }
  count++;
  }  
  // End: SJ editing
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
