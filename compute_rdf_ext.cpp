/* ----------------------------------------------------------------------
   compute rdf/ext is a child class of "compute", developed based on two 
   classes of "compute msd" and "compute rdf", provided by LAMMPS.
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
#include <string.h>
#include "compute_rdf_ext.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix_store.h"
#include "error.h"

#include <iostream>
#include "memory.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeRDFExtend::ComputeRDFExtend(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute rdf/ext command");

  scalar_flag = 1;
  //size_scalar = 4;
  extscalar = 0;
  create_attribute = 0;
  char filen[] = "rdfs.dat";
  filename = new char[strlen(filen)+1];
  strcpy(filename,filen); 

  // optional args
  binsize = 1000;     // initial number of bins
  WriteFileEvery = 1000;  // write files every
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"Nbin") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute rdf/ext command");
      binsize = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Nwrite") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute rdf/ext command");
      WriteFileEvery = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ordern command");
      delete [] filename;
      filename = new char[strlen(arg[iarg+1])+1];
      strcpy(filename,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute rdf/ext command");
  }


  MPI_Comm_rank(world, &me);            // Assigning the rank of a molecule for each core
  MPI_Comm_size(world, &nprocs);        // The number of processors 
  tmprecvcnts = new int[nprocs]; 
  recvcnts = new int[nprocs];
  displs = new int[nprocs];
  sendbuff = new double[atom->natoms*5];
  recvbuff = new double[atom->natoms*5];
  memory->create(TmpPos,atom->natoms,4,"compute/rdf/ext:TmpPos");
  memory->create(Groups,atom->natoms,3,"compute/rdf/ext:Groups");
  memory->create(Gg,binsize,MAXGROUPS,MAXGROUPS,"compute/rdf/ext:Gg");
  memory->create(Gglocal,binsize,MAXGROUPS,MAXGROUPS,"compute/rdf/ext:Gglocal");
  memory->create(PartSum,binsize,MAXGROUPS,MAXGROUPS,"compute/rdf/ext:PartSum");

}

/* ---------------------------------------------------------------------- */

ComputeRDFExtend::~ComputeRDFExtend()
{
  
  delete [] tmprecvcnts;
  delete [] recvcnts;
  delete [] displs;
  delete [] sendbuff;
  delete [] recvbuff;
  memory->destroy(TmpPos);
  memory->destroy(Groups);
  memory->destroy(Gg);
  memory->destroy(Gglocal);
  memory->destroy(PartSum);
  delete [] filename; 
}

/* ---------------------------------------------------------------------- */

void ComputeRDFExtend::init()
{
  
  tmpnumgroup = 0;
  count = 0;
  for (int ii=0; ii<binsize; ii++)	// be sure everything (0 to binsize) is zero
  {
    for (int jj = 1; jj<MAXGROUPS;jj++)
    {
      for (int kk = 1; kk<MAXGROUPS; kk++)
      {
        Gglocal[ii][jj][kk] = 0;
	      PartSum[ii][jj][kk] = 0;
      }
    }
  }
  Ggt = 0.0;
  boxsize = domain->xprd;
  hBox = 0.5 * boxsize ;
  maxdist = hBox * sqrt(2.0);
  Delta = maxdist/binsize;
  Delta_1 = 1.0 / Delta;
  Cube_deltaboxsize = CUBE(Delta/boxsize);
  Sqr_deltaboxsize = SQR(Delta/boxsize);

  
  
}

/* ---------------------------------------------------------------------- */

double ComputeRDFExtend::compute_scalar()
{  
  invoked_scalar = update->ntimestep;

  // dx,dy,dz = displacement of atom from reference position
  // reference unwrapped position is stored by fix
  // for triclinic, need to unwrap current atom coord via h matrix

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double *h = domain->h;

  int atomgroup , atomgroupbit;            // The ID of the group, its bit
  int numgroup = group->ngroup;            // The total # of available groups
  int natom = atom->natoms;                // total # of atoms in system, could be 0
  int *currentgroupbit = group->bitmask;    // The bitmask of a group

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) { 
	      sendbuff[5*i] = x[i][0];
	      sendbuff[5*i+1] = x[i][1];
	      sendbuff[5*i+2] = x[i][2];
	      sendbuff[5*i+3] = (double) ((atom->tag[i]) - 1);
	      sendbuff[5*i+4] = (double) mask[i];
      }
  } else error->all(FLERR,"Not a cubic simulation box for rdf/ext");
  
  // Sending data from all cores (0, 1, ..., n) to all cores
  for (int ii = 0; ii < nprocs; ii++) 
    tmprecvcnts[ii] = (ii == me) ?  5*nlocal : 0 ;
  MPI_Allreduce(tmprecvcnts, recvcnts, nprocs, MPI_INT , MPI_SUM, world);
  for (int ii = 0; ii < nprocs; ii++)
    for (int jj = ii; jj < nprocs ; jj++)
    {
      if (ii == 0) displs[jj] = 0;
      else displs[jj] += recvcnts[ii-1];
    }
  MPI_Allgatherv(sendbuff,5*nlocal, MPI_DOUBLE, recvbuff, recvcnts, displs, MPI_DOUBLE, world);
  // Definining the group number of each molecules 
  for (int ii = 0; ii < natom; ii++)
  {
    realatom = (int) (recvbuff[5*ii+3]+0.1);
    TmpPos[realatom][0] = recvbuff[5*ii];
    TmpPos[realatom][1] = recvbuff[5*ii+1];
    TmpPos[realatom][2] = recvbuff[5*ii+2];
    // Finding corresponding groups to each atom at the first time
    if (count == 0)   // only run during the first time step
    {
      TmpPos[realatom][3] = recvbuff[5*ii+4]+0.1;
      tmpfoundgroup = 0;
      if (tmpnumgroup > 0 ) // First try to find the group number if the group is available
      {
    	for (int jj = 1; jj <= tmpnumgroup; jj++)
        {
          if (((int) TmpPos[realatom][3]) & tmpgroup[jj][1])
          {
            Groups[realatom][0] = jj;
            Groups[realatom][1] = tmpgroup[jj][1];
            Groups[realatom][2] = tmpgroup[jj][0];
            tmpgroup[jj][2] += 1;
            tmpfoundgroup = 1;
            break;
          }
        }
      }
      if (tmpfoundgroup == 0) // If couldn't find any, tries to find the new group
      {
        for (int jj = 1; jj < MAXGROUPS; jj++)
        {
          if ( ((int) TmpPos[realatom][3]) & group->bitmask[jj])
          {
            tmpnumgroup++;
            tmpfoundgroup = 1;
            tmpgroup[tmpnumgroup][0] = jj;
            tmpgroup[tmpnumgroup][1] = group->bitmask[jj];
            tmpgroup[tmpnumgroup][2] = 1;
            Groups[realatom][0] = tmpnumgroup; 
            Groups[realatom][1] = group->bitmask[jj]; 
            Groups[realatom][2] = jj;
            break; 
          }
        }
        if (tmpfoundgroup == 0) // Finally, this atom does not belong to the available groups
        {
          Groups[realatom][0] = -1; 
          Groups[realatom][1] = -1; 
          Groups[realatom][2] = -1;
        }
      }
    }
  }
  // CALCULATING RDF
  if (count == 0)  
  {
    samplerate = (update->ntimestep);
  } else if (count == 1) {
    samplerate = (update->ntimestep) - samplerate;
    timeinterval = (double) (samplerate)*(update->dt);
  }
  Ggt+=1.0;
  for(int i=0;i<nlocal;i++)
  {
    realatom = (atom->tag[i]) - 1;
    groupatom1 = Groups[realatom][0];  
    if (groupatom1>0)
    {
      for(int j=realatom+1;j<natom;j++)
      {
        groupatom2 = Groups[j][0];      
        if ( (groupatom2>0) )
        {
          dr[0] = TmpPos[realatom][0] - TmpPos[j][0];
          dr[1] = TmpPos[realatom][1] - TmpPos[j][1];
          dr[2] = TmpPos[realatom][2] - TmpPos[j][2];

          dr[0] -= boxsize * rint(dr[0]/boxsize);
          dr[1] -= boxsize * rint(dr[1]/boxsize);
          dr[2] -= boxsize * rint(dr[2]/boxsize);

          r2 = sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );
    
          if (r2 < maxdist)
          { 
            whichbin = (int) (r2*Delta_1);
            Gglocal[whichbin][groupatom1][groupatom2] += 1;
	        }
	      }
      }
    }
  } 
  // Writting in the file
  if ( ((count > 0) && (count%WriteFileEvery == 0))  )
  {
    int sizeofGg = binsize * MAXGROUPS * MAXGROUPS ; // Everything should be passed
    // pass the address of the first element of a multidimensional array
    MPI_Reduce(&Gglocal[0][0][0] , &Gg[0][0][0] , sizeofGg , MPI_DOUBLE , MPI_SUM , 0 , world);
    if ( me == 0)
    {
      // Start RDF postprocessing (Only the first binsize data should be processed)
      for ( int i = 0; i<binsize; i++)
      {
        for ( int j = 1; j <= tmpnumgroup ; j++ )
        {
          for ( int k = j+1; k <= tmpnumgroup ; k++ )
          {
	          Gg[i][j][k] += Gg[i][k][j];
	          Gg[i][k][j] = 0;
	        }
        }
      }
      for ( int i = 0; i<binsize; i++)
      {
        for ( int j = 1; j <= tmpnumgroup ; j++)
	      {
	        for (int k = j ; k <=tmpnumgroup ; k++)
	        {
	          if (i>0)
	          {
	            PartSum[i][j][k] = PartSum[i-1][j][k] + Gg[i][j][k];
	          } else {
	            PartSum[i][j][k] = Gg[i][j][k];
	          }
	        }
	      }
      }
      // Open a file and write the header
      FilePtr = fopen(filename,"w");
      fprintf(FilePtr , "# Radius \t");
      for (int j = 1 ; j <= tmpnumgroup ; j++)
      {
        for (int k = j ; k <= tmpnumgroup ; k++)
        {
	        // PRINT g(r), g(r)*correction)
	        fprintf(FilePtr , "finiterdf__%s_%s\t",group->names[tmpgroup[j][0]],group->names[tmpgroup[k][0]]);
	        fprintf(FilePtr , "rdf__%s_%s\t",group->names[tmpgroup[j][0]],group->names[tmpgroup[k][0]]);
	      }
      }
      fprintf(FilePtr , "\n");

      // Write all the bins
      for ( int i = 0 ; i < binsize ; i++ )
      {
        r = (i+0.5)*Delta;
	      rin = 1.0*i*Delta;
	      rout = 1.0*(i+1.0)*Delta;
	      if (rin<hBox) {vin = (4.0*M_PI/3.0)* Cube_deltaboxsize * CUBE(1.0*i) ;}
	      else { vin = (-M_PI/12.0) * (3.0 - 36.0*Sqr_deltaboxsize*SQR(1.0*i) + 32.0*Cube_deltaboxsize*CUBE(1.0*i) ) ;}
	      if (rout<hBox) {vout = (4.0*M_PI/3.0)* Cube_deltaboxsize * CUBE(i+1.0) ;}
	      else { vout = (-M_PI/12.0) * (3.0 - 36.0*Sqr_deltaboxsize*SQR(i+1.0) + 32.0*Cube_deltaboxsize*CUBE(i+1.0) ) ;}
        sf = (vout - vin) * Ggt;
	      SphereVolFrac = vout;

        fprintf(FilePtr,"%g ",r);
        NumberOfParticles = 0; 
        for (int j = 1 ; j <= tmpnumgroup ; j++)
        {
          for (int k = j ; k <= tmpnumgroup ; k++)
	        {
	          if ( j != k )
	          {
	            pairs = 1.0 * sf * tmpgroup[j][2] * tmpgroup[k][2] ; 
	            gr_correction = tmpgroup[k][2] * (1.0 - SphereVolFrac) / (tmpgroup[k][2] - PartSum[i][j][k]/(tmpgroup[j][2]*Ggt));
	          } else {
	            pairs = 0.5 * sf * tmpgroup[j][2] * tmpgroup[k][2] ; 
	            gr_correction = tmpgroup[k][2] * (1.0 - SphereVolFrac) / (tmpgroup[k][2] - 2.0*PartSum[i][j][k]/(tmpgroup[j][2]*Ggt)-1.0);
	          }
	          // PRINT g(r),  g(r)*correction)
	          fprintf(FilePtr , "%g\t",Gg[i][j][k]/pairs);
	          fprintf(FilePtr , "%g\t",gr_correction*Gg[i][j][k]/pairs);
	        }
	        NumberOfParticles += tmpgroup[j][2] ;
        }
        fprintf(FilePtr , "\n");
      }
      fclose(FilePtr);
    }
  }
  
  count++;
  return 0.0;   // JUST TO RETURN SOMETHING
}

