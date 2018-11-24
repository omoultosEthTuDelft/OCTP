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
    compute_ordern_diffusion is based on the compute_msd command.
*/

#include <string.h>
#include "compute_position.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
//#include "modify.h"
//#include "fix_store.h"
#include "error.h"

using namespace LAMMPS_NS;

// SJ: start changing //
#include <iostream>
#include <string>
#include "memory.h"
// SJ: end changing //


/* ---------------------------------------------------------------------- */

ComputePosition::ComputePosition(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
//TB:  if (narg < 3) error->all(FLERR,"Illegal compute msd command");
  if (narg != 3) error->all(FLERR,"Illegal compute position command");

/*
//TB: what are these variables good for? They are not used!
//TB: maybe important for the created vector!!!
  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  create_attribute = 1;

  // optional args

//TB: no flags needed
  comflag = 0;
  avflag = 0;


  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute msd command");
      if (strcmp(arg[iarg+1],"no") == 0) comflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) comflag = 1;
      else error->all(FLERR,"Illegal compute msd command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"average") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute msd command");
      if (strcmp(arg[iarg+1],"no") == 0) avflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) avflag = 1;
      else error->all(FLERR,"Illegal compute msd command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute msd command");
  }
  
  // create a new fix STORE style for reference positions
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

//TB: we do not need to store the initial positions
  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "1";
  newarg[5] = (char *) "3";
  modify->add_fix(6,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;
*/
    double **x = atom->x;
    int *mask = atom->mask;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

/*
//TB: no initial positions no flags
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
      else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;

    // adjust for COM if requested

    if (comflag) {
      double cm[3];
      masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,cm);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          xoriginal[i][0] -= cm[0];
          xoriginal[i][1] -= cm[1];
          xoriginal[i][2] -= cm[2];
        }
    }

    // initialize counter for average positions if requested

    naverage = 0;
  }
*/
  // displacement vector
  
  vector = new double[4];  //TB: This has to be changed [4] -> natoms*5 or MAX_NUMBER_ATOMS*5

  // SJ: Start changing //
  MPI_Comm_rank(world, &me);            // Assigning the rank of a molecule for each core
  MPI_Comm_size(world, &nprocs);        // The number of processors 
  tmprecvcnts = new int[nprocs]; 
  recvcnts = new int[nprocs];
  displs = new int[nprocs];

  //TB: what is memory->create good for?
  //memory->create(BlockDATA,MAX_NUMBER_OF_BLOCKS,atom->natoms,MAX_NUMBER_OF_BLOCKELEMENTS,3,"compute/msd_norder:BlockDATA");
  // SJ: end changing  //

}

/* ---------------------------------------------------------------------- */

ComputePosition::~ComputePosition()
{
  // check nfix in case all fixes have already been deleted

  //if (modify->nfix) modify->delete_fix(id_fix);

  //delete [] id_fix;
  delete [] vector;
  
  // SJ: Start changing //
  delete [] tmprecvcnts;
  delete [] recvcnts;
  delete [] displs;
  //memory->destroy(BlockDATA);
  // SJ: end changing //
}

/* ---------------------------------------------------------------------- */

void ComputePosition::init()
{

  //TB: init() is required by Lammps, but we do not need to do anything?

  // set fix which stores reference atom coords

  //int ifix = modify->find_fix(id_fix);
  //if (ifix < 0) error->all(FLERR,"Could not find compute msd fix ID");
  //fix = (FixStore *) modify->fix[ifix];

  // nmsd = # of atoms in group

  //nmsd = group->count(igroup);
  //masstotal = group->mass(igroup);

  // SJ: start changing //
  //WriteFileEvery = 10000;  
  //tmpnumgroup = 0;
  count = 0;
  NumberOfBlocks = 1;
  //for (int ii = 0; ii< MAX_NUMBER_OF_BLOCKS; ii++)
  //{
  //    BlockLength[ii]=1;
  //}
  // SJ: end changing //

}

/* ---------------------------------------------------------------------- */

void ComputePosition::compute_vector()
{
  invoked_vector = update->ntimestep;

  // cm = current center of mass

  //double cm[3];
  //if (comflag) group->xcm(igroup,masstotal,cm);
  //else cm[0] = cm[1] = cm[2] = 0.0;

  // dx,dy,dz = displacement of atom from reference position
  // reference unwrapped position is stored by fix
  // relative to center of mass if comflag is set
  // for triclinic, need to unwrap current atom coord via h matrix

  //double **xoriginal = fix->astore;

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  //double dx,dy,dz;
  int xbox,ybox,zbox;

  //double msd[4];
  //msd[0] = msd[1] = msd[2] = msd[3] = 0.0;

  double xtmp, ytmp, ztmp;

  // update number of averages if requested

  /* 
  //TB: no flags
  double navfac;
  if (avflag) {
    naverage++;
    navfac = 1.0/(naverage+1);
  }
  */

  // SJ: start chaning //
  //TB: what is this? do we need any of this?
  int CurrentBlockLength; 
  int atomgroup , atomgroupbit;            // The ID of the group, its bit
  int numgroup = group->ngroup;            // The total # of available groups
  int natom = atom->natoms;                // total # of atoms in system, could be 0
  int *currentgroupbit = group->bitmask;    // The bitmask of a group
  //int ii = count/pow(MAX_NUMBER_OF_BLOCKELEMENTS,NumberOfBlocks);
  //while (ii != 0)
  //{
  // NumberOfBlocks++;
  //  ii /= MAX_NUMBER_OF_BLOCKELEMENTS;
  //}
  // SJ: end changing //

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {  
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
	xtmp = x[i][0] + xbox*xprd;
	ytmp = x[i][1] + ybox*yprd;
	ztmp = x[i][2] + zbox*zprd;
        // SJ: start chaning //
	realatom = (atom->tag[i]) - 1;
	sendbuff[5*i] = xtmp;
	sendbuff[5*i+1] = ytmp;
	sendbuff[5*i+2] = ztmp;
	sendbuff[5*i+3] = (double) realatom;
	sendbuff[5*i+4] = (double) mask[i];
        // SJ: end changing //

	// use running average position for reference if requested
/*
	if (avflag) {
	  xoriginal[i][0] = (xoriginal[i][0]*naverage + xtmp)*navfac;
	  xoriginal[i][1] = (xoriginal[i][1]*naverage + ytmp)*navfac;
	  xoriginal[i][2] = (xoriginal[i][2]*naverage + ztmp)*navfac;
	}

        dx = xtmp - xoriginal[i][0];
        dy = ytmp - xoriginal[i][1];
        dz = ztmp - xoriginal[i][2];
        msd[0] += dx*dx;
        msd[1] += dy*dy;
        msd[2] += dz*dz;
        msd[3] += dx*dx + dy*dy + dz*dz;
*/
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        xtmp = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
        ytmp = x[i][1] + h[1]*ybox + h[3]*zbox;
        ztmp = x[i][2] + h[2]*zbox;
        // SJ: start chaning //
	realatom = (atom->tag[i]) - 1;
	sendbuff[5*i] = xtmp;
	sendbuff[5*i+1] = ytmp;
	sendbuff[5*i+2] = ztmp;
	sendbuff[5*i+3] = (double) realatom;
	sendbuff[5*i+4] = (double) mask[i];
        // SJ: end changing //

	// use running average position for reference if requested

/*
	if (avflag) {
	  xoriginal[i][0] = (xoriginal[i][0]*naverage + xtmp)*navfac;
	  xoriginal[i][1] = (xoriginal[i][0]*naverage + xtmp)*navfac;
	  xoriginal[i][2] = (xoriginal[i][0]*naverage + xtmp)*navfac;
	}

        dx = xtmp - xoriginal[i][0];
        dy = ytmp - xoriginal[i][1];
        dz = ztmp - xoriginal[i][2];
        msd[0] += dx*dx;
        msd[1] += dy*dy;
        msd[2] += dz*dz;
        msd[3] += dx*dx + dy*dy + dz*dz;
*/
      }
  }
  
  // SJ: start changing //
  // Sending the position of all atoms from all cores (0, 1, ..., n) to all cores
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

  //TB: test start -> this has to be changed -> only to test if the program will compile
    vector[0] = 0.0;
    vector[1] = 1.0;
    vector[2] = 2.0;
    vector[3] = 3.0;
  //TB: test end
  

  // Definining the group number of each molecules on rank 0 
  //TB: maybe not needed? could also be done i the fix
  /*
  if (me == 0) 
  {
    for (int ii = 0; ii < natom; ii++)
    {
      int tmpgroups[MAX_NUMBER_OF_TYPES][2];
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
	      tmpfoundgroup = 1;
	      break;
	    }
	  }
	}
        if (tmpfoundgroup == 0) // If couldn't find any, tries to find the new group
	{
	  for (int jj = 1; jj < MAX_NUMBER_OF_TYPES; jj++)
	  {
	    if ( ((int) TmpPos[realatom][3]) & group->bitmask[jj])
	    {
	      tmpnumgroup++;
	      tmpfoundgroup = 1;
	      tmpgroup[tmpnumgroup][0] = jj;
	      tmpgroup[tmpnumgroup][1] = group->bitmask[jj];
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
  }
  */

  /*
  //TB: not need
  // Calculating the elements of Dself, Cii*, and Cij 
  if (me == 0) 
  {
    if (count>=1) 
    {
      if (count == 1)   
      {
        timeinterval = (double) (update->ntimestep - 0)*(update->dt);
        samplerate = (update->ntimestep);
      }
      //loop over all the blocks to test which blocks need sampling
      for(int CurrentBlock = 0; CurrentBlock < NumberOfBlocks; CurrentBlock++)
      {
        if ((count)%((int)pow(MAX_NUMBER_OF_BLOCKELEMENTS,CurrentBlock))==0)
        {
          //compute the current length of the block, limited to size ’MAX_NUMBER_OF_BLOCKELEMENTS’
          CurrentBlockLength = MIN(BlockLength[CurrentBlock],MAX_NUMBER_OF_BLOCKELEMENTS);
          //loop over the molecules in the system
          for(int i = MAX_NUMBER_OF_BLOCKELEMENTS-1;i >= MAX_NUMBER_OF_BLOCKELEMENTS - CurrentBlockLength;i--)
          { 
            count_samples[CurrentBlock][i] += 1.0;
	    for(int k1 = 1 ; k1 <= tmpnumgroup ; k1++ )
            {
               PosCorrSum[CurrentBlock][i][k1][0] = 0;
               PosCorrSum[CurrentBlock][i][k1][1] = 0;
               PosCorrSum[CurrentBlock][i][k1][2] = 0;
            }
            for (int realatom = 0; realatom < natom ; realatom++) // realatom = 0;realatom < natom;realatom++)
            {
              //realatom = (atom->tag[i]) - 1;
              groupatom1 = Groups[realatom][0];
              if (groupatom1 > 0)
              {
                PosC_ii[CurrentBlock][i][groupatom1] += 
          	                       ( SQR(BlockDATA[CurrentBlock][realatom][i][0] - TmpPos[realatom][0])
                                         + SQR(BlockDATA[CurrentBlock][realatom][i][1] - TmpPos[realatom][1])
                                         + SQR(BlockDATA[CurrentBlock][realatom][i][2] - TmpPos[realatom][2]) );
                PosCorrSum[CurrentBlock][i][groupatom1][0] += (BlockDATA[CurrentBlock][realatom][i][0] - TmpPos[realatom][0]);
                PosCorrSum[CurrentBlock][i][groupatom1][1] += (BlockDATA[CurrentBlock][realatom][i][1] - TmpPos[realatom][1]);
                PosCorrSum[CurrentBlock][i][groupatom1][2] += (BlockDATA[CurrentBlock][realatom][i][2] - TmpPos[realatom][2]);
              }
            }
            for ( int k1 = 1 ; k1 <= tmpnumgroup ; k1++)
            {
              for (int k2 = 1 ; k2 <= tmpnumgroup ; k2++)
              {
                PosC_ij[CurrentBlock][i][k1][k2] += 
                   + PosCorrSum[CurrentBlock][i][k1][0] * PosCorrSum[CurrentBlock][i][k2][0]
              	   + PosCorrSum[CurrentBlock][i][k1][1] * PosCorrSum[CurrentBlock][i][k2][1]
              	   + PosCorrSum[CurrentBlock][i][k1][2] * PosCorrSum[CurrentBlock][i][k2][2];
              }
            }
          }
          //increase the current blocklength
          BlockLength[CurrentBlock]++;
          //shift to the left, set last index to the correlation value
          for( realatom = 0 ; realatom < natom ; realatom++)
          { 
            for(int k = 1;k < MAX_NUMBER_OF_BLOCKELEMENTS;k++)
            {
              BlockDATA[CurrentBlock][realatom][k-1][0] = BlockDATA[CurrentBlock][realatom][k][0];
              BlockDATA[CurrentBlock][realatom][k-1][1] = BlockDATA[CurrentBlock][realatom][k][1];
              BlockDATA[CurrentBlock][realatom][k-1][2] = BlockDATA[CurrentBlock][realatom][k][2];
            }
            BlockDATA[CurrentBlock][realatom][MAX_NUMBER_OF_BLOCKELEMENTS-1][0] = TmpPos[realatom][0];
            BlockDATA[CurrentBlock][realatom][MAX_NUMBER_OF_BLOCKELEMENTS-1][1] = TmpPos[realatom][1];
            BlockDATA[CurrentBlock][realatom][MAX_NUMBER_OF_BLOCKELEMENTS-1][2] = TmpPos[realatom][2];
          }
        }
      } 
    }  else   {
      for(realatom = 0; realatom < natom ; realatom++)
      {
        for(int CurrentBlock = 0; CurrentBlock < MAX_NUMBER_OF_BLOCKS; CurrentBlock++)
        {
          BlockDATA[CurrentBlock][realatom][MAX_NUMBER_OF_BLOCKELEMENTS-1][0] = TmpPos[realatom][0];
          BlockDATA[CurrentBlock][realatom][MAX_NUMBER_OF_BLOCKELEMENTS-1][1] = TmpPos[realatom][1];
          BlockDATA[CurrentBlock][realatom][MAX_NUMBER_OF_BLOCKELEMENTS-1][2] = TmpPos[realatom][2];
          for (int k = 0; k < MAX_NUMBER_OF_BLOCKELEMENTS; k++)
          {
            count_samples[CurrentBlock][k] = 0.0;
            for (int k1 = 0; k1 < MAX_NUMBER_OF_TYPES; k1++)
            {  
              PosC_ii[CurrentBlock][k][k1] =0.0;
              for (int k2 = 0; k2 < MAX_NUMBER_OF_TYPES ; k2++)
              {
                PosC_ij[CurrentBlock][k][k1][k2] = 0.0;
              }
            }
          }
        } 
      }
    }
    // Writting in the file
    if ( ((count > 0) && (count%WriteFileEvery == 0))  || ( ((update->endstep)-(update->ntimestep)) < samplerate  )  )
    {
      FilePtrii = fopen("msd_norder_ii.dat","w");
      FilePtrij = fopen("msd_norder_ij.dat","w");
      // Writing the header
      fprintf(FilePtrii,"#  Time  \t");
      fprintf(FilePtrij,"#  Time  \t");
      for (int k = 1; k <= tmpnumgroup; k++)
      {
        //fprintf(FilePtrii,"Ds_%s\tCii*_%s\t",group->names[tmpgroup[k][0]],group->names[tmpgroup[k][0]]);
        fprintf(FilePtrii,"Ds__%s\t",group->names[tmpgroup[k][0]]);
        for (int kk = k; kk <= tmpnumgroup; kk++)
        {
          fprintf(FilePtrij,"C__%s_%s\t",group->names[tmpgroup[k][0]],group->names[tmpgroup[kk][0]]);
          //fprintf(FilePtrij,"C__%d_%d\t",tmpgroup[k][0],tmpgroup[kk][0]);
        }
      }
      fprintf(FilePtrii,"\n");
      fprintf(FilePtrij,"\n");
      
      for(int CurrentBlock = 0;CurrentBlock < MIN(MAX_NUMBER_OF_BLOCKS,NumberOfBlocks); CurrentBlock++)
      {
        CurrentBlockLength = MIN(BlockLength[CurrentBlock],MAX_NUMBER_OF_BLOCKELEMENTS);
        for(int j = 1; j <= CurrentBlockLength; j++)
        {
          fprintf(FilePtrii,"%-9g\t",(double)(j*(timeinterval)*pow(MAX_NUMBER_OF_BLOCKELEMENTS,CurrentBlock)));
          fprintf(FilePtrij,"%-9g\t",(double)(j*(timeinterval)*pow(MAX_NUMBER_OF_BLOCKELEMENTS,CurrentBlock)));
          msdcounting = count_samples[CurrentBlock][MAX_NUMBER_OF_BLOCKELEMENTS-j];	
          for(int k = 1; k <= tmpnumgroup; k++)
          {
            sum_dself = PosC_ii[CurrentBlock][MAX_NUMBER_OF_BLOCKELEMENTS-j][k]/(msdcounting);
            fprintf(FilePtrii,"%-10g\t",sum_dself);
            for(int kk = k; kk <= tmpnumgroup; kk++)
            {
              sum_cij = PosC_ij[CurrentBlock][MAX_NUMBER_OF_BLOCKELEMENTS-j][k][kk]/(msdcounting) ; 
              fprintf(FilePtrij,"%-10g\t",sum_cij); 
            }
          }
          fprintf(FilePtrii,"\n");
          fprintf(FilePtrij,"\n");
        }
      }
      fclose(FilePtrii);
      fclose(FilePtrij);
    }
  }
  count++;
  // SJ: end changing //

  MPI_Allreduce(msd,vector,4,MPI_DOUBLE,MPI_SUM,world);
  if (nmsd) {
    vector[0] /= nmsd;
    vector[1] /= nmsd;
    vector[2] /= nmsd;
    vector[3] /= nmsd;
  }
*/

}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */
/*
TB: not needed if we want to determine the current position

void ComputePosition::set_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}
*/
