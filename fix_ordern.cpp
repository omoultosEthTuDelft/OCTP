/* ----------------------------------------------------------------------
   fix ordern is a child class of "fix", developed based on two classes
   of "fix ave/time" and "fix ave/correlate/long", provided by LAMMPS.
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

#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "fix_ordern.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

#include "group.h"
#include "domain.h"
#include "atom.h"


using namespace LAMMPS_NS;
using namespace FixConst;

enum{SCALAR,VECTOR};
enum{VISCOSITY,THERMCOND,DIFFUSIVITY};

#define INVOKED_VECTOR 2


/* ---------------------------------------------------------------------- */

FixOrderN::FixOrderN(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // At least 7 arguments are needed: [0-2], MODE, nevery, nwrite, value 
  if (narg < 7) error->all(FLERR,"Illegal fix ordern command");

  MPI_Comm_rank(world,&me);

  // Initial values
  startstep = 0;
  flag_Dxyz = 0;
  tnb = 10;
  tnbe = 10;
  fp1 = NULL;
  fp2 = NULL;
  title = NULL;
  format_user = NULL;
  format = (char *) " %g";
  dynamic_group_allow = 0;  // the groups should not be modified.

  // SPECIFYING THE MAIN ARGUMENTS
  // Define the type of transport property calculation
  if (strcmp(arg[3],"diffusivity") == 0) {
    mode = DIFFUSIVITY;
    char filen1[] = "selfdiffusivity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);
    char filen2[] = "onsagercoefficient.dat";
    filename2 = new char[strlen(filen2)+1];
    strcpy(filename2,filen2);
  } else if (strcmp(arg[3],"viscosity") == 0) {
    mode = VISCOSITY;
    char filen1[] = "viscosity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);   
    filename2 = NULL;
  } else if (strcmp(arg[3],"thermalconductivity") == 0) {
    mode = THERMCOND;
    char filen1[] = "thermconductivity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);
    filename2 = NULL;
  } else 
    error->all(FLERR,"Illegal fix ordern command with no transport property");
  // rate of sampling (end_of_step())
  nevery = force->inumeric(FLERR,arg[4]); 
  // rate of writing files
  nfreq = force->inumeric(FLERR,arg[5]);  
  global_freq = nevery;  

  // OBTAINING THE ID OF COMPUTE FOR THIS FIX
  // number of input values (it must be only one compute)
  nvalues = 0;
  int iarg = 6;
  while (iarg < narg) {
    if ((strncmp(arg[iarg],"c_",2) == 0)) {
      nvalues++;
      iarg++;
    } else break;
  }
  if (nvalues == 0) error->all(FLERR,"No values in fix ordern command");
  if (nvalues > 1) error->all(FLERR,"More than 1 value in fix ordern command");
  char *suffix = new char[strlen(arg[6])];
  strcpy(suffix,&arg[6][2]);
  char *ptr = strchr(suffix,'[');
  if (ptr) error->all(FLERR,"fix ordern requires all components of the compute");
  idcompute = new char[strlen(suffix) + 1];
  strcpy(idcompute,suffix);
  delete [] suffix;
  icompute = modify->find_compute(idcompute);     // id of the compute (int)
  nrows = modify->compute[icompute]->size_vector; // get the number of rows (int)
  Compute *compute = modify->compute[icompute];   // the whole compute class 

  // PARSING OPTIONAL ARGUMENTS
  iarg = 7;
  while (iarg < narg) {
    // add more file options for mode == visocisty/diffusion/thermcond
    if (strcmp(arg[iarg],"file") == 0) {
      if (mode == DIFFUSIVITY) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ordern command");
        delete [] filename1;
        filename1 = new char[strlen(arg[iarg+1])+1];
        strcpy(filename1,arg[iarg+1]);
        delete [] filename2;
        filename2 = new char[strlen(arg[iarg+2])+1];
        strcpy(filename2,arg[iarg+2]);
        iarg += 3;
      } else if (mode == THERMCOND || mode == VISCOSITY) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix ordern command");
        delete [] filename1;
        filename1 = new char[strlen(arg[iarg+1])+1];
        strcpy(filename1,arg[iarg+1]);
        iarg += 2;
      } else error->all(FLERR,"Illegal fix ordern command");
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ordern command");
      startstep = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"nb") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ordern command");
      tnb = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"nbe") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix ordern command");
      tnbe = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ordern command");
      delete [] format_user;
      int n = strlen(arg[iarg+1]) + 2;
      format_user = new char[n];
      sprintf(format_user," %s",arg[iarg+1]);
      format = format_user;
      iarg += 2;
    } else if (strcmp(arg[iarg],"title") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ordern command");
      delete [] title;
      int n = strlen(arg[iarg+1]) + 1;
      title = new char[n];
      strcpy(title,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Dxyz") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ordern command");
      if (strcmp(arg[iarg+1],"no") == 0) flag_Dxyz = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) flag_Dxyz = 1;
      else error->all(FLERR,"Illegal fix ordern command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ordern command");
  }

  // SETUP & ERROR CHECK
  // for fix inputs, check that fix frequency is acceptable
  // set variable_length if any compute is variable length
  if (nevery <= 0 || nfreq <= 0 || startstep < 0)
    error->all(FLERR,"Illegal fix ordern command: illegal number");
  if (nfreq % (2*nevery) )  // in case of integration that it cannot write
    error->all(FLERR,"Illegal fix ordern command: nevery is not a factor of nfreq");
  if (startstep % (nevery) )
    error->all(FLERR,"Illegal fix ordern command: nevery is not a factor of start");
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix ordern does not exist");
  if (modify->compute[icompute]->vector_flag == 0)
    error->all(FLERR,"Fix ordern compute does not calculate a vector");
  if (modify->compute[icompute]->size_vector_variable)
    error->all(FLERR,"The size of the vector should be kept fixed");

  // DEFINING THE PARAMETERS AND VARILABLES
  count = -1;   // the number of samp
  cnb = 1;
  boltz = force->boltz;
  nktv2p = force->nktv2p;
	  
  // Specific variables for each mode
  if (mode == DIFFUSIVITY)  {
    deltat = (double) (nevery)*(update->dt);
    vecsize = 3*atom->natoms;   // DOUBLE CHECK
    ngroup = 0;
  } else if (mode == VISCOSITY) {
    deltat = (double) (2.0*nevery)*(update->dt);
    vecsize = 7;
    sampsize = 8;
    sumP = 0;
    numP = 0;
    
  } else if (mode == THERMCOND) {
    deltat = (double) (2.0*nevery)*(update->dt);
    vecsize = 3;
    sampsize = 3;
  }

  // Order-n algorithm-specific parameters
  if ( (mode == VISCOSITY) || (mode == THERMCOND) )
  {
    memory->create(data,vecsize,"fix/ordern:data");
    memory->create(simpf0,vecsize,"fix/ordern:simpf0");
    memory->create(simpf1,vecsize,"fix/ordern:simpf1");
    memory->create(samp,tnb,tnbe,sampsize,"fix/ordern:samp");
  } else if ( mode == DIFFUSIVITY)
  {
    memory->create(PosC_ii,tnb,tnbe,MAXGROUPS,"fix/ordern:PosC_ii");
    memory->create(PosC_iix,tnb,tnbe,MAXGROUPS,"fix/ordern:PosC_iix");
    memory->create(PosC_iiy,tnb,tnbe,MAXGROUPS,"fix/ordern:PosC_iiy");
    memory->create(PosC_iiz,tnb,tnbe,MAXGROUPS,"fix/ordern:PosC_iiz");
    memory->create(PosC_ij,tnb,tnbe,MAXGROUPS,MAXGROUPS,"fix/ordern:PosC_ij");
    memory->create(PosC_ijx,tnb,tnbe,MAXGROUPS,MAXGROUPS,"fix/ordern:PosC_ijx");
    memory->create(PosC_ijy,tnb,tnbe,MAXGROUPS,MAXGROUPS,"fix/ordern:PosC_ijy");
    memory->create(PosC_ijz,tnb,tnbe,MAXGROUPS,MAXGROUPS,"fix/ordern:PosC_ijz");
    memory->create(PosCorrSum,tnb,tnbe,MAXGROUPS,3,"fix/ordern:PosCorrSum");
    memory->create(atomingroup,atom->natoms,2,"fix/ordern:atomingroup");
    memory->create(groupinfo,MAXGROUPS,2,"fix/ordern:groupinfo");
  }
  memory->create(recdata,nrows,"fix/ordern:recdata"); // data passed from compute
  memory->create(nsamp,tnb,tnbe,"fix/ordern:nsamp");
  memory->create(nbe,tnb,"fix/ordern:nbe");
  memory->create(oldint,tnb,tnbe,vecsize,"fix/ordern:oldint");
  memory->create(rint,vecsize,"fix/ordern:rint");
  for (int i = 0; i < tnb; i++)  nbe[i]=1;

  // this fix produces a global scalar 
  // intensive/extensive flags set by compute that produces value
  // This fix produces only a SCALAR value that I don't know yet (DOUBLE CHECK)
  scalar_flag = 1;
  vector_flag = 0;
  extscalar = 0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked
  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);

  // Opening new data files to output data  
  if (me == 0)  {
    if (mode == DIFFUSIVITY) {
      fp1 = fopen(filename1,"w");
      fp2 = fopen(filename2,"w");
      if (fp1 == NULL || fp2 == NULL) {
        error->all(FLERR,"Cannot open fix ordern files");
      }
    } else if (mode == THERMCOND || mode == VISCOSITY) {
      fp1 = fopen(filename1,"w");
      if (fp1 == NULL) {
        char str[128];
        snprintf(str,128,"Cannot open fix ordern file %s",filename1);
        error->one(FLERR,str);
      }
    }
  }

  // Writing the header lines to files
  if (fp1 && me == 0) {
    clearerr(fp1);
    if (title) fprintf(fp1,"%s\n",title);
    if (mode == DIFFUSIVITY)  {
      fprintf(fp1,"#NOTE: divide self-diffusivities by the ");
      fprintf(fp1,"number of molecules of species i (N_i).\n");
    } else if (mode == VISCOSITY) {
      fprintf(fp1,"#NOTE: divide shear viscosities by the temperature (T).\n");
      fprintf(fp1,"#Time\teta_xx\teta_yy\teta_zz\teta_xy\teta_xz\teta_yz\t");
      fprintf(fp1,"eta_off\teta_all\teta_bulk\n");
    } else if (mode == THERMCOND) {
      fprintf(fp1,"#NOTE: divide thermal conductivities by the temperature^2).\n");
      fprintf(fp1,"#Time\tlambda_x\tlambda_y\tlambda_z\tlambda_total\n");
    }
    if (ferror(fp1)) error->one(FLERR,"Error writing file header");
    filepos1 = ftell(fp1);
  }
  if (fp2 && me == 0)  {
    clearerr(fp2);
    if (title) fprintf(fp2,"%s\n",title);
    if (mode == DIFFUSIVITY) {
      fprintf(fp2,"#NOTE: divide Onsager coefficients ");
      fprintf(fp2,"by the total number of molecules (N).\n");
    }
  if (ferror(fp2)) error->one(FLERR,"Error writing file header");
  filepos2 = ftell(fp2);
  }
  delete [] title;
}

/* ---------------------------------------------------------------------- */

FixOrderN::~FixOrderN()
{
  delete [] format_user;
  if (fp1 && me == 0) fclose(fp1);
  if (fp2 && me == 0) fclose(fp2);
  if ( (mode == VISCOSITY) || (mode == THERMCOND) )
  {
    memory->destroy(data);
    memory->destroy(simpf0);
    memory->destroy(simpf1);
    memory->destroy(samp);
  } else if ( mode == DIFFUSIVITY)
  {
    memory->destroy(PosC_ii);
    memory->destroy(PosC_iix);
    memory->destroy(PosC_iiy);
    memory->destroy(PosC_iiz);
    memory->destroy(PosC_ij);
    memory->destroy(PosC_ijx);
    memory->destroy(PosC_ijy);
    memory->destroy(PosC_ijz);
    memory->destroy(PosCorrSum);
    memory->destroy(atomingroup);
    memory->destroy(groupinfo);
  }
  memory->destroy(recdata);
  memory->destroy(nsamp);
  memory->destroy(oldint);
  memory->destroy(rint);
  memory->destroy(nbe);
}

/* ----------------------------------------------------------------------
   defines when fix can be called (at the end of step)
------------------------------------------------------------------------- */

int FixOrderN::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   Initializing the whole fix in Modify::init()
------------------------------------------------------------------------- */

void FixOrderN::init()
{
  // set current indices for all computes,fixes,variables
  int icompute_new = modify->find_compute(idcompute);
  if (icompute < 0 || icompute_new != icompute)
    error->all(FLERR,"Compute ID for fix ordern does not exist");    
  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed
  if (nvalid < update->ntimestep) {
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixOrderN::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixOrderN::end_of_step()
{
  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner
  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix order");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;
  invoke_scalar(ntimestep);
}

void FixOrderN::invoke_scalar(bigint ntimestep)
{
  
  int i,j,k,l;
  double scalar = ntimestep;
  // update the next timestep to call end_of_step()
  modify->addstep_compute(ntimestep+nevery);

  // invoke compute vector if not previously invoked
  // get the data from compute_vector and store it in recdata
  // the size of recdata is nrows
  modify->clearstep_compute();
  Compute *compute = modify->compute[icompute];  // the whole compute class
  if (!(compute->invoked_flag & INVOKED_VECTOR)) {
    compute->compute_vector();
    compute->invoked_flag |= INVOKED_VECTOR;
  }
  double *cvector = compute->vector;
  for (i = 0; i < nrows; i++)
    recdata[i] = cvector[i];
    
  // update the next timestep to call end_of_step()
  nvalid += nevery;
  modify->addstep_compute(nvalid);
  if (count < 0)  
  {
    count = 0;
    icount = 0;
  } else {
    icount++;
    if ( (mode == DIFFUSIVITY) || (icount%2 == 0) )
      count++;
  }

  // From now, everything is computed only on the main core
  if (me != 0)  return;

  // Preliminary calculations for each transport property
  // Fill in the vector "data" accordingly
  if (mode == DIFFUSIVITY)  // DIFFUSION
  {
    tnatom = atom->natoms;
    if (count == 0)   // only run during the first time step
    { 
      natom = 0;
      // Finding corresponding groups to each atom at the first time
      for (sortID = 0; sortID < tnatom; sortID++)
      {
        ID = (int) (recdata[5*sortID+3]+0.1);
        atommask = (int) recdata[5*sortID+4]+0.1;
        int groupfound = 0;
        if (ngroup > 0 ) // First try to match with available groups
        {
      	  for (j = 1; j <= ngroup; j++)
          {
            if ( atommask  & groupinfo[j][1] )
            {
              atomingroup[ID][0] = natom;  // ID starts from 0
              atomingroup[ID][1] = j;  // groupID starts from 1
              natom++;
              groupfound = 1;
              break;
            }
          }
        }
        if (groupfound == 0) // If no match, try to find a new group
        {
          for (k = 1 ; k < MAXGROUPS; k++)
          {
            if ( atommask & group->bitmask[k] )
            {
              ngroup++;
              groupinfo[ngroup][0] = k;
              groupinfo[ngroup][1] = group->bitmask[k];
              atomingroup[ID][0] = natom;  // ID starts from 0
              atomingroup[ID][1] = ngroup;  // groupID starts from 1
              natom++;
              groupfound = 1;
              break; 
            }
          }
        }
        // if this atom doesn't belong to any group
        if (groupfound == 0)
          atomingroup[ID][0] = atomingroup[ID][1] = -1;
      }

      // redefine the arrays for storing the positions
      memory->destroy(oldint);
      memory->create(oldint,tnb,tnbe,3*natom,"fix/ordern:oldint");
      memory->destroy(rint);
      memory->create(rint,3*natom,"fix/ordern:rint");
      // initializing the arrays
      for(sortID = 0; sortID < tnatom ; sortID++)
      {
        ID = (int) (recdata[5*sortID+3]+0.1);
        if (atomingroup[ID][0] < 0)
          continue;
        atomID = atomingroup[ID][0];
        for( i = 0; i < tnb; i++)
        {
          oldint[i][tnbe-1][3*atomID] = recdata[5*sortID];
          oldint[i][tnbe-1][3*atomID+1] = recdata[5*sortID+1];
          oldint[i][tnbe-1][3*atomID+2] = recdata[5*sortID+2];
          for ( j = 0; j < tnbe; j++)
          {
            nsamp[i][j] = 0.0;
            for ( k = 0; k < MAXGROUPS; k++)
            {  
              PosC_ii[i][j][k] = 0.0;
              PosC_iix[i][j][k] = 0.0;
              PosC_iiy[i][j][k] = 0.0;
              PosC_iiz[i][j][k] = 0.0;
              for ( l = 0; l < MAXGROUPS ; l++)
              {
                PosC_ij[i][j][k][l] = 0.0;
                PosC_ijx[i][j][k][l] = 0.0;
                PosC_ijy[i][j][k][l] = 0.0;
                PosC_ijz[i][j][k][l] = 0.0;
              }
            }
          }
        } 
      }
    }
    if (count == 0) return; // nothing to do at this timestep

    for(sortID = 0; sortID < tnatom ; sortID++)
    {
      ID = (int) (recdata[5*sortID+3]+0.1);
      if (atomingroup[ID][0] < 0)
        continue;
      atomID = atomingroup[ID][0];
      rint[3*atomID] = recdata[5*sortID];
      rint[3*atomID+1] = recdata[5*sortID+1];
      rint[3*atomID+2] = recdata[5*sortID+2];
    }
  } else if (mode == VISCOSITY) // VISCOSITY
  {
    if (count == 0) {
      for ( i = 0; i<tnb; i++ )
        for ( j = 0; j < tnbe; j++ )  {
          nsamp[i][j] = 0.0;        
          for (int k = 0; k < vecsize; k++) {
            oldint[i][j][k] = 0.0;
            samp[i][j][k] = 0.0;
          }
          samp[i][j][vecsize] = 0.0;
        }
    }
    data[6] = (recdata[0]+recdata[1]+recdata[2])/3.0;
    for (i = 0; i < 3; i++) data[i] = (recdata[i] - data[6]);
    for (i = 3; i < 6; i++) data[i] = recdata[i];
    sumP += data[6];
    numP += 1.0;
  } else if (mode == THERMCOND)  // THERMAL CONDUCTIVITY
  {
    if (count == 0) {
      for ( i = 0; i<tnb; i++ )
        for ( j = 0; j < tnbe; j++ )  {
          nsamp[i][j] = 0.0;        
          for (int k = 0; k < vecsize; k++) {
            oldint[i][j][k] = 0.0;
            samp[i][j][k] = 0.0;
          }
        }
    }
    for (i = 0; i < 3; i++) data[i] = recdata[i];
  }

  // INTEGRATION according to Simpson's rule
  if (mode == VISCOSITY || mode == THERMCOND)
  {
    integrate();
    if ((icount % 2) == 1)  return; // return for odd timesteps
  }

  // ORDER-N ALGORITHM
  // loop over all blocks
  i = count/(pow(tnbe,cnb));
  while (i != 0)
  {
    cnb++;
    i /= tnbe;
  }
  for(i = 0; i < cnb; i++)
  {
    if ((count)%((int)pow(tnbe,i))==0)
	  {
	    cnbe = MIN(nbe[i],tnbe);
	    for (j=tnbe-1; j>=tnbe-cnbe; j--)
	    {
        nsamp[i][j] += 1.0;
        if ( (mode == VISCOSITY) || (mode == THERMCOND) )
	      {
	        for ( k = 0; k < vecsize; k++)
          {
		        dist = rint[k]-oldint[i][j][k];
            // Add the coefficient later
	          samp[i][j][k] += (dist*dist);	
		        if ( (mode == VISCOSITY) && (k == vecsize-1) )
		        {
		          samp[i][j][vecsize] += (dist);
		        }
	        } 
        } else if (mode == DIFFUSIVITY)
        {
          for( k = 1 ; k <= ngroup ; k++ )
          {
            PosCorrSum[i][j][k][0] = 0;
            PosCorrSum[i][j][k][1] = 0;
            PosCorrSum[i][j][k][2] = 0;
          }
          double distx, disty, distz;
          double distx2, disty2, distz2;
          for (ID = 0; ID < tnatom ; ID++)
          {
            if (atomingroup[ID][0] < 0)
              continue;
            atomID = atomingroup[ID][0];
            atomgroup = atomingroup[ID][1];
            distx = oldint[i][j][3*atomID] - rint[3*atomID];
            disty = oldint[i][j][3*atomID+1] - rint[3*atomID+1];
            distz = oldint[i][j][3*atomID+2] - rint[3*atomID+2];
            distx2 = distx*distx;
            disty2 = disty*disty;
            distz2 = distz*distz;
            // Add coefficients (2.0 and 6.0) later
            PosC_iix[i][j][atomgroup] += distx2;
            PosC_iiy[i][j][atomgroup] += disty2;
            PosC_iiz[i][j][atomgroup] += distz2;
            PosC_ii[i][j][atomgroup] += (distx2 + disty2 + distz2);
            PosCorrSum[i][j][atomgroup][0] += distx;
            PosCorrSum[i][j][atomgroup][1] += disty;
            PosCorrSum[i][j][atomgroup][2] += distz;
          }
          for ( k = 1 ; k <= ngroup ; k++)
          {
            for ( l = 1 ; l <= ngroup ; l++)
            {
              distx2 = PosCorrSum[i][j][k][0] * PosCorrSum[i][j][l][0];
              disty2 = PosCorrSum[i][j][k][1] * PosCorrSum[i][j][l][1];
              distz2 = PosCorrSum[i][j][k][2] * PosCorrSum[i][j][l][2];
              // Add coefficients (2.0 and 6.0) later
              PosC_ijx[i][j][k][l] += distx2;
              PosC_ijy[i][j][k][l] += disty2;
              PosC_ijz[i][j][k][l] += distz2;
              PosC_ij[i][j][k][l] += ( distx2 + disty2 + distz2 ); 
            }
          }
        }
	    }
      //increase the current blocklength
    	nbe[i]++;
      //shift to the left, set last index to the correlation value
      if ( (mode == VISCOSITY) || (mode == THERMCOND) )
      {
	      for (int k=0; k < vecsize; k++)
	      {
	        for (int j=1; j < tnbe; j++)
	          oldint[i][j-1][k] = oldint[i][j][k] ;
	        oldint[i][tnbe-1][k] = rint[k]; 
	      }
      } else if (mode == DIFFUSIVITY)
      {
        for (ID = 0; ID < tnatom ; ID++)
        {
          if (atomingroup[ID][0] < 0)
            continue;
          atomID = atomingroup[ID][0];
          for(j = 1 ; j < tnbe ; j++)
          {
            oldint[i][j-1][3*atomID] = oldint[i][j][3*atomID];
            oldint[i][j-1][3*atomID+1] = oldint[i][j][3*atomID+1];
            oldint[i][j-1][3*atomID+2] = oldint[i][j][3*atomID+2];
          }
          oldint[i][tnbe-1][3*atomID] = rint[3*atomID];
          oldint[i][tnbe-1][3*atomID+1] = rint[3*atomID+1];
          oldint[i][tnbe-1][3*atomID+2] = rint[3*atomID+2];
        }
      }
	  }
  }

  // Output results to files (fp1 and fp2) if time == nfreq
  if (ntimestep % nfreq)  return;
  if (mode == DIFFUSIVITY)
    write_diffusivity();
  else if (mode == VISCOSITY)
    write_viscosity();
  else if (mode == THERMCOND)
    write_thermcond();
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   lower bound is the smallest multiple of nevery larger than startstep
   used only in the initialization
------------------------------------------------------------------------- */
bigint FixOrderN::nextvalid()
{
  bigint nvalid = update->ntimestep;
  if (startstep > nvalid) nvalid = startstep;
  if (nvalid % nevery) nvalid = (nvalid/nevery)*nevery+nevery;
  return nvalid;
}

/*-------------------------------------------------------------------------
   Integrating Dynamical Variables According to the Simpson's Rule
------------------------------------------------------------------------- */
void FixOrderN::integrate()
{
  if (icount == 0) {
    for (int i = 0; i < vecsize; i++) 
    {
      rint[i] = 0.0;
      simpf0[i] = data[i];
    }  
  } else if ((icount % 2) == 1) {
    for (int i = 0; i < vecsize; i++) 
      simpf1[i] = data[i];
  } else {
    for (int i = 0; i < vecsize; i++) 
    {
      // delta = dt*(f0+4f1+f2)/6
      rint[i] += deltat*(simpf0[i]+4*simpf1[i]+data[i])/6.0;
      simpf0[i] = data[i];
      simpf1[i] = 0.0;
    }
  }
  return;
}

/*-------------------------------------------------------------------------
   Writing Order-n Results for Diffusivity into a File
------------------------------------------------------------------------- */
void FixOrderN::write_diffusivity()
{
  fseek(fp1,filepos1,SEEK_SET);
  fseek(fp2,filepos2,SEEK_SET);
  int i, j, k, l, xyz;
  // Writing the header
  fprintf(fp1,"#Time\t");
  fprintf(fp2,"#Time\t");
  for ( k = 1; k <= ngroup; k++)
  {
    if (flag_Dxyz)
    {
      fprintf(fp1,"Ds__%s_x\t",group->names[groupinfo[k][0]]);
      fprintf(fp1,"Ds__%s_y\t",group->names[groupinfo[k][0]]);
      fprintf(fp1,"Ds__%s_z\t",group->names[groupinfo[k][0]]);
    }
    fprintf(fp1,"Ds__%s\t",group->names[groupinfo[k][0]]);
    for ( l = k; l <= ngroup; l++)
    {
      if (flag_Dxyz)
      {
        fprintf(fp2,"C__%s_%s_x\t",
          group->names[groupinfo[k][0]],group->names[groupinfo[l][0]]);
        fprintf(fp2,"C__%s_%s_y\t",
          group->names[groupinfo[k][0]],group->names[groupinfo[l][0]]);
        fprintf(fp2,"C__%s_%s_z\t",
          group->names[groupinfo[k][0]],group->names[groupinfo[l][0]]);
      }
      fprintf(fp2,"C__%s_%s\t",
        group->names[groupinfo[k][0]],group->names[groupinfo[l][0]]);
    }
  }
  fprintf(fp1,"\n");
  fprintf(fp2,"\n");
  
  for( i = 0;i < MIN(tnb,cnb); i++ )
  {
    cnbe = MIN(nbe[i],tnbe);
	  if (i == MIN(tnb,cnb)-1)
	    cnbe = cnbe-1;
    for( j = 1; j <= cnbe; j++ )
    {
	    time = (double) ((1.0*j)*(deltat)*pow(tnbe,i));
      fprintf(fp1,format,time);
      fprintf(fp2,format,time);
      for( k = 1; k <= ngroup; k++ )
      {
        if (flag_Dxyz)
        {
          fprintf(fp1,format,PosC_iix[i][tnbe-j][k]/nsamp[i][tnbe-j]/2.0);
          fprintf(fp1,format,PosC_iiy[i][tnbe-j][k]/nsamp[i][tnbe-j]/2.0);
          fprintf(fp1,format,PosC_iiz[i][tnbe-j][k]/nsamp[i][tnbe-j]/2.0);
        }
        fprintf(fp1,format,PosC_ii[i][tnbe-j][k]/nsamp[i][tnbe-j]/6.0);
        for( l = k; l <= ngroup; l++ )
        {
          if (flag_Dxyz)
          {
            fprintf(fp2,format,PosC_ijx[i][tnbe-j][k][l]/nsamp[i][tnbe-j]/2.0);
            fprintf(fp2,format,PosC_ijy[i][tnbe-j][k][l]/nsamp[i][tnbe-j]/2.0);
            fprintf(fp2,format,PosC_ijz[i][tnbe-j][k][l]/nsamp[i][tnbe-j]/2.0);
          }
          fprintf(fp2,format,PosC_ij[i][tnbe-j][k][l]/nsamp[i][tnbe-j]/6.0); 
        }
      }
      fprintf(fp1,"\n");
      fprintf(fp2,"\n");
    }
  }
  fflush(fp1);
  fflush(fp2);
  // delete all unnecessary text from the output file
  long fileend1 = ftell(fp1);
  if (fileend1 > 0) ftruncate(fileno(fp1),fileend1);
  long fileend2 = ftell(fp2);
  if (fileend2 > 0) ftruncate(fileno(fp2),fileend2); 
}

/*-------------------------------------------------------------------------
   Writing Order-n Results for Viscosity into a File
------------------------------------------------------------------------- */
void FixOrderN::write_viscosity()
{
  fseek(fp1,filepos1,SEEK_SET);
  int i, j, k;
	double totalall;
	double totaloff;
  double stresscomp;
  double volume = (domain->xprd * domain->yprd * domain->zprd);
  double coef = (volume/2.0/boltz)*(1.0/nktv2p);
	avgP = sumP / numP;
	for (i=0; i<MIN(tnb,cnb); i++)
	{
    cnbe = MIN(nbe[i],tnbe);
	  if (i == MIN(tnb,cnb)-1)
	    cnbe = cnbe-1;
	  for (j = 1; j <= cnbe; j++)	// Just neglect the first data on the right (k=2)
	  {
	    time = (double) ((1.0*j)*(deltat)*pow(tnbe,i));
      fprintf(fp1,format,time);
	    // SHEAR VISCOSITY
	    totalall = 0.0;
	    totaloff = 0.0;
	    for (k = 0; k < 6; k++)
	    {
	      stresscomp = coef*samp[i][tnbe-j][k]/nsamp[i][tnbe-j];
	      if (k<3)
	      { 
          // implicit 4/3 contribution from diagonal
	        fprintf(fp1,format,stresscomp*3.0/4.0);
	        totalall += 1.0*1.0*stresscomp/10.0;	
	      }
	      else
	      {
	        fprintf(fp1,format,stresscomp);
	        totalall += 2.0*1.0*stresscomp/10.0;
	        totaloff += stresscomp/3.0;
	      }
	    }
	    fprintf(fp1,format,totaloff);
      fprintf(fp1,format,totalall);
	    // BULK VISCOSITY
	    double intp2 = samp[i][tnbe-j][vecsize-1]/nsamp[i][tnbe-j];
	    double intp = samp[i][tnbe-j][vecsize]/nsamp[i][tnbe-j];
	    double bulkvis = coef*(intp2 - 2.0*intp*(avgP*time) + (avgP*time)*(avgP*time));
	    fprintf(fp1,format,bulkvis);
      fprintf(fp1,"\n");
	  }
	}
  fflush(fp1);
  // delete all unnecessary text from the output file
  long fileend1 = ftell(fp1);
  if (fileend1 > 0) ftruncate(fileno(fp1),fileend1);
}

/*-------------------------------------------------------------------------
   Writing Order-n Results for Thermal Conductivity into a File
------------------------------------------------------------------------- */
void FixOrderN::write_thermcond()
{
  fseek(fp1,filepos1,SEEK_SET);
  int i, j, k;
	double totalcond;
  double fluxcomp;
  double volume = (domain->xprd * domain->yprd * domain->zprd);
  double coef = (1.0/volume/2.0/boltz);
	for (i=0; i<MIN(tnb,cnb); i++)
	{
    cnbe = MIN(nbe[i],tnbe);
	  if (i == MIN(tnb,cnb)-1)
	    cnbe = cnbe-1;
	  for (j = 1; j <= cnbe; j++)	// Just neglect the first data on the right (k=2)
	  {
	    time = (double) ((1.0*j)*(deltat)*pow(tnbe,i));
      fprintf(fp1,format,time);
	    totalcond = 0.0;
	    for (k = 0; k < 3; k++)
	    {
	      fluxcomp = coef*samp[i][tnbe-j][k]/nsamp[i][tnbe-j];
        fprintf(fp1,format,fluxcomp);
	      totalcond += fluxcomp/3.0;
	    }
	    fprintf(fp1,format,totalcond);
      fprintf(fp1,"\n");
	  }
	}
  fflush(fp1);
  // delete all unnecessary text from the output file
  long fileend1 = ftell(fp1);
  if (fileend1 > 0) ftruncate(fileno(fp1),fileend1);
}