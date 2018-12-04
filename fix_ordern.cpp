/* ----------------------------------------------------------------------
   fix ordern is a child class of "fix", developed based on two classes
   of "fix ave/time" and "fix ave/correlate/long", provided in LAMMPS.
   This command is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

  //char str[128];
  //snprintf(str,128,"Cannot open fix ordern file %sQQQ",filename2);
  //error->one(FLERR,str);


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
   Contributing author: Pieter in 't Veld (SNL)
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
  tnb = 10;
  tnbe = 10;
  fp1 = NULL;
  fp2 = NULL;
  title = NULL;
  format_user = NULL;
  format = (char *) "\t%g";
  dynamic_group_allow = 0;  // the groups should not be modified.

  // SPECIFYING THE MAIN ARGUMENTS
  // Define the type of transport property calculation
  if (strcmp(arg[3],"diffusivity") == 0) {
    mode = DIFFUSIVITY;
    char filen1[] = "selfdiffusivity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);
    char filen2[] = "onsagercoefficients.dat";
    filename2 = new char[strlen(filen2)+1];
    strcpy(filename2,filen2);
  } else if (strcmp(arg[3],"viscosity") == 0) {
    mode = VISCOSITY;
    char filen1[] = "shearviscosity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);   
    char filen2[] = "bulkviscosity.dat";
    filename2 = new char[strlen(filen2)+1];
    strcpy(filename2,filen2);
  } else if (strcmp(arg[3],"thermalconductivity") == 0) {
    mode = THERMCOND;
    char filen1[] = "thermalconductivity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);
    filename2 = NULL;
  } else 
    error->all(FLERR,"Illegal fix ordern command with no transport property");
  // rate of sampling (end_of_step())
  nevery = force->inumeric(FLERR,arg[4]); 
  // rate of writing files
  nfreq = force->inumeric(FLERR,arg[5]);  
  global_freq = nfreq;  

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
        if (mode == DIFFUSIVITY || mode == VISCOSITY) {
          if (iarg+3 > narg) error->all(FLERR,"Illegal fix ordern command");
          delete [] filename1;
          filename1 = new char[strlen(arg[iarg+1])+1];
          strcpy(filename1,arg[iarg+1]);
          delete [] filename2;
          filename2 = new char[strlen(arg[iarg+2])+1];
          strcpy(filename2,arg[iarg+2]);
          iarg += 3;
        } else if (mode == THERMCOND) {
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
  // Order-n algorithm-specific parameters
  memory->create(recdata,nrows,"ordern:recdata"); // data passed from compute
  memory->create(samp,tnb,tnbe,8,"fix/ordern:samp");
  memory->create(nsamp,tnb,tnbe,8,"fix/ordern:nsamp");
  memory->create(oldint,tnb,tnbe,7,"fix/ordern:oldint");
  memory->create(nbe,tnb,"fix/ordern:nbe");

  count = -1;   // the number of samp
  cnb = 1;
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  volume = (domain->xprd * domain->yprd * domain->zprd);
	
  
  for (int i = 0; i < tnb; i++)  nbe[i]=1;
  
  
  // Specific variables for each mode
  if (mode == DIFFUSIVITY)  {
    deltat = (double) (nevery)*(update->dt);
    numgroup = group->ngroup;            // The total # of available groups

  } else if (mode == VISCOSITY) {
    deltat = (double) (2*nevery)*(update->dt);
    coef = (volume/2.0/boltz)*(1.0/nktv2p);
    
  } else if (mode == THERMCOND) {
    deltat = (double) (2*nevery)*(update->dt);
    

  }

  // this fix produces a global scalar 
  // intensive/extensive flags set by compute that produces value
  // This fix produces only a SCALAR value that I don't know yet (DOUBLE CHECK)
  scalar_flag = 1;
  extscalar = compute->extscalar;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked
  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);

// Opening new data files to output data
  if (me == 0)  {
    if (mode == DIFFUSIVITY || mode == VISCOSITY) {
      fp1 = fopen(filename1,"w");
      fp2 = fopen(filename2,"w");
      if (fp1 == NULL || fp2 == NULL) {
        error->all(FLERR,"Cannot open fix ordern files");
      }
    } else if (mode == THERMCOND) {
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
      fprintf(fp1,"#Time\t");
      for (int k = 0; k <= numgroup; k++) {
        fprintf(fp1,"Ds__%s\t",group->names[k]);      // DOUBLE CHECK  (tmpgroup)
      }
      fprintf(fp1,"\n");
    } else if (mode == VISCOSITY) {
      fprintf(fp1,"#NOTE: divide shear viscosities by the temperature (T).\n");
      fprintf(fp1,"#Time\teta_xx\teta_yy\teta_zz\t");
      fprintf(fp1,"eta_xy\teta_xz\teta_yz\teta_off\teta_total\n");
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
      fprintf(fp2,"#Time\t");
      for (int k = 0; k <= numgroup; k++) {
        for (int l = 0; l <= k; l++) {
          // DOUBLE CHECK (tmpgroup)
          fprintf(fp2,"Lambda__%s_%s\t",group->names[k],group->names[l]);
        }
      }
      fprintf(fp1,"\n");
    } else if (mode == VISCOSITY) {
      fprintf(fp2,"#NOTE: divide bulk viscosities by the temperature (T).\n");
      fprintf(fp2,"#Time\teta_b_xx\teta_b_yy\teta_b_zz\teta_b_total\n");
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
  // DOUBLE CHECK: DELETE All arrays here
  if (mode == DIFFUSIVITY)  {
    int define_arrays = 0;
  } else if (mode == VISCOSITY) {
    int define_arrays = 0;
  } else if (mode == THERMCOND) {
    int define_arrays = 0;
  }


  //delete [] which;
  //delete [] argindex;
  //delete [] value2index;
  //delete [] offcol;
  //delete [] varlen;
  //for (int i = 0; i < nvalues; i++) delete [] ids[i];
  //delete [] ids;

  //delete [] extlist;


  memory->destroy(recdata);
  memory->create(samp);
  memory->create(nsamp);
  memory->create(oldint);
  memory->create(nbe);


  //delete [] vector;
  //delete [] vector_total;
  //memory->destroy(array);
  //memory->destroy(array_total);
  //memory->destroy(array_list);
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

  // DOUBLE CHECK: Correct initialization of all variables
  if (mode == DIFFUSIVITY)  {
    
    
  } else if (mode == VISCOSITY) {
    sumP = 0;
    numP = 0;
  } else if (mode == THERMCOND) {


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
    error->all(FLERR,"Invalid timestep reset for fix ave/time");
  if (ntimestep != nvalid) return;
  
  // update the next timestep to call end_of_step()
  nvalid_last = nvalid;
  nvalid += nevery;
  modify->addstep_compute(nvalid);

  int i,j,k,l;
  // DOUBLE CHECK: HOW TO DEFINE first timestep
  if (count < 0)  
  {
    count = 0;
    icount = 0;
  } else 
  {
    icount++;
    if ((mode == VISCOSITY || mode == THERMCOND) && (icount%2 == 0))
      count++;
  }
  
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

  // From now-on, it should be only on the master core
  if (me != 0)  return;
    
  // Preliminary calculations for each transport property
  if (mode == DIFFUSIVITY)  // DIFFUSION
  {
    // NO NEED
    
  } else if (mode == VISCOSITY) // VISCOSITY
  {
    data[6] = (recdata[0]+recdata[1]+recdata[2])/3.0;
    for (i = 1; i < 4; i++) data[i] = (recdata[i] - data[6]);
    for (i = 4; i < 7; i++) data[i] = recdata[i];
    sumP += data[6];
    numP += 1.0;
  } else if (mode == THERMCOND)  // THERMAL CONDUCTIVITY
  {
    
  }
  // INTEGRATION according to Simpson's rule
  if (mode == VISCOSITY || mode == THERMCOND)
  {
    if (icount == 0) {
      for (i = 0; i < 7; i++) 
      {
        rint[i] = 0.0;
        simpf0[i] = data[i];
      }  
    } else if ((icount % 2) == 1) {
      for (i = 0; i < 7; i++) 
        simpf1[i] = data[i];
      return;   // Not continuing till next timestep
    } else {
      for (i = 0; i < 7; i++) 
      {
        rint[i] += (nevery/3.0)*(simpf0[i]+4*simpf1[i]+data[i]);
        simpf0[i] = data[i];
        simpf1[i] = 0.0;
      }
    }
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
        //if ((countint)%((j+1)*((int)pow(tnbe,i)))==0)
	      {
	        for (k = 0; k < 7; k++)
          {
	          // correction for dimensions in total, and adding kb and volume (WITHOUT TEMPERATURE)
		        dist = accint[k]-oldint[i][j][k];
	          //samp[i][j][k] += (dist*dist)*(1.0/inv_volume/2.0/boltz)*(1.0/nktv2p);
	          samp[i][j][k] += (dist*dist);	// ADD THE COEFFICIENT LATER
	          nsamp[i][j][k] += 1.0;
		        if (k == 6)
		        {
		          samp[i][j][7] += (dist);
	            nsamp[i][j][7] += 1.0;
		        }
	        }
	      }
	    }
    	nbe[i]++;
	    for (int k=0; k < 7; k++)
	    {
	      for (int j=1; j < tnbe; j++)
	        oldint[i][j-1][k] = oldint[i][j][k] ;
	      oldint[i][tnbe-1][k] = accint[k]; 
	    }
	  }
  }

  



  // OUTPUT RESULTS TO THE FILES (fp1 and fp2) IF TIME == nfreq
  if (ntimestep % nfreq)  return;
  fseek(fp1,filepos1,SEEK_SET);
	for (i=0; i<MIN(tnb,cnb); i++)
	{
	  if (i == MIN(tnb,cnb)-1)
	    cnbe = MIN(nbe[i],tnbe)-1;
	  else
	    cnbe = MIN(nbe[i],tnbe);
	  for (j = 1; j <= cnbe; j++)	// Just neglect the first data on the right (k=2)
	  {
	    time = ((1.0*j)*(deltat)*pow(tnbe,i));
      fprintf(fp1,format,time);
	    // SHEAR VISCOSITY
	    double totalstress = 0.0;
	    double totalshear = 0.0;
	    for (k = 0; k< 6; k++)
	    {
	      stresscomp = coef*samp[i][tnbe-j][k]/nsamp[i][tnbe-j][k];
	      fprintf(fp1,format,stresscomp);
	      if (k<3)
	      { 
          // implicit 4/3 contribution from diagonal
	        totalstress += 1.0*1.0*stresscomp/10.0;	
	      }
	      else
	      {
	        totalstress += 2.0*1.0*stresscomp/10.0;
	        totalshear += stresscomp/3.0;
	      }
	    }
	    fprintf(fp1,format,totalshear);
      fprintf(fp1,format,totalstress);
	    // BULK VISCOSITY
	    double intp2 = samp[i][tnbe-j][6]/nsamp[i][tnbe-j][6];
	    double intp = samp[i][tnbe-j][7]/nsamp[i][tnbe-j][7];
	    avgP = sumP / numP;
	    double bulkvis = coef*(intp2 - 2.0*intp*(avgP*time) + (avgP*time)*(avgP*time));
	    fprintf(fp1,format,bulkvis);
	  }
	}
  fflush(fp1);
  // delete all unnecessary text from the output file
  long fileend1 = ftell(fp1);
  if (fileend1 > 0) ftruncate(fileno(fp1),fileend1);

  // First output file (DOUBLE CHECK)
  /*if (fp1 && me == 0) {
    // getting the position to the end of header
    fseek(fp1,filepos1,SEEK_SET);
    // write the data 
    fprintf(fp1,BIGINT_FORMAT " %d\n",ntimestep,nrows);
    for (i = 0; i < nrows; i++) {
      //fprintf(fp1,"%d",i+1);
      fprintf(fp1,format,recdata[i]);  // user-defined format
      fprintf(fp1,"\n");
    }
    fflush(fp1);
    // delete all unnecessary text from the output file
    long fileend1 = ftell(fp1);
    if (fileend1 > 0) ftruncate(fileno(fp1),fileend1);
  }*/
  // Second output file (DOUBLE CHECK)

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
