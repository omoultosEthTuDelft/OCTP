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

//enum{COMPUTE,FIX,VARIABLE};
//enum{ONE,RUNNING,WINDOW};
enum{SCALAR,VECTOR};
enum{VISCOSITY,THERMCOND,DIFFUSIVITY};

#define INVOKED_VECTOR 2


/* ---------------------------------------------------------------------- */

FixOrderN::FixOrderN(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
  /* DELETE ALL THESE INTIALIZATIONS,
  nvalues(0), which(NULL), argindex(NULL), value2index(NULL),
  offcol(NULL), varlen(NULL), ids(NULL),
  vector(NULL), vector_total(NULL), vector_list(NULL),
  column(NULL), array(NULL), array_total(NULL), array_list(NULL)*/
{
  // At least 7 arguments are needed: [0-2], MODE, nevery, nwrite, value 
  if (narg < 7) error->all(FLERR,"Illegal fix ordern command");

  MPI_Comm_rank(world,&me);

  // Initial values
  fp1 = NULL;
  fp2 = NULL;
  startstep = 0;
  //noff = 0;
  //offlist = NULL;
  format_user = NULL;
  format = (char *) "\t%g";
  title = NULL;
  
  // Define the type of transport property calculation
  if (strcmp(arg[3],"diffusivity") == 0) {
    mode = DIFFUSIVITY;
    char filen1[] = "selfdiffusivity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);
    //delete [] filen1;
    char filen2[] = "onsagercoefficients.dat";
    filename2 = new char[strlen(filen2)+1];
    strcpy(filename2,filen2);
    //delete [] filen2;
  } else if (strcmp(arg[3],"viscosity") == 0) {
    mode = VISCOSITY;
    char filen1[] = "shearviscosity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);   
    //delete [] filen1;
    char filen2[] = "bulkviscosity.dat";
    filename2 = new char[strlen(filen2)+1];
    strcpy(filename2,filen2);
    //delete [] filen2;
  } else if (strcmp(arg[3],"thermalconductivity") == 0) {
    mode = THERMCOND;
    char filen1[] = "thermalconductivity.dat";
    filename1 = new char[strlen(filen1)+1];
    strcpy(filename1,filen1);
    //delete [] filen1;
    filename2 = NULL;
  } else 
    error->all(FLERR,"Illegal fix ordern command with no transport property");

  nevery = force->inumeric(FLERR,arg[4]); // rate of sampling (end_of_step())
  nfreq = force->inumeric(FLERR,arg[5]);  // rate of writing files
  global_freq = nfreq;   //// DOUBLE CHECK THE MEANING

  dynamic_group_allow = 0;  // the groups should not be modified.

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

  // getting the ID of compute for the transport property
  char *suffix = new char[strlen(arg[6])];
  strcpy(suffix,&arg[6][2]);
  char *ptr = strchr(suffix,'[');
  if (ptr) error->all(FLERR,"All components of the compute are required for fix ordern");
  idcompute = new char[strlen(suffix) + 1];
  strcpy(idcompute,suffix);
  delete [] suffix;
  
  // Optional arguments
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
    }  else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
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
    } else error->all(FLERR,"Illegal fix ave/time command");
  }

  /*char str[128];
  snprintf(str,128,"Cannot open fix ordern file %sQQQ",idcompute);
  error->one(FLERR,str);*/





  /* START COMPLETE DELETE */

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  //int expand = 0;
  //char **earg;
  //nvalues = input->expand_args(nvalues,&arg[6],VECTOR,earg);

  //if (earg != &arg[6]) expand = 1;
  //arg = earg;

  // parse values

  //which = new int[nvalues];           // COMPUTE, VARIABLE, FIX
  //argindex = new int[nvalues];        // Index of each property
  //value2index = new int[nvalues];     // The conversion from ID of a property to its value
  //offcol = new int[nvalues];          // Irrelevant as no array
  //varlen = new int[nvalues];        // NOT NEEDED
  //ids = new char*[nvalues];

  /*for (int i = 0; i < nvalues; i++) {
    if (arg[i][0] == 'c') which[i] = COMPUTE;
    else if (arg[i][0] == 'f') which[i] = FIX;
    else if (arg[i][0] == 'v') which[i] = VARIABLE;

    int n = strlen(arg[i]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[i][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      if (suffix[strlen(suffix)-1] != ']')
        error->all(FLERR,"Illegal fix ave/time command");
      argindex[i] = atoi(ptr+1);
      *ptr = '\0';
    } else argindex[i] = 0;

    n = strlen(suffix) + 1;
    ids[i] = new char[n];
    strcpy(ids[i],suffix);
    delete [] suffix;
  }*/

  // set off columns now that nvalues is finalized

  /*for (int i = 0; i < nvalues; i++) offcol[i] = 0;
  for (int i = 0; i < noff; i++) {
    if (offlist[i] < 1 || offlist[i] > nvalues)
      error->all(FLERR,"Invalid fix ave/time off column");
    offcol[offlist[i]-1] = 1;
  }*/

  /* FINISH COMPLETE DELETE */


  // setup and error check
  // for fix inputs, check that fix frequency is acceptable
  // set variable_length if any compute is variable length
  if (nevery <= 0 || nfreq <= 0 || startstep < 0)
    error->all(FLERR,"Illegal fix ordern command");
  if (nfreq % nevery)
    error->all(FLERR,"Illegal fix ordern command: nevery should be a factor of nfreq");
  if (startstep % nevery)
    error->all(FLERR,"Illegal fix ordern command: nevery should be a factor of start");
  icompute = modify->find_compute(idcompute);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix ordern does not exist");
  if (modify->compute[icompute]->vector_flag == 0)
    error->all(FLERR,"Fix ordern compute does not calculate a vector");
  if (modify->compute[icompute]->size_vector_variable)
    error->all(FLERR,"The size of the vector should be kept fixed");

  nrows = modify->compute[icompute]->size_vector; // get the number of rows
  Compute *compute = modify->compute[icompute];  // the whole compute class
  if (mode == DIFFUSIVITY)
    numgroup = group->ngroup;            // The total # of available groups
    // DOUBLE CHECK:  find the ID of groups to name them for diffusion








  // NOTE: None of compute commands are variable length; thus, no LOCK is needed





  /* START COMPLETE DELETE */

  /* for (int i = 0; i < nvalues; i++) {
    varlen[i] = 0;

    if (which[i] == COMPUTE && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,"Fix ave/time compute does not calculate a vector");
      if (argindex[i] && modify->compute[icompute]->array_flag == 0)
        error->all(FLERR,"Fix ave/time compute does not calculate an array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_array_cols)
        error->all(FLERR,"Fix ave/time compute array is accessed out-of-range");
      if (argindex[i] == 0 && modify->compute[icompute]->size_vector_variable)
        varlen[i] = 1;
      if (argindex[i] && modify->compute[icompute]->size_array_rows_variable)
        varlen[i] = 1;
    } else if (which[i] == FIX && mode == VECTOR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,"Fix ave/time fix does not calculate a vector");
      if (argindex[i] && modify->fix[ifix]->array_flag == 0)
        error->all(FLERR,"Fix ave/time fix does not calculate an array");
      if (argindex[i] && modify->fix[ifix]->size_array_rows_variable)
        error->all(FLERR,"Fix ave/time fix array cannot be variable length");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_array_cols)
        error->all(FLERR,"Fix ave/time fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/time not computed at compatible time");

    } else if (which[i] == VARIABLE && mode == VECTOR) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/time does not exist");
      if (argindex[i] == 0 && input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/time variable is not vector-style variable");
      if (argindex[i])
        error->all(FLERR,"Fix ave/time mode vector variable cannot be indexed");
      varlen[i] = 1;
    }
  } */

  // all_variable_length = 1 if all values are variable length
  // any_variable_length = 1 if any values are variable length

  /*
  all_variable_length = 1;
  any_variable_length = 0;
  for (int i = 0; i < nvalues; i++) {
    if (varlen[i] == 0) all_variable_length = 0;
    if (varlen[i]) any_variable_length = 1;
  }
  */

  // if VECTOR mode, check that all columns are same length
  // nrows = # of rows in output array
  // if all columns are variable length, just set nrows = 1 for now

  /*column = NULL;
  if (all_variable_length == 0) 
    error->all(FLERR,"Fix ordern does not accept arrays");
  else nrows = 1;*/


  // enable locking of row count by this fix for computes of variable length
  // only if nrepeat > 1 or ave = RUNNING/WINDOW,
  //   so that locking spans multiple timesteps

  /*if (any_variable_length &&
      (nrepeat > 1 || ave == RUNNING || ave == WINDOW)) {
    for (int i = 0; i < nvalues; i++)
      if (varlen[i] && which[i] == COMPUTE) {
        int icompute = modify->find_compute(ids[i]);
        modify->compute[icompute]->lock_enable();
      }
    lockforever = 0;
  }*/

  // print file comment lines
  // for mode = VECTOR, cannot use arg to print
  // since array args may have been expanded to multiple vectors


// Opening new files to output data
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
      fprintf(fp1,"#NOTE: divide self-diffusivities by the number of molecules of species i (N_i).\n");
      fprintf(fp1,"#Time\t");
      for (int k = 0; k <= numgroup; k++) {
        fprintf(fp1,"Ds__%s\t",group->names[k]);      // DOUBLE CHECK  (tmpgroup)
      }
      fprintf(fp1,"\n");
    } else if (mode == VISCOSITY) {
      fprintf(fp1,"#NOTE: divide shear viscosities by the temperature (T).\n");
      fprintf(fp1,"#Time\teta_xx\teta_yy\teta_zz\teta_xy\teta_xz\teta_yz\teta_off\teta_total\n");
    } else if (mode == THERMCOND) {
      fprintf(fp1,"#NOTE: divide thermal conductivities by the temperature^2 (i.e.,T^2)).\n");
      fprintf(fp1,"#Time\tlambda_x\tlambda_y\tlambda_z\tlambda_total\n");
    }
    if (ferror(fp1)) error->one(FLERR,"Error writing file header");
    filepos1 = ftell(fp1);
    if (fp2)  {
      clearerr(fp2);
      if (title) fprintf(fp2,"%s\n",title);
      if (mode == DIFFUSIVITY) {
        fprintf(fp2,"#NOTE: divide Onsager coefficients by the total number of molecules (N).\n");
        fprintf(fp2,"#Time\t");
        for (int k = 0; k <= numgroup; k++) {
          for (int l = 0; l <= k; l++) {
            fprintf(fp2,"Lambda__%s_%s\t",group->names[k],group->names[l]); // DOUBLE CHECK (tmpgroup)
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
  }
  delete [] title;

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  /*if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }*/

  // allocate memory for the order-n algorithm
  // DOUBLE CHECK: ADD All arrays here (DOUBLE CHECK)

  memory->create(recdata,nrows,"ordern:recdata");

  if (mode == DIFFUSIVITY)  {

    int define_arrays = 0;

  } else if (mode == VISCOSITY) {

    count = -1;
    numpressure = 0;
    sumpressure = 0;
    avgpressure = 0;
    countint = 1;
    NumberOfBlock = 1;
    for (int i = 0; i < NUMBER_OF_BLOCKS; i++)
      BlockLength[i]=1;
    int define_arrays = 0;

  } else if (mode == THERMCOND) {

    int define_arrays = 0;

  }


  

  /*vector = vector_total = NULL;
  vector_list = NULL;
  array = array_total = NULL;
  array_list = NULL;

  if (mode == SCALAR) {
    vector = new double[nvalues];
    vector_total = new double[nvalues];
    if (ave == WINDOW)
      memory->create(vector_list,nwindow,nvalues,"ave/time:vector_list");
  } else allocate_arrays();*/



  // this fix produces either a global scalar or vector or array
  // SCALAR mode produces either a scalar or vector
  // VECTOR mode produces either a vector or array
  // intensive/extensive flags set by compute,fix,variable that produces value
    // This fix produces only a SCALAR value that I don't know yet (DOUBLE CHECK)
  scalar_flag = 1;
  extscalar = compute->extscalar;
  
  
  //extlist = NULL;
  /*if (mode == SCALAR) {
    if (nvalues == 1) {
      scalar_flag = 1;
      if (which[0] == COMPUTE) {
        Compute *compute = modify->compute[modify->find_compute(ids[0])];
        if (argindex[0] == 0) extscalar = compute->extscalar;
        else if (compute->extvector >= 0) extscalar = compute->extvector;
        else extscalar = compute->extlist[argindex[0]-1];
      } else if (which[0] == FIX) {
        Fix *fix = modify->fix[modify->find_fix(ids[0])];
        if (argindex[0] == 0) extscalar = fix->extscalar;
        else if (fix->extvector >= 0) extscalar = fix->extvector;
        else extscalar = fix->extlist[argindex[0]-1];
      } else if (which[0] == VARIABLE) {
        extscalar = 0;
      }

    } else {
      vector_flag = 1;
      size_vector = nrows = nvalues;
      extvector = -1;
      extlist = new int[nvalues];
      for (int i = 0; i < nvalues; i++) {
        if (which[i] == COMPUTE) {
          Compute *compute = modify->compute[modify->find_compute(ids[i])];
          if (argindex[i] == 0) extlist[i] = compute->extscalar;
          else if (compute->extvector >= 0) extlist[i] = compute->extvector;
          else extlist[i] = compute->extlist[argindex[i]-1];
        } else if (which[i] == FIX) {
          Fix *fix = modify->fix[modify->find_fix(ids[i])];
          if (argindex[i] == 0) extlist[i] = fix->extscalar;
          else if (fix->extvector >= 0) extlist[i] = fix->extvector;
          else extlist[i] = fix->extlist[argindex[i]-1];
        } else if (which[i] == VARIABLE) {
          extlist[i] = 0;
        }
      }
    }

  } else {
    if (nvalues == 1) {
      vector_flag = 1;
      size_vector = nrows;
      if (all_variable_length) size_vector_variable = 1;
      if (which[0] == COMPUTE) {
        Compute *compute = modify->compute[modify->find_compute(ids[0])];
        if (argindex[0] == 0) {
          extvector = compute->extvector;
          if (extvector == -1) {
            extlist = new int[nrows];
            for (int i = 0; i < nrows; i++) extlist[i] = compute->extlist[i];
          }
        } else extvector = compute->extarray;
      } else if (which[0] == FIX) {
        Fix *fix = modify->fix[modify->find_fix(ids[0])];
        if (argindex[0] == 0) {
          extvector = fix->extvector;
          if (extvector == -1) {
            extlist = new int[nrows];
            for (int i = 0; i < nrows; i++) extlist[i] = fix->extlist[i];
          }
        } else extvector = fix->extarray;
      } else if (which[0] == VARIABLE) {
        extlist = new int[nrows];
        for (int i = 0; i < nrows; i++) extlist[i] = 0;
      }

    } else {
      array_flag = 1;
      size_array_rows = nrows;
      size_array_cols = nvalues;
      if (all_variable_length) size_array_rows_variable = 1;
      int value;
      for (int i = 0; i < nvalues; i++) {
        if (which[i] == COMPUTE) {
          Compute *compute = modify->compute[modify->find_compute(ids[i])];
          if (argindex[i] == 0) value = compute->extvector;
          else value = compute->extarray;
        } else if (which[i] == FIX) {
          Fix *fix = modify->fix[modify->find_fix(ids[i])];
          if (argindex[i] == 0) value = fix->extvector;
          else value = fix->extarray;
        } else if (which[i] == VARIABLE) {
          value = 0;
        }
        if (value == -1)
          error->all(FLERR,"Fix ave/time cannot set output array "
                     "intensive/extensive from these inputs");
        if (i == 0) extarray = value;
        else if (value != extarray)
          error->all(FLERR,"Fix ave/time cannot set output array "
                     "intensive/extensive from these inputs");
      }
    }
  }*/

  // initializations
  // set vector_total to zero since it accumulates
  // array_total already zeroed in allocate_arrays

  /*irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  if (mode == SCALAR)
    for (int i = 0; i < nvalues; i++) vector_total[i] = 0.0;
  */

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked
  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);


}

/* ---------------------------------------------------------------------- */

FixOrderN::~FixOrderN()
{
  // decrement lock counter in compute chunk/atom, it if still exists

  /*if (any_variable_length &&
      (nrepeat > 1 || ave == RUNNING || ave == WINDOW)) {
    for (int i = 0; i < nvalues; i++)
      if (varlen[i]) {
        int icompute = modify->find_compute(ids[i]);
        if (icompute >= 0) {
          if (ave == RUNNING || ave == WINDOW)
            modify->compute[icompute]->unlock(this);
          modify->compute[icompute]->lock_disable();
        }
      }
  }*/

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

  int i,j;
  // DOUBLE CHECK: HOW TO DEFINE first timestep


  /*if (irepeat == 0) {    

    // DELET: THERE IS NO VARIABLE-SIZE COMPUTE
    if (any_variable_length) {
      modify->clearstep_compute();
      int nrows_new = column_length(1);
      modify->addstep_compute(ntimestep+nevery);
      modify->addstep_compute(ntimestep+nfreq);

      if (all_variable_length && nrows_new != nrows) {
        nrows = nrows_new;
        memory->destroy(column);
        memory->create(column,nrows,"ave/time:column");
        allocate_arrays();
      }

      bigint ntimestep = update->ntimestep;
      int lockforever_flag = 0;
      for (i = 0; i < nvalues; i++) {
        if (!varlen[i] || which[i] != COMPUTE) continue;
        if (nrepeat > 1 && ave == ONE) {
          Compute *compute = modify->compute[value2index[i]];
          compute->lock(this,ntimestep,ntimestep+(nrepeat-1)*nevery);
        } else if ((ave == RUNNING || ave == WINDOW) && !lockforever) {
          Compute *compute = modify->compute[value2index[i]];
          compute->lock(this,update->ntimestep,-1);
          lockforever_flag = 1;
        }
      }
      if (lockforever_flag) lockforever = 1;
    }

    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array[i][j] = 0.0;
  }*/

  
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


  if (ntimestep % nfreq)  return;

  // OUTPUT RESULTS TO THE FILES (fp1 and fp2) IF TIME == nfreq
  // DOUBLE CHECK
  if (fp1 && me == 0) {
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
  }
  // DO THE SAME WITH fp2
  // DOUBLE CHECK

}



/* ----------------------------------------------------------------------
   reallocate arrays for mode = VECTOR of size Nrows x Nvalues
------------------------------------------------------------------------- */

/*void FixOrderN::allocate_arrays()
{
  memory->destroy(array);
  memory->destroy(array_total);
  memory->create(array,nrows,nvalues,"ave/time:array");
  memory->create(array_total,nrows,nvalues,"ave/time:array_total");
  if (ave == WINDOW) {
    memory->destroy(array_list);
    memory->create(array_list,nwindow,nrows,nvalues,"ave/time:array_list");
  }

  // reinitialize regrown array_total since it accumulates

  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < nvalues; j++) array_total[i][j] = 0.0;
}*/

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
