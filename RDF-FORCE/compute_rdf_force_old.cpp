// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL), Jeff Greathouse (SNL)
------------------------------------------------------------------------- */

#include "compute_rdf_force.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "KIM/kim_units.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;


/* ---------------------------------------------------------------------- */

ComputeRDFForce::ComputeRDFForce(LAMMPS *lmp, int narg, char **arg) :
        Compute(lmp, narg, arg),
        rdfpair(nullptr), nrdfpair(nullptr), hist(nullptr), histall(nullptr) {
    if (narg < 4) error->all(FLERR, "Illegal compute rdf command");

    array_flag = 1;
    extarray = 0;

    nbin = utils::inumeric(FLERR, arg[3], false, lmp);
    if (nbin < 1) error->all(FLERR, "Illegal compute rdf command");

    // optional args
    // nargpair = # of pairwise args, starting at iarg = 4

    cutflag = 0;

    int iarg;
    for (iarg = 4; iarg < narg; iarg++)
        if (strcmp(arg[iarg], "cutoff") == 0) break;

    int nargpair = iarg - 4;

    while (iarg < narg) {
        if (strcmp(arg[iarg], "cutoff") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal compute rdf command");
            cutoff_user = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
            if (cutoff_user <= 0.0) cutflag = 0;
            else cutflag = 1;
            iarg += 2;
        } else error->all(FLERR, "Illegal compute rdf command");
    }


    // pairwise args

    if (nargpair == 0) npairs = 1;
    else {
        if (nargpair % 2) error->all(FLERR, "Illegal compute rdf command");
        npairs = nargpair / 2;
    }

    size_array_rows = nbin;
    size_array_cols = 1 + 2 * npairs;

    int ntypes = atom->ntypes;
    memory->create(rdfpair, npairs, ntypes + 1, ntypes + 1, "rdf:rdfpair");
    memory->create(nrdfpair, ntypes + 1, ntypes + 1, "rdf:nrdfpair");
    ilo = new int[npairs];
    ihi = new int[npairs];
    jlo = new int[npairs];
    jhi = new int[npairs];

    if (nargpair == 0) {
        ilo[0] = 1;
        ihi[0] = ntypes;
        jlo[0] = 1;
        jhi[0] = ntypes;
    } else {
        iarg = 4;
        for (int ipair = 0; ipair < npairs; ipair++) {
            utils::bounds(FLERR, arg[iarg], 1, atom->ntypes, ilo[ipair], ihi[ipair], error);
            utils::bounds(FLERR, arg[iarg + 1], 1, atom->ntypes, jlo[ipair], jhi[ipair], error);
            if (ilo[ipair] > ihi[ipair] || jlo[ipair] > jhi[ipair])
                error->all(FLERR, "Illegal compute rdf command");
            iarg += 2;
        }
    }

    int i, j;
    for (i = 1; i <= ntypes; i++)
        for (j = 1; j <= ntypes; j++)
            nrdfpair[i][j] = 0;

    int ihisto;
    for (int m = 0; m < npairs; m++)
        for (i = ilo[m]; i <= ihi[m]; i++)
            for (j = jlo[m]; j <= jhi[m]; j++) {
                ihisto = nrdfpair[i][j]++;
                rdfpair[ihisto][i][j] = m;
            }

    memory->create(hist, npairs, nbin, "rdf:hist");
    memory->create(histall, npairs, nbin, "rdf:histall");
    memory->create(array, nbin, 1 + 2 * npairs, "rdf:array");
    typecount = new int[ntypes + 1];
    icount = new int[npairs];
    jcount = new int[npairs];
    duplicates = new int[npairs];

    dynamic = 0;
    natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeRDFForce::~ComputeRDFForce() {
    memory->destroy(rdfpair);
    memory->destroy(nrdfpair);
    delete[] ilo;
    delete[] ihi;
    delete[] jlo;
    delete[] jhi;
    memory->destroy(hist);
    memory->destroy(histall);
    memory->destroy(array);
    delete[] typecount;
    delete[] icount;
    delete[] jcount;
    delete[] duplicates;
}

/* ---------------------------------------------------------------------- */

void ComputeRDFForce::init() {

    if (!force->pair && !cutflag)
        error->all(FLERR, "Compute rdf requires a pair style be defined "
                          "or cutoff specified");

    if (cutflag) {
        double skin = neighbor->skin;
        mycutneigh = cutoff_user + skin;

        double cutghost;            // as computed by Neighbor and Comm
        if (force->pair)
            cutghost = MAX(force->pair->cutforce + skin, comm->cutghostuser);
        else
            cutghost = comm->cutghostuser;

        if (mycutneigh > cutghost)
            error->all(FLERR, "Compute rdf cutoff exceeds ghost atom range - "
                              "use comm_modify cutoff command");
        if (force->pair && mycutneigh < force->pair->cutforce + skin && comm->me == 0)
            error->warning(FLERR, "Compute rdf cutoff less than neighbor cutoff - "
                                  "forcing a needless neighbor list build");

        delr = cutoff_user / nbin;
    } else delr = force->pair->cutforce / nbin;

    delrinv = 1.0 / delr;

    // set 1st column of output array to bin coords

    for (int i = 0; i < nbin; i++)
        array[i][0] = (i + 0.5) * delr;

    // initialize normalization, finite size correction, and changing atom counts

    natoms_old = atom->natoms;
    dynamic = group->dynamic[igroup];
    if (dynamic_user) dynamic = 1;
    init_norm();

    // need an occasional half neighbor list
    // if user specified, request a cutoff = cutoff_user + skin
    // skin is included b/c Neighbor uses this value similar
    //   to its cutneighmax = force cutoff + skin
    // also, this NeighList may be used by this compute for multiple steps
    //   (until next reneighbor), so it needs to contain atoms further
    //   than cutoff_user apart, just like a normal neighbor list does

    auto req = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
    if (cutflag) req->set_cutoff(mycutneigh);
}

/* ---------------------------------------------------------------------- */

void ComputeRDFForce::init_list(int /*id*/, NeighList *ptr) {
    list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRDFForce::init_norm() {
    int i, j, m;

    // count atoms of each type that are also in group

    const int nlocal = atom->nlocal;
    const int ntypes = atom->ntypes;
    const int *const mask = atom->mask;
    const int *const type = atom->type;

    for (i = 1; i <= ntypes; i++) typecount[i] = 0;
    for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) typecount[type[i]]++;

    // icount = # of I atoms participating in I,J pairs for each histogram
    // jcount = # of J atoms participating in I,J pairs for each histogram
    // duplicates = # of atoms in both groups I and J for each histogram

    for (m = 0; m < npairs; m++) {
        icount[m] = 0;
        for (i = ilo[m]; i <= ihi[m]; i++) icount[m] += typecount[i];
        jcount[m] = 0;
        for (i = jlo[m]; i <= jhi[m]; i++) jcount[m] += typecount[i];
        duplicates[m] = 0;
        for (i = ilo[m]; i <= ihi[m]; i++)
            for (j = jlo[m]; j <= jhi[m]; j++)
                if (i == j) duplicates[m] += typecount[i];
    }

    int *scratch = new int[npairs];
    MPI_Allreduce(icount, scratch, npairs, MPI_INT, MPI_SUM, world);
    for (i = 0; i < npairs; i++) icount[i] = scratch[i];
    MPI_Allreduce(jcount, scratch, npairs, MPI_INT, MPI_SUM, world);
    for (i = 0; i < npairs; i++) jcount[i] = scratch[i];
    MPI_Allreduce(duplicates, scratch, npairs, MPI_INT, MPI_SUM, world);
    for (i = 0; i < npairs; i++) duplicates[i] = scratch[i];
    delete[] scratch;
}

/* ---------------------------------------------------------------------- */

void ComputeRDFForce::compute_array() {
    int i, j, m, ii, jj, inum, jnum, itype, jtype, ipair, jpair, ibin, ihisto;
    double xtmp, ytmp, ztmp, delx, dely, delz, r, Ff, t;
    int *ilist, *jlist, *numneigh, **firstneigh;
    double factor_lj, factor_coul;

    int icompute = modify->find_compute("T");
    if (icompute < 0) error->all(FLERR, "Compute rdf/force requires compute T");
    Compute *tcompute = modify->compute[icompute];
    t = tcompute->compute_scalar();

    if (natoms_old != atom->natoms) {
        dynamic = 1;
        natoms_old = atom->natoms;
    }

    // if the number of atoms has changed or we have a dynamic group
    // or dynamic updates are requested (e.g. when changing atom types)
    // we need to recompute some normalization parameters

    if (dynamic) init_norm();

    invoked_array = update->ntimestep;

    // invoke half neighbor list (will copy or build if necessary)

    neighbor->build_one(list);

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // zero the histogram counts

    for (i = 0; i < npairs; i++)
        for (j = 0; j < nbin; j++)
            hist[i][j] = 0;

    // tally the RDF
    // both atom i and j must be in fix group
    // itype,jtype must have been specified by user
    // consider I,J as one interaction even if neighbor pair is stored on 2 procs
    // tally I,J pair each time I is central atom, and each time J is central

    double **x = atom->x;
    double **f = atom->f;

    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double *special_coul = force->special_coul;
    double *special_lj = force->special_lj;

    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            factor_lj = special_lj[sbmask(j)];
            factor_coul = special_coul[sbmask(j)];
            j &= NEIGHMASK;

            // if both weighting factors are 0, skip this pair
            // could be 0 and still be in neigh list for long-range Coulombics
            // want consistency with non-charged pairs which wouldn't be in list

            if (factor_lj == 0.0 && factor_coul == 0.0) continue;

            if (!(mask[j] & groupbit)) continue;
            jtype = type[j];
            ipair = nrdfpair[itype][jtype];
            jpair = nrdfpair[jtype][itype];
            if (!ipair && !jpair) continue;

            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            // No box checks implemented
            r = sqrt(delx * delx + dely * dely + delz * delz);
            if (j < nlocal) {
                Ff = ((f[i][0]-f[j][0]) * delx + (f[i][1]-f[j][1]) * dely + 
                (f[i][2]-f[j][2]) * delz) / (r * r * r);
            }
            else {
                Ff = ((f[i][0]) * delx + (f[i][1]) * dely + (f[i][2]) * delz) / (r * r * r);
            }
            ibin = static_cast<int> (r * delrinv + 0.5);
            if (ibin >= nbin) continue;

            for (ihisto = 0; ihisto < ipair; ihisto++) {
                m = rdfpair[ihisto][itype][jtype];
                hist[m][ibin] += Ff;
            }
        }
    }

    // sum histograms across procs

    MPI_Allreduce(hist[0], histall[0], npairs * nbin, MPI_DOUBLE, MPI_SUM, world);

    // convert counts to g(r) and coord(r) and copy into output array
    // vfrac = fraction of volume in shell m
    // npairs = number of pairs, corrected for duplicates
    // duplicates = pairs in which both atoms are the same

    double constant, vfrac, gr, ncoord, rlower, rupper, normfac;

    if (domain->dimension == 3) {
        constant = 4184. * domain->xprd * domain->yprd * domain->zprd /
                (8.3144626181 * 4 * MY_PI * t);

        for (m = 0; m < npairs; m++) {
            double sum = 0.0;
            ncoord = 0.0;
            for (ibin = 0; ibin < nbin; ibin++) {
                sum += histall[m][ibin];
                gr = sum * constant / (icount[m]*jcount[m]-duplicates[m]);
                ncoord += gr;
                array[ibin][1 + 2 * m] = gr;
                array[ibin][2 + 2 * m] = ncoord;
            }
        }

    } else {
        constant = MY_PI / (domain->xprd * domain->yprd);

        for (m = 0; m < npairs; m++) {
            ncoord = 0.0;
            normfac = (icount[m] > 0) ? static_cast<double>(jcount[m])
                                        - static_cast<double>(duplicates[m]) / icount[m] : 0.0;
            for (ibin = 0; ibin < nbin; ibin++) {
                rlower = ibin * delr;
                rupper = (ibin + 1) * delr;
                vfrac = constant * (rupper * rupper - rlower * rlower);
                if (vfrac * normfac != 0.0)
                    gr = histall[m][ibin] / (vfrac * normfac * icount[m]);
                else gr = 0.0;
                if (icount[m] != 0)
                    ncoord += gr * vfrac * normfac;
                array[ibin][1 + 2 * m] = gr;
                array[ibin][2 + 2 * m] = ncoord;
            }
        }
    }
}
