# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 10:17:18 2022

@author: Jelle
"""

import numpy as np
import pandas as pd
import os
import scipy.constants as co
import uncertainties as unc
import matplotlib.pyplot as plt
plt.close('all')


class PP_OCTP:
    def __init__(self, folder, f_runs, groups, plotting=False):
        """
        This class is used to postproces octp output files and will average the
        results of multiple parallel runs if needed.

        Parameters
        ----------
        folder : string
            General folder location where the output files are located.
        f_runs : array of strings
            Intermediate folder level to separate multiple runs, like ['run_1',
            'run_2', .., 'run_n']. The results over these runs will be averaged
            and the standard deviations will be computed.
        groups : array of strings
            The array holds the names of the groups as defined for OCTP. The
            postprocessing is unsensitive of the order of these groups.
        store : string
            location where to store the output values of the postprocessing.
        plotting : boolean, optional
            State if plots of the ordern data should be shown. The default is
            False.

        Returns
        -------
        None.

        """
        self.f_folder = folder
        self.f_runs = f_runs
        self.f_file = [None]*len(self.f_runs)
        for i in range(len(self.f_runs)):
            self.f_file[i] = self.f_folder + self.f_runs[i]

        self.filenames(Default=True)
        self.changefit()
        self.groups = groups

        self.plotting = plotting
        # Preparing results dataframe
        self.results = pd.DataFrame()
        self.results['Legend'] = ['Mean', 'std', '# success']

        # Tool to keep track if mandatory function is ran already
        self.mandatory_ran = False

        # Tool to keep consistend graph colloring
        self.cmap = plt.get_cmap("tab10")

    def filenames(self, Default=False, density=False, volume=False,
                  total_E=False, temperature=False, pressure=False,
                  Diff_self=False, Diff_Onsag=False, viscosity=False,
                  T_conduc=False, rdf=False, log=False):
        """
        This function sets the filenames of the output files. The default
        filenames are set when the class is initiated, however single or all
        filennames can be changed lateron

        Parameters
        ----------
        Default : boolean, optional
            sets the default filenames. The default is False.
        density : string, optional
            Sets the name of the density file tracked during NPT initiation.
            The default is "density.dat".
        volume : string, optional
            Sets the name of the volume file, tracked during NPT initiation.
            The default is "volume.dat".
        total_E : string, optional
            Sets the name of the total energy file, tracked during the OCTP
            production. The default is "TotalE.dat".
        temperature : string, optional
            Sets the name of the temperature file, tracked during the OCTP
            production. The default is "temperature.dat".
        pressure : string, optional
            Sets the name of the pressure file as tracked during the OCTP
            production. The default is "pressure.dat".
        Diff_self : string, optional
            Sets the name of the OCTP self diffusion output file. The default
            is "diffself.dat".
        Diff_Onsag : string, optional
            Sets the name of the OCTP Onsager coefficient output file. The
            default is "diffonsag.dat".
        viscosity : string, optional
            Sets the name of the OCTP viscosity output file. The default is
            "viscosity.dat".
        T_conduc : string, optional
            Sets the name of the OCTP thermal conduction output file. The
            default is "tconductivity.dat".
        rdf : string, optional
            Sets the name of the OCTP rdf output file. The default is
            "rdf.dat".
        log : string, optional
            Sets the name of the lammps log output file. The default is
            "log.lammps".

        Returns
        -------
        None.

        """
        if Default is True:
            self.f_D = '/density.dat'  # Density in g/cm^3
            self.f_V = '/volume.dat'  # Volume in A^3 = 1e-30 m^3
            self.f_totE = '/TotalE.dat'  # Total energy in kcal
            self.f_T = '/temperature.dat'  # Temperature in K
            self.f_p = '/pressure.dat'  # Pressure in atm
            self.f_Ds = '/diffself.dat'
            self.f_DO = '/diffonsag.dat'
            self.f_visc = '/viscosity.dat'
            self.f_conT = '/tconductivity.dat'
            self.f_rdf = '/rdf.dat'
            self.f_log = '/log.lammps'

        if density is not False:
            self.f_D = '/' + density

        if volume is not False:
            self.f_V = '/' + volume

        if total_E is not False:
            self.f_totE = '/' + total_E

        if temperature is not False:
            self.f_T = '/' + temperature

        if pressure is not False:
            self.f_p = '/' + pressure

        if Diff_self is not False:
            self.f_Ds = '/' + Diff_self

        if Diff_Onsag is not False:
            self.f_DO = '/' + Diff_Onsag

        if viscosity is not False:
            self.f_visc = '/' + viscosity

        if T_conduc is not False:
            self.f_conT = '/' + T_conduc

        if rdf is not False:
            self.f_rdf = '/' + rdf

        if log is not False:
            self.f_log = '/' + log

    def changefit(self, margin=0.1, Minc=15, Mmax=60, er_max=0.05):
        """
        The changefit function adjusts the fitting parameters in the
        diff_calculator function (which computes the transport properties from
        mean square displacements tracked with the fix_ordern function of OCTP)

        Parameters
        ----------
        margin : float, optional
            The part cut of on the top and bottom side the MSD. This value
            avoids fitting in the balistic regime and in statistically relevant
            parts. The default is 0.1.
        Minc : integer, optional
            The minimum ammount of points that shall be included in a fit. The
            default is 15.
        Mmax : integer, optional
            The maximum ammount of points that shall be investigated for
            fitting. The default is 60.
        er_max : float, optional
            The maximum deviation of a slope of 1 in the log-log space. If this
            is set too tight, the fitting might return nan. The default is
            0.05.

        Returns
        -------
        None.

        """
        self.margin = margin
        self.Minc = Minc
        self.Mmax = Mmax
        self.er_max = er_max

    def check_succesfull(self, T_min):
        """
        This function removes the uncompleted (statistically untrustworthy)
        runs from the list of parallel runs. This funcion checks succesfull
        runtime with the temperature file, so make sure to set that file name
        correctly before calling this function.

        Parameters
        ----------
        T_min : float or integer
            Minimum succesfull runtime of production phase in ns.

        Returns
        -------
        None.

        """
        f_run = []
        for i in range(len(self.f_runs)):
            t = np.array(pd.read_table(self.f_folder + self.f_runs[i]+self.f_T,
                                       delimiter=' ', header=None,
                                       skiprows=2))[:, 0]
            if t[-1] > T_min*1e6:
                f_run.append(self.f_runs[i])
        print(len(f_run), 'statistically succesful runs for', self.f_folder)
        self.f_runs = f_run

        # Set the runfile list again with now only succesfull runs.
        self.f_file = [None]*len(self.f_runs)
        for i in range(len(self.f_runs)):
            self.f_file[i] = self.f_folder + self.f_runs[i]

    def mandatory_properties(self):
        """
        MANDATORY FUNCTION.This function will compute all mandatory properties.
        These properties will be used for OCTP outputs and include:

        1. The volume of the OCTP production run (average result of NPT run) in
        m^3.

        2. The box size of the OCTP production run (derived from volume) in m.

        3. The average temperature of the OCTP production run in K.

        4. The number of particles per group.

        Before this function is called, the filenames have to be set correctly.
        As this function is mandatory, it will automatically run before any
        non-Mandatory function is called for, exempting when the mandatory
        properties are already determined.

        Returns
        -------
        None.

        """
        # VOLUME SECTION
        # Here the volume and box size is computed and stored
        V = np.zeros(len(self.f_file))
        for i in range(len(self.f_file)):
            V[i] = np.array(pd.read_table(self.f_file[i]+self.f_V,
                                          delimiter=' ', header=None,
                                          skiprows=2))[:, 1]*1e-30
        L = np.power(V, 1/3)

        # Storing the results to be available per run in class
        self.V = V
        self.L = L

        # Storing the results in dataframe
        self.results['Volume/[m^3]'] = self.repacking_results(V)
        self.results['Box size/[m]'] = self.repacking_results(L)

        # TEMPERATURE SECTION
        # Here the temperature is computed and stored
        T_l = [None]*len(self.f_file)
        T = np.zeros(len(self.f_file))
        for i in range(len(self.f_file)):
            T_i = np.array(pd.read_table(self.f_file[i]+self.f_T,
                                         delimiter=' ', header=None,
                                         skiprows=2))[:, 1]
            T_l[i] = unc.ufloat(np.mean(T_i), np.std(T_i))
            T[i] = np.mean(T_i)

        # Storing the results to be available per run in class
        self.T = T

        # Storing the results in dataframe
        T_ave, T_sig = np.mean(T_l).n, np.mean(T_l).s
        T = [T_ave, T_sig, len(self.f_file)]
        self.results['Temperature/[K]'] = T

        # NUMBER OF PARTICLES PER GROUP SECTION
        # Here we will retrieve the number of particles per group from the log
        # file of lammps. This number is assumed to be constant for all
        # possible parallel runs.
        N_per_group = np.zeros(len(self.groups))
        for i in range(len(self.groups)):
            with open(self.f_file[0]+self.f_log) as fp:
                for l_no, line in enumerate(fp):
                    if 'group ' + self.groups[i] in line:
                        break
                words = fp.readline(l_no + 1).split()
                N_per_group[i] = int(words[0])

        # Store the result to be available in code
        self.N_per_group = N_per_group
        self.N_tot = np.sum(N_per_group)

        # Store the results in dataframe
        for i in range(len(self.groups)):
            self.results['Number of ' + self.groups[i]] = [N_per_group[i], 0,
                                                           len(self.f_file)]

        self.mandatory_ran = True

    def pressure(self):
        """
        OPTIONAL FUNCTION. The pressure is optimal and will not be reused. This
        function computes the average pressure and computes the mean and std
        of the system pressure during the OCTP sampling. This provides insight
        in the performance of the system. The pressure will be reported in Pa.

        Returns
        -------
        None.

        """
        # run mandatory properties run if that has not occured yet
        if self.mandatory_ran is False:
            self.mandatory_properties()

        p_l = [None]*len(self.f_file)
        for i in range(len(self.f_file)):
            # Pressure
            p_i = np.array(pd.read_table(self.f_file[i]+self.f_p,
                                         delimiter=' ', header=None,
                                         skiprows=2))[:, 1]
            p_l[i] = unc.ufloat(np.mean(p_i)*co.atm, np.std(p_i)*co.atm)

        # computing total mean and handling uncertainty correctly
        p_ave, p_sig = np.mean(p_l).n, np.mean(p_l).s
        self.results['pressure/[Pa]'] = [p_ave, p_sig, len(self.f_file)]

    def density(self):
        """
        OPTIONAL FUNCTION. The density is optional and will not be reused.
        This function computes the density per simulation, based on an average
        from a NPT ensemble. As a check, the density could be compared with the
        total mass and volume of the final ensemble. The unit of the density
        computed is kg/m^3.

        Returns
        -------
        None.

        """
        # run mandatory properties run if that has not occured yet
        if self.mandatory_ran is False:
            self.mandatory_properties()

        rho = [None]*len(self.f_file)
        for i in range(len(self.f_file)):
            rho_i = np.array(pd.read_table(self.f_file[i]+self.f_D,
                                           delimiter=' ', header=None,
                                           skiprows=2))[:, 1]
            rho[i] = np.mean(rho_i)*1000
        # Storing the result in the dataframe
        rho = [np.mean(rho), np.std(rho), len(self.f_file)]
        self.results['density/[kg/m^3]'] = rho

    def molarity(self, group):
        """
        OPTIONAL FUNCTION. This function computes the molarity of a specific
        group in a mixture, depending on group name. This is provided in
        mol/l_solution.

        Parameters
        ----------
        group : string
            The name of the group of which identifies solute.

        Returns
        -------
        None.

        """
        # run mandatory properties run if that has not occured yet
        if self.mandatory_ran is False:
            self.mandatory_properties()

        i = self.groups.index(group)  # Finding group index
        N = self.N_per_group[i]  # Getting number of particles in group
        M = (N/co.N_A)/(self.V*1000)

        # Storing the results in dataframe
        self.results['molarity/[mol/l]'] = [np.mean(M), np.std(M),
                                            len(self.f_file)]

    def molality(self, group, group_solvent, u_solvent):
        """
        OPTIONAL FUNTION. This function computes the molarity of a specific
        group in a mixture depending on the group and solvent name. The unit
        is mol/kg_solvent.

        Parameters
        ----------
        group : string
            The name of the group of identifies solute.
        group_solven : string
            The name of the group which indentifies the solvent.
        u_solvent : float
            The molecular weight of of the solvent in g/mol.

        Returns
        -------
        None.

        """
        # run mandatory properties run if that has not occured yet
        if self.mandatory_ran is False:
            self.mandatory_properties()

        i = self.groups.index(group)  # Finding group index
        N = self.N_per_group[i]  # Getting number of particles in group

        i = self.groups.index(group_solvent)  # Finding group index
        N_solv = self.N_per_group[i]  # Getting number of particles in group

        m = N/(N_solv*(u_solvent/1000))  # in mol/kg_solvent

        # Storing the result in dataframe
        self.results['molality/[mol/kg]'] = [m, 0, len(self.f_file)]

    def viscosity(self):
        """
        OPTIONAL FUNCTION: The viscosy is optional, however, it computes data
        mandatory for the finite size corrections for the self diffusion
        coefficients. The viscosity will be returned in Pa*s.

        Returns
        -------
        None.

        """
        # run mandatory properties run if that has not occured yet
        if self.mandatory_ran is False:
            self.mandatory_properties()

        visc = np.zeros(len(self.f_file))  # in Pa*s
        visc_bulk = np.zeros(len(self.f_file))  # in Pa*s

        for i in range(len(self.f_file)):
            # Getting data
            data = read(self.f_file[i]+self.f_visc, 'tmp')
            t = np.array(data['Time'])

            # Shear viscosity
            fact = 1e-3*1.01325e-7/self.T[i]
            msd_all = np.array(data['MSD_all'])
            D, t_fit, fit = self.diff_calculator(t, msd_all)
            visc[i] = D*fact
            # Plot if asked for
            if self.plotting is True:
                plt.figure('viscosity shear')
                plt.loglog(t, np.abs(msd_all)*fact, marker=".",
                           label=self.f_runs[i], color=self.cmap(i))
                plt.loglog(t_fit, fit*fact, ':', color=self.cmap(i))

            # Bulk viscosity
            msd_bulk = np.array(data['MSD_bulkvisc'])
            D, t_fit, fit = self.diff_calculator(t, msd_bulk)
            visc_bulk[i] = D*fact
            # Plot if asked for
            if self.plotting is True:
                plt.figure('viscosity bulk')
                plt.loglog(t, np.abs(msd_bulk)*fact, marker=".",
                           label=self.f_runs[i], color=self.cmap(i))
                plt.loglog(t_fit, fit*fact, ':', color=self.cmap(i))

        if self.plotting is True:
            plt.figure('viscosity shear')
            plt.title('viscosity shear')
            plt.grid('on')
            plt.ylabel('MSD pressure fluctuations/(Pa*s*fs)')
            plt.xlabel('time in fs')
            plt.legend()

            plt.figure('viscosity bulk')
            plt.title('viscosity bulk')
            plt.grid('on')
            plt.ylabel('MSD pressure fluctuations/(Pa*s*fs)')
            plt.xlabel('time in fs')
            plt.legend()
        # Storing the results to be availlable per run in class
        self.visc = visc

        # Storing results in dataframe. If nan is fit output, it is ignored
        visc = self.repacking_results(visc)
        self.results['viscosity shear/[Pa*s]'] = visc
        visc_bulk = self.repacking_results(visc_bulk)
        self.results['viscosity bulk/[Pa*s]'] = visc_bulk

    def thermal_conductivity(self):
        """
        OPTIONAL FUNCTION: The themal conductivity is optional and will be
        returned in units of W/(m*K).

        Returns
        -------
        None.

        """
        # run mandatory properties run if that has not occured yet
        if self.mandatory_ran is False:
            self.mandatory_properties()

        T_con = np.zeros(len(self.f_file))
        T_conb = np.zeros(len(self.f_file))
        for i in range(len(self.f_file)):
            data = read(self.f_file[i] + self.f_conT, 'tmp')
            t = np.array(data['Time'])
            msd_all = np.array(data['MSD_all'])
            fact = 6.9477e4/(self.T[i]**2)
            D, t_fit, fit = self.diff_calculator(t, msd_all)
            T_con[i] = D*fact

            D, t_fitb, fitb = self.diff_calculator(t, msd_all)
            T_conb[i] = D*fact

            # Plotting if asked for
            if self.plotting is True:
                plt.figure('thermal conductivity')
                plt.loglog(t, np.abs(msd_all)*fact, marker=".",
                           label=self.f_runs[i], color=self.cmap(i))
                plt.loglog(t_fit, fit*fact, ':', color=self.cmap(i))

        # Plotting if asked for
        if self.plotting is True:
            plt.figure('thermal conductivity')
            plt.title('thermal conductivity')
            plt.grid('on')
            plt.ylabel('MSD heat flux fluctuations/(W*fs/m/K)')
            plt.xlabel('time in fs')
            plt.legend()

        # Storing the results in dataframe
        T_con = self.repacking_results(T_con)
        self.results['Thermal conductivity/[W/(m*K)]'] = T_con

    def self_diffusivity(self, YH_correction=False, box_size_check=False):
        """
        OPTIONAL FUNCTION: The self-diffusion coefficients of the independent
        species are computed in m^2/s.

        Parameters
        ----------
        YH_correction : Boolean, optional
            Turns on the system size corrections by using the method proposed
            by Yeh and Hummer. For this correction, the viscosity has to be
            determined also. The default is False, True enables the correction.

        box_size_check : Boolean, optional
            Turns on the comparison between the box size and the length scales
            of the selected fits. If the lengthscales are smaller than the box
            size, a warning is raise, however the computations continue.
        Returns
        -------
        None.

        """
        # run mandatory properties run if that has not occured yet
        if self.mandatory_ran is False:
            self.mandatory_properties()

        D_YH = [0]*len(self.f_file)
        if YH_correction is True:
            # run viscosity if that has not occured yet
            try:
                visc = self.visc
            except AttributeError:
                self.viscosity()
                visc = self.visc
            xi = 2.837297
            D_YH = xi*co.k*self.T/(6*np.pi*visc*self.L)

        A_to_m = 1e-10  # ratio from angstrom to m
        fs_to_s = 1e-15  # ratio from femtosecond to second

        for j in range(len(self.groups)):
            N = self.N_per_group[j]
            D_s = np.zeros(len(self.f_file))

            for i in range(len(self.f_file)):
                data = read(self.f_file[i]+self.f_Ds, 'tmp')
                t = np.array(data['Time'])
                selfdif = np.array(data['MSD__'+self.groups[j]])   # Ansgrom^2

                if box_size_check is True:
                    length = self.L[i]/A_to_m
                else:
                    length = 0

                fact = A_to_m**2/(N*fs_to_s)
                D, t_fit, fit = self.diff_calculator(t, selfdif, m=length)
                D_s[i] = D*fact
                D_s[i] += D_YH[i]  # adds 0 if no correction is needed

                # Plotting if asked for
                if self.plotting is True:
                    plt.figure('D_self '+self.groups[j])
                    plt.loglog(t, selfdif*fact, marker=".",
                               label=self.f_runs[i], color=self.cmap(i))
                    plt.loglog(t_fit, fit*fact, ':', color=self.cmap(i))

            # Plotting if asked for
            if self.plotting is True:
                plt.figure('D_self ' + self.groups[j])
                plt.grid('on')
                plt.title('D_self ' + self.groups[j])
                plt.ylabel('msd ' + self.groups[j] + '/m^2')
                plt.xlabel('time in fs')
                plt.legend()

            # Storing self diffusion to dataframe
            words = 'Self Diffusion ' + self.groups[j] + '/[m^2/s]'
            D_self = self.repacking_results(D_s)
            self.results[words] = D_self

    def onsager_coeff(self, box_size_check=False):
        """
        OPTIONAL FUNCTION: The Onsager Coefficients of the independent species
        interactions are computed in m^2/s.

        Parameters
        ----------
        box_size_check : Boolean, optional
            Turns on the comparison between the box size and the length scales
            of the selected fits. If the lengthscales are smaller than the box
            size, a warning is raise, however the computations continue.

        Returns
        -------
        None.

        """
        # run mandatory properties run if that has not occured yet
        if self.mandatory_ran is False:
            self.mandatory_properties()

        N_tot = self.N_tot
        A_to_m = 1e-10  # ratio from angstrom to m
        fs_to_s = 1e-15  # ratio from femtosecond to second

        for i in range(len(self.groups)):
            for j in range(len(self.groups)-i):
                D_O = np.zeros(len(self.f_file))
                for k in range(len(self.f_file)):
                    # i indexes first species
                    # j indexes second species (i+j)
                    # k indexes parallel run number
                    data = read(self.f_file[k]+self.f_DO, 'tmp')
                    t = np.array(data['Time'])

                    try:
                        diff = np.array(data['MSD__' + self.groups[i]+'_' +
                                             self.groups[j+i]])  # Angstrom^2
                    except KeyError:
                        diff = np.array(data['MSD__' + self.groups[j+i]+'_' +
                                             self.groups[i]])  # Angstrom^2

                    if box_size_check is True:
                        length = self.L[k]/A_to_m
                    else:
                        length = 0

                    fact = A_to_m**2/(N_tot*fs_to_s)
                    D, t_fit, fit = self.diff_calculator(t, diff, m=length)
                    D_O[k] = D*fact

                    # Plotting if asked for
                    if self.plotting is True:
                        plt.figure('D_Onsag '+self.groups[i] + ' ' +
                                   self.groups[j+i])
                        plt.loglog(t, np.abs(diff)*fact,  marker=".",
                                   label=self.f_runs[k], color=self.cmap(k))
                        plt.loglog(t_fit, fit*fact, ':', color=self.cmap(k))

                # Plotting if asked for
                if self.plotting is True:
                    plt.figure('D_Onsag '+self.groups[i] + ' ' +
                               self.groups[j+i])
                    plt.title('D_Onsag '+self.groups[i] + ' ' +
                              self.groups[j+i])
                    plt.grid('on')
                    plt.ylabel('msd ' + self.groups[j] + ' ' +
                               self.groups[j+i] + '/m^2')
                    plt.xlabel('time in fs')
                    plt.legend()

                # Storing onsager coefficients to dataframe
                words = 'Onsager '+self.groups[i]+self.groups[j+i]+'/[m^2/s]'
                D_O = self.repacking_results(D_O)
                self.results[words] = D_O

    def store(self, location=False, name='postprocessed.csv'):
        """
        OPTIONAL FUNCTION: this function stores the dataframe in the
        intended location and under the itended name

        Parameters
        ----------
        location : string, optional
            Path to the file location where the dataframe should be stored. The
            default is self.f_folder (set automatically by using False).
        name : string, optional
            The intended filename. The default is 'postprocessed.csv'.

        Returns
        -------
        None.

        """
        if location is False:
            location = self.f_folder
        self.results.to_csv(location+name, index=False)

    def diff_calculator(self, t, MSD_in, m=False):
        """
        This function executes the fitting to a slope of 1 in the log-log
        domain. It tries all possible begin and start points for the fit and
        returns the fit whith a slope closest to 1. If no such fit can be
        found, NaN is returned for the transport property.

        Parameters
        ----------
        t : array of floats
            The timepoints returned by the OCTP algorithm.
        MSD_in : array of floats
            The mean-square-displacement values returned by the OCTP algorithm.
        m : float, optional
            Box size, scaled with the relevant number of particles. This
            provides the opportunity to be warned when the MSD for self-
            diffusivity and the Onsager coefficients are smaller than the
            box, which would be bad practice to use. The default is False.

        Returns
        -------
        D : float
            The value of the transport property which can be computed using the
            OCTP output data. This value is not yet corrected to the correct
            units.
        t_fit : array of floats
            The timepoints belonging to the fit. This can then be visualised
            when the individual runs are plotted.
        fit : array of floats
            The fit, to be visualised with the timepoints of the fit.

        """
        # Settings for the margins and fit method
        margin = self.margin  # cut away range at left and right side
        Minc = self.Minc  # minimum amount of points included in the fit
        Mmax = self.Mmax  # maximum ammount of points included in the fit
        er_max = self.er_max  # maximum allowed error

        t_log = np.log10(t)
        MSD_log_in = np.log10(np.abs(MSD_in))
        ibest = 'failed'
        jbest = 'failed'
        rebest = 'failed'
        mbest = 0

        for i in range(int(margin*len(t_log)), int((1-margin)*len(t_log))-Minc):
            for j in range(Minc, min(Mmax, int((1-margin)*len(t_log))-Minc-i)):
                if (t[i] != t[i+1]):
                    p, res, aa, aa1, aa3 = np.polyfit(t_log[i:i+j],
                                                      MSD_log_in[i:i+j], 1,
                                                      full=True)
                    mlog = p[0]
                    if (mlog > (1-er_max) and mlog < (1+er_max) and abs(mbest-1) > abs(mlog-1)):
                        mbest = mlog
                        jbest = j
                        ibest = i

        # Make sure to return nan (not included in np.nanmean() for averaging)
        if ibest == 'failed':
            D = np.nan
            t_fit = t[0]
            fit = MSD_in[0]

        else:
            D, b = np.polyfit(t[ibest:ibest+jbest],
                              MSD_in[ibest:ibest+jbest], 1)

            # Test box size to displacement comparison
            if np.abs(MSD_in[ibest+jbest]-MSD_in[ibest]) < m**2 and type(m) is not bool:
                print('MSD fit is smaller than simulation box',
                      MSD_in[ibest+jbest]-MSD_in[ibest], 'versus', m**2)

            t_fit = t[ibest:ibest+jbest]
            fit = D*t_fit + b
        return D, t_fit, fit

    def repacking_results(self, D_array):
        """
        This function repacks the outputs of multiple runs. If a list with
        fitted transport properties is inputted, it will output the mean value,
        its standard deviation and the number of non-nan type inputs. The
        output is packaged in such a way that it can be included in the results
        database immediatly.

        Parameters
        ----------
        D_array : array of floats
            A set of transport properties, for example the D_self of 5
            equivalent runs which have to be treated statistically.

        Returns
        -------
        list
            The returned list contains the mean, std and number of non-nan type
            items in the inputted D_array.

        """
        D_mean = np.nanmean(D_array)
        D_std = np.nanstd(D_array)
        n = D_array.size - np.count_nonzero(np.isnan(D_array))
        return [D_mean, D_std, n]


def read(datafile: str, export: str) -> pd.DataFrame:
    """
    This function reads the OCTP output files and reshapes it in a way that
    pandas can read the data effectively. The use of a temporary file is
    choosen as to not permanently eddid the OCTP output files for other uses.

    Parameters
    ----------
    datafile : str
        The name of the OCTP output file.
    export : str
        The name of the temporary file which will be read.

    Returns
    -------
    df : dataframe
        A pandas dataframe which makes it possible to acces the correct OCTP
        data needed for fitting.

    """
    inputFile = open(datafile, 'r')
    exportFile = open(export, 'w')
    for i, line in enumerate(inputFile):
        if line.lstrip()[:5] == '#Time':
            row = i
        new_line1 = line.replace('#', '').lstrip()
        new_line2 = new_line1.replace('\t', ' ').lstrip()
        exportFile.write(new_line2)

    inputFile.close()
    exportFile.close()

    df = pd.read_csv('tmp', sep=r'\s+', header=row)
    os.remove('tmp')
    return df
