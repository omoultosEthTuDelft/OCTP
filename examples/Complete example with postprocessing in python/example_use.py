# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 14:34:46 2022

@author: Jelle
"""

import OCTP_postprocess_CLASS as post

"""
After calling the class, it is important to prepare some initial data. It is
advisable to have the lammps-OCTP output in an separate folder and for proper
treatment of the results, multiple similar simulations (with for example
different velocity seeds, but otherwise the same inputs). The can also run with
only a sigle output set, however it will report 0 for the standard deviation.
Additionally, the group names, as set in LAMMPS should, be prepared as well.
"""
folder = 'example_lammps_output/'  # Path to the main folder
f_runs = ['run1', 'run2', 'run3', 'run4', 'run5']  # All internal runs
groups = ['wat', 'met']  # Example group files

"""
Calling the initialisation of the class file
"""
octp = post.PP_OCTP(folder, f_runs, groups, plotting=True)

"""
If needed, some small adjustments can be made to the file names and the fitting
behaviour/limits can be adjusted as well. In this case, a special name is set
for the density filename and the fitting allows for slightly larger deviations
from 1 then default. Lastly, failed, short simulations, can be excluded from
the postprocessing by implementing a check to see how many succesfull runs are
completed. All simulations that fail this minimum runtime in nanoseconds are
then excluded from the further calculations for the results.

The change in the fitting properties can also be applied
between the functions calling for fitting. If for example, the viscosity
only fits correctly with broader margins, the
OCTP_postprocess_CLASS.changefit() can be called twice, once before the
OCTP_postprocess_CLASS.viscosity() function, setting broader margins and then
after the fitting setting the standard margins for the rest of the transport
properties.
"""
octp.filenames(density='density.dat')  # example of changing a file name
octp.changefit(er_max=0.08, Minc=12)  # alowes for 0.08 deviation of slope 1
octp.check_succesfull(8)  # only include runs of at least 8 ns

"""
The next section shows the real use of the OCTP postprocessing class.
Properties are then computed per type and some functions allow for aditional
inputs. For example, the Yeh-Hummer correction for system-size effects of self-
diffusion coefficients can be activated. As this is based on the viscosity of
the system, it will automatically call the compute viscosity function if this
has not occured yet.
"""
octp.pressure()  # Compute the average pressure during the sampling
octp.density()  # Compute the system density
octp.molality('met', 'wat', 18.01528)  # Compute molality of methanol in water
octp.molarity('met')  # Compute the molarity of methanol in water
octp.viscosity()  # Compute the viscosity
octp.thermal_conductivity()  # Compute the thermal conductivity
octp.self_diffusivity(YH_correction=True)
octp.onsager_coeff()  # Compute Onsager Coefficients

"""
Storing the results in a .csv format which can be opened with pandas and Excel
for further postprocessing and plotting of the results.
"""
octp.store(location='')  # Use location of this file instead of default folder.
