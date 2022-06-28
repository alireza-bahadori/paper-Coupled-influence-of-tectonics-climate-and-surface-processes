# Coupled influence of tectonics, climate, and surface processes on landscape evolution in southwestern North America

### This repository belongs to the paper

### Bahadori et al. (2022) “Coupled influence of tectonics, climate, and surface processes on landscape evolution in southwestern North America”

### and contains the input files (*.txt, *.xml, *.h5, *.hdf5, *.xmf) to the UWGeodynamics code used to compute the model results in the paper.

### The simulations presented in the manuscript were run using UWGeodynamics version 2.10.1.

The original version of the code used to perform the simulations in this study can be accessed at: 
https://github.com/underworldcode/UWGeodynamics

A user guide including basic information about model setup for simulations with UWGeodynamics is available at:   https://uwgeodynamics.readthedocs.io/en/latest/UserGuide.html

UWGeodynamics Documentation can be accessed at: https://uwgeodynamics.readthedocs.io/_/downloads/en/latest/pdf/

We recommend visualization of the output files using the open-source software “Paraview”: www.paraview.org/

## Instructions for use

Detailed instructions for installation of UWGeodynamics can be found at:  
https://github.com/underworldcode/UWGeodynamics

This repository includes a Jupyter Notebook file “WUS_4D_simulation.ipynb” that can be run to reproduce the simulations presented in this work. 

The python file “WUS_4D_simulation.py” can be used to run the simulation inside a docker container for UWGeodynamics on a Linux system.

You need to copy the directory “Bahadori_et_al_Nature_Communications_2022” into your docker container for UWGeodynamcis and run the following command inside the directory:

python WUS_4D_simulation.py

The directories named “outbdls” and “Output_WUS_3D” include the simulation output files from 36-35.5 Ma for surface processes and geodynamic models, respectively. The python file named “WUS_4D_simulation.py” will restart the simulation from step #5 or from 35.5 Ma onward and will stop the simulations at 35.0 Ma. 

“WUS_4D_simulation.py” file uses input data (e.g., temperature at 35.5Ma, traction field at 35.5Ma, rainfall maps at 35.5Ma, sea level fluctuation data at 35.5Ma, and model parameter file for surface processes at 35.5Ma) provided in directories named:

sea_level
UWG_temperature 
UWG_traction 
maped_rain
rain_mapes

After your simulations from 35.5-35.0Ma is completed, you need to perform the following command:

cp -r Output_WUS_3D_results/. Output_WUS_3D

You can then run the simulation from step #10 onward or from 35.0-34.5 Ma by updating the name of the input files for temperature, traction, and surface processes parameters at time 35.0Ma in “WUS_4D_simulation.py”.

The singularity container provided here (uwgeodynamics_latest.sif) can be used to run the simulation on a cluster using the following command:

singularity exec uwgeodynamics_latest.sif python WUS_4D_simulation.py

### Here is an example of a successful initiation of the simulation at step #5 (or at 35.5Ma) with one CPU core using input files provided in this repository: 

![Screen Shot 2022-06-28 at 2 08 55 PM](https://user-images.githubusercontent.com/54119695/176306825-7faca797-7f6c-4d2d-a2fa-19015aa35803.png)
![Screen Shot 2022-06-28 at 2 09 53 PM](https://user-images.githubusercontent.com/54119695/176306915-3817853d-6c2f-4f50-a367-a846d9551cb9.png)
