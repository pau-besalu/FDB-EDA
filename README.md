# FDB-EDA
Script associated with the streamline generation of inputs and output analysis of FDB-EDA calculations.

This script needs to be adapted for each queue system, and version of AMS/ADF. In this particular case we mused
The main structure, however, can be addapted and utilized. For any questions regarding the implementation of this software please contact the main developer: paubesalu@gmail.com

The repository contains the following files:
  - File 1: eef.sh - the script that will be used to generate all inputs and fetch data from output
  - File 2: eef.inp - input required to run eef.sh
  - File 3: field.inp - data on which fields to use for FDB expansion, and 2D expansion
  - File 4: int.inp - ADF input semi-ready to automatically compute all the field-dependent EDA terms
  - File 5: strain.inp - ADF input semi-ready to automatically compute the FDB-AMS.
Also, there are three *.xyz files that are uploaded as an example for the user. In this case correspond to the Diels-Alder
reaction analysis. Nevertheless, they can be changed by the user.
  - File 6: cp.xyz - coordinates of the first reactant of the DA reaction.
  - File 7: ma.xyz - coordinates of the second reactant of the DA reaction.
  - File 8: TS_endo.xyz - coordinates and specifications of the supramolecule created by the two fragments.
  - File 9: TS_endo_cp.xyz - coordiantes of reactant "cp" in the supramolecule.
  - File 10: TS_endo_ma.xyz - coordinates of reactant "ma" in the supramolecule.
  - File 11: ams_oeef.f90 - once compilated, the derivatives calculator.

The workflow of the script is the following:
First, the eef.inp must be set up. In this file we have 6 lines which are the settings for the calculation.
  Line 1: what to do, among:
      (1) Perform only 1D calculations 
      (2) only 1D data process 
      (0) both (1)+(2)
      (3) 2D data generation and analysis 
      (4) 2D data analysis only
   
   Line 2: name of the supramolecule (in this example TS_endo)
   Line 3: name of reactant 1 (in this example cp)
   Line 4: name of reactant 2 (in this example ma)
   Line 5: number of electric fields used for the finite-field and romberg numerical calculations.
   Line 6: Which cartesian coordinate to analyze (only for 1D calculations)

  Then, the field.inp file must be set up. It contains 10 lines. The first 9 first lines must contain the electric field strengths in which the 1D analysis must take place (the output will be a file named ams_oeef.out that will contain the FDB-AMS/EDA analysis for the specified 9 field-strenghts). The 10th line corresponds to the grid of the 2D results plotting. In the example provided in this repository, the grid is 100 (10x10 2D grid). For publishable-quality material, grid 1000 is advised. Note that all this is computed *after* the ADF calculations analitically, and it takes virtually no time.

  As explained below, the FDB-AMS/EDA analysis can be run multiple times withouth recomputing the electronic structure with the eef.sh script. If the 9 EF strenghts are not enough, or you want to do some more testing, the file field.inp might be edited and the eef.sh run as many times as needed with no use of CPU time.

## General workflow:
The workflow of the script is summarized below:

1. **Configuration:** Set up `eef.inp` to specify the calculation mode, the supramolecule and reactant names, the number of electric fields used for the finite-field and Romberg calculations, and the Cartesian coordinate to analyze. Configure `field.inp` to define the electric-field strengths and, for 2D calculations, the grid size.

2. **Input generation:** Run `eef.sh` through the queueing system. The script creates a directory named after the supramolecule and automatically generates the input files required for the selected calculations.

3. **Electronic-structure calculations:** The generated inputs are used to perform the field-dependent interaction-energy and strain-energy calculations using ADF/AMS. The calculations are organized into the `int_f` and `strain` subdirectories.

4. **Output collection and processing:** Once the calculations are complete, the script extracts the relevant energy contributions from the output files, including the total energy, Pauli repulsion, orbital interaction, electrostatic interaction, and strain energy. The collected data are summarized in `summary.dat`.

5. **1D analysis:** For one-dimensional calculations, the collected energies are processed to obtain the field-dependent energy contributions and their derivatives. The corresponding numerical results, including the Romberg triangles, are made available for further analysis.

6. **2D analysis:** If 2D calculations are requested, additional calculations are performed to obtain the crossed-field derivatives. The `ams_oeef_2D.exe` post-processing program then uses `summary.dat` to generate the corresponding energy data files for each pair of Cartesian field directions.

7. **Data visualization:** The resulting files can be imported into Python or other data-analysis software for plotting and further interpretation. For 2D calculations, the grid size specified in `field.inp` determines the number of field combinations included in each output file.

**Note:** The queueing-system configuration, AMS/ADF commands, and output-parsing instructions may need to be adapted to the computational environment and software version being used.

Once the eef.inp file is set up, we can proceed with the calculations by sending the eef.sh script into the queue system. Please, adapt the header of eef.sh for proper alignement with your queue system.
The eef.sh script will create a new folder (with the name of the supramolecule) in which all the calculations will take place. This folder will contain two subfolders called int_f and strain, in which the field-dependent interaction energy and the strain calculations will take place. All the necessary outputs for the data analysis (after collecting data from the ADF outputs) will be stored at the base of the newly created folder. The Specifically, lines 1-4 and 81-84 of eef.sh are susceptible of change by the user, depending on the queue system (we assume slurm) and the ADF version. Note that we used ADF/2019.102, for newer versions blended with AMS package, the script might have to change in more places not just in the adverted lines of code.

Then, all inputs will be prepared automatically, according to the settings in eef.inp

The script is mainly divided into two main parts: 1D-fied type calculations and 2D-field type calculations.

At line 310 of code, the 1D module starts preparing a.ll files that are needed for the FDB-EDA run, and the calculations prepared previously are executed. Please, adapt this section to your requirements.

At line 420 of eef.sh the data recollection starts for the 1D calculations. Output files will be prepared summarizing the data fetched. This is the total bonding energies, the Pauli repulsion the orbital interaction, and electrostatics. Note that this section contains several "grep" commands. Adapt the necessary if you are using newer versions of AMS, and the data appears in a different format.

If 2D data is requested, additional crossed-term derivatives are required. In this case, after line 550 of the script the 2D-specific calculations are sent to queue, collected, and processed.

The computation of 2D crossed derivatives is done by ams_oeef_2D.exe software. The code for this is also provided in this repository as a *.f90 file. Please compile it in your cluster and route it properly before execution.

All the required energies will be summarized at a field labeled as summary.dat. This will contain either 1D or 2D data, and it will be located at the newly created folder in which the inputs and outputs are stored.

The ams_oeef_2D.exe software will collect the data from the summary.dat and will prepare several files labelled in the general form "A_F1F2.dat". A can be: Str, pau, oi, inten, toten; which are shortcuts for Strain energy, Pauli repulsion, orbital interaction, interaction energy and total energy. Then F1 and F2 are the two directions of the field in which this derivatives are computed. For example Str_XY.dat contains the data of the strain energy with electric-fields in X and Y directions.

All these ouptut files present three columns, the two first correspond to the electric-field strength in the two specified directions, and the third column corresponds to the energy term computed. For instance in Str_XY.dat first column represents the EF strenght in X direction, the second column represents the EF strenght in Y direction and the third column represents the strain energy for each combination of X and Y electric-field strenghts.

This data will be ready for eventual plotting using usual python plotting libraries like matplotlib. The lenght of the file will be the grid size for the plot. This can be controlled by the field.inp file. This field contains 9 values representing the EF strenghts to compute the 1D analysis, and the last (10th) row represents the grid size for the plotting. In this way, the grid (and the lenght of the files) can be controlled by the field.inp. In the example provided in this repository, the grid is 100 (10x10 2D grid). For publishable-quality material, grid 1000 is advised. Note that all this is computed *after* the ADF calculations analitically, and it takes virtually no time.

Also, raw romberg triangles for 1st and 2nd derivatives will be available at deriv1.dat and deriv2.dat. This way, the Romberg triangles from which the derivatives are computed can be manually inspected if necessary.

More data that might be useful (for 1D analysis) can be obtained at ams_oeef.out



