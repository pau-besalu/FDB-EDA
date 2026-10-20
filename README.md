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

Once the eef.inp file is set up, we can proceed with the calculations by sending the eef.sh script into the queue system. Please, adapt the header of eef.sh for proper alignement with your queue system.
The eef.sh script will prepare a new folder in which all the calculations will take place.





