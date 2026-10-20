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
  - File 6: cp.xyz - coordinates of the first reactant of the DA reaction
  - File 7: ma.xyz - coordinates of the second reactant of the DA reaction
  - File 8: TS_endo.xyz - coordinates and specifications of the supramolecule created by the other two fragments.



