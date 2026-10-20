#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-00:02:00

##
##      This script generates the input and sends the calculations for the obtention
##      of the data for the Taylor expansions of the field-dependent bonding energy
##
##      Version 2.0 --> New feature: 2D map data generation
##
##      Required files:
##          - supra.xyz       contains the geometry of the whole system
##          - supra_sub1.xyz  contains the geometry of the subsystem 1 in the supramolecule
##          - supra_sub2.xyz  contains the geometry of the subsystem 2 in the supramolecule
##          - strain.inp      contains a sample of strain calculation
##          - int.inp         contains a sample of interaction energy w/electric-field calculation
##
##      Coded by Pau Besalú-Sala, PhD candidate, Amsterdam, 2021.
##
##

 ################################################################
 ##                                                            ##
 ##                   Initialization module                    ##
 ##                                                            ##
 ################################################################


 #Reading of input parameters from file eef.inp
 declare -a input_array
 k=0
 while IFS= read -r line; do
   input_array[k]=$line
   k=$((k+1))
 done < eef.inp

 whatdo=${input_array[0]}       # Perform only 1D calculations (1) only data process (2) or both (0). 2D data generation and analysis (3) 2D data analysis only (4)
 name=${input_array[1]}         # Name of the supra system
 sub1=${input_array[2]}         # Name of the subsystem1
 sub2=${input_array[3]}         # Name of the subsystem2
 num_fields=${input_array[4]}   # Number of EFs applied
 if [ $whatdo -ne 3 ]; then
    dirF=${input_array[5]}      # Direction of the EF (x, y or z) - in working mode 3 this is not used.
 fi

 ref="no" #Comes from old versions. Whether to compute the references (yes) or read from xyz (no)

 #Starting information
 echo " "
 echo "     Job information"
 echo "     ---------------"
 echo " "
 echo "     Job name:                        " $name
 echo "         Subsystem 1:                 " $sub1
 echo "         Subsystem 2:                 " $sub2
 echo " "
 if [ $whatdo -ne 3 ]; then
 echo "     Direction of the EF is:          " $dirF
 else
 echo "     2D data generation requested. All directions for the EF computed. "
 fi
 echo "     Num. fields for FF-RR procedure: " $num_fields
 echo " "
 echo "     I/O directory:                   " $PWD
 echo " "
 if [ $ref == "yes" ]; then
    echo "     WARNING! THE REQUEST OF OPTIMIZATION DOES NOT SUPPORT YET "
    echo "     THE SP AND POLAR CALCULATIONS FOR THE OPTIMIZED GEOMETRY  "
 fi

 if [ $whatdo -ne 2 -a $whatdo -ne 4 ]; then # We will compute something. Need the directory empty.
   if [ -d $PWD/${name}_run ]; then
     echo "JOB FOLDER ALREADY EXISTS. EXITING"
     exit 1
   fi
 fi

 # Are all the required files there?
 if [ $whatdo -ne 2 -a $whatdo -ne 4 ]; then # We will compute something. Need the files
   #Charging modules
   cd $SLURM_SUBMIT_DIR
   module load 2021
   module load ADF/2019.102-intelmpi

   echo "     Required files:              "
   for i in $name.xyz ${name}_${sub1}.xyz ${name}_${sub2}.xyz strain.inp int.inp ; do
    if [ -e $i ]; then
      echo "        " $i "exists"
    else
      echo "        " $i " is not found in " $SLURM_SUBMIT_DIR " EXITING"
      echo " "
      exit 1
    fi
   done
   if [ $ref == "yes" ]; then #we dont need to optimize the subsystems, thus the coordinates must be present
    for i in ${sub1}.xyz ${sub2}.xyz ; do
      if [ -e $i ]; then
        echo "        " $i "exists"
      else
        echo "        " $i " is not found in " $SLURM_SUBMIT_DIR " EXITING"
        echo " "
        exit 1
      fi
    done
   fi

 echo " "
   echo "Generating directory tree"
   jobdir=${name}_run
   mkdir $jobdir
   strain=${SLURM_SUBMIT_DIR}/$jobdir/strain
   int=${SLURM_SUBMIT_DIR}/$jobdir/int_f$dirF
   mkdir $strain
   mkdir $int

   echo "Generating inputs "
   cp strain.inp $strain
   cp int.inp $int
   echo " "
 else # we only need to define some paths for the data process
   jobdir=${name}_run
   strain=$PWD/$jobdir/strain
   int=$PWD/$jobdir/int_f$dirF
 fi

 ################## END OF INITIALIZATION MODULE ################


 ################################################################
 ##                                                            ##
 ##                   1D data generation module                ##
 ##                                                            ##
 ################################################################

 if [ $whatdo -eq 0 ] || [ $whatdo -eq 1 ]; then # We must prepare inputs for 1D calculations
 ################### INPUTS FOR STRAIN CALCULATIONS ###################
 ###################         1D calculation         ###################

 # Prepare strain calculation
 cd $strain

 # Transform general strain.inp input into the proper input
 for i in ${name}_$sub1 ${name}_$sub2; do
   cp strain.inp ${i}_sp.inp
   cp strain.inp ${i}_polar.inp
   sed -i "s,filename,${i}_sp,g" ${i}_sp.inp
   sed -i "s,filename,${i}_polar,g" ${i}_polar.inp
   sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_sp.inp
   sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_polar.inp
   if [ $ref == "yes" ]; then #we need to generate the optimization inputs
      cp strain.inp ${i}_opt.inp
      sed -i "s,filename,${i}_opt,g" ${i}_opt.inp
      sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_opt.inp
   fi
 done
 if [ $ref == "no" ]; then #the sp and polar calculations for the isolated reactants are prepared
   for i in $sub1 $sub2 ; do
     cp strain.inp ${i}_sp.inp
     cp strain.inp ${i}_polar.inp
     sed -i "s,filename,${i}_sp,g" ${i}_sp.inp
     sed -i "s,filename,${i}_polar,g" ${i}_polar.inp
     sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_sp.inp
     sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_polar.inp
   done
 fi

 # If BETA wants to be computed, this could lead to an inconsistency
 # and code must be checked
 #First, save the line where the Response keyword is found
 #Then, delate polar calculation request for the single point calculations
 p=`grep -in "AOResponse" ${i}_polar.inp | sed "s,:, ,g" | awk '{print $1}' `
 sed -i "${p}d" *_sp.inp # supress such line and the following two
 sed -i "${p}d" *_sp.inp

 if [ $ref == "yes" ]; then
   #Same than before for the optimizations
   sed -i "${p}d" *_opt.inp # supress such line and the following two
   sed -i "${p}d" *_opt.inp

   ## Now, we request that the *opt.inp files have the optimization keyword + freq calculation ENDINPUT is used.
   sed -i '/^ENDINPUT/i Geometry'  *opt.inp
   sed -i '/^Geometry/a Iterations 299' *opt.inp
   sed -i '/^Iterations/a Optim Delocalized' *opt.inp
   sed -i '/^Optim Delocalized/a Converge E=1e-5 grad=1e-5' *opt.inp
   sed -i '/^Converge/a End' *opt.inp
   sed -i '/^ENDINPUT/i AnalyticalFreq' *opt.inp
   sed -i '/^AnalyticalFreq/a End' *opt.inp
   sed -i '/^ENDINPUT/i ScanFreq -1000 200' *opt.inp
   sed -i '/^ENDINPUT/i \ ' *opt.inp

   sed -i 's,Iterations 299,    Iterations 299,g' *opt.inp
   sed -i 's,Optim Delocalized,    Optim Delocalized,g' *opt.inp
   sed -i 's,Converge E=1e-5 grad=1e-5,    Converge E=1e-5 grad=1e-5,g' *opt.inp
 fi

 rm strain.inp
 cd ../ #Exiting strain directory

 #################### INPUTS FOR INTERACTION ENERGY ########################
 ####################        Fragment inputs        ########################

 # Prepare interaction calculation with field
 cd $int
 if [ $dirF == "x" ]; then
  sed -i "s,EField xx yy zz,EField xxchange 0.0 0.0,g" int.inp
 elif [ $dirF == "y" ]; then
  sed -i "s,EField xx yy zz,EField 0.0 xxchange 0.0,g" int.inp
 else  # If "x" or "y" not specified, use "z" by default
  sed -i "s,EField xx yy zz,EField 0.0 0.0 xxchange,g" int.inp
 fi

 #Zero-field calculation for fragments
 cp int.inp int_f${dirF}0.inp
 sed -i "/EField/d" int_f${dirF}0.inp #Remove EF keyword
 #Prepare filename and coordinates
 cp int_f${dirF}0.inp ${name}_${sub1}_f${dirF}0.inp
 cp int_f${dirF}0.inp ${name}_${sub2}_f${dirF}0.inp
 sed -i "s,filename,${name}_${sub1},g" ${name}_${sub1}_f${dirF}0.inp
 sed -i "s,filename,${name}_${sub2},g" ${name}_${sub2}_f${dirF}0.inp
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub1}.xyz,g" ${name}_${sub1}_f${dirF}0.inp
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub2}.xyz,g" ${name}_${sub2}_f${dirF}0.inp

 #Field-dependent calculations for fragments
 current_field=1
 for i in $( seq 1 $num_fields ); do
    intname1=${name}_${sub1}_f${dirF}${current_field}.inp
    intname2=${name}_${sub2}_f${dirF}${current_field}.inp
    intname1_neg=${name}_${sub1}_f${dirF}-${current_field}.inp
    intname2_neg=${name}_${sub2}_f${dirF}-${current_field}.inp

    for j in $intname1 $intname2 $intname1_neg $intname2_neg; do cp int.inp $j ; done

    sed -i "s,filename,${name}_${sub1}_f${dirF}${current_field},g"  $intname1
    sed -i "s,filename,${name}_${sub2}_f${dirF}${current_field},g"  $intname2
    sed -i "s,filename,${name}_${sub1}_f${dirF}-${current_field},g" $intname1_neg
    sed -i "s,filename,${name}_${sub2}_f${dirF}-${current_field},g" $intname2_neg

    #Change EF keyword
    for j in $intname1 $intname2 $intname1_neg $intname2_neg; do
       display_field=$(printf "%04d" $current_field)
       if [ $j == $intname1 ]; then #Positive field
         sed -i "s,xxchange,0.${display_field},g" $j
       elif [ $j == $intname2 ]; then #Also positive field
         sed -i "s,xxchange,0.${display_field},g" $j
       else #Negative field
         sed -i "s,xxchange,-0.${display_field},g" $j
       fi
    done

    current_field=$((2*current_field))
 done

 #Change filename and coordinates
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub1}.xyz,g" ${name}_${sub1}*
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub2}.xyz,g" ${name}_${sub2}*

 #################### INPUTS FOR INTERACTION ENERGY ########################
 ####################         Supramolecule         ########################

 mkdir interaction && cd interaction && mv ../int.inp .

 #Uncomment the fragment option
 sed -i "s,!Fragments,Fragments,g" int.inp
 sed -i "s,! f, f,g" int.inp
 sed -i "s,!End,End,g" int.inp
 sed -i "s,directory,$int,g" int.inp

 #Change geometry
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}.xyz,g" int.inp

 #Free-field calculation preparation
 cp int.inp ${name}_f${dirF}0.inp
 sed -i "s,filename,${name}_f${dirF}0,g" ${name}_f${dirF}0.inp
 sed -i "/EField/d" ${name}_f${dirF}0.inp
 sed -i "s,sub1.,interaction/${name}_${sub1}_f${dirF}0.,g" ${name}_f${dirF}0.inp
 sed -i "s,sub2.,interaction/${name}_${sub2}_f${dirF}0.,g" ${name}_f${dirF}0.inp

 #Field-dependent calculations for the whole molecule
 current_field=1
 for i in $( seq 1 $num_fields ); do
    intname1=${name}_f${dirF}${current_field}.inp
    intname1_neg=${name}_f${dirF}-${current_field}.inp
    display_field=$(printf "%04d" $current_field)

    for j in $intname1 $intname1_neg ; do cp int.inp $j ; done

    #Change filename
    sed -i "s,filename,${name}_f${dirF}${current_field},g"  $intname1
    sed -i "s,filename,${name}_f${dirF}-${current_field},g" $intname1_neg

    #Change EF keyword --> needs to be modified
    sed -i "s,xxchange,0.${display_field},g"  $intname1
    sed -i "s,xxchange,-0.${display_field},g" $intname1_neg

    #Change tape21 filename
    sed -i "s,sub1.,interaction/${name}_${sub1}_f${dirF}${current_field}.,g"  $intname1
    sed -i "s,sub2.,interaction/${name}_${sub2}_f${dirF}${current_field}.,g"  $intname1
    sed -i "s,sub1.,interaction/${name}_${sub1}_f${dirF}-${current_field}.,g" $intname1_neg
    sed -i "s,sub2.,interaction/${name}_${sub2}_f${dirF}-${current_field}.,g" $intname1_neg

    current_field=$((2*current_field))
 done

 cd ../../ #We are back to $SLURM_SUBMIT_DIR

 ####################      SENDING CALCULATIONS     ########################

 echo "Sending calculations: strain energy"
 echo " "

 srun mkdir -p $TMPDIR/$USER/$SLURM_JOBID
 #Copy all the files to the temporal directory
 cp -r $SLURM_SUBMIT_DIR/$jobdir $TMPDIR/$USER/$SLURM_JOBID/.
 cd $TMPDIR/$USER/$SLURM_JOBID/$jobdir
 echo $SLURM_SUBMIT_DIR > DIR

 #Iterative sending of calculations
 cd strain/

   ns1=${name}_${sub1}
   ns2=${name}_${sub2}
   out_path=$SLURM_SUBMIT_DIR/$jobdir/strain
   for job in ${ns1}_sp ${ns2}_sp ${ns1}_polar ${ns2}_polar ; do
     echo "Sending " $job $PWD
     JOBNAME=$job
     $ADFBIN/adf < ${job}.inp > $out_path/$JOBNAME.out
     #Copy output data
     mv TAPE21 $out_path/$JOBNAME.t21
     if [ -e TAPE13 ]; then
       mv logfile $out_path/$JOBNAME.logfile
     fi
   done
   if [ $ref == "yes" ]; then #also send opt calculations
     for job in ${ns1}_opt ${ns2}_opt; do
       echo "Sending " $job $PWD
       JOBNAME=$job
       $ADFBIN/adf < ${job}.inp > $out_path/$JOBNAME.out
       #Copy output data
       mv TAPE21 $out_path/$JOBNAME.t21
       if [ -e TAPE13 ]; then
         mv logfile $out_path/$JOBNAME.logfile
       fi
     done
   else # also send polar and sp for the reference (reactants)
     for job in ${sub1}_sp ${sub2}_sp ${sub1}_polar ${sub2}_polar ; do
       echo "Sending " $job $PWD
       JOBNAME=$job
       $ADFBIN/adf < ${job}.inp > $out_path/$JOBNAME.out
       #Copy output data
       mv TAPE21 $out_path/$JOBNAME.t21
       if [ -e TAPE13 ]; then
         mv logfile $out_path/$JOBNAME.logfile
       fi
     done
   fi
 cd ../

 echo " "
 echo "Sending calculations: fragments"
 echo " "

 cd int_f$dirF
   #Prepare the sending file
   ls *inp > send1
   sed -i "s,.inp,,g" send1 #Removing the tail (will be added automatically afterwards)
   sed -i "/int_fz0/d" send1 #Remove from the list the example input

   #Running the jobs
   out_path=$SLURM_SUBMIT_DIR/$jobdir/int_f$dirF
   while read job; do
     echo "Sending " $job $PWD
     JOBNAME=$job
     $ADFBIN/adf < ${job}.inp > $out_path/$JOBNAME.out
     #Copy output data
     cp TAPE21 $out_path/interaction/$JOBNAME.t21 # For the further step of the job
     mv TAPE21 $out_path/$JOBNAME.t21
     if [ -e TAPE13 ]; then
       mv logfile $out_path/$JOBNAME.logfile
     fi
   done < send1
 rm send1

 echo " "
 echo "Sending calculations: Int. energy"
 echo " "

 cd interaction
   #Prepares the sending file
   ls *.inp > send2
   sed -i "/int.inp/d" send2 #Removing the example input from the sending list
   sed -i "s,.inp,,g" send2 #Removing the tail (will be added automatically)

   out_path=$SLURM_SUBMIT_DIR/$jobdir/int_f$dirF/interaction
   while read job; do
     echo "Sending " $job $PWD
     JOBNAME=$job
     $ADFBIN/adf < $job.inp > $out_path/$JOBNAME.out
     #Copy output data
     mv TAPE21 $out_path/$JOBNAME.t21
     if [ -e TAPE13 ]; then
       mv logfile $out_path/$JOBNAME.logfile
     fi
   done < send2
 rm send2

 echo " "
 echo "Calculations finished."
 echo " "
 rm -rf $TMPDIR/$USER*

 ######################## END OF CALCULATIONS ##################
 ########################    1D calculation   ##################
 fi #closes the if that requests wheter to send calculations or not

 if [ $whatdo -eq 0 ] || [ $whatdo -eq 2 ] ; then #We do analyze data
 echo "Starting the analysis."
 ###################### STARTING DATA ANALYSIS #################
 ########################    1D calculation   ##################

 cd $strain
 echo "Number of electric fields used for Int energy calculation: " > ../summary.dat
 echo $num_fields >> ../summary.dat
 echo "Direction of the field" >> ../summary.dat
 if [ $dirF == "x" ]; then
  echo "1" >> ../summary.dat
 elif [ $dirF == "y" ]; then
  echo "2" >> ../summary.dat
 else
  echo "3" >> ../summary.dat
 fi

 echo "" >> ../summary.dat
 echo "Strain energy components" >> ../summary.dat
 echo "------------------------" >> ../summary.dat
 echo " " >> ../summary.dat
 # Collect data for the strain energy calculation
 be1=`grep -in1 "^Total Bonding Energy" ${name}_${sub1}_sp.out | awk '{print $4}' `
 be2=`grep -in1 "^Total Bonding Energy" ${name}_${sub2}_sp.out | awk '{print $4}' `
 be3=`grep -in1 "^Total Bonding Energy" ${sub1}_sp.out | awk '{print $4}' `
 be4=`grep -in1 "^Total Bonding Energy" ${sub2}_sp.out | awk '{print $4}' `
 echo "Bonding enery of" ${name}_${sub1} $be1 >> ../summary.dat
 echo "Bonding enery of" ${name}_${sub2} $be2 >> ../summary.dat
 echo "Bonding enery of" $sub1 $be3 >> ../summary.dat
 echo "Bonding enery of" $sub2 $be4 >> ../summary.dat
 echo " " >> ../summary.dat

 for i in ${name}_$sub1 ${name}_$sub2 $sub1 $sub2; do
   echo "Dipole moment of" $i "(Debye)" >> ../summary.dat
   dx1=`grep -A4 "Dipole Moment " ${i}_polar.out | perl -ne 'print if $. == 4' | awk '{print $3}'`
   dy1=`grep -A4 "Dipole Moment " ${i}_polar.out | perl -ne 'print if $. == 4' | awk '{print $4}'`
   dz1=`grep -A4 "Dipole Moment " ${i}_polar.out | perl -ne 'print if $. == 4' | awk '{print $5}'`
   echo "x" $dx1 >> ../summary.dat
   echo "y" $dy1 >> ../summary.dat
   echo "z" $dz1 >> ../summary.dat
   echo " " >> ../summary.dat
 done

 for i in ${name}_$sub1 ${name}_$sub2 $sub1 $sub2; do
   echo "Polarizability of" $i "(a.u.)" >> ../summary.dat
   grep -A4 "Polarizability tensor:" ${i}_polar.out >> ../summary.dat
   echo " " >> ../summary.dat
 done
 cd ../

 #Refine the document before proceeding
 sed -i "/Polarizability tensor:/d" summary.dat

 cd $int/interaction
 # Collect data for the interaction energy calculation

 echo "" >> ../../summary.dat
 echo "Electrostatic Interaction" >> ../../summary.dat
 echo "-------------------------" >> ../../summary.dat
 echo "   F    EI (a.u.)" >> ../../summary.dat
 pr=`grep --max-count=1 "Electrostatic Interaction:" ${name}_f${dirF}0.out | awk '{print $3}' `
 echo " 0000" $pr >> ../../summary.dat
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Electrostatic Interaction:" ${name}_f${dirF}${current_field}.out | awk '{print $3}' `
   display_field=$(printf "%04d" $current_field)
   echo " "$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Electrostatic Interaction:" ${name}_f${dirF}-${current_field}.out | awk '{print $3}' `
   display_field=$(printf "%04d" $current_field)
   echo -$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done

 echo "" >> ../../summary.dat
 echo "Pauli Repulsion" >> ../../summary.dat
 echo "---------------" >> ../../summary.dat
 echo "   F    PR (a.u.)" >> ../../summary.dat
 pr=`grep --max-count=1 "Total Pauli Repulsion:" ${name}_f${dirF}0.out | awk '{print $4}' `
 echo " 0000" $pr >> ../../summary.dat
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Total Pauli Repulsion:" ${name}_f${dirF}${current_field}.out | awk '{print $4}' `
   display_field=$(printf "%04d" $current_field)
   echo " "$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Total Pauli Repulsion:" ${name}_f${dirF}-${current_field}.out | awk '{print $4}' `
   display_field=$(printf "%04d" $current_field)
   echo -$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done


 echo "" >> ../../summary.dat
 echo "Orbital interaction" >> ../../summary.dat
 echo "-------------------" >> ../../summary.dat
 echo "   F    OI (a.u.)" >> ../../summary.dat
 pr=`grep --max-count=1 "Total Orbital Interactions:" ${name}_f${dirF}0.out | awk '{print $4}' `
 echo " 0000" $pr >> ../../summary.dat
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Total Orbital Interactions:" ${name}_f${dirF}${current_field}.out | awk '{print $4}' `
   display_field=$(printf "%04d" $current_field)
   echo " "$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Total Orbital Interactions:" ${name}_f${dirF}-${current_field}.out | awk '{print $4}' `
   display_field=$(printf "%04d" $current_field)
   echo -$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done

 cd ../../ #We are where the summary.dat is located
 cp ../field.inp .
 #Execute fortran code which computes the numerical derivatives and the AMS Taylor expansion
 /home/pbesalu/bin/ams_oeef.exe < field.inp > ams_oeef.out #there field.inp is required

 #Display ending message.
 echo " "
 echo "End of data analysis. Job finished."
 echo "Have a nice day!"
 echo " "
 fi #closes the if for the data analysis

 ####################   END OF DATA ANALYSIS   #################
 ####################  END OF 1D CALCULATIONS  #################


 ################################################################
 ##                                                            ##
 ##                   2D data generation module                ##
 ##                                                            ##
 ################################################################
 if [ $whatdo -eq 3 ]; then #enter de 2d Module for data generation

 # Prepare strain calculation
 cd $strain

 # Transform general strain.inp input into the proper input
 for i in ${name}_$sub1 ${name}_$sub2; do
   cp strain.inp ${i}_sp.inp
   cp strain.inp ${i}_polar.inp
   sed -i "s,filename,${i}_sp,g" ${i}_sp.inp
   sed -i "s,filename,${i}_polar,g" ${i}_polar.inp
   sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_sp.inp
   sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_polar.inp
   if [ $ref == "yes" ]; then #we need to generate the optimization inputs
      cp strain.inp ${i}_opt.inp
      sed -i "s,filename,${i}_opt,g" ${i}_opt.inp
      sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_opt.inp
   fi
 done
 if [ $ref == "no" ]; then #the sp and polar calculations for the isolated reactants are prepared
   for i in $sub1 $sub2 ; do
     cp strain.inp ${i}_sp.inp
     cp strain.inp ${i}_polar.inp
     sed -i "s,filename,${i}_sp,g" ${i}_sp.inp
     sed -i "s,filename,${i}_polar,g" ${i}_polar.inp
     sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_sp.inp
     sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${i}.xyz,g" ${i}_polar.inp
   done
 fi

 # If BETA wants to be computed, this could lead to an inconsistency
 # and code must be checked
 #First, save the line where the Response keyword is found
 #Then, delate polar calculation request for the single point calculations
 p=`grep -in "AOResponse" ${i}_polar.inp | sed "s,:, ,g" | awk '{print $1}' `
 sed -i "${p}d" *_sp.inp # supress such line and the following two
 sed -i "${p}d" *_sp.inp

 if [ $ref == "yes" ]; then
   #Same than before for the optimizations
   sed -i "${p}d" *_opt.inp # supress such line and the following two
   sed -i "${p}d" *_opt.inp

   ## Now, we request that the *opt.inp files have the optimization keyword + freq calculation ENDINPUT is used.
   sed -i '/^ENDINPUT/i Geometry'  *opt.inp
   sed -i '/^Geometry/a Iterations 299' *opt.inp
   sed -i '/^Iterations/a Optim Delocalized' *opt.inp
   sed -i '/^Optim Delocalized/a Converge E=1e-5 grad=1e-5' *opt.inp
   sed -i '/^Converge/a End' *opt.inp
   sed -i '/^ENDINPUT/i AnalyticalFreq' *opt.inp
   sed -i '/^AnalyticalFreq/a End' *opt.inp
   sed -i '/^ENDINPUT/i ScanFreq -1000 200' *opt.inp
   sed -i '/^ENDINPUT/i \ ' *opt.inp

   sed -i 's,Iterations 299,    Iterations 299,g' *opt.inp
   sed -i 's,Optim Delocalized,    Optim Delocalized,g' *opt.inp
   sed -i 's,Converge E=1e-5 grad=1e-5,    Converge E=1e-5 grad=1e-5,g' *opt.inp
 fi

 rm strain.inp
 cd ../ #Exiting strain directory

 #################### INPUTS FOR INTERACTION ENERGY ########################
 ####################        Fragment inputs        ########################

 # Prepare interaction calculation with field
 cd $int

 #Zero-field calculation for fragments
 cp int.inp int_f0.inp
 sed -i "/EField/d" int_f0.inp #Remove EF keyword
 #Prepare filename and coordinates
 cp int_f0.inp ${name}_${sub1}_f0.inp
 cp int_f0.inp ${name}_${sub2}_f0.inp
 sed -i "s,filename,${name}_${sub1},g" ${name}_${sub1}_f0.inp
 sed -i "s,filename,${name}_${sub2},g" ${name}_${sub2}_f0.inp
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub1}.xyz,g" ${name}_${sub1}_f0.inp
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub2}.xyz,g" ${name}_${sub2}_f0.inp

 for dirF in x y z; do  # Go over all three directions. Only parallel fields. Similar procedure as an 1D Calculation
 cp int.inp int_f$dirF.inp

 #Field-dependent calculations for fragments
 if [ $dirF == "x" ]; then
  sed -i "s,EField xx yy zz,EField xxchange 0.0 0.0,g" int_f$dirF.inp
 elif [ $dirF == "y" ]; then
  sed -i "s,EField xx yy zz,EField 0.0 xxchange 0.0,g" int_f$dirF.inp
 else  # If "x" or "y" not specified, use "z" by default
  sed -i "s,EField xx yy zz,EField 0.0 0.0 xxchange,g" int_f$dirF.inp
 fi

 #Field-dependent calculations for fragments
 current_field=1
 for i in $( seq 1 $num_fields ); do
    intname1=${name}_${sub1}_f${dirF}${current_field}.inp
    intname2=${name}_${sub2}_f${dirF}${current_field}.inp
    intname1_neg=${name}_${sub1}_f${dirF}-${current_field}.inp
    intname2_neg=${name}_${sub2}_f${dirF}-${current_field}.inp

    for j in $intname1 $intname2 $intname1_neg $intname2_neg; do cp int_f$dirF.inp $j ; done

    sed -i "s,filename,${name}_${sub1}_f${dirF}${current_field},g"  $intname1
    sed -i "s,filename,${name}_${sub2}_f${dirF}${current_field},g"  $intname2
    sed -i "s,filename,${name}_${sub1}_f${dirF}-${current_field},g" $intname1_neg
    sed -i "s,filename,${name}_${sub2}_f${dirF}-${current_field},g" $intname2_neg

    #Change EF keyword
    for j in $intname1 $intname2 $intname1_neg $intname2_neg; do
       display_field=$(printf "%04d" $current_field)
       if [ $j == $intname1 ]; then #Positive field
         sed -i "s,xxchange,0.${display_field},g" $j
       elif [ $j == $intname2 ]; then #Also positive field
         sed -i "s,xxchange,0.${display_field},g" $j
       else #Negative field
         sed -i "s,xxchange,-0.${display_field},g" $j
       fi
    done

    current_field=$((2*current_field))
 done

 #Change filename and coordinates
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub1}.xyz,g" ${name}_${sub1}*
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub2}.xyz,g" ${name}_${sub2}*
 done # close the cycle for all the three fields (this is the same as in 1D calculation)

 ############ CROSSED TERMS: specific for 2D calculations #################
 for dirF in xy xz yz; do #Prepare inputs for the remaining (crossed) terms.
 cp int.inp int_f$dirF.inp
 if [ $dirF == "xy" ]; then
  sed -i "s,EField xx yy zz,EField xxchange xxchange 0.0,g" int_f$dirF.inp
 elif [ $dirF == "xz" ]; then
  sed -i "s,EField xx yy zz,EField xxchange 0.0 xxchange,g" int_f$dirF.inp
 elif [ $dirF == "yz" ]; then
  sed -i "s,EField xx yy zz,EField 0.0 xxchange xxchange,g" int_f$dirF.inp
 fi

 #Field-dependent calculations for fragments
 current_field=1
 for i in $( seq 1 $num_fields ); do
    intname1=${name}_${sub1}_f${dirF}${current_field}.inp
    intname2=${name}_${sub2}_f${dirF}${current_field}.inp
    intname1_neg=${name}_${sub1}_f${dirF}-${current_field}.inp
    intname2_neg=${name}_${sub2}_f${dirF}-${current_field}.inp

    for j in $intname1 $intname2 $intname1_neg $intname2_neg; do cp int_f$dirF.inp $j ; done

    sed -i "s,filename,${name}_${sub1}_f${dirF}${current_field},g"  $intname1
    sed -i "s,filename,${name}_${sub2}_f${dirF}${current_field},g"  $intname2
    sed -i "s,filename,${name}_${sub1}_f${dirF}-${current_field},g" $intname1_neg
    sed -i "s,filename,${name}_${sub2}_f${dirF}-${current_field},g" $intname2_neg

    #Change EF keyword
    for j in $intname1 $intname2 $intname1_neg $intname2_neg; do
       display_field=$(printf "%04d" $current_field)
       if [ $j == $intname1 ]; then #Positive field
         sed -i "s,xxchange,0.${display_field},g" $j
       elif [ $j == $intname2 ]; then #Also positive field
         sed -i "s,xxchange,0.${display_field},g" $j
       else #Negative field
         sed -i "s,xxchange,-0.${display_field},g" $j
       fi
    done

    current_field=$((2*current_field))
 done
 done #This closes the loop for the calculation of the crossed terms

 #Change filename and coordinates
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub1}.xyz,g" ${name}_${sub1}*
 sed -i "s,coordinates,$SLURM_SUBMIT_DIR/${name}_${sub2}.xyz,g" ${name}_${sub2}*

 rm int*.inp

 #################### INPUTS FOR INTERACTION ENERGY ########################
 ####################         Supramolecule         ########################

 mkdir interaction

 #Copy the files that are required (using $sub1 as the parent for the input)
 cp ${name}_${sub1}_f0.inp interaction/${name}_f0.inp #copy the F=0 input
 for i in x y z xy xz yz; do
   current_field=1
   for j in $( seq 1 $num_fields ); do
     cp ${name}_${sub1}_f${i}${current_field}.inp interaction/${name}_f${i}${current_field}.inp
     cp ${name}_${sub1}_f${i}-${current_field}.inp interaction/${name}_f${i}-${current_field}.inp
     current_field=$((2*current_field))
   done
 done

 cd interaction

 # Set up the proper title and geometry
 sed -i "s,TITLE ${name}_${sub1},TITLE ${name},g" *.inp
 sed -i "s,${name}_${sub1}.xyz,${name}.xyz,g" *.inp

 #Uncomment the fragment option
 sed -i "s,!Fragments,Fragments,g" *.inp
 sed -i "s,! f, f,g" *.inp
 sed -i "s,!End,End,g" *.inp
 sed -i "s,directory,$int,g" *.inp

 #Free-field calculation preparation
 #sed -i "s,filename,${name}_f0,g" ${name}_f0.inp
 #sed -i "/EField/d" ${name}_f0.inp
 sed -i "s,sub1.,interaction/${name}_${sub1}_f0.,g" ${name}_f0.inp
 sed -i "s,sub2.,interaction/${name}_${sub2}_f0.,g" ${name}_f0.inp

 #Field-dependent calculations for the whole molecule
 for dirF in x y z xy xz yz; do
  current_field=1
  for i in $( seq 1 $num_fields ); do
    intname1=${name}_f${dirF}${current_field}.inp
    intname1_neg=${name}_f${dirF}-${current_field}.inp
    display_field=$(printf "%04d" $current_field)

    #Change filename
    sed -i "s,filename,${name}_f${dirF}${current_field},g"  $intname1
    sed -i "s,filename,${name}_f${dirF}-${current_field},g" $intname1_neg

    #Change tape21 filename
    sed -i "s,sub1.,interaction/${name}_${sub1}_f${dirF}${current_field}.,g"  $intname1
    sed -i "s,sub2.,interaction/${name}_${sub2}_f${dirF}${current_field}.,g"  $intname1
    sed -i "s,sub1.,interaction/${name}_${sub1}_f${dirF}-${current_field}.,g" $intname1_neg
    sed -i "s,sub2.,interaction/${name}_${sub2}_f${dirF}-${current_field}.,g" $intname1_neg
    current_field=$((2*current_field))
  done #Closes all the fields preparation
 done #Closes all the directions (including crossed terms) preparation


 cd ../../ #We are back to $SLURM_SUBMIT_DIR

 ####################      SENDING CALCULATIONS     ########################

 echo "Sending calculations: strain energy"
 echo " "

 srun mkdir -p $TMPDIR/$USER/$SLURM_JOBID
 #Copy all the files to the temporal directory
 cp -r $SLURM_SUBMIT_DIR/$jobdir $TMPDIR/$USER/$SLURM_JOBID/.
 cd $TMPDIR/$USER/$SLURM_JOBID/$jobdir
 echo $SLURM_SUBMIT_DIR > DIR

 #Iterative sending of calculations
 cd strain/

   ns1=${name}_${sub1}
   ns2=${name}_${sub2}
   out_path=$SLURM_SUBMIT_DIR/$jobdir/strain
   for job in ${ns1}_sp ${ns2}_sp ${ns1}_polar ${ns2}_polar ; do
     echo "Sending " $job $PWD
     JOBNAME=$job
     $ADFBIN/adf < ${job}.inp > $out_path/$JOBNAME.out
     #Copy output data
     mv TAPE21 $out_path/$JOBNAME.t21
     if [ -e TAPE13 ]; then
       mv logfile $out_path/$JOBNAME.logfile
     fi
   done
   if [ $ref == "yes" ]; then #also send opt calculations
     for job in ${ns1}_opt ${ns2}_opt; do
       echo "Sending " $job $PWD
       JOBNAME=$job
       $ADFBIN/adf < ${job}.inp > $out_path/$JOBNAME.out
       #Copy output data
       mv TAPE21 $out_path/$JOBNAME.t21
       if [ -e TAPE13 ]; then
         mv logfile $out_path/$JOBNAME.logfile
       fi
     done
   else # also send polar and sp for the reference (reactants)
     for job in ${sub1}_sp ${sub2}_sp ${sub1}_polar ${sub2}_polar ; do
       echo "Sending " $job $PWD
       JOBNAME=$job
       $ADFBIN/adf < ${job}.inp > $out_path/$JOBNAME.out
       #Copy output data
       mv TAPE21 $out_path/$JOBNAME.t21
       if [ -e TAPE13 ]; then
         mv logfile $out_path/$JOBNAME.logfile
       fi
     done
   fi
 cd ../

 echo " "
 echo "Sending calculations: fragments"
 echo " "

 cd int_f
   #Prepare the sending file
   ls *inp > send1
   sed -i "s,.inp,,g" send1 #Removing the tail (will be added automatically afterwards)
   sed -i "/int_fz0/d" send1 #Remove from the list the example sample input

   #Running the jobs
   out_path=$SLURM_SUBMIT_DIR/$jobdir/int_f
   while read job; do
     echo "Sending " $job $PWD
     JOBNAME=$job
     $ADFBIN/adf < ${job}.inp > $out_path/$JOBNAME.out
     #Copy output data
     cp TAPE21 $out_path/interaction/$JOBNAME.t21 # For the further step of the job
     mv TAPE21 $out_path/$JOBNAME.t21
     if [ -e TAPE13 ]; then
       mv logfile $out_path/$JOBNAME.logfile
     fi
   done < send1
 rm send1

 echo " "
 echo "Sending calculations: Int. energy"
 echo " "

 cd interaction
   #Prepares the sending file
   ls *.inp > send2
   sed -i "s,.inp,,g" send2 #Removing the tail (will be added automatically)

   out_path=$SLURM_SUBMIT_DIR/$jobdir/int_f/interaction
   while read job; do
     echo "Sending " $job $PWD
     JOBNAME=$job
     $ADFBIN/adf < $job.inp > $out_path/$JOBNAME.out
     #Copy output data
     mv TAPE21 $out_path/$JOBNAME.t21
     if [ -e TAPE13 ]; then
       mv logfile $out_path/$JOBNAME.logfile
     fi
   done < send2
 rm send2

 echo " "
 echo "Calculations finished."
 echo " "
 rm -rf $TMPDIR/$USER*

 fi # End of 2D data generation module

 if [ $whatdo -eq 3 ] || [ $whatdo -eq 4 ] ; then #we should analyze data
 #################### END OF 2D DATA GENERATION MODULE #################
 #######################################################################
 #                                                                     #
 #                        2D DATA ANALYSIS MODULE                      #
 #                                                                     #
 #######################################################################
 # this module first generates the three replicas of the 1D (for x,y and z)
 # which are equivalent to the ones obtained through 1D module
 # and then generates the data to input to pyplot for the 2D plots generation

 cd $strain
 echo "Number of electric fields used for Int energy calculation: " > ../summary.dat
 echo $num_fields >> ../summary.dat

 echo "" >> ../summary.dat
 echo "Strain energy components" >> ../summary.dat
 echo "------------------------" >> ../summary.dat
 echo " " >> ../summary.dat
 # Collect data for the strain energy calculation
 be1=`grep -in1 "^Total Bonding Energy" ${name}_${sub1}_sp.out | awk '{print $4}' `
 be2=`grep -in1 "^Total Bonding Energy" ${name}_${sub2}_sp.out | awk '{print $4}' `
 be3=`grep -in1 "^Total Bonding Energy" ${sub1}_sp.out | awk '{print $4}' `
 be4=`grep -in1 "^Total Bonding Energy" ${sub2}_sp.out | awk '{print $4}' `
 echo "Bonding enery of" ${name}_${sub1} $be1 >> ../summary.dat
 echo "Bonding enery of" ${name}_${sub2} $be2 >> ../summary.dat
 echo "Bonding enery of" $sub1 $be3 >> ../summary.dat
 echo "Bonding enery of" $sub2 $be4 >> ../summary.dat
 echo " " >> ../summary.dat

 for i in ${name}_$sub1 ${name}_$sub2 $sub1 $sub2; do
   echo "Dipole moment of" $i "(Debye)" >> ../summary.dat
   dx1=`grep -A4 "Dipole Moment " ${i}_polar.out | perl -ne 'print if $. == 4' | awk '{print $3}'`
   dy1=`grep -A4 "Dipole Moment " ${i}_polar.out | perl -ne 'print if $. == 4' | awk '{print $4}'`
   dz1=`grep -A4 "Dipole Moment " ${i}_polar.out | perl -ne 'print if $. == 4' | awk '{print $5}'`
   echo "x" $dx1 >> ../summary.dat
   echo "y" $dy1 >> ../summary.dat
   echo "z" $dz1 >> ../summary.dat
   echo " " >> ../summary.dat
 done

 for i in ${name}_$sub1 ${name}_$sub2 $sub1 $sub2; do
   echo "Polarizability of" $i "(a.u.)" >> ../summary.dat
   grep -A4 "Polarizability tensor:" ${i}_polar.out >> ../summary.dat
   echo " " >> ../summary.dat
 done
 cd ../

 #Refine the document before proceeding
 sed -i "/Polarizability tensor:/d" summary.dat

 cd int_f/interaction
 # Collect data for the interaction energy calculation

 for dirF in x y z xy xz yz; do

 echo "" >> ../../summary.dat
 echo "Electrostatic Interaction" $dirF >> ../../summary.dat
 echo "-------------------------" >> ../../summary.dat
 echo "   F    EI (a.u.)" >> ../../summary.dat
 pr=`grep --max-count=1 "Electrostatic Interaction:" ${name}_f0.out | awk '{print $3}' `
 echo " 0000" $pr >> ../../summary.dat
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Electrostatic Interaction:" ${name}_f${dirF}${current_field}.out | awk '{print $3}' `
   display_field=$(printf "%04d" $current_field)
   echo " "$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Electrostatic Interaction:" ${name}_f${dirF}-${current_field}.out | awk '{print $3}' `
   display_field=$(printf "%04d" $current_field)
   echo -$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done

 echo "" >> ../../summary.dat
 echo "Pauli Repulsion" $dirF >> ../../summary.dat
 echo "---------------" >> ../../summary.dat
 echo "   F    PR (a.u.)" >> ../../summary.dat
 pr=`grep --max-count=1 "Total Pauli Repulsion:" ${name}_f0.out | awk '{print $4}' `
 echo " 0000" $pr >> ../../summary.dat
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Total Pauli Repulsion:" ${name}_f${dirF}${current_field}.out | awk '{print $4}' `
   display_field=$(printf "%04d" $current_field)
   echo " "$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Total Pauli Repulsion:" ${name}_f${dirF}-${current_field}.out | awk '{print $4}' `
   display_field=$(printf "%04d" $current_field)
   echo -$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done


 echo "" >> ../../summary.dat
 echo "Orbital interaction" $dirF >> ../../summary.dat
 echo "-------------------" >> ../../summary.dat
 echo "   F    OI (a.u.)" >> ../../summary.dat
 pr=`grep --max-count=1 "Total Orbital Interactions:" ${name}_f0.out | awk '{print $4}' `
 echo " 0000" $pr >> ../../summary.dat
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Total Orbital Interactions:" ${name}_f${dirF}${current_field}.out | awk '{print $4}' `
   display_field=$(printf "%04d" $current_field)
   echo " "$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done
 current_field=1
 for i in $( seq 1 $num_fields ); do
   pr=`grep --max-count=1 "Total Orbital Interactions:" ${name}_f${dirF}-${current_field}.out | awk '{print $4}' `
   display_field=$(printf "%04d" $current_field)
   echo -$display_field $pr >> ../../summary.dat
   current_field=$((2*current_field))
 done

 done #closing the iterations over the several EF

 cd ../../ #We are where the summary.dat is located
 cp ../field.inp .
 #Execute fortran code which computes the numerical derivatives and the AMS Taylor expansion
 /home/pbesalu/bin/ams_oeef_2D.exe < field.inp > ams_oeef.out #there field.inp is required

 #Display ending message.
 echo " "
 echo "End of data analysis. Job finished."
 echo " "
 fi #closes the if for the data analysis
