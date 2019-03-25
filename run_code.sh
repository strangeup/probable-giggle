#!/bin/sh
# Quick script to batch the data aquisition
# The Directory
# Did we get a exe name?
if [ $# -eq 0 ]
  then
    echo "No arguments supplied running rectangular_sheet_kpb"
    driver_name="rectangular_sheet_kpb"
else 
  driver_name=$1
fi

# rectangle initial_area 0.033 and circle clamped/resting 0.2
initial_area=0.033 # 0.2 #0.033 #0.05 #0.1 #33 #0.032592
arguments="--case 2" #"--case 2"
# rectangle base 0.65049 and circle clamped resting 0.5 
base=0.65049 #0.5 #0.65049 #$(echo "scale =10 ; 1/sqrt(2)" | bc) #0.70711" #"0.5" #0.65049 

# Loop over 7 plot points
for it in `seq 0 6`;
do
elementarea=$(echo "scale =10 ; $initial_area * ($base)^($it)" | bc)
echo "##########################################"
echo Running:
echo "./$driver_name  $arguments --element_area $elementarea >> output"
./$driver_name  $arguments --element_area $elementarea >> output
done

# Grab the data
grep L2 output | grep -o "[0-9]\+\.[0-9e+-]*" > errors.dat
grep area output | grep -o "[0-9]\+\.[0-9e+-]*" > areas.dat
