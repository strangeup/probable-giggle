#!/bin/sh
# Quick script to batch the data aquisition
# The Directory
driver_name="rectangular_sheet_kpb"
initial_area=0.033 #0.032592
arguments="--case 2"

# Loop over 7 plot points
for it in `seq 0 6`;
do
elementarea=$(echo "scale =10 ; $initial_area * (0.65049)^($it)" | bc)
echo "##########################################"
echo Running:
echo "./$driver_name  $arguments --element_area $elementarea >> output"
 ./$driver_name  $arguments --element_area $elementarea >> output
done

# Grab the data
grep L2 output | grep -o "[0-9]\+\.[0-9e+-]*" > errors.dat
grep area output | grep -o "[0-9]\+\.[0-9e+-]*" > areas.dat
