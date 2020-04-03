#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=4


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Bernadou Basis p3
#----------------------------
cd Validation

echo "Running linear bending validation"
mkdir RESLT

echo "Running rectangular sheet validation #1"
../rectangular_sheet_kpb --validation > OUTPUT_rectangular_sheet_kpb
echo "done"
echo " " >> validation.log
echo "Linear Bending validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > rectangular_sheet_kpb_results_1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/rectangular_sheet_kpb_result_1.dat.gz   \
  rectangular_sheet_kpb_results_1.dat  >> validation.log
fi

echo "Running rectangular sheet validation #2"
../run_code.sh ../rectangular_sheet_kpb | tee OUTPUT_rectangular_sheet_kpb
echo "done"
echo " " >> validation.log
echo "Linear Bending validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat errors.dat > rectangular_sheet_kpb_results_2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/rectangular_sheet_kpb_result_2.dat.gz   \
  rectangular_sheet_kpb_results_2.dat  >> validation.log
fi

echo "Running circular sheet validation"
../circular_disc_kpb --validation > OUTPUT_circular_disc_kpb
echo "done"
echo " " >> validation.log
echo "Linear Bending validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > circular_disc_kpb_results_1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/circular_disc_kpb_result_1.dat.gz   \
  circular_disc_kpb_results_1.dat  >> validation.log
fi

echo "Running rectangular sheet demo validation"
../rectangular_sheet_kpb_demo --validate > OUTPUT_rectangular_sheet_kpb_demo
echo "done"
echo " " >> validation.log
echo "Linear Bending validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > rectangular_sheet_kpb_demo_result_1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/rectangular_sheet_kpb_demo_result_1.dat.gz   \
  rectangular_sheet_kpb_demo_result_1.dat  >> validation.log
fi
# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
 . $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
 exit 10
