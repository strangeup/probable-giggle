# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS =  rectangular_sheet_kpb circular_disc_three_pins_kpb semi_annular_kpb \
 circular_disc_kpb
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Local sources that rectangular_sheet_kpb code depends on:
rectangular_sheet_kpb_SOURCES = \
 rectangular_sheet_kirchhoff_plate_bending.cc

circular_disc_three_pins_kpb_SOURCES = \
 circular_disc_uniform_pressure_three_pins.cc

semi_annular_kpb_SOURCES = \
 semi_annular_plate_pinned.cc

circular_disc_kpb_SOURCES = \
 circular_disc_uniform_pressure.cc

#---------------------------------------------------------------------------
rectangular_sheet_kpb_LDADD = -L@libdir@ -lC1_linear_plate_bending \
 -lC1_basis -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

circular_disc_three_pins_kpb_LDADD = -L@libdir@ -lC1_linear_plate_bending \
 -lC1_basis -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

semi_annular_kpb_LDADD = -L@libdir@ -lC1_linear_plate_bending \
 -lC1_basis -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

circular_disc_kpb_LDADD = -L@libdir@ -lC1_linear_plate_bending \
 -lC1_basis -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Local sources that Jack's own code depends on: This code also uses
# objects from Jack's own library. The source code for this library
# lives in user_src/jack_cool and gets compiled in to the jack_cool 
# library. If we use it, we have to link against it...
#jacks_own_code_SOURCES = jacks_own_code.cc

# Required libraries: Jack's own library ("jack_cool") and oomph-lib's 
# and "generic" and "poisson" libraries. Following the logic
# explained above, Jack's own library should go first since it may
# (or may not) depend on objects in the other libraries, whereas they will
# certainly not depend on objects in his. 
# jacks_own_code_LDADD = -L@libdir@ -ljack_cool -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.

AM_CPPFLAGS += -I@includedir@ -I$(pathtoeigen) -I$(pathtospectra) -Wall  -g -DRANGE_CHECK -DPARANOID
