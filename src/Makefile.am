########################################################################
##- Author@mail: Nicolas.Grima@univ-brest.fr
##- Date       : November 2015
########################################################################
FORTRAN_TOOLS=.

SUFFIXES = .F90 .f90 .o .mod

.f90.mod:
	$(FC) $(FCFLAGS) -c $<

.f90.o:
	$(FC) $(FCFLAGS) -c -o $@ $<

.c.o:
	$(CC) $(CFLAGS) -c  $@ $<

clean-local:
	-rm -f *.mod *.o *.log

##
SUBDIRS	     = .
bin_PROGRAMS = mgroms
##bin_SCRIPTS  = namelist
##EXTRA_DIST   = namelist
##
##
mgroms_SOURCES =         \
	mg_grids.f90         \
	mg_mpi.f90           \
	mg_optimized.f90     \
	mg_test.f90
#	mg_variables.f90     \
#	mg_smoother.f90      \
#	mg_interpolation.f90 \
#	mg_simpleop.f90      \
#	mg_main.f90
##
