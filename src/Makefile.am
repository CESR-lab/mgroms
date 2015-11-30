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
	mg_mpi.f90           \
	mg_grids.f90         \
	mg_define_rhs.f90    \
	mg_define_matrix.f90 \
	mg_relax_line.f90    \
	mg_testhalo.f90
#	mg_testrelax.f90     \
#	mg_variables.f90     \
#	mg_smoother.f90      \
#	mg_interpolation.f90 \
#	mg_simpleop.f90      \
#	mg_main.f90
##
