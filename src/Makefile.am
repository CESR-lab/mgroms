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
SUBDIRS	       = .
bin_PROGRAMS   = testhalo testrelax testcoarsening testgather
check_PROGRAMS = $(bin_PROGRAMS)
TESTS          = $(check_PROGRAMS)
##
##
testhalo_SOURCES =       \
	mg_mpi.f90           \
	mg_grids.f90         \
	mg_define_rhs.f90    \
	mg_testhalo.f90
##
testrelax_SOURCES =      \
	mg_mpi.f90           \
	mg_grids.f90         \
	mg_define_rhs.f90    \
	mg_define_matrix.f90 \
	mg_relax_line.f90    \
	mg_testrelax.f90
##
testcoarsening_SOURCES = \
	mg_mpi.f90           \
	mg_grids.f90         \
	mg_define_rhs.f90    \
	mg_define_matrix.f90 \
	mg_relax_line.f90    \
	mg_restrict.f90      \
	mg_testcoarsening.f90
##
testgather_SOURCES =       \
	mg_mpi.f90           \
	mg_grids.f90         \
	mg_gather.f90        \
	mg_testgather.f90

