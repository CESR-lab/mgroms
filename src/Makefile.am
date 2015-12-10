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
	-rm -f *.mod *.o *.log *__genmod*

##
SUBDIRS	       = .
bin_PROGRAMS   = mainmodel testhalo testrelax testcoarsening testgather testintergrids
check_PROGRAMS = $(bin_PROGRAMS)
TESTS          = $(check_PROGRAMS)
##
##
mainmodel_SOURCES=       \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_grids.f90         \
	mg_mpi_exchange.f90  \
	mg_gather.f90        \
	mg_define_matrix.f90 \
	mg_relax.f90         \
	mg_intergrids.f90    \
	mg_solvers.f90       \
	nhydro.f90           \
	main_model.f90
##
testhalo_SOURCES =       \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_grids.f90         \
	mg_mpi_exchange.f90  \
	mg_define_rhs.f90    \
	mg_testhalo.f90
##
testrelax_SOURCES =      \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_grids.f90         \
	mg_mpi_exchange.f90  \
	mg_gather.f90        \
	mg_define_rhs.f90    \
	mg_define_matrix.f90 \
	mg_relax.f90         \
	mg_testrelax.f90
##
testcoarsening_SOURCES = \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_grids.f90         \
	mg_mpi_exchange.f90  \
	mg_gather.f90        \
	mg_define_rhs.f90    \
	mg_define_matrix.f90 \
	mg_relax.f90         \
	mg_intergrids.f90    \
	mg_testcoarsening.f90
##
testgather_SOURCES =     \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_grids.f90         \
	mg_gather.f90        \
	mg_testgather.f90
##
testintergrids_SOURCES =     \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_grids.f90         \
	mg_mpi_exchange.f90  \
	mg_gather.f90        \
	mg_intergrids.f90    \
	mg_testintergrids.f90

