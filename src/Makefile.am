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
	-rm -f *.mod *.o

cleanf:
	-rm -f *.nc fort.* *.dump *.xml gmon.out *.log *__genmod*

##
SUBDIRS	       = .
bin_PROGRAMS   =   \
	mainmodel      \
	testseamount   \
	testgalerkin   
#	testhalo       \
#	testrelax      \
#	testcoarsening \
#	testgather     \
#	testintergrids \
#	testnetcdf     \
#check_PROGRAMS = $(bin_PROGRAMS)
#TESTS          = $(check_PROGRAMS)
##
##
mainmodel_SOURCES=       \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_netcdf_out.f90    \
	mg_grids.f90         \
	mg_mpi_exchange.f90  \
	mg_gather.f90        \
	mg_define_matrix.f90 \
	mg_relax.f90         \
	mg_intergrids.f90    \
	mg_solvers.f90       \
	nhydro.f90           \
	main_model.f90
testseamount_SOURCES=    \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_netcdf_out.f90    \
	mg_grids.f90         \
	mg_mpi_exchange.f90  \
	mg_gather.f90        \
	mg_define_matrix.f90 \
	mg_relax.f90         \
	mg_intergrids.f90    \
	mg_solvers.f90       \
	nhydro.f90           \
	mg_testseamount.f90
##
testgalerkin_SOURCES =   \
	mg_mpi.f90           \
	mg_tictoc.f90        \
	mg_namelist.f90      \
	mg_netcdf_out.f90    \
	mg_grids.f90         \
	mg_mpi_exchange.f90  \
	mg_gather.f90        \
	mg_define_matrix.f90 \
	mg_intergrids.f90    \
	mg_relax.f90         \
	mg_solvers.f90       \
	mg_testgalerkin.f90
##
#testhalo_SOURCES =       \
#	mg_mpi.f90           \
#	mg_tictoc.f90        \
#	mg_namelist.f90      \
#	mg_grids.f90         \
#	mg_mpi_exchange.f90  \
#	mg_define_rhs.f90    \
#	mg_testhalo.f90
##
#testrelax_SOURCES =      \
#	mg_mpi.f90           \
#	mg_tictoc.f90        \
#	mg_namelist.f90      \
#	mg_netcdf_out.f90    \
#	mg_grids.f90         \
#	mg_mpi_exchange.f90  \
#	mg_gather.f90        \
#	mg_define_rhs.f90    \
#	mg_define_matrix.f90 \
#	mg_relax.f90         \
#	mg_testrelax.f90
##
#testcoarsening_SOURCES = \
#	mg_mpi.f90           \
#	mg_tictoc.f90        \
#	mg_namelist.f90      \
#	mg_netcdf_out.f90    \
#	mg_grids.f90         \
#	mg_mpi_exchange.f90  \
#	mg_gather.f90        \
#	mg_define_rhs.f90    \
#	mg_define_matrix.f90 \
#	mg_relax.f90         \
#	mg_intergrids.f90    \
#	mg_testcoarsening.f90
##
#testgather_SOURCES =     \
#	mg_mpi.f90           \
#	mg_tictoc.f90        \
#	mg_namelist.f90      \
#	mg_netcdf_out.f90    \
#	mg_grids.f90         \
#	mg_gather.f90        \
#	mg_testgather.f90
##
#testintergrids_SOURCES = \
#	mg_mpi.f90           \
#	mg_tictoc.f90        \
#	mg_namelist.f90      \
#	mg_netcdf_out.f90    \
#	mg_grids.f90         \
#	mg_mpi_exchange.f90  \
#	mg_gather.f90        \
#	mg_intergrids.f90    \
#	mg_testintergrids.f90
##
#testnetcdf_SOURCES =     \
#	mg_mpi.f90           \
#	mg_tictoc.f90        \
#	mg_namelist.f90      \
#	mg_netcdf_out.f90    \
#	mg_grids.f90         \
#	mg_mpi_exchange.f90  \
#	mg_gather.f90        \
#	mg_define_rhs.f90    \
#	mg_define_matrix.f90 \
#	mg_relax.f90         \
#	mg_testnetcdf.f90
##

