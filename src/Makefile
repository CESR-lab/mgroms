SRCS = mg_grids.F mg_mpi.F mg_optimized.F mg_test.F

OBJS = $(SRCS:.F=.o)

CPP = cpp
CPPFLG =

FC90 = mpif90
FFLG = -c -O3  -extend-source
OMP_FLAG = 
LIBS = 

.F.o:	$(SRCS)
	$(FC90) $(FFLG) $(OMP_FLAG) -c $*.F -o$*.o

all:	$(OBJS)
	$(FC90) $(LFLG) $(OMP_FLAG) $(OBJS) $(LIBS) -o mg_test

clean:	$(SRCS)
	$(RM) -f $(OBJS) *.mod