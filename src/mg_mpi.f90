module mg_mpi

  implicit none
!  use mpi
  include 'mpif.h'


  integer(kind=4), parameter:: is=4, rl=8

  integer(kind=is) :: myrank
  integer(kind=is) :: nprocs

  ! local dimensions which will come from the outside (ROMS)
  integer(kind=is) :: nxo, nyo, nzo

  ! http://www.open-mpi.org/doc/v1.8/man3/MPI_Recv_init.3.php

contains

  !----------------------------------------
  subroutine init_mpi(nxg, nyg, nzg, npxg, npyg)

    integer(kind=is), intent(in) :: nxg, nyg, nzg
    integer(kind=is), intent(in) :: npxg, npyg

    integer(kind=is) :: ierr

    call mpi_init(ierr)

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    call mpi_comm_size(mpi_comm_world, nprocs, ierr)

    if (nprocs /= (npxg*npyg)) then
       if(myrank.eq.0)then
          write(*,*) "Error: in number of processes !"
          write(*,'(A,I3,A)') "call with mpirun -np ",npxg*npyg," ..."
       endif
       stop
    endif

    ! WARNING, non divisibility issues !
    nxo = nxg / npxg 
    nyo = nyg / npyg
    nzo = nzg

  end subroutine init_mpi

  !----------------------------------------
  SUBROUTINE mg_mpi_finalize

    integer(kind=is) :: ierr

    ! Desactivation de MPI
    CALL MPI_FINALIZE(ierr)

  END SUBROUTINE mg_mpi_finalize

  !----------------------------------------
  subroutine global_max(maxloc,maxglo)
    ! return the global max: maxglo
    ! using the local max on each subdomain
    real(kind=rl),intent(in) :: maxloc
    real(kind=rl),intent(out) :: maxglo

    integer(kind=is) :: ierr

    ! note: the global comm using MPI_COMM_WORLD is over-kill for levels 
    ! where subdomains are gathered
    ! this is not optimal, but not wrong
    call MPI_ALLREDUCE(maxloc,maxglo,1,MPI_DOUBLE_PRECISION,MPI_max,MPI_COMM_WORLD,ierr)   

  end subroutine global_max

  !----------------------------------------
  subroutine global_sum(lev,sumloc,sumglo)
    ! return the global sum: sumglo
    ! using the local sum on each subdomain
    integer(kind=is),intent(in) :: lev
    real(kind=rl),intent(in) :: sumloc
    real(kind=rl),intent(out) :: sumglo

    integer(kind=is) :: ierr

    ! note: the global comm using MPI_COMM_WORLD is over-kill for levels 
    ! where subdomains are gathered
    call MPI_ALLREDUCE(sumloc,sumglo,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)   
    ! therefore we need to rescale the global sum
    sumglo = sumglo * (grid(lev)%npx*grid(lev)%npy)/(grid(1)%npx*grid(1)%npy)

  end subroutine global_sum

end module mg_mpi
