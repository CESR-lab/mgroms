program mg_testcoarsening

  use mpi
  use mg_mpi
  use mg_tictoc
  use mg_grids
  use mg_define_rhs
  use mg_define_matrix
  use mg_relax
  use mg_intergrids

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension
  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y
  integer(kind=ip):: it     ! iteration loop number
  integer(kind=ip):: nit    ! number of iterations

  integer(kind=ip):: nsweeps

  integer(kind=ip):: lev, ierr, np
  real(kind=rp)    :: res,res0,conv
  integer(kind=ip) :: nx, ny, nz  ! local dimensions

  !- timing
  call tic(1,'mg_testcoarsening')

  nxg   = 128
  nyg   = 128
  nzg   = 128

  npxg  = 2
  npyg  = 2

  nit     = 10
  nsweeps = 1

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  if (np /= (npxg*npyg)) then
     write(*,*) "Error: in number of processes !"
     stop -1
  endif

  nx = nxg / npxg
  ny = nyg / npyg
  nz = nzg

  !- Enter in nhydro -!
  call mg_mpi_init()

  call define_grids(npxg, npyg, nx, ny, nz)

  call define_neighbours()

  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  if (myrank.eq.0)then
    do lev=1,nlevs
       if (grid(lev)%gather.eq.0)then
          write(*,100)"lev=",lev,": ", &
                     grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                     " on ",grid(lev)%npx,' x',grid(lev)%npy," procs"
       else
          write(*,100)"lev=",lev,": ", &
                     grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                     " on ",grid(lev)%npx,' x',grid(lev)%npy," procs / gather"
       endif
    enddo
 endif
100 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)

  call define_rhs(nxg, nyg, npxg)

  write(*,*)"rhs  done"

  call define_matrices()

  write(*,*)"define matrices done"

  do lev=1,nlevs
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)
     if (myrank.eq.0)then
        write(*,'(A,I2,A,8F6.3)')"stencil on lev =",lev," / ",&
             grid(lev)%cA(:,grid(lev).nz/2,grid(lev).ny/2,grid(lev).nx/2)
!        write(*,1010)grid(lev)%cA(:,grid(lev).nz/2,grid(lev).ny/2,grid(lev).nx/2)
     endif
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  end do
1010 format(8F6.3)

  ! coarsen RHS on all grids
  do lev=1,nlevs-1
     grid(lev)%r = grid(lev)%b
     call coarse2fine(lev)
  enddo

  ! check smoothing on all grids
  do lev=1,nlevs
     if (myrank.eq.0)then
        write(*,'(A,I2,A,I3,A,I3,A,I3,A)')"check relaxation on lev = ",lev,&
             " / (",grid(lev)%nx,",",grid(lev)%ny,",",grid(lev)%nz,")"
     endif
     res0=0.
     do it=1, nit
        call relax(lev,nsweeps)
        call compute_residual(lev,res)
        conv = log(res0/res)/log(10.)
        if (myrank.eq.0) then
           write(*,'(A,I2,A,F6.3,A,F6.3)')"  ite=",it," - res=",res," - conv rate =",conv
        endif
        res0=res
     enddo
  enddo

!  call check_solution(lev)

  call mpi_finalize(ierr)

  !- timing
  call toc(1,'mg_testcoarsening')
  call print_tictoc(myrank)

end program mg_testcoarsening
