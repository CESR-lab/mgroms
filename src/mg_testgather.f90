program mg_testgather

  use mg_mpi ! everything will come from the outside !!!

  use mg_grids

  implicit none

  integer(kind=is):: nxg    ! global x dimension
  integer(kind=is):: nyg    ! global y dimension
  integer(kind=is):: nzg    ! z dimension
  !!integer(kind=is):: nhalo  ! number of halo points
  integer(kind=is):: npxg   ! number of processes in x
  integer(kind=is):: npyg   ! number of processes in y
  integer(kind=is):: it     ! iteration loop number
  integer(kind=is):: nit    ! number of iterations

  integer(kind=is):: nsweeps

  integer(kind=is):: lev, ierr
  real(kind=8)    :: res,res0,conv

  nxg   = 128
  nyg   = 128
  nzg   = 128
  nhalo = 1

  npxg  = 4
  npyg  = 4

  nit     = 10
  nsweeps = 1

  call init_mpi(nxg, nyg, nzg, npxg, npyg)

  call define_grids(nhalo,npxg,npyg)


  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  if (myrank.eq.0)write(*,*)"---------- grid hierarchy ----------"
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
 call MPI_Barrier( MPI_COMM_WORLD ,ierr)
 if (myrank.eq.0)write(*,*)"---------- check families ----------"
 call MPI_Barrier( MPI_COMM_WORLD ,ierr)
 do lev=1,nlevs-1
    if (grid(lev)%gather.eq.1)then
       call MPI_Barrier( MPI_COMM_WORLD ,ierr)
       if (myrank.eq.0)then
          write(*,'(A,I2,A)')"    ---  Level =",lev,' ---'
       endif
       call MPI_Barrier( MPI_COMM_WORLD ,ierr)
       write(*,'(A,I,A,I,A,I)')"rank=",myrank," /family=",grid(lev)%family," /color=",grid(lev)%color
    endif
    call MPI_Barrier( MPI_COMM_WORLD ,ierr)
 enddo
100 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)




  call mg_mpi_finalize()

end program mg_testgather
