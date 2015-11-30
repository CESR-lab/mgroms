program mg_testcoarsening

  use mg_mpi ! everything will come from the outside !!!

  use mg_grids
  use mg_define_rhs
  use mg_define_matrix
  use mg_relax
  use mg_restrict

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
  real(kind=8)    :: res

  nxg   = 128
  nyg   = 128
  nzg   = 128
  nhalo = 1

  npxg  = 2
  npyg  = 2

  nit     = 10
  nsweeps = 1

  call init_mpi(nxg, nyg, nzg, npxg, npyg)

  call define_grids(nhalo,npxg,npyg)

  call define_rhs(nxg, nyg, npxg)

  call define_matrix_simple()

  lev = 1

  call coarsen_matrix(lev)

  lev = lev+1
  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  if (myrank.eq.0)then
     write(*,*)grid(lev)%cA(:,grid(lev).nz/2,grid(lev).ny/2,grid(lev).nx/2)
  endif
  call MPI_Barrier( MPI_COMM_WORLD ,ierr)

  call restrict_xyz(1,2,grid(1)%b,grid(2)%b)

!!$  call compute_residual(lev,res)
!!$
  do it=1, nit
     call relax_line(lev,nsweeps)
     call compute_residual(lev,res)
     if (myrank.eq.0)then
        write(*,1000)"ite=",it," - res=",res
     endif
  enddo
1000 format(A,I2,A,F6.3)

!  call check_solution(lev)

  stop
  call mg_mpi_finalize()

end program mg_testcoarsening
