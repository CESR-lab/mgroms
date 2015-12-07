program mg_testrelax

  use mg_mpi 
  use mg_tictoc
  use mg_grids
  use mg_define_rhs
  use mg_define_matrix
  use mg_relax

  implicit none

  integer(kind=is):: nxg    ! global x dimension
  integer(kind=is):: nyg    ! global y dimension
  integer(kind=is):: nzg    ! z dimension
  integer(kind=is):: npxg   ! number of processes in x
  integer(kind=is):: npyg   ! number of processes in y
  integer(kind=is):: it     ! iteration loop number
  integer(kind=is):: nit    ! number of iterations
  integer(kind=is) :: nx, ny, nz ! local dimensions

  integer(kind=is):: nsweeps

  integer(kind=is):: lev,ierr, np
  real(kind=8)    :: res

  call tic(1,'mg_testrelax')

  !---------------!
  !- Ocean model -!
  !---------------!
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

  !-------------------!
  !- Enter in nhydro -!
  !-------------------!
  call mg_mpi_init()
  call define_grids(npxg, npyg, nx, ny, nz)
  call define_neighbours()
  !!call define_rhs(nxg, nyg, npxg)
  
  grid(1)%b = 0._8

  call random_number(grid(1)%p)

  lev = 1

  call define_matrix_simple(lev)

  call compute_residual(lev,res)
    if (myrank.eq.0) write(*,*)"ite=0 - res=",res

  do it=1, nit
     call relax(lev,nsweeps)
     call compute_residual(lev,res)
     if (myrank.eq.0)then
        write(*,1000)"ite=",it," - res=",res
     endif
  enddo
1000 format(A,I5,A,F6.3)

!  call check_solution(lev)

  call mpi_finalize(ierr)

  call toc(1,'mg_testrelax')
  call print_tictoc(myrank)

end program mg_testrelax
