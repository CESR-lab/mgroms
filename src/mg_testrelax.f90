program mg_testrelax

  use mg_mpi ! everything will come from the outside !!!

  use mg_grids
  use mg_define_rhs
  use mg_define_matrix
  use mg_relax

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

  integer(kind=is):: lev

  nxg   = 128
  nyg   = 128
  nzg   = 128
  nhalo = 1
  npxg  = 2
  npyg  = 2

  nit     = 1
  nsweeps = 2

  call init_mpi(nxg, nyg, nzg, npxg, npyg)

  call define_grids(nhalo,npxg,npyg)

  call define_rhs(nxg, nyg, nzg, npxg, npyg)

  call define_matrix_simple()

  lev = 1
  call check_solution(lev)

  do it=1, nit
     call relax_line(lev,nsweeps)
     call check_solution(lev)
  enddo

!  call check_solution(lev)

  call mg_mpi_finalize()

end program mg_testrelax
