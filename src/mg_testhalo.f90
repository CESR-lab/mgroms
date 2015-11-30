program mg_testhalo

  use mg_mpi ! everything will come from the outside !!!

  use mg_grids
  use mg_define_rhs

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
  integer(kind=is) :: south, east, north, west
  integer(kind=is) :: southwest, southeast, northeast, northwest
  integer(kind=is) :: nx, ny, nz
  integer(kind=is) :: nh

  integer(kind=is):: lev,ierr
  real(kind=8)    :: z
  logical :: test

  real(kind=8),dimension(:,:,:), pointer:: p

  nxg   = 128
  nyg   = 128
  nzg   = 128
  nhalo = 1

  npxg  = 4
  npyg  = 2

  nit     = 10
  nsweeps = 1

  call init_mpi(nxg, nyg, nzg, npxg, npyg)

  call define_grids(nhalo,npxg,npyg)

  call define_rhs(nxg, nyg, npxg)

  lev = 1

  p => grid(lev)%p
  nx = grid(lev)%nx
  ny = grid(lev)%ny
  nz = grid(lev)%nz
  nh = grid(lev)%nh
  west      = grid(lev)%neighb(4)

  p = 1._8*myrank

  call fill_halo(lev,p)

  ! test if the west halo has been properly updated
  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  if (myrank.eq.0) write(*,*)"---------- WEST ----------"
  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  if (west.ne.MPI_PROC_NULL) then
     z = sum(p(:,1:ny,1-nh:0)-west)
     test = (int(z).eq.0)
     write(*,1000)"rank = ",myrank," / west is = ",west," / test=",test
  endif


1000 format(A,I3,A,I3,A,L)

  call mg_mpi_finalize()

end program mg_testhalo
