program mg_testhalo

  use mg_mpi
  use mg_grids
  use mg_define_rhs

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension

  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y
  integer(kind=ip):: it     ! iteration loop number
  integer(kind=ip):: nit    ! number of iterations

  integer(kind=ip):: nsweeps
  integer(kind=ip) :: south, east, north, west
  integer(kind=ip) :: southwest, southeast, northeast, northwest
  integer(kind=ip) :: nx, ny, nz ! local dimensions
  integer(kind=ip) :: nh

  integer(kind=ip):: lev,ierr, np, rank
  real(kind=8)    :: z
  logical :: test

  real(kind=8),dimension(:,:,:), pointer:: p

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
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  if (np /= (npxg*npyg)) then
     write(*,*) "Error: in number of processes !"
     stop -1
  endif

  nx = nxg / npxg
  ny = nyg / npyg
  nz = nzg

  !- read the NonHydro namelist file if it is present 
  !- else default values and print them (or not).
  call read_nhnamelist(vbrank=rank)

  !-------------------!
  !- Enter in nhydro -!
  !-------------------!
  call mg_mpi_init()
  call define_grids(npxg, npyg, nx, ny, nz)
  call define_neighbours()
  call define_rhs(nxg, nyg, npxg)

  do lev = 1,nlevs

     p    => grid(lev)%p
     nx   = grid(lev)%nx
     ny   = grid(lev)%ny
     nz   = grid(lev)%nz
     nh   = grid(lev)%nh
     west = grid(lev)%neighb(4)

     p = 1._8*myrank

     call fill_halo(lev,p)

     ! test if the west halo has been properly updated
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)
     if (myrank.eq.0) write(*,"(A,I2,A)")"---------- WEST lev = ",lev," ----------"
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)

     if (west.ne.MPI_PROC_NULL) then
        z = sum(p(:,1:ny,1-nh:0)-west)
        test = (int(z).eq.0)
        write(*,1000)"rank = ",myrank," / west is = ",west," / test=",test
     else
        write(*,1000)"rank = ",myrank," no west neighbour", west
     endif
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)

  enddo

1000 format(A,I3,A,I3,A,L)

  call mpi_finalize(ierr)

end program mg_testhalo
