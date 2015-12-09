program mg_testrelax

  use mg_mpi 
  use mg_tictoc
  use mg_grids
  use mg_define_rhs
  use mg_define_matrix
  use mg_relax

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension
  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y
  integer(kind=ip):: it     ! iteration loop number
  integer(kind=ip):: nit    ! number of iterations
  integer(kind=ip) :: nx, ny, nz ! local dimensions

  integer(kind=ip):: nsweeps

  integer(kind=ip):: lev,ierr, np, nh
  real(kind=rp)    :: res

  real(kind=rp),dimension(:,:,:),allocatable  :: p0

  call tic(1,'mg_testrelax')

  !---------------!
  !- Ocean model -!
  !---------------!
  nxg   = 128
  nyg   = 128
  nzg   = 128

  npxg  = 2
  npyg  = 2

  nit     = 100
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
  
  
  nh = grid(1)%nh
  allocate(p0(nz,1-nh:ny+nh,1-nh:nx+nh))
  p0 = 0._8

  grid(1)%b = 0._8

  call random_number(grid(1)%p)

  call random_number(grid(1)%b)
  call fill_halo(1,grid(1)%b)

  lev = 1

  call define_matrix_simple(lev)

  call compute_residual(lev,res)
    if (myrank.eq.0) write(*,1000)"ite=",0," - res=",res

  do it=1, nit
!     grid(1)%p=0._8
     call relax(lev,nsweeps)
!     p0 = p0 + grid(1)%p
!     grid(1)%p = p0
     call compute_residual(lev,res)
     if (myrank.eq.0)then
        write(*,1000)"ite=",it," - res=",res
     endif
  enddo
1000 format(A,I5,A,F8.3)

!  call check_solution(lev)

  call mpi_finalize(ierr)

  call toc(1,'mg_testrelax')
  call print_tictoc(myrank)

end program mg_testrelax
