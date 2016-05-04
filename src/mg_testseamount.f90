program mg_testseamount

  use mg_mpi
  use mg_tictoc
  use mg_define_rhs
  use nhydro

  implicit none

  integer(kind=4):: nxg    ! global x dimension
  integer(kind=4):: nyg    ! global y dimension
  integer(kind=4):: nzg    ! z dimension

  integer(kind=4):: npxg   ! number of processes in x
  integer(kind=4):: npyg   ! number of processes in y

  integer(kind=4) :: nx, ny, nz  ! local dimensions

  real(kind=8), dimension(:,:,:), pointer :: u,v,w

  integer(kind=4) :: np, ierr, rank

  call tic(1,'mg_bench_seamount')

  ! global domain dimensions
  nxg   = 128
  nyg   = 128
  nzg   = 128

  Lx   =  32d3
  Ly   =  32d3
  Htot =   4d2

  npxg  = 8
  npyg  = 9

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

  !-------------------------------------!
  !- U,V,W initialisation (model vars) -!
  !-------------------------------------!
  allocate(u(nz,ny,1:nx+1))
  allocate(v(nz,1:ny+1,nx))
  allocate(w(1:nz+1,ny,nx))

  u(:,:,:)    =  0._8
  v(:,:,:)    =  0._8
  w(2:nz,:,:) = -1._8
  w(1,:,:)    =  0._8
  w(nz+1,:,:) =  0._8

  !---------------------!
  !- Initialise nhydro -!
  !---------------------!
  if (rank == 0) write(*,*)'Initialise NHydro (grids, cA, params, etc) '
  call nhydro_init(nx, ny, nz, npxg, npyg, test='seamount')

  !--------------------------!
  !- RHS initialisation (b) -!
  !--------------------------!
  if (rank == 0) write(*,*)'RHS initialisation'
  call rhs_seamount()

  !----------------------!
  !- Call nhydro solver -!
  !----------------------!
  if (rank == 0) write(*,*)'Call nhydro solver'
  call nhydro_solve(u,v,w)

  !----------------------!
  !- End Bench-seamount -!
  !----------------------!
  call mpi_finalize(ierr)

  call toc(1,'mg_bench_seamount')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testseamount

