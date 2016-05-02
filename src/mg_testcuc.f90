program mg_testcuc

  use mg_mpi 
  use mg_tictoc
  use mg_define_rhs
  use nhydro

  implicit none

  integer(kind=4):: nxg        ! global x dimension
  integer(kind=4):: nyg        ! global y dimension
  integer(kind=4):: nzg        ! z dimension
  integer(kind=4):: npxg       ! number of processes in x
  integer(kind=4):: npyg       ! number of processes in y
  integer(kind=4):: nx, ny, nz ! local dimensions

  integer(kind=4):: ierr, np, rank,inc

  real(kind=8), dimension(:,:,:), pointer :: u,v,w

  call tic(1,'mg_testcuc')

  !---------------!
  !- Ocean model -!
  !---------------!
  inc  = 1
  nxg  = 1024/inc
  nyg  = 1024/inc
  nzg  =   64/inc

  Lx   =  200d3
  Ly   =  200d3
  Htot = 4d3

  ! global variables define in mg_grids !?!
  ! Should be nhydro arguments !
  hlim    = 250._8
  theta_b =   6._8
  theta_s =   6._8

  npxg  = 2
  npyg  = 2

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
  call nhydro_init(nx, ny, nz, npxg, npyg, test='cuc')

  !--------------------------!
  !- RHS initialisation (b) -!
  !--------------------------!
  if (rank == 0) write(*,*)'RHS initialisation'
  !! call setup_random_patches()
  call rhs_random()
  
  !----------------------!
  !- Call nhydro solver -!
  !----------------------!
  if (rank == 0) write(*,*)'Call nhydro solver'
  call  nhydro_solve(u,v,w)

  !------------------!
  !- End test-model -!
  !------------------!
  call mpi_finalize(ierr)

  call toc(1,'mg_testcuc')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testcuc
