program mg_testcuc

  use mg_mpi 
  use mg_tictoc
  use mg_setup_tests
  use nhydro

  implicit none

  integer(kind=4):: nxg        ! global x dimension
  integer(kind=4):: nyg        ! global y dimension
  integer(kind=4):: nzg        ! z dimension
  integer(kind=4):: npxg       ! number of processes in x
  integer(kind=4):: npyg       ! number of processes in y
  integer(kind=4):: nx, ny, nz ! local dimensions

  integer(kind=4):: ierr, np, rank

  real(kind=8), dimension(:,:,:), pointer :: u,v,w

  real(kind=rp) :: Lx, Ly, Htot
  real(kind=rp) :: hc, theta_b, theta_s
  real(kind=rp), dimension(:,:), pointer :: dx, dy, h

  call tic(1,'mg_testcuc')

  !---------------!
  !- Ocean model -!
  !---------------!
  nxg  = 1024
  nyg  = 1024
  nzg  =   64

  Lx   = 200d3
  Ly   = 200d3
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


 !-------------------------------------------------!
  !- dx,dy,h and U,V,W initialisation (model vars) -!
  !-------------------------------------------------!
  allocate( h(0:ny+1,0:nx+1))
  allocate(dx(0:ny+1,0:nx+1))
  allocate(dy(0:ny+1,0:nx+1))

  call setup_cuc(nx, ny, nz, npxg, npyg, dx, dy, h)

  !---------------------!
  !- Initialise nhydro -!
  !---------------------!
  if (rank == 0) write(*,*)'Initialise NHydro (grids, cA, params, etc) '
  call nhydro_init(&
       nx, ny, nz, &
       npxg, npyg, &
       dx, dy, h,  &
       hc, theta_b, theta_s, &
       test='cuc' )

  !-------------------------------------!
  !- U,V,W initialisation (model vars) -!
  !-------------------------------------!
  allocate(u(0:nx+1,0:ny+1,  nz))
  allocate(v(0:nx+1,0:ny+1,  nz))
  allocate(w(0:nx+1,0:ny+1,0:nz))

  u(:,:,:)      =  0._8
  v(:,:,:)      =  0._8

  w(:,:,0)      =  0._8
  w(:,:,1:nz-1) = -1._8
  w(:,:,nz)     =  0._8

  !----------------------!
  !- Call nhydro solver -!
  !----------------------!
  if (rank == 0) write(*,*)'Call nhydro solver'
  call  nhydro_solve(nx,ny,nz,u,v,w)

  !------------------------------------------------------------!
  !- Call nhydro correct to check if nh correction is correct -!
  !------------------------------------------------------------!
  if (rank == 0) write(*,*)'Call nhydro correct'
  call check_correction(nx,ny,nz,u,v,w)

  if (netcdf_output) then
     call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='check_correction.nc',rank=myrank)
  endif

  !------------------!
  !- End test-model -!
  !------------------!
  call mpi_finalize(ierr)

  call toc(1,'mg_testcuc')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testcuc
