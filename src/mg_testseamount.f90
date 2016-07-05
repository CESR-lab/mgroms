program mg_testseamount

  use mg_mpi
  use mg_tictoc
  use mg_setup_tests
  use mg_mpi_exchange_ijk
  use nhydro

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension

  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y

  integer(kind=ip) :: nx, ny, nz  ! local dimensions

  real(kind=rp), dimension(:,:,:), allocatable, target :: u,v,w
  real(kind=rp), dimension(:,:,:), allocatable, target :: tmp_rnd
  real(kind=rp), dimension(:,:,:), pointer :: up,vp,wp
  real(kind=rp) :: Lx, Ly, Htot

  integer(kind=ip) :: np, ierr, rank

  real(kind=rp), dimension(:,:), pointer :: dx, dy, h
  real(kind=rp) :: hc, theta_b, theta_s

  integer(kind=ip) :: pi, pj
  integer(kind=ip) :: ib, ie, jb, je, kb, ke

  call tic(1,'mg_bench_seamount')

  ! global domain dimensions
  nxg   = 64
  nyg   = 64
  nzg   = 64

  npxg  = 2
  npyg  = 2

  Lx   =  1.e4_rp
  Ly   =  1.e4_rp
  Htot =  4.e3_rp 

  hc      = 0._8
  theta_b = 0._8
  theta_s = 0._8

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

  allocate( h(0:ny+1,0:nx+1))
  allocate(dx(0:ny+1,0:nx+1))
  allocate(dy(0:ny+1,0:nx+1))

  call setup_seamount(nx, ny, nz, npxg, npyg, Lx, Ly, Htot, dx, dy, h)

  !---------------------!
  !- Initialise nhydro -!
  !---------------------!
  if (rank == 0) write(*,*)'Initialise NHydro (grids, cA, params, etc) '
  call nhydro_init(nx, ny, nz, npxg, npyg, dx, dy, h, hc, theta_b, theta_s, test='seamount')

  !-------------------------------------!
  !- U,V,W initialisation (model vars) -!
  !-------------------------------------!
  if (rank == 0) write(*,*)'u, v, w initialisation...'

  allocate(u(1:nx+1,0:ny+1,1:nz))
  allocate(v(0:nx+1,1:ny+1,1:nz))
  allocate(w(0:nx+1,0:ny+1,0:nz))

  pj = myrank/npxg
  pi = mod(myrank,npxg)

  ib = 1 + pi * nx
  jb = 1 + pj * ny

  ie = ib + nx - 1
  je = jb + ny - 1

  kb = 1
  ke = nz

  allocate(tmp_rnd(1:nxg,1:nyg,1:nzg))

  call random_number(tmp_rnd)
  tmp_rnd = 2._8 * tmp_rnd - 1._8
  u(1:nx,1:ny,1:nz) = tmp_rnd(ib:ie,jb:je,kb:ke)
  up => u
  call fill_halo_ijk(nx,ny,up,'u') ! depend of mg_grids for MPI neighbours !

  call random_number(tmp_rnd)
  tmp_rnd = 2._8 * tmp_rnd - 1._8
  v(1:nx,1:ny,1:nz) = tmp_rnd(ib:ie,jb:je,kb:ke)
  vp => v
  call fill_halo_ijk(nx,ny,vp,'v') ! depend of mg_grids for MPI neighbours !

  deallocate(tmp_rnd)
  allocate(tmp_rnd(1:nxg,1:nyg,0:nzg))

  kb = 0
  ke = nz

  call random_number(tmp_rnd)
  tmp_rnd = 2._8 * tmp_rnd - 1._8
  w(1:nx,1:ny,0:nz) = tmp_rnd(ib:ie,jb:je,kb:ke)
  wp => w
  call fill_halo_ijk(nx,ny,wp,'w') ! depend of mg_grids for MPI neighbours !

  deallocate(tmp_rnd)

  if (netcdf_output) then
     call write_netcdf(up,vname='u',netcdf_file_name='u.nc',rank=myrank,iter=0)
     call write_netcdf(vp,vname='v',netcdf_file_name='v.nc',rank=myrank,iter=0)
     call write_netcdf(wp,vname='w',netcdf_file_name='w.nc',rank=myrank,iter=0)
  endif

  !----------------------!
  !- Call nhydro solver -!
  !----------------------!
  if (rank == 0) write(*,*)'Call nhydro solver'
  call nhydro_solve(nx,ny,nz,u,v,w)

  if (netcdf_output) then
     call write_netcdf(up,vname='u',netcdf_file_name='u.nc',rank=myrank,iter=2)
     call write_netcdf(vp,vname='v',netcdf_file_name='v.nc',rank=myrank,iter=2)
     call write_netcdf(wp,vname='w',netcdf_file_name='w.nc',rank=myrank,iter=2)
  endif

  !------------------------------------------------------------!
  !- Call nhydro correct to check if nh correction is correct -!
  !------------------------------------------------------------!
  if (rank == 0) write(*,*)'Call nhydro correct'
  call check_correction(nx,ny,nz,u,v,w)

  if (netcdf_output) then
     call write_netcdf(up,vname='u',netcdf_file_name='u.nc',rank=myrank,iter=3)
     call write_netcdf(vp,vname='v',netcdf_file_name='v.nc',rank=myrank,iter=3)
     call write_netcdf(wp,vname='w',netcdf_file_name='w.nc',rank=myrank,iter=3)
  endif

  if (netcdf_output) then
     call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='check_correction.nc',rank=myrank)
  endif

  !---------------------!
  !- Deallocate memory -!
  !---------------------!
  if (rank == 0) write(*,*)'Cleaning memory before to finish the program.'
  call nhydro_clean()

  !----------------------!
  !- End Bench-seamount -!
  !----------------------!
  call mpi_finalize(ierr)

  call toc(1,'mg_bench_seamount')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testseamount

