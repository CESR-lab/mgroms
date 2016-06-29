program mg_testseamount

  use mg_mpi
  use mg_tictoc
  use mg_define_rhs
  use nhydro

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension

  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y

  integer(kind=ip) :: nx, ny, nz  ! local dimensions

  real(kind=rp), dimension(:,:,:), pointer :: u,v,w

  integer(kind=ip) :: np, ierr, rank

  call tic(1,'mg_bench_seamount')

  ! global domain dimensions
  nxg   = 64
  nyg   = 64
  nzg   = 64

  npxg  = 1
  npyg  = 1

  Lx   =  1.e4_rp
  Ly   =  1.e4_rp
  Htot =  4.e3_rp 

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

  !---------------------!
  !- Initialise nhydro -!
  !---------------------!
  if (rank == 0) write(*,*)'Initialise NHydro (grids, cA, params, etc) '
  call nhydro_init(nx, ny, nz, npxg, npyg, test='seamount')

  !-------------------------------------!
  !- U,V,W initialisation (model vars) -!
  !-------------------------------------!

!!$  u(:,:,:)    =  0._8
!!$  v(:,:,:)    =  0._8
!!$  w(2:nz,:,:) = -1._8
!!$  w(1,:,:)    =  0._8
!!$  w(nz+1,:,:) =  0._8
!!$  !--------------------------!
!!$  !- RHS initialisation (b) -!
!!$  !--------------------------!
!!$  if (rank == 0) write(*,*)'RHS initialisation'
!!$  call rhs_seamount()

  allocate(u(nz  ,0:ny+1,0:nx+1))
  allocate(v(nz  ,0:ny+1,0:nx+1))
  allocate(w(nz+1,0:ny+1,0:nx+1))

  call random_number(u)
  u = 2._8 * u - 1._8
  call fill_halo(1,u)

  call random_number(v)
  v = 2._8 * v - 1._8
  call fill_halo(1,v)

  call random_number(w)
  w = 2._8 * w - 1._8
  call fill_halo(1,w)

  !----------------------!
  !- Call nhydro solver -!
  !----------------------!
  if (rank == 0) write(*,*)'Call nhydro solver'
  call nhydro_solve(u,v,w)

  !------------------------------------------------------------!
  !- Call nhydro correct to check if nh correction is correct -!
  !------------------------------------------------------------!
  if (rank == 0) write(*,*)'Call nhydro correct'
  call nhydro_correct(u,v,w)

  !----------------------!
  !- End Bench-seamount -!
  !----------------------!
  call mpi_finalize(ierr)

  call toc(1,'mg_bench_seamount')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testseamount

