program mg_testnetcdf

  use mg_mpi 
  use mg_tictoc
  use mg_netcdf_out
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

  call tic(1,'mg_testnetcdf')

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
  
  
  nh = grid(1)%nh
  allocate(p0(nz,1-nh:ny+nh,1-nh:nx+nh))
  p0 = 0._8

  grid(1)%b = 0._8

  call random_number(grid(1)%p)

  call random_number(grid(1)%b)
  call fill_halo(1,grid(1)%b)

  lev = 1

  call define_matrix_simple(lev)

  !-------------------------------------------------------
  !- Writing a 4D array in a netcdf file
  call write_netcdf(grid(lev)%cA,vname='cA1',rank=myrank)
  !-------------------------------------------------------

  call compute_residual(lev,res)
    if (myrank.eq.0) write(*,1000)"ite=",0," - res=",res

  do it=1, nit
     call relax(lev,nsweeps)
     call compute_residual(lev,res)

     !-------------------------------------------------------------------------------
     !- Writing a 3D array in a netcdf file
     call write_netcdf(grid(lev)%r,vname='r',rank=myrank,iter=it)
     !- NCO commands to sublit after the run !
     !- ncecat netcdf_file_r_000_*.nc netcdf_file_r_000.nc; \rm netcdf_file_r_000_*.nc
     !- ncecat netcdf_file_r_001_*.nc netcdf_file_r_001.nc; \rm netcdf_file_r_001_*.nc
     !- ncecat netcdf_file_r_002_*.nc netcdf_file_r_002.nc; \rm netcdf_file_r_002_*.nc
     !- ncecat netcdf_file_r_003_*.nc netcdf_file_r_003.nc; \rm netcdf_file_r_003_*.nc
     !--------------------------------------------------------------------------------

     if (myrank.eq.0)then
        write(*,1000)"ite=",it," - res=",res
     endif
  enddo
1000 format(A,I5,A,F8.3)

!  call check_solution(lev)

  call mpi_finalize(ierr)

  call toc(1,'mg_testnetcdf')
  call print_tictoc(myrank)

end program mg_testnetcdf
