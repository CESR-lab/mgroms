module nhydro

  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_tictoc
  use mg_mpi_exchange
  use mg_netcdf_out
  use mg_compute_rhs
  use mg_correct_uvw
  use mg_solvers

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(nx,ny,nz,npxg,npyg)
      
    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg

    call mg_mpi_init()

    call read_nhnamelist(vbrank=myrank)

    call define_grids(npxg,npyg,nx,ny,nz)

    call define_neighbours()

    call print_grids()

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_matrices(dx,dy,zeta,h,hc,theta_b,theta_s)

    real(kind=rp), dimension(:,:), intent(in) :: dx, dy, zeta, h
    real(kind=rp),                 intent(in) :: hc, theta_b, theta_s

    if (myrank==0) write(*,*)' nhydro_matrices:'

    hlim      = hc
    nhtheta_b = theta_b
    nhtheta_s = theta_s

    call define_matrices(dx,dy,zeta,h)

  end subroutine nhydro_matrices

  !--------------------------------------------------------------
  subroutine nhydro_solve(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz
    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz), target, intent(inout) :: ua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz), target, intent(inout) :: va
    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), target, intent(inout) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u, v, w

    real(kind=rp)    :: tol
    integer(kind=ip) :: maxite

    call tic(1,'nhydro_solve')

    tol    = solver_prec    ! solver_prec    is defined in the namelist file
    maxite = solver_maxiter ! solver_maxiter is defined in the namelist file

    u => ua
    v => va
    w => wa

    if (myrank==0) write(*,*)' nhydro_solve:'

    !- step 1 - 
    call tic(1,'compute_rhs')
    call compute_rhs(u, v, w)
    call toc(1,'compute_rhs')

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='b.nc',rank=myrank,iter=1)
    endif

    !- step 2 -
    call solve_p(tol,maxite)

    if (netcdf_output) then
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='p_end.nc',rank=myrank)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='r_end.nc',rank=myrank)
    endif

    !- step 3 -
    call correct_uvw(u,v,w)

    call toc(1,'nhydro_solve')

  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_check_nondivergence(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz), target, intent(inout) :: ua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz), target, intent(inout) :: va
    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), target, intent(inout) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u, v, w

    integer(kind=ip), save :: iter=0
    iter = iter + 1

    u => ua
    v => va
    w => wa

    if (myrank==0) write(*,*)'- check non-divergence:'

    call compute_rhs(u,v,w)

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='check.nc',rank=myrank,iter=iter)
    endif

  end subroutine nhydro_check_nondivergence

  !--------------------------------------------------------------
  subroutine nhydro_clean()

    call grids_dealloc()

  end subroutine nhydro_clean

end module nhydro
