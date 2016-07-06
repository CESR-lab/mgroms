module nhydro

  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_mpi_exchange
  use mg_netcdf_out
  use mg_compute_rhs
  use mg_correct_uvw
  use mg_solvers


  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(  &
       nx, ny, nz,         &
       npxg, npyg,         &
       dx,dy,h,            &
       hc,theta_b,theta_s, &
       test)

    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg
    real(kind=rp), dimension(:,:), intent(inout) :: dx, dy, h
    real(kind=rp),  intent(in) :: hc, theta_b, theta_s
    character(len=*), optional :: test

    call mg_mpi_init()

    hlim      = hc
    nhtheta_b = theta_b
    nhtheta_s = theta_s

    !- read the NonHydro namelist file if it is present 
    !- else default values and print them (or not).
    call read_nhnamelist(vbrank=myrank)

    call define_grids(npxg, npyg, nx, ny, nz)

    call define_neighbours()

    call print_grids()

    if (present(test)) then
       if (trim(test)=='seamount') then
          bench = 'seamount'
       endif
    endif

    call define_matrices(dx, dy, h)

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz
    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz), target, intent(inout) :: ua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz), target, intent(inout) :: va
    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), target, intent(inout) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u, v, w

    real(kind=rp)    :: tol
    integer(kind=ip) :: maxite

    tol    = solver_prec    ! solver_prec    is defined in the namelist file
    maxite = solver_maxiter ! solver_maxiter is defined in the namelist file

    u => ua
    v => va
    w => wa

    !- Step 1 - 
    call compute_rhs(u, v, w)

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='b.nc',rank=myrank,iter=1)
    endif

    if (netcdf_output) then
       call write_netcdf(u,vname='u',netcdf_file_name='u.nc',rank=myrank,iter=1)
       call write_netcdf(v,vname='v',netcdf_file_name='v.nc',rank=myrank,iter=1)
       call write_netcdf(w,vname='w',netcdf_file_name='w.nc',rank=myrank,iter=1)
    endif

    grid(1)%p(:,:,:) = 0._rp

    !- step 2 -
    call solve(tol,maxite)

    if (netcdf_output) then
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='p_end.nc',rank=myrank)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='r_end.nc',rank=myrank)
    endif

    !- Step 3 -
    call correct_uvw(u,v,w)

  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_clean()

    call grids_dealloc()

  end subroutine nhydro_clean

end module nhydro
