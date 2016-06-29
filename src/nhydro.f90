module nhydro

  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_compute_rhs
  use mg_correct_uvw
  use mg_solvers
  use mg_mpi_exchange
  use mg_netcdf_out

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
    real(kind=rp), dimension(:,:), intent(in) :: dx, dy, h
    real(kind=rp),  intent(in) :: hc, theta_b, theta_s
    character(len=*), optional :: test

    integer(kind=ip) :: inc=1

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
!!$    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz), target, intent(inout) :: ua
!!$    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz), target, intent(inout) :: va
!!$    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), target, intent(inout) :: wa

    real(kind=rp), dimension(nz  ,0:ny+1,0:nx+1), target, intent(inout) :: ua
    real(kind=rp), dimension(nz  ,0:ny+1,0:nx+1), target, intent(inout) :: va
    real(kind=rp), dimension(nz+1,0:ny+1,0:nx+1), target, intent(inout) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u, v, w

    real(kind=rp)    :: tol = 1.e-12
    integer(kind=ip) :: maxite = 50

    u => ua
    v => va
    w => wa

    write(*,*)'nhydro -> Lbound(v):',lbound(v)

    !- Step 1 - 
    call compute_rhs(u, v, w)

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='first_rhs.nc',rank=myrank)
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
