module nhydro

  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_setup_tests
  use mg_solvers
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(nx, ny, nz, npxg, npyg, test)

    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg
    character(len=*), optional :: test

    integer(kind=ip) :: inc=1

    call mg_mpi_init()

    !- read the NonHydro namelist file if it is present 
    !- else default values and print them (or not).
    call read_nhnamelist(vbrank=myrank)

    call define_grids(npxg, npyg, nx, ny, nz)

    call define_neighbours()

    call print_grids()

    if (present(test)) then

       if (trim(test)=='cuc') then
          call setup_cuc(inc)                  ! call define matrices !
       elseif (trim(test)=='seamount') then
          bench = 'seamount'
          call setup_seamount()                ! call define matrices !
       else
          write(*,*)'Error: nhydro: please enter a valid test name!', trim(test)
          stop
       endif

    endif

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    real(kind=rp)    :: tol = 1.e-12
    integer(kind=ip) :: maxite = 50

    grid(1)%p(:,:,:) = 0._rp

    call solve(tol,maxite)

    if (netcdf_output) then
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='p_end.nc',rank=myrank)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='r_end.nc',rank=myrank)
    endif

  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_clean()

    call grids_dealloc()

  end subroutine nhydro_clean

end module nhydro
