module mg_setup_tests

  use mg_mpi 
  use mg_tictoc
  use mg_grids
  use mg_zr_zw
  use mg_mpi_exchange
  use mg_define_matrix
  use mg_netcdf_out

  implicit none

contains

  !-------------------------------------------------------------------------     
  subroutine setup_cuc(inc)

    integer(kind=4):: npxg,npyg
    integer(kind=4):: nx,ny,nz,nh
    integer(kind=4):: is_err,nc_id,varid
    integer(kind=4):: i,j,i0,j0,pi,pj,inc
    real(kind=8), dimension(:,:), allocatable   :: dummy2d
    integer(kind=4), dimension(2)   :: starts,counts

    character*80 :: file,varname

    npxg = grid(1)%npx
    npyg = grid(1)%npy

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh

    pj = myrank/npxg
    pi = mod(myrank,npxg)

    ! grid definition
    allocate(h(0:ny+1,0:nx+1))
    allocate(dx(0:ny+1,0:nx+1))
    allocate(dy(0:ny+1,0:nx+1))
    ! dummy 2D to read from netcdf
    allocate(dummy2d(1:nx*inc,1:ny*inc))

    i0=2+pi*nx*inc
    j0=2+pj*ny*inc

    starts(1)=i0
    starts(2)=j0
    counts(1)=nx*inc
    counts(2)=ny*inc

    !!------------------------------------------------!!
    !! Here enter the directory where cuc_nhgrd.nc is !!
    !!------------------------------------------------!!
    file='/net/xncd1/local/tmp/1/grima/MGRoms/DATA/cuc_nhgrd.nc'

    is_err = nf90_open( trim(file), NF90_NOWRITE ,nc_id  )

    if (is_err /= nf90_noerr) then
       write(*,*)'Error: problem to open file: ', trim(file)
       stop
    end if

    ! --- read h ---
    varname='h'
    is_err = nf90_inq_varid(nc_id,trim(varname),varid)
    is_err = nf90_get_var(nc_id,varid,dummy2d,starts,counts)

    do i=1,nx
       do j=1,ny
          h(j,i)=dummy2d((i-1)*inc+1,(j-1)*inc+1)
       enddo
    enddo

    ! --- read pn ---
    varname='pn'
    is_err = nf90_inq_varid(nc_id,trim(varname),varid)
    is_err = nf90_get_var(nc_id,varid,dummy2d,starts,counts)

    do i=1,nx
       do j=1,ny
          dy(j,i)=1._8/dummy2d((i-1)*inc+1,(j-1)*inc+1)
       enddo
    enddo

    ! --- read pm ---
    varname='pm'
    is_err = nf90_inq_varid(nc_id,trim(varname),varid)
    is_err = nf90_get_var(nc_id,varid,dummy2d,starts,counts)

    do i=1,nx
       do j=1,ny
          dx(j,i)=1._8/dummy2d((i-1)*inc+1,(j-1)*inc+1)
       enddo
    enddo

    is_err = nf90_close(nc_id)

    ! unpleasant cooking to fill_halo (no 2D wrapper)
    grid(1)%p(1,:,:)=h
    grid(1)%p(2,:,:)=dx
    grid(1)%p(3,:,:)=dy

    call fill_halo(1,grid(1)%p)

    h     = grid(1)%p(1,:,:)
    dx    = grid(1)%p(2,:,:)
    dy    = grid(1)%p(3,:,:)

    do i=0,nx+1
       do j=0,ny+1
          dx(j,i)=max(1.,dx(j,i))
          dy(j,i)=max(1.,dy(j,i))
       enddo
    enddo

    grid(1)%p(:,:,:)=0._8

    call define_matrices(dx, dy, h)

  end subroutine setup_cuc

  !-------------------------------------------------------------------------     
  subroutine setup_seamount()

    integer(kind=4), parameter :: ip=4, rp=8
    integer(kind=ip):: npxg,npyg ! nb procs
    integer(kind=ip):: nxg, nyg  ! global dims
    integer(kind=ip):: nx,ny,nz,nh ! local dims
    integer(kind=ip):: pi, pj
    integer(kind=ip):: i,j

    real(kind=rp) :: x, y
    real(kind=rp) :: x0, y0

    npxg = grid(1)%npx
    npyg = grid(1)%npy

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh

    nxg = npxg * nx
    nyg = npyg * ny

    ! grid definition
    allocate( h(0:ny+1,0:nx+1))
    allocate(dx(0:ny+1,0:nx+1))
    allocate(dy(0:ny+1,0:nx+1))

    dx(:,:) = Lx/real(nxg,kind=rp)
    dy(:,:) = Ly/real(nyg,kind=rp)

    pj = myrank/npxg
    pi = myrank-pj*npxg

    x0 = Lx * 0.5_rp
    y0 = Ly * 0.5_rp
    do i = 0,nx+1 !!!  I need to know my global index range
       do j = 0,ny+1
          x = (real(i+(pi*nx),kind=rp)-0.5_rp) * dx(j,i)
          y = (real(j+(pj*ny),kind=rp)-0.5_rp) * dy(j,i)
          h(j,i) = Htot * (1._rp - 0.5_rp * exp(-(x-x0)**2._rp/(Lx/5._rp)**2._rp -(y-y0)**2._rp/(Ly/5._rp)**2._rp))
       enddo
    enddo

    if (netcdf_output) then
       call write_netcdf(h, vname= 'h',netcdf_file_name= 'h.nc',rank=myrank)
       call write_netcdf(dx,vname='dx',netcdf_file_name='dx.nc',rank=myrank)
       call write_netcdf(dy,vname='dy',netcdf_file_name='dy.nc',rank=myrank)
    endif

    call define_matrices(dx, dy, h)

  end subroutine setup_seamount

end module mg_setup_tests
