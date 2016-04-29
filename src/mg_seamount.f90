module mg_seamount

  use mg_mpi
  use mg_grids
  use mg_mpi_exchange
  use mg_define_matrix
  use mg_intergrids
  use mg_relax
  use mg_s_coord
!  use mg_netcdf_in

  implicit none

  real(kind=8), dimension(:,:), pointer   :: dx, dy
  real(kind=8), dimension(:,:), pointer   :: umask, vmask, rmask
  real(kind=8), dimension(:,:), pointer   :: h
  real(kind=8), dimension(:,:,:), pointer :: zr, zw

  ! MUST BE DEFINED BY THE USER OUTSIDE OF THIS MODULE
  ! before calling setup_realistic_matrix()

  real(kind=8) :: Lx,Ly,Htot
  !!NG real(kind=8) :: hc,theta_s,theta_b

contains

  !-------------------------------------------------------------------------     
  subroutine setup_realistic_matrix()

    integer(kind=4):: npxg,npyg
    integer(kind=4):: nx,ny,nz,nh

    npxg = grid(1)%npx
    npyg = grid(1)%npy

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh

    ! grid definition
    allocate(h(0:ny+1,0:nx+1))
    allocate(dx(0:ny+1,0:nx+1))
    allocate(dy(0:ny+1,0:nx+1))
    allocate(umask(0:ny+1,0:nx+1))
    allocate(vmask(0:ny+1,0:nx+1))
    allocate(rmask(0:ny+1,0:nx+1))
    allocate(zr(nz,0:ny+1,0:nx+1))
    allocate(zw(nz+1,0:ny+1,0:nx+1))

    dx(:,:) = Lx/real(nx*npxg,kind=8)
    dy(:,:) = Ly/real(ny*npyg,kind=8)


    call setup_fine_mask(nx,ny,nz,nh)

    call setup_fine_depth(nx,ny,nz,nh)

    !call setup_scoord(nx,ny,nz,nh)
    call setup_scoord_gene           ( &  ! Compuet zr and zw
            hlim,theta_b,theta_s,h   , &  ! input args
            zr,zw                    , &  ! output args
            coord_type='new_s_coord' )  ! optional

    call define_matrices(dx, dy, zr, zw)

  end subroutine setup_realistic_matrix

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
    allocate(umask(0:ny+1,0:nx+1))
    allocate(vmask(0:ny+1,0:nx+1))
    allocate(rmask(0:ny+1,0:nx+1))
    allocate(zr(nz,0:ny+1,0:nx+1))
    allocate(zw(nz+1,0:ny+1,0:nx+1))
    ! dummy 2D to read from netcdf
    allocate(dummy2d(1:nx*inc,1:ny*inc))

    i0=2+pi*nx*inc
    j0=2+pj*ny*inc

    starts(1)=i0
    starts(2)=j0
    counts(1)=nx*inc
    counts(2)=ny*inc

    file='./cuc_nhgrd.nc'

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

    ! --- read pm ---
    varname='mask_rho'
    is_err = nf90_inq_varid(nc_id,trim(varname),varid)
    is_err = nf90_get_var(nc_id,varid,dummy2d,starts,counts)

    do i=1,nx
       do j=1,ny
          if(h(j,i)<=21._8)then
             rmask(j,i)=0._8
          else
             rmask(j,i)=1._8
          endif
          !             rmask(j,i)=dummy2d(i,j)
       enddo
    enddo

    is_err = nf90_close(nc_id)

!!$    if(pi==0)rmask(:,0:4)=0._8
!!$    if(pi==npxg-1)rmask(:,nx-3:nx+1)=0._8
!!$    if(pj==0)rmask(0:4,:)=0._8
!!$    if(pj==npyg-1)rmask(ny-3:ny+1,:)=0._8

    ! unpleasant cooking to fill_halo umask and vmask (no 2D wrapper)
    grid(1)%p(1,:,:)=h
    grid(1)%p(2,:,:)=dx
    grid(1)%p(3,:,:)=dy
    grid(1)%p(4,:,:)=rmask

    call fill_halo(1,grid(1)%p)

    h     = grid(1)%p(1,:,:)
    dx    = grid(1)%p(2,:,:)
    dy    = grid(1)%p(3,:,:)
    rmask = grid(1)%p(4,:,:)

    do i=0,nx+1
       do j=0,ny+1
          dx(j,i)=max(1.,dx(j,i))
          dy(j,i)=max(1.,dy(j,i))
       enddo
    enddo

    !    dx = 146.4772*inc
    !    dy = dx

    !    if(pi==0)rmask(:,0:2)=0._8
    !    if(pi==npxg-1)rmask(:,nx-1:nx+1)=0._8
    !    if(pj==0)rmask(0:2,:)=0._8
    !    if(pj==npyg-1)rmask(ny-1:ny+1,:)=0._8

    do i=1,nx     
       do j=1,ny
          umask(j,i)=rmask(j,i)*rmask(j,i-1)
          vmask(j,i)=rmask(j,i)*rmask(j-1,i)
       enddo
    enddo

    ! unpleasant cooking to fill_halo umask and vmask (no 2D wrapper)
    grid(1)%p(1,:,:)=umask
    grid(1)%p(2,:,:)=vmask
    call fill_halo(1,grid(1)%p)
    umask = grid(1)%p(1,:,:)
    vmask = grid(1)%p(2,:,:)

    grid(1)%p(:,:,:)=0._8

    if (trim(vgrid) == 'topo') then
       call define_matrices(dx, dy, h)
    else
       !call setup_scoord(nx,ny,nz,nh)
       call setup_scoord_gene        ( &  ! Compuet zr and zw
            hlim,theta_b,theta_s,h   , &  ! input args
            zr,zw                    , &  ! output args
            coord_type='new_s_coord' )  ! optional

       call define_matrices(dx, dy, zr, zw)
    endif

    if (netcdf_output) then
       !     call write_netcdf(h,vname='h',netcdf_file_name='h.nc',rank=myrank)
       !     call write_netcdf(dx,vname='dx',netcdf_file_name='dx.nc',rank=myrank)
       call write_netcdf(rmask,vname='rmask',netcdf_file_name='rmask.nc',rank=myrank)
       call write_netcdf(umask,vname='umask',netcdf_file_name='umask.nc',rank=myrank)
       call write_netcdf(    h,vname=    'h',netcdf_file_name=    'h.nc',rank=myrank)
       call write_netcdf(   dx,vname=   'dx',netcdf_file_name=   'dx.nc',rank=myrank)
       call write_netcdf(   dy,vname=   'dy',netcdf_file_name=   'dy.nc',rank=myrank)
    endif

  end subroutine setup_cuc
  !-------------------------------------------------------------------------     
  subroutine setup_fine_depth(nx,ny,nz,nh)
    integer(kind=4):: nx,ny,nz,nh

    real(kind=8) :: x,y,x0,y0,dist
    integer(kind=4):: i,j,pi,pj,npxg,npyg

    npxg = grid(1)%npx
    npyg = grid(1)%npy
    pj = myrank/npxg
    pi = mod(myrank,npxg)

    x0 = Lx * 0.5_8
    y0 = Ly * 0.5_8
    do i = 0,nx+1 !!!  I need to know my global index range
       do j = 0,ny+1
          x = (real(i+(pi*nx),kind=8)-0.5_8) * dx(j,i)
          y = (real(j+(pj*ny),kind=8)-0.5_8) * dy(j,i)

          h(j,i) = Htot
          !h(j,i) = Htot * (1._8 - 0.5_8 * exp(-(x-x0)**2._8/(Lx/5._8)**2._8))
          h(j,i) = Htot * (1._8 - 0.5_8 * exp(-(x-x0)**2._8/(Lx/5._8)**2._8 -(y-y0)**2._8/(Ly/5._8)**2._8))
          ! circular slope
          dist = sqrt((x-x0)**2._8 + (y-y0)**2._8) / Lx
          !h(j,i) = Htot - Htot* ( (tanh( (dist -0.45)*10. )+1.)*0.5 *0.9 )

       enddo
    enddo

    if (netcdf_output) then
       call write_netcdf(h,vname='h',netcdf_file_name='h.nc',rank=myrank)
    endif

  end subroutine setup_fine_depth

  !-------------------------------------------------------------------------     
  subroutine setup_fine_mask(nx,ny,nz,nh)
    integer(kind=4):: nx,ny,nz,nh

    real(kind=8) :: dist,x,y
    integer(kind=4):: i,j,pi,pj,npxg,npyg

    npxg = grid(1)%npx
    npyg = grid(1)%npy
    pj = myrank/npxg
    pi = mod(myrank,npxg)

    rmask(:,:) = 0.
    do i=0,nx+1
       x=(1._8*i-0.5+pi*nx)/(npxg*nx) -0.5
       do j=0,ny+1
          y=(1._8*j-0.5+pj*ny)/(npyg*ny)-0.5
          dist = sqrt(x*x+y*y) - 0.48_8

          ! for a circular domain
          if (dist<0.) then

          ! for a square domain
          !if ((x>-0.5).and.(x<0.5).and.(y>-0.5).and.(y<0.5)) then
             rmask(j,i)  = 1.
          endif
       enddo
    enddo

    do i=1,nx+1     
       do j=1,ny+1
          umask(j,i)=rmask(j,i)*rmask(j,i-1)
          vmask(j,i)=rmask(j,i)*rmask(j-1,i)
       enddo
    enddo

    ! unpleasant cooking to fill_halo umask and vmask (no 2D wrapper)
    grid(1)%p(1,:,:)=umask
    grid(1)%p(2,:,:)=vmask
    call fill_halo(1,grid(1)%p)
    umask = grid(1)%p(1,:,:)
    vmask = grid(1)%p(2,:,:)
    grid(1)%p(:,:,:)=0._8

    if (netcdf_output) then
       call write_netcdf(umask,vname='umask',netcdf_file_name='umask.nc',rank=myrank)
       call write_netcdf(vmask,vname='vmask',netcdf_file_name='vmask.nc',rank=myrank)
    endif

  end subroutine setup_fine_mask

  !-------------------------------------------------------------------------     
  subroutine setup_rhs(nx,ny,nz,nh,rhs)
    real(kind=8), dimension(:,:,:), pointer,intent(inout)   ::  rhs
    integer(kind=4):: nx,ny,nz,nh

    integer(kind=4):: i,j,k,pi,pj,npxg,npyg
    real(kind=8) :: x,y,bet,x1,z1,x2,z2,y1,y2

    npxg = grid(1)%npx
    npyg = grid(1)%npy
    pj = myrank/npxg
    pi = mod(myrank,npxg)

    ! rhs definition
    bet = 1200._8 / (Lx*Lx)
    x1 = Lx * 0.45_8
    y1 = Ly * 0.45_8
    z1 = Htot * (0.75_8 - 1._8)
    x2 = Lx * 0.75_8
    y2 = Ly * 0.55_8
    z2 = Htot * (0.65_8 - 1._8)

    do i = 1,nx !!!  I need to know my global index range
       do j = 1,ny 
          x = (real(i+(pi*nx),kind=8)-0.5_8) * dx(j,i)
          y = (real(j+(pj*ny),kind=8)-0.5_8) * dy(j,i)
          do k = 1,nz
             rhs(k,j,i) = dx(j,i)*dy(j,i)*(zw(k+1,j,i)-zw(k,j,i)) * &
                  
                                ! pseudo 2D rhs
                  (exp(-bet * ((x-x1)**2 + (zr(k,j,i)-z1)**2)) - &
                  exp(-bet * ((x-x2)**2 + (zr(k,j,i)-z2)**2)))

             ! 3D rhs
             !                (exp(-bet * ((x-x1)**2 + (y-y1)**2 + (zr(k,j,i)-z1)**2)) - &
             !                 exp(-bet * ((x-x2)**2 + (y-y2)**2 + (zr(k,j,i)-z2)**2)))


             !           rhs(k,j,i) = (exp(-bet * ((x-x1)**2 + (zr(k,j,i)-z1)**2)) - &
             !                         exp(-bet * ((x-x2)**2 + (zr(k,j,i)-z2)**2)))
             !           rhs(k,j,i) = rhs(k,j,i) * rmask(j,i)           
          enddo
       enddo
    enddo

!    call random_number(rhs)
    rhs = 2*rhs-1.
    do i = 1,nx 
       do j = 1,ny 
          do k = 1,nz
             rhs(k,j,i) = rhs(k,j,i) * rmask(j,i)&
                  * dx(j,i)*dy(j,i)*(zw(k+1,j,i)-zw(k,j,i))
          enddo
       enddo
    enddo
    call fill_halo(1,rhs)

  end subroutine setup_rhs

  !-------------------------------------------------------------------------     
  subroutine setup_smooth_rhs(lev)
    ! generate white noise on level = lev
    ! then interpolate it back to lev=1
    ! you get a smooth random field, at the scale
    ! of 2*dx(lev)
    !
    integer(kind=4):: l,lev

    call random_number(grid(lev)%p)
    grid(lev)%p = 2.*grid(lev)%p-1.
    call fill_halo(lev,grid(lev)%p)
    

    do l=lev-1,1,-1
       grid(l)%p=0.
       grid(l)%b=0.
       call coarse2fine(l)
       call relax(l,4)
       if(myrank==0)write(*,*)l,grid(l)%p(1,5,5)
    enddo
    grid(1)%b=grid(1)%p
    grid(1)%p=0.

!     do i = 0,grid(1)%nx+1
!       do j = 0,grid(1)%ny+1 
!             grid(1)%b(:,j,i) = grid(1)%b(:,j,i) * rmask(j,i)
!          enddo
!       enddo
   
    
  end subroutine setup_smooth_rhs

  !-------------------------------------------------------------------------     
  subroutine setup_random_patches()
    
    ! define the r.h.s. as a sum of random gaussian patches
    ! sign, width, location are random

    integer(kind=4):: nx,ny,nz
    integer(kind=4):: i,j,k,l
    integer(kind=4):: npxg,npyg,pi,pj
    integer(kind=4):: nbpatch,i0,j0,k0,sign
    real(kind=8) :: sigh2,sigv2,dh2,dv2,x0,y0,z0

    nbpatch = 100

    sigh2=10.**2
    sigv2=3.**2
       
    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz


    npxg = grid(1)%npx
    npyg = grid(1)%npy
    pj = myrank/npxg
    pi = mod(myrank,npxg)

    grid(1)%b = 0.
    do l=1,nbpatch

       call random_number(x0)
       sign=floor(x0*2)*2-1

       call random_number(x0)
       sigh2=(x0*50)**2
       call random_number(x0)
       sigv2=(x0*20)**2

       call random_number(x0)
       call random_number(y0)
       call random_number(z0)
       ! we use the (i,j,k) global coordinates
       ! instead of (x,y,z)
       i0=1+floor(x0*npxg*nx)
       j0=1+floor(y0*npyg*ny)
       k0 = 1+floor(z0*nz)
       !!if(myrank==0)write(*,10)i0,j0,k0,sign,sqrt(sigh2),sqrt(sigv2)
       do i = 0,nx+1
          do j = 0,ny+1 
             do k = 1,nz
                dh2 = 1.*(i+pi*nx-i0)**2+1.*(j+pj*ny-j0)**2
                dv2 = 1.*(k-k0)**2
                grid(1)%b(k,j,i) = grid(1)%b(k,j,i) +&
                     sign*exp( -dh2 / (2.*sigh2) -dv2/(2.*sigv2))*rmask(j,i)
             enddo
          enddo
       enddo       
    enddo
 
10  format("(i0,j0,k0) =",I5,I5,I4," - sign,sigh,sigv = ",I3,F6.1,F6.1)

    call fill_halo(1,grid(1)%b)

  end subroutine setup_random_patches

!!$  !-------------------------------------------------------------------------     
!!$  subroutine setup_scoord(nx,ny,nz,nh)
!!$    ! compute zr and zw from h(j,i)
!!$
!!$    integer(kind=4):: nx,ny,nz,nh
!!$
!!$    integer(kind=4):: i,j,k
!!$    real(kind=8) :: cff,cff1,cff2,hinv,sc_w,sc_r,cff_w,cff_r,z_w0,z_r0,cs_r,cs_w
!!$    real(kind=8) :: csrf, cswf
!!$
!!$    if (trim(coord) == 'new_s_coord') then
!!$
!!$       ! vertical coordinate
!!$       do i = 0,nx+1
!!$          do j = 0,ny+1
!!$             ! new s-coord
!!$             cff=1./float(nz)
!!$             do k = 1,nz
!!$                sc_r = cff*(float(k-nz)-0.5_8)
!!$                if (theta_s.gt.0._8) then
!!$                   csrf=(1._8-cosh(theta_s*sc_r))/(cosh(theta_s)-1._8)
!!$                else
!!$                   csrf=-sc_r**2
!!$                endif
!!$                if (theta_b.gt.0._8) then
!!$                   cs_r=(exp(theta_b*csrf)-1._8)/(1._8-exp(-theta_b))
!!$                else
!!$                   cs_r=csrf
!!$                endif
!!$                sc_w = cff*float(k-1-nz)
!!$                if (theta_s.gt.0._8) then
!!$                   cswf=(1.D0-cosh(theta_s*sc_w))/(cosh(theta_s)-1._8)
!!$                else
!!$                   cswf=-sc_w**2
!!$                endif
!!$                if (theta_b.gt.0._8) then
!!$                   cs_w=(exp(theta_b*cswf)-1._8)/(1._8-exp(-theta_b))
!!$                else
!!$                   cs_w=cswf
!!$                endif
!!$                cff_w = hlim * sc_w
!!$                cff_r = hlim * sc_r
!!$                z_w0 = cff_w + cs_w*h(j,i)
!!$                z_r0 = cff_r + cs_r*h(j,i)
!!$                hinv = 1._8 / (h(j,i)+hlim)
!!$                zw(k,j,i) = z_w0 * (h(j,i)*hinv)
!!$                zr(k,j,i) = z_r0 * (h(j,i)*hinv)
!!$             enddo
!!$             zw(nz+1,j,i) = 0._8
!!$          enddo
!!$       enddo
!!$
!!$    else
!!$
!!$       cff  = 1./float(nz)
!!$       cff1 = 1./sinh(theta_s)
!!$       cff2 = 0.5/tanh(0.5*theta_s)
!!$
!!$       do i = 0,nx+1
!!$          do j = 0,ny+1
!!$             do k = 1,nz
!!$                sc_w = cff*float(k-1-nz)   !!sc_w = cff*float(k-nz)
!!$                sc_r = cff*(float(k-nz)-0.5)
!!$
!!$                cs_r = (1.-theta_b)*cff1*sinh(theta_s*sc_r) &
!!$                     +theta_b*(cff2*tanh(theta_s*(sc_r+0.5))-0.5)
!!$
!!$                cs_w = (1.-theta_b)*cff1*sinh(theta_s*sc_w) &
!!$                     +theta_b*(cff2*tanh(theta_s*(sc_w+0.5))-0.5)
!!$
!!$                cff_w = hlim * sc_w
!!$                cff_r = hlim * sc_r
!!$
!!$                z_w0 = cff_w + cs_w*h(j,i) 
!!$                z_r0 = cff_r + cs_r*h(j,i)  
!!$
!!$                hinv = 1. / (h(j,i)+hlim)
!!$
!!$                ! roms sigma coordinates
!!$                zw(k,j,i) = z_w0 * (h(j,i)*hinv)
!!$                zr(k,j,i) = z_r0 * (h(j,i)*hinv)
!!$
!!$                ! basic linear sigma coordinates
!!$                !             zr(k,j,i) = (real(k,kind=8)-0.5_8)*h(j,i)/real(nz,kind=8) - h(j,i)
!!$                !             zw(k,j,i) = (real(k,kind=8)-1.0_8)*h(j,i)/real(nz,kind=8) - h(j,i)
!!$             enddo
!!$             zw(nz+1,j,i) = 0.0_8
!!$          enddo
!!$       enddo
!!$
!!$    endif
!!$
!!$    !    if (netcdf_output) then    
!!$    !       call write_netcdf(zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank)
!!$    !    endif
!!$
!!$  end subroutine setup_scoord

end module mg_seamount
