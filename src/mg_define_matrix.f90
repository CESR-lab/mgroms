module mg_define_matrix

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_zr_zw
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

  interface define_matrices
     module procedure           &
          define_matrices_topo
  end interface define_matrices

  !NG comment: constants in a mg_cst.f90 file ?
  real(kind=rp), parameter :: one  = 1._rp
  real(kind=rp), parameter :: eigh = one/8._rp
  real(kind=rp), parameter :: qrt  = 0.25_rp
  real(kind=rp), parameter :: hlf  = 0.5_rp

contains

  !-------------------------------------------------------------------------  
  subroutine define_matrices_topo(dx, dy, zeta, h, rmask)

    real(kind=rp), dimension(:,:), intent(in) :: dx
    real(kind=rp), dimension(:,:), intent(in) :: dy
    real(kind=rp), dimension(:,:), intent(in) :: zeta
    real(kind=rp), dimension(:,:), intent(in) :: h
    real(kind=rp), dimension(:,:), pointer, intent(in) :: rmask

    integer(kind=ip)::  lev

    real(kind=rp), dimension(:,:), pointer :: dxf
    real(kind=rp), dimension(:,:), pointer :: dxc

    real(kind=rp), dimension(:,:), pointer :: dyf
    real(kind=rp), dimension(:,:), pointer :: dyc

    real(kind=rp), dimension(:,:), pointer :: zetaf
    real(kind=rp), dimension(:,:), pointer :: zetac

    real(kind=rp), dimension(:,:), pointer :: hf
    real(kind=rp), dimension(:,:), pointer :: hc

    real(kind=rp), dimension(:,:,:), pointer :: zrc
    real(kind=rp), dimension(:,:,:), pointer :: zwc

    real(kind=rp), dimension(:,:), pointer :: rmaski

    integer(kind=ip) :: nz,ny,nx
    integer(kind=ip) :: nzf,nyf,nxf
    integer(kind=ip) :: nyc,nxc

    if (myrank==0) write(*,*)'- define matrices from topography (h):'

    do lev = 1, nlevs

       if (myrank==0) write(*,*)'   lev=',lev

       nx=grid(lev)%nx
       ny=grid(lev)%ny
       nz=grid(lev)%nz

       if (lev == 1) then

          !NG: WARNING dx, dy and h model have to be defined 
          !NG: on (ny,nx) and not on (nx,ny) !!
          grid(lev)%dx(0:ny+1,0:nx+1)   = dx
          grid(lev)%dy(0:ny+1,0:nx+1)   = dy
          grid(lev)%zeta(0:ny+1,0:nx+1) = zeta
          grid(lev)%h (0:ny+1,0:nx+1)   = h

          allocate(rmaski(0:ny+1,0:nx+1))
          rmaski = rmask

       else
          nxf =grid(lev-1)%nx
          nyf =grid(lev-1)%ny
          nzf =grid(lev-1)%nz

          dxf => grid(lev-1)%dx
          dyf => grid(lev-1)%dy
          zetaf  => grid(lev-1)%zeta
          hf  => grid(lev-1)%h

          if (grid(lev)%gather == 1) then
             nxc= nx/grid(lev)%ngx
             nyc= ny/grid(lev)%ngy

             allocate(dxc(0:nyc+1,0:nxc+1))
             allocate(dyc(0:nyc+1,0:nxc+1))
             allocate(zetac(0:nyc+1,0:nxc+1))
             allocate( hc(0:nyc+1,0:nxc+1))
          else
             nxc = nx
             nyc = ny
             dxc => grid(lev)%dx
             dyc => grid(lev)%dy
             zetac  => grid(lev)%zeta
             hc  => grid(lev)%h
          endif

          if ((aggressive).and.(lev==1)) then

             write(*,*) ' define matrices (aggressive).and.(lev==1) not available !'
             STOP

          else

             ! Coarsen dx, dy and h
             dxc(1:nyc,1:nxc) = hlf      * ( &
                  dxf(1:nyf  :2,1:nxf  :2) + &
                  dxf(2:nyf+1:2,1:nxf  :2) + &
                  dxf(1:nyf  :2,2:nxf+1:2) + &
                  dxf(2:nyf+1:2,2:nxf+1:2) )

             dyc(1:nyc,1:nxc) = hlf      * ( &
                  dyf(1:nyf  :2,1:nxf  :2) + &
                  dyf(2:nyf+1:2,1:nxf  :2) + &
                  dyf(1:nyf  :2,2:nxf+1:2) + &
                  dyf(2:nyf+1:2,2:nxf+1:2) )

             zetac(1:nyc,1:nxc)  = qrt      * ( &
                  zetaf(1:nyf  :2,1:nxf  :2)  + &
                  zetaf(2:nyf+1:2,1:nxf  :2)  + &
                  zetaf(1:nyf  :2,2:nxf+1:2)  + &
                  zetaf(2:nyf+1:2,2:nxf+1:2)  )

             hc(1:nyc,1:nxc)  = qrt      * ( &
                  hf(1:nyf  :2,1:nxf  :2)  + &
                  hf(2:nyf+1:2,1:nxf  :2)  + &
                  hf(1:nyf  :2,2:nxf+1:2)  + &
                  hf(2:nyf+1:2,2:nxf+1:2)  )

          endif ! aggressive

          if (grid(lev)%gather == 1) then

             call gather(lev,dxc,grid(lev)%dx)
             call gather(lev,dyc,grid(lev)%dy)
             call gather(lev,zetac,grid(lev)%zeta)
             call gather(lev, hc,grid(lev)%h)

             deallocate(dxc)
             deallocate(dyc)
             deallocate(zetac)
             deallocate( hc)
          endif

          !- boundary mask (gather not tested!)
          if (associated(rmaski)) deallocate(rmaski)
          allocate(rmaski(0:ny+1,0:nx+1))
          rmaski(:,:) = 1._rp
          if (bmask) then
             call fill_halo_2D_bmask(lev,rmaski)
          endif

       endif ! lev == 1

       call fill_halo(lev, grid(lev)%dx)
       call fill_halo(lev, grid(lev)%dy)
       call fill_halo(lev, grid(lev)%zeta)
       call fill_halo(lev, grid(lev)%h)

       zrc => grid(lev)%zr
       zwc => grid(lev)%zw

       ! Compute zr and zw
       call setup_zr_zw                (  & 
            hlim,nhtheta_b,nhtheta_s,     &
            grid(lev)%zeta,grid(lev)%h,   &  ! input args
            grid(lev)%zr, grid(lev)%zw,   &  ! output args
            coord_type='new_s_coord'      )  ! optional

       call fill_halo(lev,grid(lev)%zr) ! Special fill_halo nh = 2
       call fill_halo(lev,grid(lev)%zw) ! Special fill_halo nh = 2

       if (netcdf_output) then
          call write_netcdf(grid(lev)%dx,vname='dx',netcdf_file_name='dx.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%dy,vname='dy',netcdf_file_name='dy.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zeta,vname='zeta',netcdf_file_name='zeta.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%h ,vname='h' ,netcdf_file_name='h.nc' ,rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zw,vname='zw',netcdf_file_name='zw.nc',rank=myrank,iter=lev)
       endif

!!$       if (myrank==0) then
!!$          write(*,*)'2 Min, max dx:', lev, minval(grid(lev)%dx),maxval(grid(lev)%dx)
!!$          write(*,*)'2 Min, max dy:', lev, minval(grid(lev)%dy),maxval(grid(lev)%dy)
!!$       endif

       ! Define matrix coefficients from dx, dy, zr and zw coarsened
       call define_matrix(lev, grid(lev)%dx, grid(lev)%dy, grid(lev)%zr, grid(lev)%zw, rmaski)

    enddo ! lev

    if (associated(rmaski)) deallocate(rmaski)

  end subroutine define_matrices_topo

  !-----------------------------------------------------------------------------------
  subroutine define_matrix(lev, dx, dy, zr, zw, rmask)

    integer(kind=ip),intent(in):: lev
    real(kind=rp), dimension(:,:)  , pointer, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: zr, zw
    real(kind=rp), dimension(:,:)  , pointer, intent(in) :: rmask

    ! Define matrix coefficients cA
    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)xd
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  ,   pointer :: umask, vmask

    real(kind=rp) :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: dzw
    real(kind=rp), dimension(:,:,:),   pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:),   pointer :: cw
    real(kind=rp), dimension(:,:,:,:), pointer :: cA

    real(kind=rp), dimension(:,:),   pointer :: zx ! dbg: to remove !!!

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz

!!$    if (myrank==0) then
!!$       write(*,*)'3 Min, max dx:', lev, minval(dx),maxval(dx)
!!$       write(*,*)'3 Min, max dy:', lev, minval(dy),maxval(dy)
!!$    endif

    cA => grid(lev)%cA 

    allocate(zx(0:ny+1,0:nx+1))

    ! Umask and Vmask from rmask
    allocate(umask(0:ny+1,0:nx+1))
    allocate(vmask(0:ny+1,0:nx+1))

    if (bmask) then
       umask(:,:)=0._rp
       vmask(:,:)=0._rp
       do i = 1,nx+1
          do j = 0,ny+1
             umask(j,i) = rmask(j,i-1) * rmask(j,i)
          enddo
       enddo
       do i = 0,nx+1
          do j = 1,ny+1
             vmask(j,i) = rmask(j-1,i) * rmask(j,i)
          enddo
       enddo
    else
       umask(:,:)=1._rp
       vmask(:,:)=1._rp
    endif

!!$    if (myrank==0) then
!!$       write(*,*)' Min max umask:', minval(umask), maxval(umask)
!!$       write(*,*)' Min max vmask:', minval(vmask), maxval(vmask)
!!$    endif

    !- Used in compute_rhs -!
    if (lev == 1) then
       dzw => grid(lev)%dzw
       do i = 0,nx+1
          do j = 0,ny+1
             dzw(1,j,i) = zr(1,j,i) - zw(1,j,i) !!
             do k = 2,nz
                dzw(k,j,i) = zr(k,j,i) - zr(k-1,j,i) !!  cell height at w-points
             enddo
             dzw(nz+1,j,i) = zw(nz+1,j,i) - zr(nz,j,i) !!
          enddo
       enddo

       !! Slopes in x- and y-direction defined at rho-points
       zxdy => grid(lev)%zxdy
       zydx => grid(lev)%zydx
       do i = 0,nx+1
          do j = 0,ny+1
             do k = 1, nz
                zydx(k,j,i) = hlf * (( zr(k,j+1,i  ) - zr(k,j-1,i  ) ) / dy(j,i) ) * dx(j,i)
                zxdy(k,j,i) = hlf * (( zr(k,j  ,i+1) - zr(k,j  ,i-1) ) / dx(j,i) ) * dy(j,i)
             enddo
          enddo
       enddo
    endif

    !- Used also in compute_rhs -!
    cw => grid(lev)%cw
    do i = 0,nx+1
       do j = 0,ny+1

          Arz = dx(j,i)*dy(j,i)

          k=1
          cw(k,j,i) = ( Arz / (zr(k,j,i)-zw(k,j,i)) ) * &
               (one + &
               ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
               ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

          do k = 2,nz
             cw(k,j,i) = ( Arz / (zr(k,j,i)-zr(k-1,j,i)) ) * &
                  (one + &
                  ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
                  ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

          enddo

          k=nz+1
          cw(k,j,i) = ( Arz / (zw(k,j,i)-zr(k-1,j,i)) ) * &
               (one + &
               ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
               ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

       enddo
    enddo

!!$    if (myrank==0) then
!!$       if (lev==4) then
!!$          write(*,*)'---> CW lev=4 test vals <----'
!!$          write(*,*) cw(1,0,0)
!!$          write(*,*) dx(0,0), dy(0,0)
!!$          write(*,*)'---> --------------- <----'
!!$       end if
!!$    endif

    !! interaction coeff with neighbours
    !XX
    !---------!
    !- K = 1 -! lower level
    !---------!
    k = 1

    do i = 1,nx
       do j = 1,ny+1

          cA(3,k,j,i) = qrt * ( &          ! couples with k+1 j-1    
               ( hlf * (zr(k+1,j+1,i)-zr(k+1,j-1,i)) / dy(j  ,i) ) * dx(j  ,i) + &
               ( hlf * (zr(k  ,j  ,i)-zr(k  ,j-2,i)) / dy(j-1,i) ) * dx(j-1,i) ) * vmask(j,i)

          cA(4,k,j,i) = & ! couples with j-1 
               ( qrt                                                                      * &  ! 
               ( zw(k+1,j,i)-zw(k,j,i)+zw(k+1,j-1,i)-zw(k,j-1,i)) * (dx(j,i)+dx(j-1,i)) ) / &  ! Ary(k,j,i) /
               ( hlf * (dy(j,i)+dy(j-1,i)))                                                 &  ! dyv(j,i)
                                ! topo terms 
               - ( (( hlf * (zr(k,j+1,i)-zr(k,j-1,i)) / dy(j,i) ) * dx(j  ,i) )**2 /        & ! zydx(k,j,i)*zydx(k,j,i) /
               ( cw(k,j,i) + cw(k+1,j,i) )                                                  & ! (cw(k,j,i)+cw(k+1,j,i)
               +   (( hlf * (zr(k,j,i)-zr(k,j-2,i)) / dy(j-1,i) ) * dx(j-1,i) )**2 /        & ! zydx(k,j-1,i)*zydx(k,j-1,i)/
               ( cw(k,j-1,i) + cw(k+1,j-1,i) ) )                                            & ! (cw(k,j-1,i)+cw(k+1,j-1,i)
                                ! from j,k cross terms
               - qrt *  ( &
               ( hlf * (zr(k,j  ,i) - zr(k,j-2,i)) / dy(j-1,i) ) * dx(j-1,i)              - & ! zydx(k,j-1,i)
               ( hlf * (zr(k,j+1,i) - zr(k,j-1,i)) / dy(j  ,i) ) * dx(j  ,i) )                ! zydx(k,j,i)

          if (bmask) then
             cA(4,k,j,i) = ( cA(4,k,j,i)  &
                                ! from i,j cross terms if lbc                                         
                  -( hlf                                                            * & ! - 0.5
                  (( hlf * (zr(k,j-1,i+1)-zr(k,j-1,i-1)) / dx(j-1,i) ) * dy(j-1,i)) * & ! * zxdy(k,j-1,i)
                  (( hlf * (zr(k,j  ,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i)) / & ! * zydx(k,j-1,i)
                  (cw(k,j-1,i)+cw(k+1,j-1,i))                                       * & ! / cw(k,j-1,i)+cw(k+1,j-1,i)
                  (umask(j-1,i+1) - umask(j-1,i))                                     & ! * umask
                  - hlf                                                             * & ! - 0.5
                  (( hlf * (zr(k,j  ,i+1)-zr(k,j  ,i-1)) / dy(j,i) ) * dx(j,i))     * & ! * zxdy(k,j,i)
                  (( hlf * (zr(k,j+1,i  )-zr(k,j-1,i  )) / dy(j,i) ) * dx(j,i))     / & ! * zydx(k,j,i)
                  (cw(k,j,i)+cw(k+1,j,i))                                           * & ! / cw(k,j,i)+cw(k+1,j,i)
                  (umask(j,i+1) - umask(j,i)) ) )                                   * & ! * umask
                  vmask(j,i)                                                            ! * vmask all Ca4
          endif

       enddo
    enddo

    do i = 1,nx+1
       do j = 1,ny

          cA(6,k,j,i) = qrt * ( &
               ( hlf * (zr(k+1,j,i+1)-zr(k+1,j,i-1)) / dx(j,i  ) ) * dy(j,i  )   + &
               ( hlf * (zr(k  ,j,i  )-zr(k  ,j,i-2)) / dx(j,i-1) ) * dy(j,i-1) ) * & ! couples with k+1 i-1
               umask(j,i)

          cA(7,k,j,i) =                                                                    &
                                ! couples with i-1
               ( qrt                                                                     * & ! 0.25
               (zw(k+1,j,i)-zw(k,j,i)+zw(k+1,j,i-1)-zw(k,j,i-1)) * (dy(j,i)+dy(j,i-1)) ) / & ! Arx(k,j,i)
               ( hlf * (dx(j,i)+dx(j,i-1)) )                                               & ! dxu(j,i)
                                ! topo terms
               - ( (( hlf * (zr(k,j,i+1)-zr(k,j,i-1)) / dx(j,i) ) * dy(j,i))**2          / & ! zxdy(k,j,i)*zxdy(k,j,i)
               ( cw(k,j,i) + cw(k+1,j,i) )                                                 & ! cw(k,j,i-1)+cw(k+1,j,i-1)
               +   (( hlf * (zr(k,j,i)-zr(k,j,i-2)) / dx(j,i-1) ) * dy(j,i-1))**2        / & ! zxdy(k,j,i-1)*zxdy(k,j,i-1)
               ( cw(k,j,i-1) + cw(k+1,j,i-1) ) )                                           & ! cw(k,j,i-1)+cw(k+1,j,i-1)
                                ! from i,k cross terms
               - qrt                                                                   * ( & ! 0.25
               ( hlf * (zr(k,j,i  )-zr(k,j,i-2)) / dx(j,i-1) ) * dy(j,i-1)               - & ! zxdy(k,j,i-1)
               ( hlf * (zr(k,j,i+1)-zr(k,j,i-1)) / dx(j,i  ) ) * dy(j,i)                 )   ! zxdy(k,j,i)

          if (bmask) then

             cA(7,k,j,i) = ( cA(7,k,j,i) &
                                ! from i,j cross terms if lbc                                         
                  -( hlf                                                            * & ! - 0.5
                  (( hlf * (zr(k,j  ,i  )-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1)) * & ! * zxdy(k,j,i-1)
                  (( hlf * (zr(k,j+1,i-1)-zr(k,j-1,i-1)) / dy(j,i-1) ) * dx(j,i-1)) / & ! * zydx(k,j,i-1)
                  (cw(k,j,i-1)+cw(k+1,j,i-1) )                                        & ! / cw(k,j,i-1)+cw(k+1,j,i-1)
                  * (vmask(j+1,i-1) - vmask(j,i-1))                                   & ! * vmask
                  - hlf                                                             * & ! - 0.5 
                  (( hlf * (zr(k,j  ,i+1)-zr(k,j  ,i-1)) / dy(j,i) ) * dx(j,i))     * & ! * zxdy(k,j,i)
                  (( hlf * (zr(k,j+1,i  )-zr(k,j-1,i  )) / dy(j,i) ) * dx(j,i))     / & ! * zydx(k,j,i)
                  (cw(k,j,i)+cw(k+1,j,i))                                             & ! / cw(k,j,i)+cw(k+1,j,i)
                  * (vmask(j+1,i) - vmask(j,i))) )                                  * & ! * vmask
                  umask(j,i)                                                            ! * umask

          endif

       enddo
    enddo

    do i = 1,nx+1
       do j = 0,ny
          ! only for k==1, couples with j+1,i-1
          cA(5,k,j,i) = &
               + hlf                                                             * &
               (( hlf * (zr(k,j+1,i+1)-zr(k,j+1,i-1)) / dx(j+1,i) ) * dy(j+1,i)) * &
               (( hlf * (zr(k,j+2,i  )-zr(k,j  ,i  )) / dy(j+1,i) ) * dx(j+1,i)) / &
               ( cw(k,j+1,i ) + cw(k+1,j+1,i))                                * &
               umask(j+1,i) * vmask(j+1,i)                                         & ! mask
               + hlf                                                             * &
               (( hlf * (zr(k,j  ,i  )-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1)) * &
               (( hlf * (zr(k,j+1,i-1)-zr(k,j-1,i-1)) / dy(j,i-1) ) * dx(j,i-1)) / &
               ( cw(k,j,i-1) + cw(k+1,j,i-1))                                    * &
               umask(j,i) * vmask(j+1,i-1)                                           ! mask
       enddo
    enddo

!!$    if (myrank==0) then
!!$       if (lev==4) then
!!$          !!write(*,*) cA(5,k,:,:)
!!$          write(*,*)'---> lev=4 test vals <----'
!!$          write(*,*) cw(k,1,1)
!!$          write(*,*) cw(k+1,1,1)
!!$          write(*,*) umask(1,1), vmask(1,1)
!!$          write(*,*) cw(k,0,0)
!!$          write(*,*) cw(k+1,0,0)
!!$          write(*,*) umask(0,1), vmask(1,0)
!!$          write(*,*) cA(5,k,0,1)
!!$          write(*,*)'---> --------------- <----'
!!$       end if
!!$    endif

    do i = 1,nx+1
       do j = 1,ny+1
          ! only for k==1, couples with j-1,i-1
          cA(8,k,j,i) = &
               - hlf                                                             * &
               (( hlf * (zr(k,j-1,i+1)-zr(k,j-1,i-1)) / dx(j-1,i) ) * dy(j-1,i)) * &
               (( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i))   / &
               (cw(k,j-1,i  ) + cw(k+1,j-1,i  ))                                 * & 
               umask(j-1,i) * vmask(j,i)                                           & ! mask
               - hlf                                                             * &
               (( hlf * (zr(k,j,i)-zr(k,j,i-2)) / dx(j,i-1) ) * dy(j,i-1))       * &
               (( hlf * (zr(k,j+1,i-1)-zr(k,j-1,i-1)) / dy(j,i-1) ) * dx(j,i-1)) / &
               (cw(k,j  ,i-1) + cw(k+1,j  ,i-1))                                 * &
               umask(j,i) * vmask(j,i-1)                                             ! mask
       enddo
    enddo

    !XX
    !---------------!
    !- K = 2, nz-1 -! interior levels
    !---------------!

    do i = 1,nx
       do j = 1,ny
          do k = 2,nz-1 
             cA(2,k,j,i) =  cw(k,j,i)                                  ! couples with k-1

             if (bmask) then
                cA(2,k,j,i) = cA(2,k,j,i)                    &
                                ! from i,k  cross terms if lbc
                     - qrt * ( &
                     ( hlf * (zr(k-1,j,i+1)-zr(k-1,j,i-1)) / dx(j,i) ) * dy(j,i)   - & ! zxdy(k-1,j,i)
                     ( hlf * (zr(k  ,j,i+1)-zr(k  ,j,i-1)) / dx(j,i) ) * dy(j,i) ) * & ! zxdy(k,j,i)
                     (umask(j,i+1) - umask(j,i))             &
                                ! from j,k  cross terms if lbc
                     - qrt * ( &
                     ( hlf * (zr(k-1,j+1,i) - zr(k-1,j-1,i)) / dy(j,i) ) * dx(j,i)     - & ! zydx(k-1,j,i)
                     ( hlf * (zr(k  ,j+1,i) - zr(k  ,j-1,i)) / dy(j,i) ) * dx(j,i)   ) * & ! zydx(k,j,i)
                     (vmask(j+1,i) - vmask(j,i))
             endif

          enddo
       enddo
    enddo

    do i = 1,nx
       do j = 1,ny+1
          do k = 2,nz-1 

             cA(3,k,j,i) =  qrt * ( &
                  ( hlf * (zr(k+1,j+1,i  )-zr(k+1,j-1,i  )) / dy(j,i) ) * dx(j,i) + &
                  ( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i) ) * & ! couples with k+1 j-1
                  vmask(j,i)

             cA(4,k,j,i) =  ( qrt * &
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  (dx(j,i)+dx(j-1,i)) ) / ( hlf * (dy(j,i)+dy(j-1,i)) )     * & ! couples with j-1
                  vmask(j,i)

             cA(5,k,j,i) =- qrt * ( &
                  (( hlf * (zr(k-1,j+1,i  )-zr(k-1,j-1,i  )) / dy(j,i) ) * dx(j,i)) + &
                  (( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i)) ) * & ! couples with k-1 j-1
                  vmask(j,i)

          enddo
       enddo
    enddo

    do i = 1,nx+1
       do j = 1,ny 
          do k = 2,nz-1 

             cA(6,k,j,i) =  qrt * ( &
                  (( hlf * (zr(k+1,j  ,i+1)-zr(k+1,j  ,i-1)) / dx(j,i) ) * dy(j,i)) + &
                  (( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1)) ) * &    ! Couples with k+1 i-1
                  umask(j,i)

             cA(7,k,j,i) =   (qrt * &
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  (dy(j,i)+dy(j,i-1)) ) / ( hlf * (dx(j,i)+dx(j,i-1)) )     * & ! Couples with i-1
                  umask(j,i)

             cA(8,k,j,i) =- qrt * ( &
                  (( hlf * (zr(k-1,j  ,i+1)-zr(k-1,j  ,i-1)) / dx(j,i) ) * dy(j,i)) + &
                  (( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1)) ) * & ! Couples with k-1 i-1
                  umask(j,i)

          enddo
       enddo
    enddo

    !XX
    !----------!
    !- K = nz -! upper level
    !----------!
    k = nz

    do i = 1,nx    
       do j = 1,ny 
          cA(2,k,j,i) = cw(k,j,i)                                    ! couples with k-1
       enddo
    enddo

    do i = 1,nx
       do j = 1,ny+1

          cA(4,k,j,i) = & ! couples with j-1
               ( qrt * ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * (dx(j,i)+dx(j-1,i)) / & ! Ary(k,j,i)
               ( hlf * (dy(j,i)+dy(j-1,i)) )                                                             & ! dyv(j,i)
               + qrt * (                                                                                 &
               - (( hlf * (zr(k,j  ,i)-zr(k,j-2,i)) / dy(j-1,i) ) * dx(j-1,i))                           & ! zydx(k,j-1,i)
               + (( hlf * (zr(k,j+1,i)-zr(k,j-1,i)) / dy(j  ,i) ) * dx(j  ,i)) ) )                     * & ! zydx(k,j,i)
               vmask(j,i)

          cA(5,k,j,i) =- qrt * ( &
               (( hlf * (zr(k-1,j+1,i  )-zr(k-1,j-1,i  )) / dy(j,i) ) * dx(j,i)) + &
               (( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i)) ) * & ! couples with k-1 j-1
               vmask(j,i)

       enddo
    enddo

    do i = 1,nx+1
       do j = 1,ny 

          cA(7,k,j,i) =  & ! Couples with i-1
               (qrt * ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * (dy(j,i)+dy(j,i-1))  / & ! Arx(k,j,i)
               ( hlf * (dx(j,i)+dx(j,i-1)) )                                                            & ! dxu(j,i)
               + qrt * (                                                                                 &
               - (( hlf * (zr(k,j,i  )-zr(k,j,i-2)) / dx(j,i-1)) * dy(j,i-1))                            & ! zxdy(k,j,i-1)
               + (( hlf * (zr(k,j,i+1)-zr(k,j,i-1)) / dx(j,i  )) * dy(j,i  )) ) )                      * & ! zxdy(k,j,i)
               umask(j,i)

          cA(8,k,j,i) =- qrt * ( &
               (( hlf * (zr(k-1,j  ,i+1)-zr(k-1,j  ,i-1)) / dx(j,i) ) * dy(j,i)) + &
               (( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1)) ) * & ! Couples with k-1 i-1
               umask(j,i)

       enddo
    enddo

    if (bmask) then
       call fill_halo(lev,cA)
    endif

    !! interaction coeff with itself
    do i = 1,nx
       do j = 1,ny

          k = 1 !lower level
          cA(1,k,j,i) =                     &
               -cA(2,k+1,j,i)               &
               -cA(4,k,j,i)-cA(4,k,j+1,i)   &
               -cA(7,k,j,i)-cA(7,k,j,i+1)   &
               -cA(6,k,j,i)-cA(8,k+1,j,i+1) &
               -cA(3,k,j,i)-cA(5,k+1,j+1,i) &
               -cA(5,k,j,i)-cA(5,k,j-1,i+1) &
               -cA(8,k,j,i)-cA(8,k,j+1,i+1)

          do k = 2,nz-1 !interior levels
             cA(1,k,j,i) = &
                  -cA(2,k,j,i)-cA(2,k+1,j,i)   &
                  -cA(4,k,j,i)-cA(4,k,j+1,i)   &
                  -cA(7,k,j,i)-cA(7,k,j,i+1)   &
                  -cA(6,k,j,i)-cA(6,k-1,j,i+1) &
                  -cA(8,k,j,i)-cA(8,k+1,j,i+1) & 
                  -cA(3,k,j,i)-cA(3,k-1,j+1,i) &
                  -cA(5,k,j,i)-cA(5,k+1,j+1,i)

          enddo

          k = nz !upper level
          cA(1,k,j,i) =                    &
               - cA(2,k,j,i)               &
               - cw(k+1,j,i)               &
               + hlf * (hlf * (zr(k,j  ,i+2) - zr(k,j  ,i  )) / dx(j  ,i+1)) * dy(j  ,i+1) & ! + 2 * 0.25 * zxdy(nz,j,i+1)
               - hlf * (hlf * (zr(k,j  ,i  ) - zr(k,j  ,i-2)) / dx(j  ,i-1)) * dy(j  ,i-1) & ! - 2 * 0.25 * zxdy(nz,j,i-1)
               + hlf * (hlf * (zr(k,j+2,i  ) - zr(k,j  ,i  )) / dy(j+1,i  )) * dx(j+1,i  ) & ! + 2 * 0.25 * zydx(nz,j+1,i)
               - hlf * (hlf * (zr(k,j  ,i  ) - zr(k,j-2,i  )) / dy(j-1,i  )) * dx(j-1,i  ) & ! - 2 * 0.25 * zydx(nz,j-1,i)
               - cA(4,k  ,j,i  ) - cA(4,k,j+1,i  ) &
               - cA(7,k  ,j,i  ) - cA(7,k,j  ,i+1) &
               - cA(6,k-1,j,i+1)                   &
               - cA(8,k  ,j,i  )                   &
               - cA(3,k-1,j+1,i)                   &
               - cA(5,k  ,j,i  )

       enddo ! j
    enddo ! i

    k = nz
    do i = 0,nx+1
       do j = 0, ny+1
          zx(j,i) = hlf * (zr(k,j,i+1) - zr(k,j,i-1)) / dx(j ,i)
      enddo ! j
    enddo ! i

    

    if (netcdf_output) then
       if (myrank==0) write(*,*)'       write cA in a netcdf file'
       call write_netcdf(grid(lev)%cA,vname='ca',netcdf_file_name='cA.nc',rank=myrank,iter=lev)
       call write_netcdf(zx,vname='zx',netcdf_file_name='zx.nc',rank=myrank,iter=lev)
       call write_netcdf(zr(nz,0:ny+1,0:nx+1),vname='zr2',netcdf_file_name='zr2.nc',rank=myrank,iter=lev)
    endif

    deallocate(zx)

    deallocate(umask)
    deallocate(vmask)

  end subroutine define_matrix

end module mg_define_matrix
