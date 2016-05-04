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
  real(kind=rp), parameter :: eigh = one/8._8
  real(kind=rp), parameter :: qrt  = 0.25_rp
  real(kind=rp), parameter :: hlf  = 0.5_rp

contains

  !-------------------------------------------------------------------------  
  subroutine define_matrices_topo(dx, dy, h)
    real(kind=rp), dimension(:,:), pointer, intent(in) :: dx
    real(kind=rp), dimension(:,:), pointer, intent(in) :: dy
    real(kind=rp), dimension(:,:), pointer, intent(in) :: h

    integer(kind=ip)::  lev

    real(kind=rp), dimension(:,:), pointer :: dxf
    real(kind=rp), dimension(:,:), pointer :: dxc

    real(kind=rp), dimension(:,:), pointer :: dyf
    real(kind=rp), dimension(:,:), pointer :: dyc

    real(kind=rp), dimension(:,:), pointer :: hf
    real(kind=rp), dimension(:,:), pointer :: hc

    real(kind=rp), dimension(:,:,:), pointer :: zrc
    real(kind=rp), dimension(:,:,:), pointer :: zwc

    integer(kind=ip) :: nz,ny,nx,nh
    integer(kind=ip) :: nzf,nyf,nxf,nhf
    integer(kind=ip) :: nyc,nxc

    if (myrank==0) write(*,*)'- define matrices from topography (h):'

    do lev = 1, nlevs

       if (myrank==0) write(*,*)'   lev=',lev

       nx=grid(lev)%nx
       ny=grid(lev)%ny
       nz=grid(lev)%nz
       nh=grid(lev)%nh

       if (lev == 1) then

          grid(lev)%dx(1-nh:ny+nh,1-nh:nx+nh) = dx
          grid(lev)%dy(1-nh:ny+nh,1-nh:nx+nh) = dy
          grid(lev)%h (1-nh:ny+nh,1-nh:nx+nh) =  h

       else
          nxf =grid(lev-1)%nx
          nyf =grid(lev-1)%ny
          nzf =grid(lev-1)%nz
          nhf =grid(lev-1)%nh

          dxf => grid(lev-1)%dx
          dyf => grid(lev-1)%dy
          hf  => grid(lev-1)%h

          if (grid(lev)%gather == 1) then
             nxc= nx/grid(lev)%ngx
             nyc= ny/grid(lev)%ngy

             allocate(dxc(1-nh:nyc+nh,1-nh:nxc+nh))
             allocate(dyc(1-nh:nyc+nh,1-nh:nxc+nh))
             allocate( hc(1-nh:nyc+nh,1-nh:nxc+nh))
          else
             nxc = nx
             nyc = ny
             dxc => grid(lev)%dx
             dyc => grid(lev)%dy
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

             hc(1:nyc,1:nxc)  = qrt      * ( &
                  hf(1:nyf  :2,1:nxf  :2)  + &
                  hf(2:nyf+1:2,1:nxf  :2)  + &
                  hf(1:nyf  :2,2:nxf+1:2)  + &
                  hf(2:nyf+1:2,2:nxf+1:2)  )

          endif ! aggressive

          if (grid(lev)%gather == 1) then

             call gather(lev,dxc,grid(lev)%dx)
             call gather(lev,dyc,grid(lev)%dy)
             call gather(lev, hc,grid(lev)%h)

             deallocate(dxc)
             deallocate(dyc)
             deallocate( hc)
          endif

       endif ! lev == 1

       call fill_halo(lev, grid(lev)%dx)
       call fill_halo(lev, grid(lev)%dy)
       call fill_halo(lev, grid(lev)%h )

       zrc => grid(lev)%zr
       zwc => grid(lev)%zw

       ! Compute zr and zw
       if (trim(bench) == 'seamount') then
          call setup_zr_zw(grid(lev)%h,grid(lev)%zr,grid(lev)%zw)
          !!call setup_zr_zw                     ( & 
          !!     hlim,theta_b,theta_s,grid(lev)%h, &  ! input args
          !!     grid(lev)%zr, grid(lev)%zw      )    ! output args
       else
          call setup_zr_zw                     ( & 
               hlim,theta_b,theta_s,grid(lev)%h, &  ! input args
               grid(lev)%zr, grid(lev)%zw      , &  ! output args
               coord_type='new_s_coord'        )    ! optional
       endif

       if (netcdf_output) then
          call write_netcdf(grid(lev)%dx,vname='dx',netcdf_file_name='dx.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%dy,vname='dy',netcdf_file_name='dy.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%h ,vname='h' ,netcdf_file_name='h.nc' ,rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zw,vname='zw',netcdf_file_name='zw.nc',rank=myrank,iter=lev)
       endif

       ! Define matrix coefficients from dx, dy, zr and zw coarsened
       call define_matrix(lev, grid(lev)%dx, grid(lev)%dy, grid(lev)%zr, grid(lev)%zw)

    enddo ! lev

  end subroutine define_matrices_topo

  !-----------------------------------------------------------------------------------
  subroutine define_matrix(lev, dx, dy, zr, zw)

    integer(kind=ip),intent(in):: lev
    real(kind=rp), dimension(:,:),   pointer, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: zr, zw

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
    integer(kind=ip):: nh

    real(kind=rp), dimension(:,:,:,:), pointer :: cA
    real(kind=rp), dimension(:,:)  ,   pointer :: dxu, dyv
    real(kind=rp), dimension(:,:,:),   pointer :: Arx, Ary
    real(kind=rp), dimension(:,:)  ,   pointer :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: dz
    real(kind=rp), dimension(:,:,:),   pointer :: dzw
    real(kind=rp), dimension(:,:,:),   pointer :: zy,zx
    real(kind=rp), dimension(:,:,:),   pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:),   pointer :: zxw,zyw
    real(kind=rp), dimension(:,:,:),   pointer :: cw

    ! TODO NG
    ! zw,zr can change in time
    ! dx,dy constant in time

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    cA => grid(lev)%cA 

    !NG comment: perf ? allocate-deallocate -> mg_grid ?

    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i) !!  cell height at rho-points
          enddo
       enddo
    enddo

    allocate(dzw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          dzw(1,j,i) = zr(1,j,i)-zw(1,j,i) !!
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i) !!  cell height at w-points
          enddo
          dzw(nz+1,j,i) = zw(nz+1,j,i)-zr(nz,j,i) !!
       enddo
    enddo

    !! Cell widths
    allocate(dxu(ny,nx+1))
    do i = 1,nx+1
       do j = 1,ny
          dxu(j,i) = hlf * (dx(j,i)+dx(j,i-1))
       enddo
    enddo

    allocate(dyv(ny+1,nx))
    do i = 1,nx
       do j = 1,ny+1
          dyv(j,i) = hlf * (dy(j,i)+dy(j-1,i))
       enddo
    enddo

    !!  Areas
    allocate(Arx(nz,ny,nx+1))
    do i = 1,nx+1
       do j = 1,ny
          do k = 1,nz
             Arx(k,j,i) = qrt * (dz(k,j,i)+dz(k,j,i-1)) * (dy(j,i)+dy(j,i-1))
          enddo
       enddo
    enddo

    allocate(Ary(nz,ny+1,nx))
    do i = 1,nx
       do j = 1,ny+1
          do k = 1,nz
             Ary(k,j,i) = qrt * (dz(k,j,i)+dz(k,j-1,i)) * (dx(j,i)+dx(j-1,i)) 
          enddo
       enddo
    enddo
    allocate(Arz(0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          Arz(j,i) = dx(j,i)*dy(j,i)
       enddo
    enddo

    !! Slopes in x- and y-direction defined at rho-points
    allocate(zx(nz,0:ny+1,0:nx+1))
    allocate(zy(nz,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             zy(k,j,i) = hlf * (zr(k,j+1,i  )-zr(k,j-1,i  )) / dy(j,i)
             zx(k,j,i) = hlf * (zr(k,j  ,i+1)-zr(k,j  ,i-1)) / dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(lev,zy)
    call fill_halo(lev,zx) 

    allocate(zxdy(nz,0:ny+1,0:nx+1))
    allocate(zydx(nz,0:ny+1,0:nx+1))
    do k = 1,nz
       zydx(k,:,:) = zy(k,:,:)*dx(:,:)
       zxdy(k,:,:) = zx(k,:,:)*dy(:,:)
    enddo

    allocate(zyw(nz+1,0:ny+1,0:nx+1))
    allocate(zxw(nz+1,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz+1
             zyw(k,j,i) = hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i)
             zxw(k,j,i) = hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(lev,zyw)
    call fill_halo(lev,zxw)

    allocate(cw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             cw(k,j,i) = Arz(j,i)/dzw(k,j,i) * (one + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i))
          enddo
       enddo
    enddo

    !! interaction coeff with neighbours
    do i = 1,nx
       do j = 1,ny

          !---------!
          !- K = 1 -! lower level
          !---------!
          k = 1
          cA(3,k,j,i) = qrt * ( zydx(k+1,j,i) + zydx(k,j-1,i) ) ! couples with k+1 j-1
          cA(4,k,j,i) =                                                                &
               ! couples with j-1
               Ary(k,j,i)/dyv(j,i)                                                     & 
               ! topo terms 
               - ( zydx(k,j  ,i) * zydx(k,j  ,i) / ( cw(k,j  ,i) + cw(k+1,j  ,i) )     &
               +   zydx(k,j-1,i) * zydx(k,j-1,i) / ( cw(k,j-1,i) + cw(k+1,j-1,i) )   ) &
               ! from j,k cross terms
               - ( qrt * zydx(k,j-1,i) - qrt * zydx(k,j,i)                           ) &  
               ! from i,j cross terms if lbc  
               - ( hlf * zxdy(k,j-1,i) * zydx(k,j-1,i) / (cw(k,j-1,i) + cw(k+1,j-1,i)) &
               -   hlf * zxdy(k,j  ,i) * zydx(k,j  ,i) / (cw(k,j  ,i) + cw(k+1,j  ,i)) )
          cA(6,k,j,i) = qrt * ( zxdy(k+1,j,i) + zxdy(k,j,i-1) )  ! couples with k+1 i-1
          cA(7,k,j,i) =                                                                &
               ! couples with i-1
               Arx(k,j,i)/dxu(j,i)                                                     &  
               ! topo terms
               - ( zxdy(k,j,i  ) * zxdy(k,j,i  ) / ( cw(k,j,i  ) + cw(k+1,j,i  ) )     &
               +   zxdy(k,j,i-1) * zxdy(k,j,i-1) / ( cw(k,j,i-1) + cw(k+1,j,i-1) )   ) &
               ! from i,k cross terms
               - ( qrt*zxdy(k,j,i-1) - qrt*zxdy(k,j,i)                               ) &
               ! from i,j cross terms if lbc    
               - ( hlf * zxdy(k,j,i-1) * zydx(k,j,i-1) / (cw(k,j,i-1) + cw(k+1,j,i-1)) &
               -   hlf * zxdy(k,j,i  ) * zydx(k,j,i  ) / (cw(k,j,i  ) + cw(k+1,j,i  )) )
          ! only for k==1, couples with j+1,i-1
          cA(5,k,j,i) = &
               + hlf * zxdy(k,j+1,i  ) * zydx(k,j+1,i  ) / ( cw(k,j+1,i  ) + cw(k+1,j+1,i  ))  &
               + hlf * zxdy(k,j  ,i-1) * zydx(k,j  ,i-1) / ( cw(k,j  ,i-1) + cw(k+1,j  ,i-1))   
          ! only for k==1, couples wit
          cA(8,k,j,i) = &
               - hlf * zxdy(k,j-1,i  ) * zydx(k,j-1,i  ) / (cw(k,j-1,i  ) + cw(k+1,j-1,i  )) & 
               - hlf * zxdy(k,j  ,i-1) * zydx(k,j  ,i-1) / (cw(k,j  ,i-1) + cw(k+1,j  ,i-1)) 

          !---------------!
          !- K = 2, nz-1 -! interior levels
          !---------------!
          do k = 2,nz-1 
             cA(2,k,j,i) =  cw(k,j,i)                    &          ! couples with k-1
                  - qrt * ( zxdy(k-1,j,i) - zxdy(k,j,i)) &
                  - qrt * ( zydx(k-1,j,i) - zydx(k,j,i)) 
             cA(3,k,j,i) =  qrt * ( zydx(k+1,j,i) + zydx(k,j-1,i))  ! couples with k+1 j-1
             cA(4,k,j,i) =  Ary(k,j,i) / dyv(j,i)                   ! couples with j-1
             cA(5,k,j,i) =- qrt * ( zydx(k-1,j,i) + zydx(k,j-1,i))  ! couples with k-1 j-1
             cA(6,k,j,i) =  qrt * ( zxdy(k+1,j,i) + zxdy(k,j,i-1))  ! Couples with k+1 i-1
             cA(7,k,j,i) =  Arx(k,j,i) / dxu(j,i)                   ! Couples with i-1
             cA(8,k,j,i) =- qrt * ( zxdy(k-1,j,i) + zxdy(k,j,i-1))  ! Couples with k-1 i-1
          enddo

          !----------!
          !- K = nz -! upper level
          !----------!
          k = nz
          cA(2,k,j,i) = cw(k,j,i)                                    ! couples with k-1
          cA(4,k,j,i) = Ary(k,j,i) / dyv(j,i) &                      ! couples with j-1
               + qrt * ( zydx(k,j-1,i) - zydx(k,j,i))
          cA(5,k,j,i) =- qrt * ( zydx(k-1,j,i) + zydx(k,j-1,i) )     ! couples with k-1 j-1
          cA(7,k,j,i) =  Arx(k,j,i)/dxu(j,i) &                       ! Couples with i-1
               - qrt * (- zxdy(k,j,i-1) + zxdy(k,j,i))
          cA(8,k,j,i) =- qrt * ( zxdy(k-1,j,i) + zxdy(k,j,i-1) )     ! Couples with k-1 i-1

       enddo ! j
    enddo ! i

    call fill_halo(lev,cA)

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
          cA(1,k,j,i) =                   &
               -cA(2,k,j,i)               &
               -cw(k+1,j,i)               &
               -cA(4,k,j,i)-cA(4,k,j+1,i) &
               -cA(7,k,j,i)-cA(7,k,j,i+1) &
               -cA(6,k-1,j,i+1)           &
               -cA(8,k,j,i)               &
               -cA(3,k-1,j+1,i)           &
               -cA(5,k,j,i)

       enddo ! j
    enddo ! i

    deallocate(Arx)
    deallocate(Ary)
    deallocate(Arz)
    deallocate(zxdy)
    deallocate(zydx)
    deallocate(cw)

    if (netcdf_output) then
       if (myrank==0) write(*,*)'       write cA in a netcdf file'
       call write_netcdf(grid(lev)%cA,vname='ca',netcdf_file_name='cA.nc',rank=myrank,iter=lev)
    endif

  end subroutine define_matrix

end module mg_define_matrix
