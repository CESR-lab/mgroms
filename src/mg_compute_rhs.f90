module mg_compute_rhs

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine compute_rhs(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw

    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:)  , pointer :: Arz

    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: dz,dzw
    real(kind=rp), dimension(:,:,:), pointer :: zy,zx
    real(kind=rp), dimension(:,:,:), pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:), pointer :: zxw,zyw
    real(kind=rp), dimension(:,:,:), pointer :: cw
    real(kind=rp), dimension(:,:,:), pointer :: uf,vf,wf,rhs

    !NG comment: constants in a mg_cst.f90 file ?
    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- compute rhs:'

    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))

    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             dz(k,j,i) = zw(k+1,j,i) - zw(k,j,i) !!  cell height at rho-points
          enddo
       enddo
    enddo

    allocate(dzw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          dzw(1,j,i) = zr(1,j,i) - zw(1,j,i) !!
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i) - zr(k-1,j,i) !!  cell height at w-points
          enddo
          dzw(nz+1,j,i) = zw(nz+1,j,i) - zr(nz,j,i) !!
       enddo
    enddo

    !! Cell widths
    allocate(dxu(0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          dxu(j,i) = hlf * ( dx(j,i) + dx(j,i-1) )
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = hlf * ( dy(j,i) + dy(j-1,i) )
       enddo
    enddo

    !!  Areas
    allocate(Arx(nz,0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          do k = 1,nz
             Arx(k,j,i) = qrt * (dz(k,j,i)+dz(k,j,i-1)) * ( dy(j,i)+dy(j,i-1) ) 
          enddo
       enddo
    enddo

    allocate(Ary(nz,ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             Ary(k,j,i) = qrt *( dz(k,j,i)+dz(k,j-1,i)) * (dx(j,i)+dx(j-1,i)) 
          enddo
       enddo
    enddo

    allocate(Arz(0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          Arz(j,i) = dx(j,i) * dy(j,i)
       enddo
    enddo

    !! Slopes in x- and y-direction defined at rho-points
    allocate(zx(nz,0:ny+1,0:nx+1))
    allocate(zy(nz,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             zy(k,j,i) = hlf * (zr(k,j+1,i)-zr(k,j-1,i)) / dy(j,i)
             zx(k,j,i) = hlf * (zr(k,j,i+1)-zr(k,j,i-1)) / dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(1,zy)
    call fill_halo(1,zx) 

    allocate(zxdy(nz,0:ny+1,0:nx+1))
    allocate(zydx(nz,0:ny+1,0:nx+1))
    do k = 1, nz
       zydx(k,:,:) = zy(k,:,:) * dx(:,:)
       zxdy(k,:,:) = zx(k,:,:) * dy(:,:)
    enddo

    allocate(zyw(nz+1,0:ny+1,0:nx+1))
    allocate(zxw(nz+1,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz+1
             zyw(k,j,i) = hlf * (zw(k,j+1,i)-zw(k,j-1,i)) / dy(j,i)
             zxw(k,j,i) = hlf * (zw(k,j,i+1)-zw(k,j,i-1)) / dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(1,zyw)
    call fill_halo(1,zxw)

    allocate(cw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             cw(k,j,i) = Arz(j,i) / dzw(k,j,i) * (one + zxw(k,j,i)*zxw(k,j,i) + zyw(k,j,i)*zyw(k,j,i))
          enddo
       enddo
    enddo


    !! Fluxes
    allocate(uf(nz  ,0:ny+1,0:nx+1))
    allocate(vf(nz  ,0:ny+1,0:nx+1))
    allocate(wf(nz+1,0:ny+1,0:nx+1))

   if (myrank==0) write(*,*)'- compute rhs: uf'

    do i = 1,nx+1  ! West  to East
       do j = 1,ny ! South to North

          k = 1 ! lower level

          uf(k,j,i) =  Arx(k,j,i)/dxu(j,i) * dxu(j,i) * u(i,j,k) &
               - qrt * ( &
               + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(i  ,j,k+1-1)   &
               + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(i-1,j,k+1-1) ) &
               
               -( &
               + zxdy(k,j,i  )*zxdy(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  )) &
               + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) ) * dxu(j,i)*u(i,j,k) &
               
               -( & 
               + zxdy(k,j,i  )*zydx(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  ))   &
               * hlf * (dyv(j,i  )*v(i,j,k ) + dyv(j+1,i  )*v(i,j+1,k ))   & 
               + zxdy(k,j,i-1)*zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))   &
               * hlf * (dyv(j,i-1)*v(i-1,j,k) + dyv(j+1,i-1)*v(i-1,j+1,k)) )

          do k = 2,nz-1 !interior levels

             uf(k,j,i) =  Arx(k,j,i)/dxu(j,i) * dxu(j,i) * u(i,j,k) &
                  - qrt * ( &
                  + zxdy(k,j,i  )*dzw(k  ,j,i  )*w(i  ,j,k  -1) &
                  + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(i  ,j,k+1-1) &
                  + zxdy(k,j,i-1)*dzw(k  ,j,i-1)*w(i-1,j,k  -1) &
                  + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(i-1,j,k+1-1) &
                  )  ! umask
          enddo

          k = nz !upper level

          uf(k,j,i) = Arx(k,j,i)/dxu(j,i) * dxu(j,i) * u(i,j,k) &
               - qrt * ( &
               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(i  ,j,k  -1) &
               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(i  ,j,k+1-1) &
               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(i-1,j,k  -1) &
               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(i-1,j,k+1-1) &
               )  ! umask

       enddo
    enddo

    call fill_halo(1,uf,lbc_null='u')

   if (myrank==0) write(*,*)'- compute rhs: vf'
    do i = 1,nx
       do j = 1,ny+1

          k = 1 !lower level

          vf(k,j,i) = Ary(k,j,i)/dyv(j,i) * dyv(j,i) * v(i,j,k) &
               - qrt * ( &
               + zydx(k,j  ,i)*dzw(k+1,j  ,i)*w(i,j  ,k+1-1) &
               + zydx(k,j-1,i)*dzw(k+1,j-1,i)*w(i,j-1,k+1-1) &
               ) &
               
               -( &
               + zydx(k,j  ,i)*zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i)) &
               + zydx(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) ) * dyv(j,i)*v(i,j,k) &
               
               - ( &
               + zxdy(k,j  ,i)*zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i))  &
               * hlf * (dxu(j  ,i)*u(i,j,k) + dxu(j  ,i+1)*u(i+1,j,k)) &
               + zxdy(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i))   &
               * hlf * (dxu(j-1,i)*u(i,j-1,k) + dxu(j-1,i+1)*u(i+1,j-1,k)) &
               ) 

          do k = 2,nz !interior levels

             vf(k,j,i) = Ary(k,j,i)/dyv(j,i) * dyv(j,i) * v(i,j,k) &
                  - qrt * ( &
                  + zydx(k,j  ,i)*dzw(k  ,j  ,i)*w(i,j  ,k  -1) &
                  + zydx(k,j  ,i)*dzw(k+1,j  ,i)*w(i,j  ,k+1-1) &
                  + zydx(k,j-1,i)*dzw(k  ,j-1,i)*w(i,j-1,k  -1) &
                  + zydx(k,j-1,i)*dzw(k+1,j-1,i)*w(i,j-1,k+1-1) &
                  )  !* vmask(j,i)

          enddo

          k = nz !upper level

          vf(k,j,i) =  Ary(k,j,i)/dyv(j,i) * dyv(j,i) * v(i,j,k) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i)*w(i,j  ,k  -1) &
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*w(i,j  ,k+1-1) &
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*w(i,j-1,k  -1) &
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*w(i,j-1,k+1-1) &
               ) !* vmask(j,i)

       enddo
    enddo

    call fill_halo(1,vf,lbc_null='v')

   if (myrank==0) write(*,*)'- compute rhs: wf'

    do i = 1,nx
       do j = 1,ny

          k = 1 !bottom

          wf(k,j,i) = zero

          do k = 2,nz !interior levels

             wf(k,j,i) = cw(k,j,i) * dzw(k,j,i) * w(i,j,k-1) &
                  
                  - qrt * ( &
                  + zxdy(k  ,j,i)*dxu(j,i  )*u(i  ,j,k  ) &
                  + zxdy(k  ,j,i)*dxu(j,i+1)*u(i+1,j,k  ) &
                  + zxdy(k-1,j,i)*dxu(j,i  )*u(i  ,j,k-1) &
                  + zxdy(k-1,j,i)*dxu(j,i+1)*u(i+1,j,k-1) ) &
                  
                  - qrt * ( &
                  + zydx(k  ,j,i)*dyv(j  ,i)*v(i,j  ,k  ) &
                  + zydx(k  ,j,i)*dyv(j+1,i)*v(i,j+1,k  ) &
                  + zydx(k-1,j,i)*dyv(j  ,i)*v(i,j  ,k-1) &
                  + zydx(k-1,j,i)*dyv(j+1,i)*v(i,j+1,k-1) )
          enddo

          k = nz+1 !surface

          wf(k,j,i) = cw(k,j,i) * dzw(k,j,i) * w(i,j,k-1)&
               
               - hlf *( &
               + zxdy(k-1,j,i)*dxu(j,i  )*u(i  ,j,k-1) &
               + zxdy(k-1,j,i)*dxu(j,i+1)*u(i+1,j,k-1) ) &
               
               - hlf *( &
               + zydx(k-1,j,i)*dyv(j  ,i)*v(i,j  ,k-1) &
               + zydx(k-1,j,i)*dyv(j+1,i)*v(i,j+1,k-1) )

       enddo
    enddo

    if (netcdf_output) then
       call write_netcdf(uf,vname='uf',netcdf_file_name='uf.nc',rank=myrank)
       call write_netcdf(vf,vname='vf',netcdf_file_name='vf.nc',rank=myrank)
       call write_netcdf(wf,vname='wf',netcdf_file_name='wf.nc',rank=myrank)
    endif

    !! Divergence
    rhs => grid(1)%b

    rhs(:,:,:) = zero

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz

             rhs(k,j,i) =                       &
                  + uf(k  ,j  ,i+1) - uf(k,j,i) &
                  + vf(k  ,j+1,i  ) - vf(k,j,i) &
                  + wf(k+1,j  ,i  ) - wf(k,j,i)

          enddo
       enddo
    enddo

    deallocate(uf)
    deallocate(vf)
    deallocate(wf)

    deallocate(Arx)
    deallocate(Ary)
    deallocate(Arz)
    deallocate(zxdy)
    deallocate(zydx)
    deallocate(cw)

    if (myrank==0) write(*,*)'- compute rhs (finish)'

  end subroutine compute_rhs

end module mg_compute_rhs
