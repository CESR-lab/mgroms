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
  subroutine compute_rhs(dx,dy,rmask,zr,zw,u,v,w)

    real(kind=rp), dimension(:,:)  , pointer, intent(in) :: dx,dy,rmask
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz, nh

    real(kind=rp), dimension(:,:)  ,   pointer :: umask,vmask
    real(kind=rp), dimension(:,:)  ,   pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:),   pointer :: Arx,Ary
    real(kind=rp), dimension(:,:)  ,   pointer :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: dz,dzw
    real(kind=rp), dimension(:,:,:),   pointer :: zy,zx
    real(kind=rp), dimension(:,:,:),   pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:),   pointer :: zxw,zyw
    real(kind=rp), dimension(:,:,:),   pointer :: cw
    real(kind=rp), dimension(:,:,:),   pointer :: uf,vf,wf,rhs

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh

    if (myrank==0) write(*,*)'- compute rhs:'

    !! umask and vmask
    allocate(umask(0:ny+1,0:nx+1))
    allocate(vmask(0:ny+1,0:nx+1))
!!$    becare
    do i = 1,nx+1
       do j = 1,ny+1
          umask(j,i) = rmask(j,i-1)*rmask(j,i)
          vmask(j,i) = rmask(j-1,i)*rmask(j,i)
       enddo
    enddo
!!$      call fill_halo(lev,umask)
!!$      call fill_halo(lev,vmask) 

  if (netcdf_output) then
     call write_netcdf(umask,vname='umask',netcdf_file_name='umask.nc',rank=myrank)
     call write_netcdf(vmask,vname='vmask',netcdf_file_name='vmask.nc',rank=myrank)
  endif
   
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
    allocate(dxu(0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          dxu(j,i) = 0.5_8*(dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = 0.5_8*(dy(j,i)+dy(j-1,i))
       enddo
    enddo

    !!  Areas
    allocate(Arx(nz,0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          do k = 1,nz
             Arx(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j,i-1))*(dy(j,i)+dy(j,i-1)) 
          enddo
       enddo
    enddo
    allocate(Ary(nz,ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             Ary(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j-1,i))*(dx(j,i)+dx(j-1,i)) 
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
             zy(k,j,i) = 0.5_8*(zr(k,j+1,i)-zr(k,j-1,i))/dy(j,i)
             zx(k,j,i) = 0.5_8*(zr(k,j,i+1)-zr(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(1,zy)
    call fill_halo(1,zx) 

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
             zyw(k,j,i) = 0.5_8*(zw(k,j+1,i)-zw(k,j-1,i))/dy(j,i)
             zxw(k,j,i) = 0.5_8*(zw(k,j,i+1)-zw(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(1,zyw)
    call fill_halo(1,zxw)

    allocate(cw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             cw(k,j,i) = Arz(j,i)/dzw(k,j,i) * (1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i))
          enddo
       enddo
    enddo


    !! Fluxes
    allocate(uf(nz,0:ny+1,0:nx+1))
    allocate(vf(nz,0:ny+1,0:nx+1))
    allocate(wf(1:nz+1,0:ny+1,0:nx+1))

    do i = 1,nx+1
       do j = 1,ny

          k = 1 !lower level

          uf(k,j,i) = ( Arx(k,j,i)*u(k,j,i) &

               - 0.25*( zxdy(k,j,i  )*dzw(k+1,j,i  )*w(k+1,j,i  )*rmask(j,i  ) &
                      + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(k+1,j,i-1)*rmask(j,i-1) ) &

               -( zxdy(k,j,i  )*zxdy(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  )) &
                + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) ) * dxu(j,i)*u(k,j,i) &

               -( zxdy(k,j,i  )*zydx(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  )) &
                * 0.5_8*(dyv(j,i  )*v(k,j,i  )*vmask(j,i  )+dyv(j+1,i  )*v(k,j+1,i  )*vmask(j+1,i  )) &
                + zxdy(k,j,i-1)*zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) &
                * 0.5_8*(dyv(j,i-1)*v(k,j,i-1)*vmask(j,i-1)+dyv(j+1,i-1)*v(k,j+1,i-1)*vmask(j+1,i-1)) ) ) * umask(j,i)

          do k = 2,nz-1 !interior levels

             uf(k,j,i) = ( Arx(k,j,i)*u(k,j,i) &

                  - 0.25*( zxdy(k,j,i  )*dzw(k  ,j,i  )*w(k  ,j,i  )*rmask(j,i  ) &
                         + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(k+1,j,i  )*rmask(j,i  ) &
                         + zxdy(k,j,i-1)*dzw(k  ,j,i-1)*w(k  ,j,i-1)*rmask(j,i-1) &
                         + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(k+1,j,i-1)*rmask(j,i-1) ) ) * umask(j,i)
          enddo

          k = nz !upper level

             uf(k,j,i) = ( Arx(k,j,i)*u(k,j,i) &

                  - 0.25*( zxdy(k,j,i  )*   dzw(k  ,j,i  )*w(k  ,j,i  )*rmask(j,i  ) &
                         + zxdy(k,j,i  )*2.*dzw(k+1,j,i  )*w(k+1,j,i  )*rmask(j,i  ) &
                         + zxdy(k,j,i-1)*   dzw(k  ,j,i-1)*w(k  ,j,i-1)*rmask(j,i-1) &
                         + zxdy(k,j,i-1)*2.*dzw(k+1,j,i-1)*w(k+1,j,i-1)*rmask(j,i-1) ) ) * umask(j,i)

       enddo
    enddo

    do i = 1,nx
       do j = 1,ny+1

          k = 1 !lower level

          vf(k,j,i) = ( Ary(k,j,i)*v(k,j,i) &

               - 0.25*( zydx(k,j  ,i)*dzw(k+1,j  ,i)*w(k+1,j  ,i)*rmask(j  ,i) &
                      + zydx(k,j-1,i)*dzw(k+1,j-1,i)*w(k+1,j-1,i)*rmask(j-1,i) ) &
               
               -( zydx(k,j  ,i)*zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i)) &
                + zydx(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) ) * dyv(j,i)*v(k,j,i) &

               -( zxdy(k,j  ,i)*zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i)) &
                * 0.5_8*(dxu(j  ,i)*u(k,j  ,i)*umask(j  ,i)+dxu(j  ,i+1)*u(k,j  ,i+1)*umask(j  ,i+1)) &
                + zxdy(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) &
                * 0.5_8*(dxu(j-1,i)*u(k,j-1,i)*umask(j-1,i)+dxu(j-1,i+1)*u(k,j-1,i+1)*umask(j-1,i+1)) ) ) * vmask(j,i)

          do k = 2,nz !interior levels

             vf(k,j,i) = ( Ary(k,j,i)*v(k,j,i) &

                  - 0.25*( zydx(k,j  ,i)*dzw(k  ,j  ,i)*w(k  ,j  ,i)*rmask(j  ,i) &
                         + zydx(k,j  ,i)*dzw(k+1,j  ,i)*w(k+1,j  ,i)*rmask(j  ,i) &
                         + zydx(k,j-1,i)*dzw(k  ,j-1,i)*w(k  ,j-1,i)*rmask(j-1,i) &
                         + zydx(k,j-1,i)*dzw(k+1,j-1,i)*w(k+1,j-1,i)*rmask(j-1,i) ) ) * vmask(j,i)

          enddo

          k = nz !upper level

             vf(k,j,i) = ( Ary(k,j,i)*v(k,j,i) &

                  - 0.25*( zydx(k,j  ,i)*   dzw(k  ,j  ,i)*w(k  ,j  ,i)*rmask(j  ,i) &
                         + zydx(k,j  ,i)*2.*dzw(k+1,j  ,i)*w(k+1,j  ,i)*rmask(j  ,i) &
                         + zydx(k,j-1,i)*   dzw(k  ,j-1,i)*w(k  ,j-1,i)*rmask(j-1,i) &
                         + zydx(k,j-1,i)*2.*dzw(k+1,j-1,i)*w(k+1,j-1,i)*rmask(j-1,i) ) ) * vmask(j,i)

       enddo
    enddo

    do i = 1,nx
       do j = 1,ny

          k = 1 !bottom

          wf(k,j,i) = 0._8

          do k = 2,nz !interior levels

             wf(k,j,i) = Arz(j,i)*w(k,j,i)*(1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i)) &

                  - 0.25*( zxdy(k  ,j,i)*dxu(j,i  )*u(k  ,j,i  )*umask(j,i  ) &
                         + zxdy(k  ,j,i)*dxu(j,i+1)*u(k  ,j,i+1)*umask(j,i+1) &
                         + zxdy(k-1,j,i)*dxu(j,i  )*u(k-1,j,i  )*umask(j,i  ) &
                         + zxdy(k-1,j,i)*dxu(j,i+1)*u(k-1,j,i+1)*umask(j,i+1) ) &

                  - 0.25*( zydx(k  ,j,i)*dyv(j  ,i)*v(k  ,j  ,i)*vmask(j  ,i) &
                         + zydx(k  ,j,i)*dyv(j+1,i)*v(k  ,j+1,i)*vmask(j+1,i) &
                         + zydx(k-1,j,i)*dyv(j  ,i)*v(k-1,j  ,i)*vmask(j  ,i) &
                         + zydx(k-1,j,i)*dyv(j+1,i)*v(k-1,j+1,i)*vmask(j+1,i) )
          enddo

          k = nz+1 !surface

          wf(k,j,i) = Arz(j,i)*w(k,j,i)*(1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i)) &

                  - 0.5 *( zxdy(k-1,j,i)*dxu(j,i  )*u(k-1,j,i  )*umask(j,i  ) &
                         + zxdy(k-1,j,i)*dxu(j,i+1)*u(k-1,j,i+1)*umask(j,i+1) ) &

                  - 0.5 *( zydx(k-1,j,i)*dyv(j  ,i)*v(k-1,j  ,i)*vmask(j  ,i) &
                         + zydx(k-1,j,i)*dyv(j+1,i)*v(k-1,j+1,i)*vmask(j+1,i) )

       enddo
    enddo


    !! Divergence
    rhs => grid(1)%b

    rhs(:,:,:) = 0._8

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz

             rhs(k,j,i) = uf(k,j,i+1)-uf(k,j,i) &
                        + vf(k,j+1,i)-vf(k,j,i) &
                        + wf(k+1,j,i)-wf(k,j,i)

          enddo
       enddo
    enddo

    deallocate(umask)
    deallocate(vmask)
    deallocate(Arx)
    deallocate(Ary)
    deallocate(Arz)
    deallocate(zxdy)
    deallocate(zydx)
    deallocate(cw)

  end subroutine compute_rhs

end module mg_compute_rhs
