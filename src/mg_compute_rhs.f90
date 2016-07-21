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

    real(kind=rp) :: Arx, Ary
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:), pointer :: cw
    real(kind=rp), dimension(:,:,:), pointer :: uf,vf,wf,rhs

    !NG comment: constants in a mg_cst.f90 file ?
    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

!    integer(kind=ip), save :: iter=0
!    iter = iter + 1

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- compute rhs:'

    !- Cell heights (computed in define matices)
    dzw => grid(1)%dzw

    !- Slopes in x- and y-direction defined at rho-points (computed in define matices)
    zxdy => grid(1)%zxdy
    zydx => grid(1)%zydx

    !- (computed in define matices)
    cw  => grid(1)%cw

    !--------------!
    !- DIVERGENCE -!
    !--------------!
    !- Divergence is stored in rhs array pointer
    !- It is calculated progressively using first uf, 
    !- then using vf and at the e wf
    rhs => grid(1)%b
    rhs(:,:,:) = zero

    !------!
    !- UF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: uf'

    uf => grid(1)%dummy3Dnz

    do i = 1,nx+1  ! West  to East
       do j = 1,ny ! South to North

          k = 1 ! lower level

          Arx = qrt                                                      * & 
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( dy(j,i) + dy(j,i-1) )

          uf(k,j,i) =  Arx * u(i,j,k) &
               - qrt * ( &
               + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(i  ,j,k+1-1)   &
               + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(i-1,j,k+1-1) ) &
               
               -( &
               + zxdy(k,j,i  )*zxdy(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  )) &
               + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) ) * &
               (hlf * ( dx(j,i) + dx(j,i-1) )) * u(i,j,k) &
               
               -( & 
               + zxdy(k,j,i) * zydx(k,j,i)/(cw(k,j,i  )+cw(k+1,j,i  ))   &
               * hlf * ( &
               hlf * ( dy(j  ,i) + dy(j-1,i) )*v(i,j,k ) + &
               hlf * ( dy(j+1,i) + dy(j  ,i) )*v(i,j+1,k ))   & 
               + zxdy(k,j,i-1) * zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))   &
               * hlf * ( &
               hlf * ( dy(j  ,i-1) + dy(j-1,i-1) )*v(i-1,j,k) + &
               hlf * ( dy(j+1,i-1) + dy(j  ,i-1) )*v(i-1,j+1,k)) )

          do k = 2,nz-1 !interior levels

             Arx = qrt                                                      * & 
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  ( dy(j,i) + dy(j,i-1) )

             uf(k,j,i) =  Arx * u(i,j,k) &
                  - qrt * ( &
                  + zxdy(k,j,i  )*dzw(k  ,j,i  )*w(i  ,j,k  -1) &
                  + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(i  ,j,k+1-1) &
                  + zxdy(k,j,i-1)*dzw(k  ,j,i-1)*w(i-1,j,k  -1) &
                  + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(i-1,j,k+1-1) &
                  )  ! umask
          enddo

          k = nz !upper level

          Arx = qrt                                                      * & 
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( dy(j,i) + dy(j,i-1) )

          uf(k,j,i) = Arx * u(i,j,k) &
               - qrt * ( &
               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(i  ,j,k  -1) &
               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(i  ,j,k+1-1) &
               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(i-1,j,k  -1) &
               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(i-1,j,k+1-1) &
               )  ! umask

       enddo
    enddo

    call fill_halo(1,uf,lbc_null='u')

!!$    if (netcdf_output) then
!!$       call write_netcdf(uf,vname='uf',netcdf_file_name='uf.nc',rank=myrank,iter=iter)
!!$    endif

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz

             rhs(k,j,i) = uf(k  ,j  ,i+1) - uf(k,j,i) 

          enddo
       enddo
    enddo

    uf => null()

    !------!
    !- VF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: vf'

    vf => grid(1)%dummy3Dnz

    do i = 1,nx
       do j = 1,ny+1

          k = 1 !lower level

          Ary = qrt * &
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( dx(j,i) + dx(j-1,i) )

          vf(k,j,i) = Ary * v(i,j,k) &
               - qrt * ( &
               + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(i,j  ,k+1-1) &
               + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(i,j-1,k+1-1) &
               ) &
               
               -( &
               + zydx(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i)) &
               + zydx(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) ) * &
               hlf * ( dy(j,i) + dy(j-1,i) )*v(i,j,k) &
               
               - ( &
               + zxdy(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i))  &
               * hlf * ( &
               hlf * ( dx(j,i  ) + dx(j,i-1) ) * u(i,j,k) + &
               hlf * ( dx(j,i+1) + dx(j,i  ) ) * u(i+1,j,k)) &
               + zxdy(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i))   &
               * hlf * ( &
               hlf * ( dx(j-1,i  ) + dx(j-1,i-1) ) * u(i,j-1,k) + &
               hlf * ( dx(j-1,i+1) + dx(j-1,i  ) ) * u(i+1,j-1,k)) &
               ) 

          do k = 2,nz-1 !interior levels

             Ary = qrt * &
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  ( dx(j,i) + dx(j-1,i) )

             vf(k,j,i) = Ary * v(i,j,k) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i)*w(i,j  ,k  -1) &
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(i,j  ,k+1-1) &
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i)*w(i,j-1,k  -1) &
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(i,j-1,k+1-1) &
                  )  !* vmask(j,i)

          enddo

          k = nz !upper level

          Ary = qrt * &
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( dx(j,i) + dx(j-1,i) )

          vf(k,j,i) =  Ary * v(i,j,k) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i)*w(i,j  ,k  -1) &
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*w(i,j  ,k+1-1) &
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*w(i,j-1,k  -1) &
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*w(i,j-1,k+1-1) &
               ) !* vmask(j,i)

       enddo
    enddo

    call fill_halo(1,vf,lbc_null='v')

!!$    if (netcdf_output) then
!!$       call write_netcdf(vf,vname='vf',netcdf_file_name='vf.nc',rank=myrank,iter=iter)
!!$    endif

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz

             rhs(k,j,i) =  rhs(k,j,i) + vf(k  ,j+1,i  ) - vf(k,j,i)

          enddo
       enddo
    enddo

    vf => null()

    !------!
    !- WF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: wf'

    wf => grid(1)%dummy3Dnzp

    do i = 1,nx
       do j = 1,ny

          k = 1 !bottom

          wf(k,j,i) = zero

          do k = 2,nz !interior levels

             wf(k,j,i) = cw(k,j,i) * dzw(k,j,i) * w(i,j,k-1) &
                  
                  - qrt * ( &
                  + zxdy(k  ,j,i)*hlf * ( dx(j,i  ) + dx(j,i-1) )*u(i  ,j,k  ) &
                  + zxdy(k  ,j,i)*hlf * ( dx(j,i+1) + dx(j,i  ) )*u(i+1,j,k  ) &
                  + zxdy(k-1,j,i)*hlf * ( dx(j,i  ) + dx(j,i-1) )*u(i  ,j,k-1) &
                  + zxdy(k-1,j,i)*hlf * ( dx(j,i+1) + dx(j,i  ) )*u(i+1,j,k-1) ) &
                  
                  - qrt * ( &
                  + zydx(k  ,j,i)*hlf * ( dy(j  ,i) + dy(j-1,i) )*v(i,j  ,k  ) &
                  + zydx(k  ,j,i)*hlf * ( dy(j+1,i) + dy(j  ,i) )*v(i,j+1,k  ) &
                  + zydx(k-1,j,i)*hlf * ( dy(j  ,i) + dy(j-1,i) )*v(i,j  ,k-1) &
                  + zydx(k-1,j,i)*hlf * ( dy(j+1,i) + dy(j  ,i) )*v(i,j+1,k-1) )
          enddo

          k = nz+1 !surface

          wf(k,j,i) = cw(k,j,i) * dzw(k,j,i) * w(i,j,k-1)&
               
               - hlf *( &
               + zxdy(k-1,j,i)*hlf * ( dx(j,i) + dx(j,i-1) )*u(i  ,j,k-1) &
               + zxdy(k-1,j,i)*hlf * ( dx(j,i+1) + dx(j,i) )*u(i+1,j,k-1) ) &
               
               - hlf *( &
               + zydx(k-1,j,i)*hlf * ( dy(j,i) + dy(j-1,i) )*v(i,j  ,k-1) &
               + zydx(k-1,j,i)*hlf * ( dy(j+1,i) + dy(j,i) )*v(i,j+1,k-1) )

       enddo
    enddo

!!$    if (netcdf_output) then
!!$       call write_netcdf(wf,vname='wf',netcdf_file_name='wf.nc',rank=myrank,iter=iter)
!!$    endif

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz

             rhs(k,j,i) = rhs(k,j,i) + wf(k+1,j  ,i  ) - wf(k,j,i)

          enddo
       enddo
    enddo

    wf => null()

    !- if (myrank==0) write(*,*)'- compute rhs (finish)'

  end subroutine compute_rhs

end module mg_compute_rhs
