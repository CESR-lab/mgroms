module mg_define_matrix

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

    real(kind=rp), parameter :: one=1._rp
    real(kind=rp), parameter :: eigh=one/8._8
    real(kind=rp), parameter :: qrt=0.25_rp
    real(kind=rp), parameter :: hlf=0.5_rp

contains
  !-------------------------------------------------------------------------     
  subroutine define_matrices(dx, dy, zr, zw)

    real(kind=rp), dimension(:,:)  , pointer, optional, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:,:), pointer, optional, intent(in) :: zr, zw

    integer(kind=ip)::  lev
    character(len = 16) :: filen
    real(kind=rp) :: cff

    real(kind=rp), dimension(:,:,:), pointer :: zrf
    real(kind=rp), dimension(:,:,:), pointer :: zwf

    real(kind=rp), dimension(:,:,:), pointer :: zrc
    real(kind=rp), dimension(:,:,:), pointer :: zwc

    integer(kind=ip) :: nz,ny,nx
    integer(kind=ip) :: nzf,nyf,nxf
    integer(kind=ip) :: nzc,nyc,nxc


    if (myrank==0) write(*,*)'- define matrix:'

    do lev = 1, nlevs

       nx=grid(lev)%nx
       ny=grid(lev)%ny
       nz=grid(lev)%nz

       if (lev == 1) then
          grid(lev)%zr(:,0:ny+1,0:nx+1) = zr
          grid(lev)%zw(:,0:ny+1,0:nx+1) = zw
       else
          nxf=grid(lev-1)%nx
          nyf=grid(lev-1)%ny
          nzf=grid(lev-1)%nz

          zrf => grid(lev-1)%zr
          zwf => grid(lev-1)%zw

          if (grid(lev)%gather == 1) then
             nxc= nx/grid(lev)%ngx
             nyc= ny/grid(lev)%ngy
             allocate(zrc(nz,-1:nyc+2,-1:nxc+2))
             allocate(zwc(nz+1,-1:nyc+2,-1:nxc+2))
          else
             nxc=nx
             nyc=ny
             zrc => grid(lev)%zr
             zwc => grid(lev)%zw
          endif

          if ((aggressive).and.(lev==1)) then

             write(*,*) ' define matrices (aggressive).and.(lev==1) not available !'
             STOP

          elseif (grid(lev)%nz == 1) then

             zrc(1,1:nyc+1,1:nxc+1) = qrt         * ( &
                  zrf(1,0:nyf:2,0:nxf:2)          + &
                  zrf(1,1:nyf+1:2,0:nxf:2)        + &
                  zrf(1,0:nyf:2,1:nxf+1:2)        + &
                  zrf(1,1:nyf+1:2,1:nxf+1:2)      )

             zwc(:,1:nyc+1,1:nxc+1) = qrt            * ( &
                  zwf(:,0:nyf:2,0:nxf:2)             + &
                  zwf(:,1:nyf+1:2,0:nxf:2)           + &
                  zwf(:,0:nyf:2,1:nxf+1:2)           + &
                  zwf(:,1:nyf+1:2,1:nxf+1:2)         )

          else

             zrc(:,1:nyf+1,1:nxf+1) = eigh * ( &
                  zrf(1:nz-1:2,0:nyf:2,0:nxf:2)     + zrf(2:nz:2,0:nyf:2,0:nxf:2)     + &
                  zrf(1:nz-1:2,1:nyf+1:2,0:nxf:2)   + zrf(2:nz:2,1:nyf+1:2,0:nxf:2)   + &
                  zrf(1:nz-1:2,0:nyf:2,1:nxf+1:2)   + zrf(2:nz:2,0:nyf:2,1:nxf+1:2)   + &
                  zrf(1:nz-1:2,1:nyf+1:2,1:nxf+1:2) + zrf(2:nz:2,1:nyf+1:2,1:nxf+1:2) )

             zwc(:,1:ny+1,1:nx+1) = qrt * ( &
                  zwf(0:nz:2,0:nyf:2,0:nxf:2)       + &
                  zwf(0:nz:2,1:nyf+1:2,0:nxf:2)     + &
                  zwf(0:nz:2,0:nyf:2,1:nxf+1:2)     + &
                  zwf(0:nz:2,1:nyf+1:2,1:nxf+1:2)   )

          end if

          if (grid(lev)%gather == 1) then
             call gather(lev,zrc,grid(lev)%zr,nzi=nz)
             call gather(lev,zwc,grid(lev)%zw,nzi=nz+1)

             deallocate(zrc)
             deallocate(zwc)
          endif

       endif


       call fill_halo_3D_nb(lev,grid(lev)%zr,nhi=2)
       call fill_halo_3D_nb(lev,grid(lev)%zw,nhi=2)

       call define_matrix(lev, dx, dy, grid(lev)%zr, grid(lev)%zw)

    enddo

    if (netcdf_output) then
       do lev=1,nlevs
          !          write(filen,'("cA_",i1,".nc")') lev
          !          call write_netcdf(grid(lev)%cA,vname='cA',netcdf_file_name=filen,rank=myrank)
          write(filen,'("msk_",i1,".nc")') lev
          call write_netcdf(grid(lev)%rmask*one,vname='msk',netcdf_file_name=filen,rank=myrank)
       enddo
    endif

  end subroutine define_matrices

!----------------------------------------
  subroutine define_matrix(lev, dx, dy, zr, zw)

    integer(kind=ip),intent(in):: lev
    real(kind=rp), dimension(:,:), pointer, intent(in) :: dx, dy
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

    integer(kind=ip):: l, k, j, i
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

    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i)    !!  cell height at rho-points
          enddo
       enddo
    enddo
    allocate(dzw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          dzw(1,j,i) = zr(1,j,i)-zw(1,j,i)        !!
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i)   !!  cell height at w-points
          enddo
          dzw(nz+1,j,i) = zw(nz+1,j,i)-zr(nz,j,i) !!
       enddo
    enddo

    !! Cell widths
    allocate(dxu(ny,nx+1))
    do i = 1,nx+1
       do j = 1,ny
          dxu(j,i) = hlf*(dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,nx))
    do i = 1,nx
       do j = 1,ny+1
          dyv(j,i) = hlf*(dy(j,i)+dy(j-1,i))
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
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             zy(k,j,i) = hlf*(zr(k,j+1,i)-zr(k,j-1,i))/dy(j,i)
             zx(k,j,i) = hlf*(zr(k,j,i+1)-zr(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo

    allocate(zxdy(nz,0:ny+1,0:nx+1))
    allocate(zydx(nz,0:ny+1,0:nx+1))
    do k = 1,nz
       zydx(k,:,:) = zy(k,:,:)*dx(:,:)
       zxdy(k,:,:) = zx(k,:,:)*dy(:,:)
    enddo

    allocate(zyw(nz+1,0:ny+1,0:nx+1))
    allocate(zxw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             zyw(k,j,i) = hlf*(zw(k,j+1,i)-zw(k,j-1,i))/dy(j,i)
             zxw(k,j,i) = hlf*(zw(k,j,i+1)-zw(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo

    allocate(cw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             cw(k,j,i) = Arz(j,i)/dzw(k,j,i) * (one + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i))
          enddo
       enddo
    enddo

    cA(:,:,:,:)=0._8

    !! interaction coeff with neighbours
    do i = 1,nx
       do j = 1,ny

          k = 1 !lower level
          cA(3,k,j,i) = qrt*(zydx(k+1,j,i) + zydx(k,j-1,i))                     !! couples with k+1 j-1
          cA(4,k,j,i) = ( Ary(k,j,i)/dyv(j,i) &                                 !! couples with j-1 
               -(zydx(k,j,i)*zydx(k,j,i)/(cw(k,j,i)+cw(k+1,j,i)) &
               + zydx(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i))) &   ! from j,k cross terms
               -(qrt*zydx(k,j-1,i) - qrt*zydx(k,j,i)) &                       ! from i,j cross terms if lbc   
               -(hlf*zxdy(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) & 
               - hlf*zxdy(k,j,i)*zydx(k,j,i)/(cw(k,j,i)+cw(k+1,j,i))))
          cA(6,k,j,i) = ( qrt*zxdy(k+1,j,i) + qrt*zxdy(k,j,i-1) )             !! couples with k+1 i-1
          cA(7,k,j,i) = ( Arx(k,j,i)/dxu(j,i) &                               !! couples with i-1
               -(zxdy(k,j,i)*zxdy(k,j,i)/(cw(k,j,i)+cw(k+1,j,i)) &
               + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))) &
               -(qrt*zxdy(k,j,i-1) - qrt*zxdy(k,j,i)) &                        ! from i,j cross terms if lbc  
               -(hlf*zxdy(k,j,i-1)*zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) & 
               - hlf*zxdy(k,j,i)*zydx(k,j,i)/(cw(k,j,i)+cw(k+1,j,i)) ))
          cA(5,k,j,i) = +hlf*zxdy(k,j+1,i)*zydx(k,j+1,i)/(cw(k,j+1,i)+cw(k+1,j+1,i)) &  !! only for k==1, couples with j+1,i-1
               +hlf*zxdy(k,j,i-1)*zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))               
          cA(8,k,j,i) =-hlf*zxdy(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) &   !! only for k==1, couples with j-1,i-1
               -hlf*zxdy(k,j,i-1)*zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))                                        

          do k = 2,nz-1 !interior levels
             cA(2,k,j,i) = cw(k,j,i)                              !! couples with k-1
             cA(3,k,j,i) = qrt*(zydx(k+1,j,i) + zydx(k,j-1,i))    !! couples with k+1 j-1
             cA(4,k,j,i) = Ary(k,j,i)/dyv(j,i)                    !! couples with j-1
             cA(5,k,j,i) =-qrt*(zydx(k-1,j,i) + zydx(k,j-1,i))    !! couples with k-1 j-1
             cA(6,k,j,i) = qrt*(zxdy(k+1,j,i) + zxdy(k,j,i-1))    !! Couples with k+1 i-1
             cA(7,k,j,i) = Arx(k,j,i)/dxu(j,i)                    !! Couples with i-1
             cA(8,k,j,i) =-qrt*(zxdy(k-1,j,i) + zxdy(k,j,i-1))    !! Couples with k-1 i-1
          enddo

          k = nz !upper level
          cA(2,k,j,i) = cw(k,j,i)                                  !! couples with k-1
          cA(4,k,j,i) = &                                          !! couples with j-1
               ( Ary(k,j,i)/dyv(j,i) - (-qrt*zydx(k,j-1,i) + qrt*zydx(k,j,i)) )
          cA(5,k,j,i) =-( qrt*zydx(k-1,j,i) + qrt*zydx(k,j-1,i) )  !! couples with k-1 j-1
          cA(7,k,j,i) =  &                                         !! Couples with i-1
               ( Arx(k,j,i)/dxu(j,i)-(-qrt*zxdy(k,j,i-1) + qrt*zxdy(k,j,i)) )
          cA(8,k,j,i) =-( qrt*zxdy(k-1,j,i) + qrt*zxdy(k,j,i-1) )  !! Couples with k-1 i-1
       enddo
    enddo

    !JM: we want to make it obsolete !
    call fill_halo(lev,cA)

    !! interaction coeff with itself
    do i = 1,nx
       do j = 1,ny

          k = 1 !lower level
          cA(1,k,j,i) = -cA(2,k+1,j,i)      &
               -cA(4,k,j,i)-cA(4,k,j+1,i)   &
               -cA(7,k,j,i)-cA(7,k,j,i+1)   &
               -cA(6,k,j,i)-cA(8,k+1,j,i+1) &
               -cA(3,k,j,i)-cA(5,k+1,j+1,i) &
               -cA(5,k,j,i)-cA(5,k,j-1,i+1) &
               -cA(8,k,j,i)-cA(8,k,j+1,i+1)

          do k = 2,nz-1 !interior levels
             !JM: check if diagonal elements => zero
             cA(1,k,j,i) =                     &
                  -cA(2,k,j,i)-cA(2,k+1,j,i)   &
                  -cA(4,k,j,i)-cA(4,k,j+1,i)   &
                  -cA(7,k,j,i)-cA(7,k,j,i+1)   &
                  -cA(6,k,j,i)-cA(6,k-1,j,i+1) &
                  -cA(8,k,j,i)-cA(8,k+1,j,i+1) & 
                  -cA(3,k,j,i)-cA(3,k-1,j+1,i) &
                  -cA(5,k,j,i)-cA(5,k+1,j+1,i)   
          enddo

          k = nz !upper level
          cA(1,k,j,i) = -cA(2,k,j,i)-cw(k+1,j,i) &
               -cA(4,k,j,i)-cA(4,k,j+1,i)-cA(7,k,j,i)-cA(7,k,j,i+1) &
               -cA(6,k-1,j,i+1)-cA(8,k,j,i)-cA(3,k-1,j+1,i)-cA(5,k,j,i)
       enddo
    enddo

    deallocate(Arx)
    deallocate(Ary)
    deallocate(Arz)
    deallocate(zxdy)
    deallocate(zydx)
    deallocate(cw)

  end subroutine define_matrix

!-------------------------------------------------------------------

end module mg_define_matrix
