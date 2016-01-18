module mg_define_matrix

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine define_matrices(dx, dy, zr, zw)

    real(kind=rp), dimension(:,:)  , allocatable, optional, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:,:), allocatable, optional, intent(in) :: zr, zw

    integer(kind=ip)::  lev

    lev = 1
    if (cmatrix == 'simple') then
       call define_matrix_simple(lev)
    elseif (cmatrix == 'real') then
       call define_matrix_real(lev, dx, dy, zr, zw)
    else
       stop
    endif

    do lev = 1,nlevs-1
       call coarsen_matrix(lev)
    enddo

  end subroutine define_matrices

  !----------------------------------------
  subroutine define_matrix_simple(lev)

    integer(kind=ip),intent(in):: lev

    character(len=15) :: matrixname

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

    real(kind=rp), dimension(:,:,:,:), pointer :: cA
    integer(kind=ip):: k, j, i
    real(kind=rp):: dxi, dyi, dzi
    integer(kind=ip):: nx, ny, nz
    integer(kind=ip):: nh

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    cA => grid(1)%cA ! check the syntax / lighten the writing

    dxi=1._8   !/dx
    dyi=1._8   !/dy
    dzi=1._8*16   !/dz

    matrixname='7points' ! or '15points'

    !extended loops will be a pain for the real matrix
    do i = 1-nh,nx+nh
       do j = 1-nh,ny+nh
          do k = 1,nz
             if(trim(matrixname)=='7points')then
! --- regular 7 points Laplacian ---
                cA(1,k,j,i) = 2._8*(-dxi*dxi-dyi*dyi-dzi*dzi)
                cA(2,k,j,i) = dzi*dzi
                cA(3,k,j,i) = 0.0_8
                cA(4,k,j,i) = dyi*dyi
                cA(5,k,j,i) = 0.0_8
                cA(6,k,j,i) = 0.0_8
                cA(7,k,j,i) = dxi*dxi
                cA(8,k,j,i) = 0.0_8
             endif
! --- extended stencil with diagonal coupling: better convergence rate ---
             if(trim(matrixname)=='15points')then
                cA(1,k,j,i) = 2._8*(-dxi*dxi-dyi*dyi-dzi*dzi)-4*(dxi*dzi+dyi*dzi)
                cA(2,k,j,i) = dzi*dzi
                cA(3,k,j,i) = 0.5*dyi*dzi
                cA(4,k,j,i) = dyi*dyi
                cA(5,k,j,i) = 0.5*dyi*dzi
                cA(6,k,j,i) = 0.5*dxi*dzi
                cA(7,k,j,i) = dxi*dxi
                cA(8,k,j,i) = 0.5*dxi*dzi
             endif
          enddo
          cA(1,nz,j,i) = cA(1,nz,j,i) - dzi*dzi 
          cA(1,1,j,i)  = cA(1,1,j,i)  + dzi*dzi 
       enddo
    enddo

  end subroutine define_matrix_simple

!----------------------------------------
  subroutine define_matrix_real(lev, dx, dy, zr, zw)

    integer(kind=ip),intent(in):: lev
    real(kind=rp), dimension(:,:)  , allocatable, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:,:), allocatable, intent(in) :: zr, zw

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

    real(kind=rp), dimension(:,:,:,:), pointer :: cA
    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz
    integer(kind=ip):: nh

    real(kind=rp) ::  zxu, zyu, zxv, zyv

    real(kind=rp), dimension(:,:)  , allocatable :: dxu, dyv
    real(kind=rp), dimension(:,:,:), allocatable :: Arx, Ary
    real(kind=rp), dimension(:,:)  , allocatable :: Arz
    real(kind=rp), dimension(:,:,:), allocatable :: dz
    real(kind=rp), dimension(:,:,:), allocatable :: dzw
    real(kind=rp), dimension(:,:,:), allocatable :: zy,zx
    real(kind=rp), dimension(:,:,:), allocatable :: zydx,zxdy
    real(kind=rp), dimension(:,:,:), allocatable :: zxw, zyw

    ! TODO NG
    ! zw,zr can change in time
    ! dx,dy constant in time
!!! I'm assuming that I'm getting zw,zr,dx,dy from outside this routine

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh
    cA => grid(1)%cA 

!ND
    call write_netcdf(zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank)
    call write_netcdf(zw,vname='zw',netcdf_file_name='zw.nc',rank=myrank)

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
          dzw(1,j,i) = zr(1,j,i) - zw(1,j,i)
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i) !!  cell height at w-points
          enddo
          dzw(nz+1,j,i) = zw(nz+1,j,i) - zr(nz,j,i)
       enddo
    enddo
    allocate(dxu(ny,nx+1))
    do i = 1,nx+1
       do j = 1,ny
          dxu(j,i) = 0.5*(dx(j,i) + dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,nx))
    do i = 1,nx
       do j = 1,ny+1
          dyv(j,i) = 0.5*(dy(j,i) + dy(j-1,i))
       enddo
    enddo

    !!  Areas
    allocate(Arx(nz,ny,nx+1))
    do i = 1,nx+1
       do j = 1,ny
          do k = 1,nz
             Arx(k,j,i) = 0.25*(dz(k,j,i)+dz(k,j,i-1))*(dy(j,i)+dy(j,i-1)) 
          enddo
       enddo
    enddo
    allocate(Ary(nz,ny+1,nx))
    do i = 1,nx
       do j = 1,ny+1
          do k = 1,nz
             Ary(k,j,i) = 0.25*(dz(k,j,i)+dz(k,j-1,i))*(dx(j,i)+dx(j-1,i)) 
          enddo
       enddo
    enddo
    allocate(Arz(ny,nx))
    do i = 1,nx
       do j = 1,ny
          Arz(j,i) = dx(j,i)*dy(j,i)
       enddo
    enddo

    !! Slopes in x- and y-direction defined at rho-points
    allocate(zx(nz+1,0:ny+1,0:nx+1))
    allocate(zy(nz+1,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             zy(k,j,i) = 0.5*(zr(k,j+1,i)-zr(k,j-1,i))/dy(j,i)
             zx(k,j,i) = 0.5*(zr(k,j,i+1)-zr(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo
    zy(nz+1,:,:) = 0._8
    zy(:,0,:)    = 0._8
    zy(:,ny+1,:) = 0._8
    zy(:,:,0)    = 0._8
    zy(:,:,nx+1) = 0._8

    zx(nz+1,:,:) = 0._8
    zx(:,0,:)    = 0._8
    zx(:,ny+1,:) = 0._8
    zx(:,:,0)    = 0._8
    zx(:,:,nx+1) = 0._8

    allocate(zxdy(nz+1,0:ny+1,0:nx+1))
    allocate(zydx(nz+1,0:ny+1,0:nx+1))
    do k = 1,nz
       zydx(k,:,:) = zy(k,:,:) *dx(:,:)
       zxdy(k,:,:) = zx(k,:,:) *dy(:,:)
    enddo

    allocate(zyw(nz+1,ny,nx))
    allocate(zxw(nz+1,ny,nx))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz+1
             zyw(k,j,i) = 0.5*(zw(k,j+1,i)-zw(k,j-1,i))/dy(j,i)
             zxw(k,j,i) = 0.5*(zw(k,j,i+1)-zw(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo

    !! cA : levels 2 to nz
    do i = 1,nx
       do j = 1,ny
          do k = 2,nz
             cA(2,k,j,i) = Arz(j,i)/dzw(k,j,i)*(1+zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i)) !! couples with k-1
             cA(3,k,j,i) = 0.25*zydx(k+1,j,i)+0.25*zydx(k,j-1,i)                               !! couples with k+1 j-1
             cA(4,k,j,i) = Ary(k,j,i)/dyv(j,i)                                                 !! couples with j-1
             cA(5,k,j,i) =-0.25*zydx(k-1,j,i)-0.25*zydx(k,j-1,i)                               !! couples with k-1 j-1
             cA(6,k,j,i) = 0.25*zxdy(k+1,j,i)+0.25*zxdy(k,j,i-1)                               !! Couples with k+1 i-1
             cA(7,k,j,i) = Arx(k,j,i)/dxu(j,i)                                                 !! Couples with i-1
             cA(8,k,j,i) =-0.25*zxdy(k-1,j,i)-0.25*zxdy(k,j,i-1)                               !! Couples with k-1 i-1
          enddo
          k = nz
          cA(4,k,j,i) = cA(4,k,j,i) - 0.25*zydx(k,j-1,i) + 0.25*zydx(k,j,i)
          cA(7,k,j,i) = cA(7,k,j,i) - 0.25*zxdy(k,j,i-1) + 0.25*zxdy(k,j,i)
       enddo
    enddo
!ND
    call fill_halo(1,cA) !fill_halo_4d in the future!
!
    do i = 1,nx
       do j = 1,ny
          do k = 2,nz-1
!ND
             cA(1,k,j,i) = -cA(2,k,j,i)-cA(2,k+1,j,i)-cA(4,k,j,i)-cA(4,k,j+1,i)-cA(7,k,j,i)-cA(7,k,j,i+1)
!             cA(1,k,j,i) = -cA(2,k,j,i)-cA(2,k+1,j,i)-cA(7,k,j,i)-cA(7,k,j,i+1) ! for comparing with matlab code 
          enddo
          k = nz
!ND
          cA(1,k,j,i) = -3*cA(2,k,j,i)-cA(4,k,j,i)-cA(4,k,j+1,i)-cA(7,k,j,i)-cA(7,k,j,i+1)
!          cA(1,k,j,i) = -3*cA(2,k,j,i)-cA(7,k,j,i)-cA(7,k,j,i+1) ! for comparing with matlab code 
       enddo
    enddo

    !! cA : bottom level
    k = 1
    do i = 1,nx
       do j = 1,ny
          zxv = 0.5_8 * (zx(k,j,i) + zx(k,j-1,i))
          zyv = 0.5_8 * (zy(k,j,i) + zy(k,j-1,i))
          zxu = 0.5_8 * (zx(k,j,i) + zx(k,j,i-1))
          zyu = 0.5_8 * (zy(k,j,i) + zy(k,j,i-1))
          cA(3,k,j,i) = 0.25*zydx(k+1,j,i)+0.25*zydx(k,j-1,i)                              !! couples with k+1 j-1
          cA(4,k,j,i) = Ary(k,j,i)/dyv(j,i)*(1 - 0.5*zyv*zyv/(1+zxv*zxv+zyv*zyv))          !! couples with j-1
          cA(6,k,j,i) = 0.25*zxdy(k+1,j,i)+0.25*zxdy(k,j,i-1)                              !! Couples with k+1 i-1
          cA(7,k,j,i) = Arx(k,j,i)/dxu(j,i)*(1 - 0.5*zxu*zxu/(1+zxu*zxu+zyu*zyu))          !! Couples with i-1
          cA(5,k,j,i) = 0.125*zx(k,j+1,i)*zy(k,j+1,i) + 0.125*zx(k,j,i-1)*zy(k,j,i-1)      !! only for k==1, couples with j+1,i-1
          cA(8,k,j,i) =-0.125*zx(k,j-1,i)*zy(k,j-1,i) - 0.125*zx(k,j,i-1)*zy(k,j,i-1)      !! only for k==1, couples with j-1,i-1
       enddo
    enddo

    call fill_halo(1,cA) !fill_halo_4d in the future!

    do i = 1,nx
       do j = 1,ny
!ND
          cA(1,k,j,i) = -cA(2,k+1,j,i)-cA(4,k,j,i)-cA(4,k,j+1,i)-cA(7,k,j,i)-cA(7,k,j,i+1)
!          cA(1,k,j,i) = -cA(2,k+1,j,i)-cA(7,k,j,i)-cA(7,k,j,i+1) ! for comparing with matlab code 
       enddo
    enddo

!ND
    call write_netcdf(cA,vname='cA',netcdf_file_name='cA.nc',rank=myrank)

    !!
    !!! For the moment, we will implement side bc by means of a buffer. The next phase will include a mask, whereupon we must
    !! implement horizontal condition in the matrix coefficients

    !! West Boundary 
    i = 1
    do j = 1,ny
       do k = 1,nz
          !         cA(2,k,j,i) = Arz(j,i)/dzw(k,j,i)*(1+0.5*zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i)) !! couples with k-1
       enddo
    enddo
    do j = 1,ny
       do k = 1,nz
          !         cA(1,k,j,i) = -cA(2,k,j,i)- cA(2,k+1,j,i)-cA(4,k,j,i)- cA(4,k,j+1,i)-cA(7,k,j,i+1)
          !         cA(2,k,j,i) = cA(2,k,j,i) - 0.25*zxdy(k-1,j,i) + 0.25*zxdy(k,j,i)
       enddo
    enddo
    !! East Boundary 
    i = nx
    do j = 1,ny
       do k = 1,nz
          !         cA(2,k,j,i) = Arz(j,i)/dzw(k,j,i)*(1+0.5*zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i)) !! couples with k-1
          ! &                   + 0.25*zxdy(k-1,j,i) - 0.25*zxdy(k,j,i)
       enddo
    enddo
    !! South Boundary 
    j = 1
    do i = 1,nx
       do k = 1,nz
          !         cA(2,k,j,i) = Arz(j,i)/dzw(k,j,i)*(1+zxw(k,j,i)*zxw(k,j,i)+0.5*zyw(k,j,i)*zyw(k,j,i)) !! couples with k-1
          ! &                   - 0.25*zydx(k-1,j,i) + 0.25*zydx(k,j,i)
       enddo
    enddo
    !! North Boundary 
    j = ny
    do i = 1,nx
       do k = 1,nz
          !         cA(2,k,j,i) = Arz(j,i)/dzw(k,j,i)*(1+zxw(k,j,i)*zxw(k,j,i)+0.5*zyw(k,j,i)*zyw(k,j,i)) !! couples with k-1
          ! &                   + 0.25*zydx(k-1,j,i) - 0.25*zydx(k,j,i)
       enddo
    enddo

    deallocate(Arx)
    deallocate(Ary)
    deallocate(Arz)
    deallocate(zxdy)
    deallocate(zydx)

  end subroutine define_matrix_real

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix(lev)
    integer(kind=ip),intent(in):: lev

    real(kind=rp),dimension(:,:,:,:),pointer :: Ac
    real(kind=rp),dimension(:,:,:,:),pointer :: Af

    integer(kind=ip) :: nx, ny, nz
    integer(kind=ip) :: l, nd

    nx = grid(lev+1)%nx
    ny = grid(lev+1)%ny
    nz = grid(lev+1)%nz

    ! the matrix on the fine grid
    Af => grid(lev)%cA

    if (grid(lev+1)%gather == 1) then
       Ac => grid(lev+1)%cAdummy
       nx = grid(lev+1)%nx / grid(lev+1)%ngx
       ny = grid(lev+1)%ny / grid(lev+1)%ngy
!       if(myrank == 0) write(*,*)"gather lev=",lev+1,"nx,ny,nz=",nx,ny,nz
    else
       Ac => grid(lev+1)%cA
!       if(myrank == 0) write(*,*)"F2C   lev=",lev+1,"nx,ny,nz=",nx,ny,nz
    endif


    if ((aggressive).and.(lev==1)) then
!       call coarsen_matrix_aggressive(lev)

    elseif (grid(lev+1)%nz == 1) then
       call coarsen_matrix_2D(Af,Ac,nx,ny,nz)
       ! fill the halo
       do l=1,3
          grid(lev+1)%r = grid(lev+1)%cA(l,:,:,:)
          call fill_halo(lev+1,grid(lev+1)%r)
          grid(lev+1)%cA(l,:,:,:) = grid(lev+1)%r
       enddo

    else
       call tic(lev,'coarsen_matrix_3D')
       call coarsen_matrix_3D(Af,Ac,nx,ny,nz)
       call toc(lev,'coarsen_matrix_3D')
    end if

    if (grid(lev+1)%gather == 1) then
       nd = size(Ac,1)       
       do l=1,nd
          grid(lev+1)%dummy3 = Ac(l,:,:,:)
          call gather(lev+1,grid(lev+1)%dummy3,grid(lev+1)%r)
          ! fill the halo
          call fill_halo(lev+1,grid(lev+1)%r)
          grid(lev+1)%cA(l,:,:,:) = grid(lev+1)%r
       enddo
    else
       ! fill the halo
       nd = size(Ac,1) 
       do l=1,nd
          grid(lev+1)%r = grid(lev+1)%cA(l,:,:,:)
          call fill_halo(lev+1,grid(lev+1)%r)
          grid(lev+1)%cA(l,:,:,:) = grid(lev+1)%r
       enddo
    endif


  end subroutine coarsen_matrix

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_aggressive(lev)
    integer(kind=ip),intent(in):: lev

    integer(kind=ip) :: idum
    idum = lev

    write(*,*)'Error: coarsen matrix aggressive not available yet !'
    stop -1
  end subroutine coarsen_matrix_aggressive


  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_2D(cA,cA2,nx2,ny2,nz2) ! from lev to lev+1

    integer(kind=ip):: nx2, ny2, nz2! on lev+1
    real(kind=rp), dimension(:,:,:,:), pointer :: cA,cA2

    integer(kind=ip):: k, j, i
    integer(kind=ip):: km, jm, im
    integer(kind=ip):: k2, j2, i2
    integer(kind=ip):: d

    real(kind=rp)   :: diag,cff


    k = 1
    km= 2
    k2= 1
    ! how many diagonal in the fine matrix? 3 or 8 ?
    d = size(cA,1) 

    ! I'm pretty sure it depends on whether d==3 or d==8

    if (d ==8) then
       cff = 1._8/16._8 ! check this value!!!
        ! fine matrix was 3D
       do i2 = 1,nx2
          i = 2*i2-1
          im = i+1
          do j2 = 1,ny2
             j = 2*j2-1
             jm = j+1     
             ! cA2(2,:,:,:) plays the role of cA(4,:,:,:)
             ! cA2(3,:,:,:) plays the role of cA(7,:,:,:)

             ! TODO: CHECK THESE FORMULA, I'm not completely sure

             cA2(2,k2,j2,i2) = cff*(cA(4,k,j,i)+cA(4,km,j,i)+cA(4,k,j,im)+cA(4,km,j,im))
             cA2(3,k2,j2,i2) = cff*(cA(7,k,j,i)+cA(7,km,j,i)+cA(7,k,jm,i)+cA(7,km,jm,i))
             ! 
             diag = cA(4,k,jm,i)+cA(4,km,jm,i)+cA(4,k,jm,im)+cA(4,km,jm,im)
             diag = cA(7,k,j,im)+cA(7,km,j,im)+cA(7,k,jm,im)+cA(7,km,jm,im) + diag

             diag = diag + diag
             !
             diag = cA(1,k,j,i) +cA(1,km,j,i) +cA(1,k,jm,i) +cA(1,km,jm,i) &
                  +cA(1,k,j,im)+cA(1,km,j,im)+cA(1,k,jm,im)+cA(1,km,jm,im) + diag

             cA2(1,k2,j2,i2) = cff*diag
          enddo
       enddo
    else
       cff = 1._8/4._8 ! check this value!!!
       ! fine matrix was already 2D
       do i2 = 1,nx2
          i = 2*i2-1
          im = i+1
          do j2 = 1,ny2
             j = 2*j2-1
             jm = j+1     
             cA2(2,k2,j2,i2) = cff*(cA(2,k,j,i)+cA(2,k,j,im))
             cA2(3,k2,j2,i2) = cff*(cA(3,k,j,i)+cA(3,k,jm,i))
             ! 
             diag = cA(2,k,jm,i)+cA(2,k,jm,im)
             diag = cA(3,k,jm,i)+cA(3,k,jm,im) + diag
             diag = diag + diag
             diag = cA(1,k,j,i) + cA(1,k,jm,i) + cA(1,k,j,im) + cA(1,k,jm,im) + diag
             cA2(1,k2,j2,i2) = cff*diag
          enddo
       enddo
    endif


  end subroutine coarsen_matrix_2D

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_3D(cA,cA2,nx2,ny2,nz2) ! from lev to lev+1

!    integer(kind=ip),intent(in):: lev
    integer(kind=ip):: nx2, ny2, nz2! on lev+1
    real(kind=rp), dimension(:,:,:,:), pointer :: cA,cA2
!    real(kind=rp), dimension(:,:,:) , pointer :: dummy3D

    integer(kind=ip):: k, j, i
    integer(kind=ip):: km, jm, im
    integer(kind=ip):: k2, j2, i2

    real(kind=rp)   :: diag,cff


!    cA  => grid(lev)%cA
!    cA2 => grid(lev+1)%cA
!    nx2 = grid(lev+1)%nx
!    ny2 = grid(lev+1)%ny
!    nz2 = grid(lev+1)%nz
!    nh  = grid(lev+1)%nh

    ! the coefficients should be rescaled with 1/16
    cff = 1._8/16._8

    do i2 = 1,nx2
       i = 2*i2-1
       im = i+1
       do j2 = 1,ny2
          j = 2*j2-1
          jm = j+1     
          do k2 = 1,nz2
             k = 2*k2-1
             km = k+1
             ! cA(2,:,:,:)      -> p(k-1,j,i)
             cA2(2,k2,j2,i2) = cff*(cA(2,k,j,i)+cA(2,k,jm,i)+cA(2,k,j,im)+cA(2,k,jm,im))
             ! cA(3,:,:,:)      -> p(k+1,j-1,i)
             cA2(3,k2,j2,i2) = cff*(cA(3,k,j,i)+cA(3,k,j,im))
             ! cA(4,:,:,:)      -> p(k,j-1,i)
             cA2(4,k2,j2,i2) = cff*(cA(4,k,j,i)+cA(4,km,j,i)+cA(4,k,j,im)+cA(4,km,j,im))
             ! cA(5,:,:,:)      -> p(k-1,j-1,i)
             cA2(5,k2,j2,i2) = cff*(cA(5,k,j,i)+cA(5,k,j,im))
             ! cA(6,:,:,:)      -> p(k+1,j,i-1)
             cA2(6,k2,j2,i2) = cff*(cA(6,k,j,i)+cA(6,k,jm,i))
             ! cA(7,:,:,:)      -> p(k,j,i-1)
             cA2(7,k2,j2,i2) = cff*(cA(7,k,j,i)+cA(7,km,j,i)+cA(7,k,jm,i)+cA(7,km,jm,i))
             ! cA(8,:,:,:)      -> p(k-1,j,i-1)
             cA2(8,k2,j2,i2) =cff*( cA(8,k,j,i)+cA(8,k,jm,i))

             ! the diagonal term is the sum of 48 terms ...             
             ! why?
             ! easy to see: the coarse cell (call it the box) is made of 8 fine cells
             ! take one fine cell, lay the 15 points on this cell
             ! count how many fine cells within this box are connected to it
             ! you should find 6
             ! multifly that by the number of fine cells

             ! here is the first 20
             diag = cA(2,km,j,i)+cA(2,km,jm,i)+cA(2,km,j,im)+cA(2,km,jm,im)
             diag = cA(3,k,jm,i)+cA(3,k,jm,im)                              + diag
             diag = cA(4,k,jm,i)+cA(4,km,jm,i)+cA(4,k,jm,im)+cA(4,km,jm,im) + diag
             diag = cA(5,km,j,im)+cA(5,km,jm,im)                            + diag
             diag = cA(6,k,j,im)+cA(6,k,jm,im)                              + diag
             diag = cA(7,k,j,im)+cA(7,km,j,im)+cA(7,k,jm,im)+cA(7,km,jm,im) + diag
             diag = cA(8,km,j,im)+cA(8,km,jm,im)                            + diag

             ! double that to account for symmetry of connections, we've now 40 terms
             diag = diag+diag

             ! add the 8 self-interacting terms
             diag = cA(1,k,j,i) +cA(1,km,j,i) +cA(1,k,jm,i) +cA(1,km,jm,i) &
                  +cA(1,k,j,im)+cA(1,km,j,im)+cA(1,k,jm,im)+cA(1,km,jm,im) + diag

             ! here we go!
             cA2(1,k2,j2,i2) = cff*diag
          enddo
       enddo
    enddo

    return

!!$!    if (myrank.eq.0)write(*,*)"coefficients computed"
!!$
!!$    ! fill the halo
!!$    ! the data should be contiguous in memory to use fill_halo 
!!$    ! no need to allocate an extra buffer
!!$    ! use the residual as a dummy variable
!!$    dummy3D => grid(lev+1)%r 
!!$
!!$    !- we should consider a specific fill_halo(4D) -!
!!$
!!$
!!$
!!$    do d = 1,8       
!!$       !!if (myrank.eq.0)write(*,*)"updating halo of coef(",d,",:,:,:)"
!!$       do i2 = 1,nx2
!!$          do j2 = 1,ny2
!!$             do k2 = 1,nz2
!!$                dummy3D(k2,j2,i2) = cA2(d,k2,j2,i2)
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       call fill_halo(lev+1,dummy3D)
!!$
!!$       do i2 = 1-nh,nx2+nh
!!$          do j2 = 1-nh,ny2+nh
!!$             do k2 = 1,nz2
!!$                cA2(d,k2,j2,i2) = dummy3D(k2,j2,i2) * cff
!!$             enddo
!!$          enddo
!!$       enddo
!!$       ! a way of improvement (only if it impacts perfs):
!!$       ! copy from cA2 to dummy3 only the interior ring used to fill the halo
!!$       ! copy from dummy to cA2 only the halo
!!$    enddo

!!$    do d = 1,8       
!!$       if (myrank.eq.0)write(*,*)"updating halo of coef(",d,",:,:,:)"
!!$
!!$       dummy3D(1:nx2,1:ny2,1:nz2) = cA2(d,1:nx2,1:ny2,1:nz2)
!!$
!!$       call fill_halo(lev+1, dummy3D)
!!$
!!$       cA2(d,:,:,:) =  dummy3D(:,:,:) * cff
!!$
!!$       ! a way of improvement (only if it impacts perfs):
!!$       ! copy from cA2 to dummy3 only the interior ring used to fill the halo
!!$       ! copy from dummy to cA2 only the halo
!!$    enddo

!    if (myrank.eq.0) write(*,*)"coarsening done"



  end subroutine coarsen_matrix_3D

end module mg_define_matrix
