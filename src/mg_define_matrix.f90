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
  subroutine define_matrices(dx, dy, zr, zw, umask, vmask, rmask)

    real(kind=rp), dimension(:,:)  , pointer, optional, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:)  , pointer, optional, intent(in) :: umask, vmask, rmask
    real(kind=rp), dimension(:,:,:), pointer, optional, intent(in) :: zr, zw

    integer(kind=ip)::  lev,ndf,ndc,i,j,nx,ny,nh
    character(len = 16) :: filen

    if (myrank==0) write(*,*)'- define matrix:'

    lev = 1
    if (cmatrix == 'simple') then
       call define_matrix_simple(lev, dx, dy, zr, zw, umask, vmask)
    elseif (cmatrix == 'real') then
       call define_matrix_real(lev, dx, dy, zr, zw, umask, vmask, rmask)
    else
       stop
    endif

    call define_loc()

    do lev = 1,nlevs-1
       ndf = size(grid(lev)%cA,1)
       ndc = size(grid(lev+1)%cA,1)
       if(myrank==0)write(*,'(A,I2,A,I2,A,I2)')'coarsening matrix from lev=',lev,' / nd fine=',ndf,' => nd coarse=',ndc

       nh=grid(lev)%nh
       nx=grid(lev)%nx
       ny=grid(lev)%ny
       do i=1-nh,nx+nh
          do j=1-nh,ny+nh
             if ( abs(grid(lev)%cA(1,1,j,i)) < 1d-12 ) then 
                grid(lev)%rmask(j,i)=0
                grid(lev)%cA(:,:,j,i)=0.
             else
                grid(lev)%rmask(j,i)=1
             endif
          enddo
       enddo
       call coarsen_matrix(lev)
!       if(myrank==0)write(*,*)"  / done"
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    enddo
    lev=nlevs
    nh=grid(lev)%nh
    nx=grid(lev)%nx
    ny=grid(lev)%ny
    do i=1-nh,nx+nh
       do j=1-nh,ny+nh
          if ( abs(grid(lev)%cA(1,1,j,i)) < 1d-12 ) then 
             grid(lev)%rmask(j,i)=0
             grid(lev)%cA(:,:,j,i)=0.
          else
             grid(lev)%rmask(j,i)=1
          endif
       enddo
    enddo


!!$    do lev=1,nlevs
!!$       nh=grid(lev)%nh
!!$       nx=grid(lev)%nx
!!$       ny=grid(lev)%ny
!!$       do i=1-nh,nx+nh
!!$          do j=1-nh,ny+nh
!!$             if ( abs(grid(lev)%cA(1,1,j,i)) < 1d-12 ) then 
!!$                grid(lev)%rmask(j,i)=0
!!$                grid(lev)%cA(:,:1,j,i)=0.
!!$             else
!!$                grid(lev)%rmask(j,i)=1
!!$             endif
!!$          enddo
!!$       enddo
!!$    enddo

    if (netcdf_output) then
       do lev=1,nlevs
          write(filen,'("cA_",i1,".nc")') lev
          call write_netcdf(grid(lev)%cA,vname='cA',netcdf_file_name=filen,rank=myrank)
          write(filen,'("msk_",i1,".nc")') lev
          call write_netcdf(grid(lev)%rmask*1._8,vname='msk',netcdf_file_name=filen,rank=myrank)

       enddo
    endif

  end subroutine define_matrices

  !----------------------------------------
  subroutine define_loc()
    !loc relates each matrix coefficient cA(l,:,:,:) to a 3D direction


!    cA(01)(k,j,i) -> p(k+0,j+0,i+0)  ***
!    cA(02)(k,j,i) -> p(k+1,j+0,i+0)  ***
!    cA(03)(k,j,i) -> p(k-1,j+0,i+0)  *** ! cA3(k+1,j,i) == cA2(k,j,i) ?

!    cA(04)(k,j,i) -> p(k+1,j+1,i+1)
!    cA(05)(k,j,i) -> p(k+1,j+1,i+0)
!    cA(06)(k,j,i) -> p(k+1,j+1,i-1)
!    cA(07)(k,j,i) -> p(k+1,j+0,i+1)
!    cA(08)(k,j,i) -> p(k+1,j+0,i-1)
!    cA(09)(k,j,i) -> p(k+1,j-1,i+1)
!    cA(10)(k,j,i) -> p(k+1,j-1,i+0)
!    cA(11)(k,j,i) -> p(k+1,j-1,i-1)

!    cA(12)(k,j,i) -> p(k+0,j+1,i+1)
!    cA(13)(k,j,i) -> p(k+0,j+1,i+0)
!    cA(14)(k,j,i) -> p(k+0,j+1,i-1)
!    cA(15)(k,j,i) -> p(k+0,j+0,i+1)
!    cA(16)(k,j,i) -> p(k+0,j+0,i-1)
!    cA(17)(k,j,i) -> p(k+0,j-1,i+1)
!    cA(18)(k,j,i) -> p(k+0,j-1,i+0)
!    cA(19)(k,j,i) -> p(k+0,j-1,i-1)

!    cA(20)(k,j,i) -> p(k-1,j+1,i+1)
!    cA(21)(k,j,i) -> p(k-1,j+1,i+0)
!    cA(22)(k,j,i) -> p(k-1,j+1,i-1)
!    cA(23)(k,j,i) -> p(k-1,j+0,i+1)
!    cA(24)(k,j,i) -> p(k-1,j+0,i-1)
!    cA(25)(k,j,i) -> p(k-1,j-1,i+1)
!    cA(26)(k,j,i) -> p(k-1,j-1,i+0)
!    cA(27)(k,j,i) -> p(k-1,j-1,i-1)

    integer(kind=ip):: di,dj,dk,l
    loc(:,:)=0
    loc(2,1)=1
    loc(3,1)=-1
    l=3
    do dk=-1,1
       do dj=-1,1
          do di=-1,1
             if (not( (di==0).and.(dj==0) ))then
                l = l+1
                loc(l,1)=-dk
                loc(l,2)=-dj
                loc(l,3)=-di
             endif
          enddo
       enddo
    enddo

  end subroutine define_loc

  !----------------------------------------
  subroutine define_matrix_simple(lev, dx, dy, zr, zw, umask, vmask)

    integer(kind=ip),intent(in):: lev
    real(kind=rp), dimension(:,:), pointer, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:), pointer, intent(in) :: umask, vmask
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: zr, zw

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
    real(kind=rp), dimension(:,:,:), pointer :: dz

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    cA => grid(1)%cA ! check the syntax / lighten the writing

    !ND
    if (netcdf_output) then
       call write_netcdf(zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank)
       call write_netcdf(zw,vname='zw',netcdf_file_name='zw.nc',rank=myrank)
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

    dxi=1._8/dx(1,1)
    dyi=1._8/dy(1,1)
    dzi=1._8/dz(1,1,1)

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

    if (netcdf_output) then
       call write_netcdf(cA,vname='cA',netcdf_file_name='cA.nc',rank=myrank)
    endif

  end subroutine define_matrix_simple

!----------------------------------------
  subroutine define_matrix_real(lev, dx, dy, zr, zw, umask, vmask, rmask)

    integer(kind=ip),intent(in):: lev
    real(kind=rp), dimension(:,:), pointer, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:), pointer, intent(in) :: umask, vmask, rmask
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

    cA => grid(1)%cA 

    !ND
    if (netcdf_output) then
       call write_netcdf(zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank)
       call write_netcdf(zw,vname='zw',netcdf_file_name='zw.nc',rank=myrank)
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
    allocate(dxu(ny,nx+1))
    do i = 1,nx+1
       do j = 1,ny
          dxu(j,i) = 0.5_8*(dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,nx))
    do i = 1,nx
       do j = 1,ny+1
          dyv(j,i) = 0.5_8*(dy(j,i)+dy(j-1,i))
       enddo
    enddo

    !!  Areas
    allocate(Arx(nz,ny,nx+1))
    do i = 1,nx+1
       do j = 1,ny
          do k = 1,nz
             Arx(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j,i-1))*(dy(j,i)+dy(j,i-1)) 
          enddo
       enddo
    enddo
    allocate(Ary(nz,ny+1,nx))
    do i = 1,nx
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

    cA(:,:,:,:)=0.
    !! interaction coeff with neighbours
    do i = 1,nx
       do j = 1,ny

          k = 1 !lower level
          cA(3,k,j,i) = ( 0.25_8*zydx(k+1,j,i) + 0.25_8*zydx(k,j-1,i) ) * vmask(j,i)      !! couples with k+1 j-1
          cA(4,k,j,i) = ( Ary(k,j,i)/dyv(j,i) &                                           !! couples with j-1
                                ! topo terms                                           
               -(zydx(k,j,i)*zydx(k,j,i)/(cw(k,j,i)+cw(k+1,j,i)) &
               + zydx(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i))) & 
                                ! from j,k cross terms
               -(0.25_8*zydx(k,j-1,i) - 0.25_8*zydx(k,j,i)) & 
                                ! from i,j cross terms if lbc                                         
               -(0.5_8*zxdy(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) &
               * (umask(j-1,i+1) - umask(j-1,i)) & 
               - 0.5_8*zxdy(k,j,i)*zydx(k,j,i)/(cw(k,j,i)+cw(k+1,j,i)) &
               * (umask(j,i+1) - umask(j,i))) ) * vmask(j,i)
          cA(6,k,j,i) = ( 0.25_8*zxdy(k+1,j,i) + 0.25_8*zxdy(k,j,i-1) ) * umask(j,i)      !! couples with k+1 i-1
          cA(7,k,j,i) = ( Arx(k,j,i)/dxu(j,i) &                                           !! couples with i-1
                                ! topo terms                                                   
               -(zxdy(k,j,i)*zxdy(k,j,i)/(cw(k,j,i)+cw(k+1,j,i)) &
               + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))) &
                                ! from i,k cross terms
               -(0.25_8*zxdy(k,j,i-1) - 0.25_8*zxdy(k,j,i)) &
                                ! from i,j cross terms if lbc                                         
               -(0.5_8*zxdy(k,j,i-1)*zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) &
               * (vmask(j+1,i-1) - vmask(j,i-1)) & 
               - 0.5_8*zxdy(k,j,i)*zydx(k,j,i)/(cw(k,j,i)+cw(k+1,j,i)) &
               * (vmask(j+1,i) - vmask(j,i))) ) * umask(j,i) 
          cA(5,k,j,i) = +0.5_8*zxdy(k,j+1,i)*zydx(k,j+1,i)/(cw(k,j+1,i)+cw(k+1,j+1,i)) &  !! only for k==1, couples with j+1,i-1
               * umask(j+1,i) * vmask(j+1,i) &
               +0.5_8*zxdy(k,j,i-1)*zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) &
               * umask(j,i) * vmask(j+1,i-1)              
          cA(8,k,j,i) =-0.5_8*zxdy(k,j-1,i)*zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) &   !! only for k==1, couples with j-1,i-1
               * umask(j-1,i) * vmask(j,i) &
               -0.5_8*zxdy(k,j,i-1)*zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) &
               * umask(j,i) * vmask(j,i-1)                                        

          do k = 2,nz-1 !interior levels
             cA(2,k,j,i) = cw(k,j,i) &                                                    !! couples with k-1
                                ! from i,k  cross terms if lbc
                  -(0.25_8*zxdy(k-1,j,i) - 0.25_8*zxdy(k,j,i)) &
                  * (umask(j,i+1) - umask(j,i)) &
                                ! from j,k  cross terms if lbc
                  -(0.25_8*zydx(k-1,j,i) - 0.25_8*zydx(k,j,i)) &
                  * (vmask(j+1,i) - vmask(j,i))
             cA(3,k,j,i) = (0.25_8*zydx(k+1,j,i) + 0.25_8*zydx(k,j-1,i)) * vmask(j,i)     !! couples with k+1 j-1
             cA(4,k,j,i) = Ary(k,j,i)/dyv(j,i) * vmask(j,i)                               !! couples with j-1
             cA(5,k,j,i) =-(0.25_8*zydx(k-1,j,i) + 0.25_8*zydx(k,j-1,i)) * vmask(j,i)     !! couples with k-1 j-1
             cA(6,k,j,i) = (0.25_8*zxdy(k+1,j,i) + 0.25_8*zxdy(k,j,i-1)) * umask(j,i)     !! Couples with k+1 i-1
             cA(7,k,j,i) = Arx(k,j,i)/dxu(j,i) * umask(j,i)                               !! Couples with i-1
             cA(8,k,j,i) =-(0.25_8*zxdy(k-1,j,i) + 0.25_8*zxdy(k,j,i-1)) * umask(j,i)     !! Couples with k-1 i-1
          enddo

          k = nz !upper level
          cA(2,k,j,i) = cw(k,j,i)                                                         !! couples with k-1
          cA(4,k,j,i) = ( Ary(k,j,i)/dyv(j,i) &                                           !! couples with j-1
               -(-0.25_8*zydx(k,j-1,i) + 0.25_8*zydx(k,j,i)) ) * vmask(j,i)
          cA(5,k,j,i) =-( 0.25_8*zydx(k-1,j,i) + 0.25_8*zydx(k,j-1,i) )  * vmask(j,i)     !! couples with k-1 j-1
          cA(7,k,j,i) = ( Arx(k,j,i)/dxu(j,i) &                                           !! Couples with i-1
               -(-0.25_8*zxdy(k,j,i-1) + 0.25_8*zxdy(k,j,i)) ) * umask(j,i)
          cA(8,k,j,i) =-( 0.25_8*zxdy(k-1,j,i) + 0.25_8*zxdy(k,j,i-1) )  * umask(j,i)     !! Couples with k-1 i-1
       enddo
    enddo

    call fill_halo(1,cA)

    !! interaction coeff with itself
    do i = 1,nx
       do j = 1,ny

          k = 1 !lower level
          cA(1,k,j,i) = -cA(2,k+1,j,i) &
               -cA(4,k,j,i)-cA(4,k,j+1,i)-cA(7,k,j,i)-cA(7,k,j,i+1) &
               -cA(6,k,j,i)-cA(8,k+1,j,i+1)-cA(3,k,j,i)-cA(5,k+1,j+1,i) &
               -cA(5,k,j,i)-cA(5,k,j-1,i+1)-cA(8,k,j,i)-cA(8,k,j+1,i+1)
          !for comparing with matlab 2d code 
          !          cA(1,k,j,i) = -cA(2,k+1,j,i) &
          !                        -cA(7,k,j,i)-cA(7,k,j,i+1) &
          !                        -cA(6,k,j,i)-cA(8,k+1,j,i+1) 

          do k = 2,nz-1 !interior levels
             cA(1,k,j,i) = -cA(2,k,j,i)-cA(2,k+1,j,i) &
                  -cA(4,k,j,i)-cA(4,k,j+1,i)-cA(7,k,j,i)-cA(7,k,j,i+1) &
                  -cA(6,k,j,i)-cA(6,k-1,j,i+1)-cA(8,k,j,i)-cA(8,k+1,j,i+1) & 
                  -cA(3,k,j,i)-cA(3,k-1,j+1,i)-cA(5,k,j,i)-cA(5,k+1,j+1,i)   
             !for comparing with matlab 2d code
             !             cA(1,k,j,i) = -cA(2,k,j,i)-cA(2,k+1,j,i) &
             !                           -cA(7,k,j,i)-cA(7,k,j,i+1) &
             !                           -cA(6,k,j,i)-cA(6,k-1,j,i+1)-cA(8,k,j,i)-cA(8,k+1,j,i+1) 
          enddo

          k = nz !upper level
          cA(1,k,j,i) = -cA(2,k,j,i)-cw(k+1,j,i) &
               -cA(4,k,j,i)-cA(4,k,j+1,i)-cA(7,k,j,i)-cA(7,k,j,i+1) &
               -cA(6,k-1,j,i+1)-cA(8,k,j,i)-cA(3,k-1,j+1,i)-cA(5,k,j,i)
          !for comparing with matlab 2d code 
          !          cA(1,k,j,i) = -cA(2,k,j,i)-cw(k+1,j,i) &
          !                        -cA(7,k,j,i)-cA(7,k,j,i+1) &
          !                        -cA(6,k-1,j,i+1)-cA(8,k,j,i)
       enddo
    enddo

    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             do l = 1,8
                cA(l,k,j,i)=cA(l,k,j,i)*rmask(j,i)
             enddo
          enddo
       enddo
    enddo

    call fill_halo(lev,cA)

    if (netcdf_output) then
       call write_netcdf(cA,vname='cA',netcdf_file_name='cA.nc',rank=myrank)
    endif
    

    deallocate(Arx)
    deallocate(Ary)
    deallocate(Arz)
    deallocate(zxdy)
    deallocate(zydx)
    deallocate(cw)

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
    nd = size(grid(lev)%cA,1) 

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

    if ((trim(interp_type)=='nearest') .and. (trim(restrict_type)=='avg')) then

       if ((aggressive).and.(lev==1)) then
          !       call coarsen_matrix_aggressive(lev)

       elseif (grid(lev+1)%nz == 1) then
          call coarsen_matrix_2D_nearest_avg(lev,Af,Ac,nx,ny,nz)

       else
          call tic(lev,'coarsen_matrix_3D')
          call coarsen_matrix_3D_nearest_avg(lev,Af,Ac,nx,ny,nz)
          call toc(lev,'coarsen_matrix_3D')
       end if

    elseif (( trim(interp_type)=='linear') .and. (trim(restrict_type)=='avg')) then

       if ((aggressive).and.(lev==1)) then
          !       call coarsen_matrix_aggressive(lev)

       elseif (grid(lev+1)%nz == 1) then
          call coarsen_matrix_2D_linear_avg(Af,Ac,nx,ny,nz)

       else
          call tic(lev,'coarsen_matrix_3D')
          if (lev==1) then
                call coarsen_matrix_3D_linear_avg_8(lev,Af,Ac,nx,ny,nz)
          else
                call coarsen_matrix_3D_linear_avg_27(lev,Af,Ac,nx,ny,nz)
          endif
!          Ac = Ac /2
          call toc(lev,'coarsen_matrix_3D')
       end if

    elseif (( trim(interp_type)=='nearest') .and. (trim(restrict_type)=='linear')) then
       ! todo if we think it's important

    endif


    if (grid(lev+1)%gather == 1) then
       !ND
!       nd = size(Ac,1)       
       do l=1,nd
          grid(lev+1)%dummy3 = Ac(l,:,:,:)
          call gather(lev+1,grid(lev+1)%dummy3,grid(lev+1)%r)
          !          ! fill the halo
          !          call fill_halo(lev+1,grid(lev+1)%r)
          grid(lev+1)%cA(l,:,:,:) = grid(lev+1)%r
       enddo
    endif

    call fill_halo(lev+1,grid(lev+1)%cA)

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
  subroutine coarsen_matrix_2D_nearest_avg(lev,cA,cA2,nx2,ny2,nz2) ! from lev to lev+1

    integer(kind=ip),intent(in):: lev
    integer(kind=ip):: nx2, ny2, nz2! on lev+1
    real(kind=rp), dimension(:,:,:,:), pointer :: cA,cA2

    integer(kind=ip):: k, j, i
    integer(kind=ip):: km, jm, im
    integer(kind=ip):: k2, j2, i2
    integer(kind=ip):: nd

    real(kind=rp)   :: diag,cff,cffm

    k = 1
    km= 2
    k2= 1
    ! how many diagonal in the fine matrix? 3 or 8 ?
    nd = size(cA,1) 

!    if(myrank>=0)write(*,*)"nd=",nd
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if (nd ==8) then

! -------- 2D ------------ <=  -------- 3D ------------ 
! cA(1,k,j,i)*p(k,j  ,i  ) <=  cA(1,k,j,i)*p(k,j  ,i  )
! cA(2,k,j,i)*p(k,j-1,i  ) <=  cA(4,k,j,i)*p(k,j-1,i  )
! cA(3,k,j,i)*p(k,j  ,i-1) <=  cA(7,k,j,i)*p(k,j  ,i-1)
! cA(4,k,j,i)*p(k,j-1,i-1) <=  cA(8,k,j,i)*p(k,j-1,i-1)
! cA(5,k,j,i)*p(k,j+1,i-1) <=  cA(5,k,j,i)*p(k,j+1,i-1)

       cff = 1._8/16._8 
       cffm= 1.
        ! fine matrix was 3D
       do i2 = 1,nx2
          i = 2*i2-1
          im = i+1
          do j2 = 1,ny2
             j = 2*j2-1
             jm = j+1     
             ! ex cA(4,k2,j2,i2)
             cA2(2,k2,j2,i2) = cff*(cA(4,k,j,i) +cA(4,km,j,i)*cffm &
                                   +cA(4,k,j,im)+cA(4,km,j,im)*cffm &
                                   +cA(5,km,j,i)*cffm+cA(5,km,j,im)*cffm & 
                                   +cA(3,k ,j,i)+cA(3,k ,j,im) ) 
             ! ex cA2(7,k2,j2,i2)
             cA2(3,k2,j2,i2) = cff*(cA(7,k,j,i)+cA(7,km,j,i)*cffm &
                                   +cA(7,k,jm,i)+cA(7,km,jm,i)*cffm &
                                   +cA(8,km,j,i)*cffm+cA(8,km,jm,i)*cffm &
                                   +cA(6,k ,j,i)+cA(6,k ,jm,i) )
             ! ex cA2(8,k2,j2,i2)
             cA2(4,k2,j2,i2) = cff*cA(8,k,j,i) 
             ! ex cA2(5,k2,j2,i2)
             cA2(5,k2,j2,i2) = cff*cA(5,k,jm,i) 
             !
             diag = 0._8
             
             diag = (cA(2,km,j,i)+cA(2,km,jm,i)+cA(2,km,j,im)+cA(2,km,jm,im))*cffm + diag
             diag = (cA(5,km,jm,i)+cA(5,km,jm,im))*cffm                            + diag
             diag = (cA(8,km,j,im)+cA(8,km,jm,im))*cffm                            + diag
             diag = cA(3,k,jm,i)+cA(3,k,jm,im)                              + diag
             diag = cA(6,k,j,im)+cA(6,k,jm,im)                              + diag
             diag = cA(4,k,jm,i)+cA(4,km,jm,i)*cffm+cA(4,k,jm,im)+cA(4,km,jm,im)*cffm + diag
             diag = cA(7,k,j,im)+cA(7,km,j,im)*cffm+cA(7,k,jm,im)+cA(7,km,jm,im)*cffm + diag

!             if (k2 == 1) then
             ! add the horizontal diagonal connections in the bottom level
             diag = cA(5,k,j,im)  + diag
             diag = cA(8,k,jm,im) + diag
!             endif

             ! double that to account for symmetry of connections, we've now 40 terms
             diag = diag+diag

             ! add the 8 self-interacting terms
             diag = cA(1,k,j,i) +cA(1,km,j,i)*cffm +cA(1,k,jm,i) +cA(1,km,jm,i)*cffm &
                   +cA(1,k,j,im)+cA(1,km,j,im)*cffm+cA(1,k,jm,im)+cA(1,km,jm,im)*cffm + diag

             ! here we go!
             cA2(1,k2,j2,i2) = cff*diag
             
          enddo
       enddo
    elseif(nd==5)then
!       write(*,*)"coarsening 5 diags => 5 diags"
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       cff = 1._8/8._8 ! 
       ! fine matrix was already 2D
       do i2 = 1,nx2
          i = 2*i2-1
          im = i+1
          do j2 = 1,ny2
             j = 2*j2-1
             jm = j+1     
             cA2(2,k2,j2,i2) = cff*(cA(2,k,j,i)+cA(2,k,j,im)+cA(4,k,j,im)+cA(5,k,j-1,im))
             cA2(3,k2,j2,i2) = cff*(cA(3,k,j,i)+cA(3,k,jm,i)+cA(4,k,jm,i)+cA(5,k,j,i))
             cA2(4,k2,j2,i2) = cff*(cA(4,k,j,i))
             cA2(5,k2,j2,i2) = cff*(cA(5,k,jm,i))
             ! 
             diag = cA(2,k,jm,i)+cA(2,k,jm,im)
             diag = cA(3,k,jm,i)+cA(3,k,jm,im) + diag
             diag = cA(4,k,jm,im)+cA(5,k,j,im) + diag
             diag = diag + diag
             diag = cA(1,k,j,i) + cA(1,k,jm,i) + cA(1,k,j,im) + cA(1,k,jm,im) + diag
             cA2(1,k2,j2,i2) = cff*diag
          enddo
       enddo
    else
       write(*,*)"2D matrix is supposed to have 5 diagonals, not 3"
       stop -1
    endif

  end subroutine coarsen_matrix_2D_nearest_avg

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_2D_linear_avg(cA,cA2,nx2,ny2,nz2) ! from lev to lev+1

    integer(kind=ip):: nx2, ny2, nz2! on lev+1
    real(kind=rp), dimension(:,:,:,:), pointer :: cA,cA2

  end subroutine coarsen_matrix_2D_linear_avg

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_3D_nearest_avg(lev,cA,cA2,nx2,ny2,nz2) ! from lev to lev+1

    integer(kind=ip),intent(in):: lev
    integer(kind=ip):: nx2, ny2, nz2! on lev+1
    real(kind=rp), dimension(:,:,:,:), pointer :: cA,cA2
    !    real(kind=rp), dimension(:,:,:) , pointer :: dummy3D

    integer(kind=ip):: k, j, i
    integer(kind=ip):: km, jm, im, kp
    integer(kind=ip):: k2, j2, i2, sm

    real(kind=rp)   :: diag,cff,cffm,cffp
    !!real(kind=rp),dimension(0:4)::c3

    integer(kind=1),dimension(:,:),pointer::msk

    !!data c3/0.,1.,0.5,0.3333333333333333333333333333,0.25/
    msk => grid(lev)%rmask
    !    cA  => grid(lev)%cA
    !    cA2 => grid(lev+1)%cA
    !    nx2 = grid(lev+1)%nx
    !    ny2 = grid(lev+1)%ny
    !    nz2 = grid(lev+1)%nz
    !    nh  = grid(lev+1)%nh

    !    write(*,*)'coarsening the matrix'

    ! the coefficients should be rescaled with 1/16
    cff = 1._8/16._8
    cffp = 1._8/16._8
    do i2 = 1,nx2
       i = 2*i2-1
       im = i+1
       do j2 = 1,ny2
          j = 2*j2-1
          jm = j+1   
          do k2 = 1,nz2
             k = 2*k2-1
             km = k+1
             kp = k-1
             cffm=1.
             !if(k2==nz2)cffm=0.5!cffm=(1.-4./8**lev)!0.5_8
             if((lev==2).and.(k2==nz2))cffm=2.

             !             cffp = 1._8/16._8
             !             if(k2==nz2)cffp=1._8/12._8


             !gr: we may be tempted to not define cA2(2,...) cA2(5,...) and cA2(8,...) at k2=1 
             !    because their fine grid correspond are not defined (they point downward at the bottom level)
             !    it's actually not a pb because these coefficients (at k2=1) are never used
             !
             !    second thought: we can't define cA2(2,1,j2,i2) because cA(3,...) and cA(6,...) ain't defined

             if(k2.gt.1)then
                ! cA(2,:,:,:)      -> p(k-1,j,i)
                sm = msk(j,i)+msk(jm,i)+msk(j,im)+msk(jm,im)
                !                cff = c3(sm)/4._8

                cA2(2,k2,j2,i2) = cff*(cA(2,k ,j ,i )+cA(2,k ,jm, i) &
                     +cA(2,k ,j ,im)+cA(2,k ,jm,im) &
                     +cA(5,k ,jm,i )+cA(5,k ,jm,im) &
                     +cA(8,k ,j ,im)+cA(8,k ,jm,im) &
                     +cA(3,kp,jm,i )+cA(3,kp,jm,im) &
                     +cA(6,kp,j ,im)+cA(6,kp,jm,im) )

                ! cA(5,:,:,:)      -> p(k-1,j-1,i)

                sm = msk(j,i)+msk(j,im)
                !                cff = c3(sm)/8._8

                cA2(5,k2,j2,i2) = cff*(cA(5,k,j,i)+cA(5,k,j,im)) 
                ! cA(8,:,:,:)      -> p(k-1,j,i-1)

                sm = msk(j,i)+msk(jm,i)
                !                cff = c3(sm)/8._8
                cA2(8,k2,j2,i2) = cff*(cA(8,k,j,i)+cA(8,k,jm,i)) 
             else ! bottom level = horizontal diagonal cross terms
                ! cA(5,:,:,:)      -> p(k,j+1,i-1)
                ! cA(8,:,:,:)      -> p(k,j-1,i-1)
                cA2(2,k2,j2,i2) = 0._8
                !             cff = c3(msk(jm,i))/16._8
                cA2(5,k2,j2,i2) = cff*cA(5,k,jm,i) 
                !             cff = c3(msk(j,i))/16._8
                cA2(8,k2,j2,i2) = cff*cA(8,k,j,i) 
             endif

             if (k2 < nz2) then
                ! cA(3,:,:,:)      -> p(k+1,j-1,i)

                sm = msk(j,i)+msk(j,im)
                !                cff = c3(sm)/8._8
                cA2(3,k2,j2,i2) = cffp*(cA(3,km,j,i)+cA(3,km,j,im))*cffm
                ! cA(6,:,:,:)      -> p(k+1,j,i-1)

                sm = msk(j,i)+msk(jm,i)
                !                cff = c3(sm)/8._8             
                cA2(6,k2,j2,i2) = cffp*(cA(6,km,j,i)+cA(6,km,jm,i)) *cffm
             else
                cA2(3,k2,j2,i2) = 0._8
                cA2(6,k2,j2,i2) = 0._8
             endif
             ! cA(4,:,:,:)      -> p(k,j-1,i)

             sm = msk(j,i)+msk(j,im)
             !                cff = c3(sm)/8._8
             cA2(4,k2,j2,i2) = cff*(cA(4,k,j,i) +cA(4,km,j,i)*cffm &
                  +cA(4,k,j,im)+cA(4,km,j,im)*cffm &
                  +cA(5,km,j,i)*cffm+cA(5,km,j,im)*cffm & 
                  +cA(3,k ,j,i)+cA(3,k ,j,im) )              

             ! cA(7,:,:,:)      -> p(k,j,i-1)

             sm = msk(j,i)+msk(jm,i)
             !                cff = c3(sm)/8._8             
             cA2(7,k2,j2,i2) = cff*(cA(7,k,j,i)+cA(7,km,j,i)*cffm &
                  +cA(7,k,jm,i)+cA(7,km,jm,i)*cffm &
                  +cA(8,km,j,i)*cffm+cA(8,km,jm,i)*cffm &
                  +cA(6,k ,j,i)+cA(6,k ,jm,i) )

             if(k2==1)then
                cA2(4,k2,j2,i2) = cA2(4,k2,j2,i2) &
                     +cff*( cA(5,k,jm,i) +cA(8,k,j,i) )
                cA2(7,k2,j2,i2) = cA2(7,k2,j2,i2)&
                     +cff*( cA(5,k,j-1,im) +cA(8,k,j,im) )
             endif


             ! the diagonal term is the sum of 48 terms ...             
             ! why?
             ! easy to see: the coarse cell (call it the box) is made of 8 fine cells
             ! take one fine cell, lay the 15 points on this cell
             ! count how many fine cells within this box are connected to it
             ! you should find 6
             ! multifly that by the number of fine cells

             ! here is the first 20

             diag = 0._8

             diag = (cA(2,km,j,i)+cA(2,km,jm,i)+cA(2,km,j,im)+cA(2,km,jm,im))*cffm + diag
             diag = (cA(5,km,jm,i)+cA(5,km,jm,im))*cffm                            + diag!bug fixed in cA5
             diag = (cA(8,km,j,im)+cA(8,km,jm,im))*cffm                            + diag
             diag = cA(3,k,jm,i)+cA(3,k,jm,im)                              + diag
             diag = cA(6,k,j,im)+cA(6,k,jm,im)                              + diag
             diag = cA(4,k,jm,i)+cA(4,km,jm,i)*cffm+cA(4,k,jm,im)+cA(4,km,jm,im)*cffm + diag
             diag = cA(7,k,j,im)+cA(7,km,j,im)*cffm+cA(7,k,jm,im)+cA(7,km,jm,im)*cffm + diag

             if (k2 == 1) then
                ! add the horizontal diagonal connections in the bottom level
                diag = cA(5,k,j,im)  + diag
                diag = cA(8,k,jm,im) + diag
             endif

             ! double that to account for symmetry of connections, we've now 40 terms
             diag = diag+diag


             sm = msk(j,i)+msk(jm,i)+msk(j,im)+msk(jm,im)
             !             cff = c3(sm)/4._8
             ! add the 8 self-interacting terms
             diag = cA(1,k,j,i) +cA(1,km,j,i)*cffm +cA(1,k,jm,i) +cA(1,km,jm,i)*cffm &
                  +cA(1,k,j,im)+cA(1,km,j,im)*cffm+cA(1,k,jm,im)+cA(1,km,jm,im)*cffm + diag

             ! here we go!
             cA2(1,k2,j2,i2) = cff*diag
          enddo
       enddo
    enddo

  end subroutine coarsen_matrix_3D_nearest_avg

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_3D_linear_avg_8(lev,cA,cA2,nx2,ny2,nz2) ! from lev to lev+1

    integer(kind=ip),intent(in):: lev
    integer(kind=ip):: nx2, ny2, nz2! on lev+1
    real(kind=rp), dimension(:,:,:,:), pointer :: cA,cA2

    !local
    integer(kind=ip):: nx, ny, nz
    integer(kind=ip):: l, k, j, i
    integer(kind=ip):: k2, j2, i2
    integer(kind=ip):: kx, jx, ix
    integer(kind=ip):: ky, jy, iy
    integer(kind=ip):: dk2, dj2, di2
    real(kind=rp)   :: cx, cy, cz, z1, z2
    real(kind=rp), dimension(4,4,4) :: x
    real(kind=rp), dimension(2,2,2) :: y
    real(kind=rp), dimension(-1:1,1:4) :: cI
    real(kind=rp) :: yy
    real(kind=rp),dimension(0:4)::c3
    real(kind=rp) :: cff
    integer(kind=ip):: sm

    integer(kind=1),dimension(:,:),pointer::msk

    nx = nx2*2
    ny = ny2*2
    nz = nz2*2

    data c3/0.,1.,0.5,0.3333333333333333333333333333,0.25/
      
    msk => grid(lev)%rmask

    ! weights for linear interpolation, cI has 2 entries 
    !
    ! cI(di,ix), di={-1,0,1}, ix={1,2,3,4}
    !
    ! cI(di,ix) is the weight to put on point x(ix) for coarse grid point
    ! located at di
    !
    cI(:,:)=0.

    cI(-1,1)=0.75_8
    cI(-1,2)=0.25_8

    cI(0,1)=0.25_8
    cI(0,2)=0.75_8
    cI(0,3)=0.75_8
    cI(0,4)=0.25_8

    cI(+1,3)=0.25_8
    cI(+1,4)=0.75_8

    ! i1,j1,k1 are index on the fine grid
    ! i2,j2,k2 are index on the coarse grid

    ! l direction

    ! di2,dj2,dk2 index shift on the coarse grid (corresponding to 'l')

    ! map relating 'l' to directions
    ! l=1 is the central point (di=dj=dk=0)
    !
    !
    ! map for di = 0
    ! +---+---+---+---+---+---+
    ! | o | o | o | o | o | o |
    ! +--10---+---2---+---5---+
    ! | o | x | x | x | x | o |   ^
    ! +---+---+---+---+---+---+   |
    ! | o | x | y | y | x | o |   |
    ! +--18---+---1---+--13---+  k axis
    ! | o | x | y | y | x | o |   |
    ! +---+---+---+---+---+---+   |
    ! | o | x | x | x | x | o |
    ! +--26---+---3---+--21---+
    ! | o | o | o | o | o | o |
    ! +---+---+---+---+---+---+
    !
    !      -- j axis --> 

    !    cA(01)(k,j,i) -> p(k+0,j+0,i+0)  ***
    !    cA(02)(k,j,i) -> p(k+1,j+0,i+0)  ***
    !    cA(03)(k,j,i) -> p(k-1,j+0,i+0)  *** ! cA3(k+1,j,i) == cA2(k,j,i) ?

    !    cA(04)(k,j,i) -> p(k+1,j+1,i+1)
    !    cA(05)(k,j,i) -> p(k+1,j+1,i+0)
    !    cA(06)(k,j,i) -> p(k+1,j+1,i-1)
    !    cA(07)(k,j,i) -> p(k+1,j+0,i+1)
    !    cA(08)(k,j,i) -> p(k+1,j+0,i-1)
    !    cA(09)(k,j,i) -> p(k+1,j-1,i+1)
    !    cA(10)(k,j,i) -> p(k+1,j-1,i+0)
    !    cA(11)(k,j,i) -> p(k+1,j-1,i-1)

    !    cA(12)(k,j,i) -> p(k+0,j+1,i+1)
    !    cA(13)(k,j,i) -> p(k+0,j+1,i+0)
    !    cA(14)(k,j,i) -> p(k+0,j+1,i-1)
    !    cA(15)(k,j,i) -> p(k+0,j+0,i+1)
    !    cA(16)(k,j,i) -> p(k+0,j+0,i-1)
    !    cA(17)(k,j,i) -> p(k+0,j-1,i+1)
    !    cA(18)(k,j,i) -> p(k+0,j-1,i+0)
    !    cA(19)(k,j,i) -> p(k+0,j-1,i-1)

    !    cA(20)(k,j,i) -> p(k-1,j+1,i+1)
    !    cA(21)(k,j,i) -> p(k-1,j+1,i+0)
    !    cA(22)(k,j,i) -> p(k-1,j+1,i-1)
    !    cA(23)(k,j,i) -> p(k-1,j+0,i+1)
    !    cA(24)(k,j,i) -> p(k-1,j+0,i-1)
    !    cA(25)(k,j,i) -> p(k-1,j-1,i+1)
    !    cA(26)(k,j,i) -> p(k-1,j-1,i+0)
    !    cA(27)(k,j,i) -> p(k-1,j-1,i-1)


    do i2 = 1,nx2
       do j2 = 1,ny2
!          if(grid(lev+1)%rmask(j2,i2)==1)then
          do k2 = 1,nz2
             !
             do l=1,27! loop over coarse points=each connection
                ! loc relates the index l to each direction
                ! direction takes value in [-1,0,+1]
                di2 = loc(l,3)
                dj2 = loc(l,2)
                dk2 = loc(l,1)
                ! compute cA2 only if this is a valid direction in z
                if ((k2+dk2>=1).and.(k2+dk2<=nz2)) then


                   ! step 1: INTERPOLATE COARSE POINT 'l' on x's
                   !
                   ! x is a 4x4x4 array of fine points surrounding (k2,j2,i2)
                   !
                   do ix=1,4                    
                      cx = cI(di2,ix)
                      i = (i2+di2)*2+ix-3
                      do jx=1,4
                         cy = cI(dj2,jx)  
                         j = (j2+dj2)*2+jx-3
                         do kx=1,4
                            cz = cI(dk2,kx)
                            if(k2==1)then
                               ! forces x to 0 because x(k) is below the bottom boundary
                               if ((dk2==0).and.(kx==1)) cz = 0._8 
                               if ((dk2==0).and.(kx==2)) cz = 1._8 ! nearest interpolation
                            endif
                            if(k2==nz2)then
                               ! forces x to 0 because x(k) is above the top boundary
                               if ((dk2==0).and.(kx==4)) cz = 0._8 
                               if ((dk2==0).and.(kx==3)) cz = 0.5_8 ! nearest interpolation
                            endif
                            x(kx,jx,ix) = cx*cy*cz*grid(lev)%rmask(j,i)
                         enddo
                      enddo
                   enddo

                   ! step 2: COMPUTE FINE GRID LAPLACIAN on y's
                   !
                   ! y is a 2x2x2 array of fine points surrounding (k2,j2,i2)
                   !
                   do iy=1,2        ! iy,jy,ky are coordinate of the "y" / the 2x2x2 cube
                      i = i2*2+iy-2 ! i,j,k    are absolute position on the fine grid
                      ix = iy+1     ! ix,jx,kx are coordinate of the "x" / the 4x4x4 cube
                      do jy=1,2
                         j = j2*2+jy-2
                         jx=jy+1
                         do ky=1,2
                            k = k2*2+ky-2
                            kx=ky+1
                            if (k==1) then ! lower level

                               yy = cA(1,k,j,i)*x(kx  ,jx  ,ix)                      &
                                                + cA(2,k+1,j  ,i)  *x(kx+1,jx  ,ix)  &
                  + cA(3,k,j,i)*x(kx+1,jx-1,ix)                                      &
                  + cA(4,k,j,i)*x(kx  ,jx-1,ix) + cA(4,k  ,j+1,i)  *x(kx  ,jx+1,ix)&
                                                + cA(5,k+1,j+1,i)  *x(kx+1,jx+1,ix)&
                  + cA(6,k,j,i)*x(kx+1,jx,ix-1)                                      &
                  + cA(7,k,j,i)*x(kx  ,jx,ix-1) + cA(7,k  ,j  ,i+1)*x(kx  ,jx  ,ix+1)&
                                                + cA(8,k+1,j  ,i+1)*x(kx+1,jx  ,ix+1)
                               if (cmatrix == 'real') then
          !- Exception for the redefinition of the coef for the bottom level
                                  yy = yy &
               + cA(5,k,j,i)*x(kx,jx+1,ix-1) + cA(5,k,j-1,i+1)*x(kx,jx-1,ix+1) &
               + cA(8,k,j,i)*x(kx,jx-1,ix-1) + cA(8,k,j+1,i+1)*x(kx,jx+1,ix+1)
                               endif


                            elseif (k==nz) then

                               yy = cA(1,k,j,i)*x(kx  ,jx  ,ix)                      &
                  + cA(2,k,j,i)*x(kx-1,jx  ,ix)                                      &
                                                + cA(3,k-1,j+1,i)  *x(kx-1,jx+1,ix)  &
                  + cA(4,k,j,i)*x(kx  ,jx-1,ix) + cA(4,k  ,j+1,i)  *x(kx  ,jx+1,ix)&
                  + cA(5,k,j,i)*x(kx-1,jx-1,ix)                                    &
                                                + cA(6,k-1,j  ,i+1)*x(kx-1,jx  ,ix+1)&
                  + cA(7,k,j,i)*x(kx  ,jx,ix-1) + cA(7,k  ,j  ,i+1)*x(kx  ,jx  ,ix+1)&
                  + cA(8,k,j,i)*x(kx-1,jx,ix-1)                                      

                            else

                               yy = cA(1,k,j,i)*x(kx  ,jx  ,ix)                      &
                  + cA(2,k,j,i)*x(kx-1,jx  ,ix) + cA(2,k+1,j  ,i)  *x(kx+1,jx  ,ix)  &
                  + cA(3,k,j,i)*x(kx+1,jx-1,ix) + cA(3,k-1,j+1,i)  *x(kx-1,jx+1,ix)  &
                  + cA(4,k,j,i)*x(kx  ,jx-1,ix) + cA(4,k  ,j+1,i)  *x(kx  ,jx+1,ix)&
                  + cA(5,k,j,i)*x(kx-1,jx-1,ix) + cA(5,k+1,j+1,i)  *x(kx+1,jx+1,ix)&
                  + cA(6,k,j,i)*x(kx+1,jx,ix-1) + cA(6,k-1,j  ,i+1)*x(kx-1,jx  ,ix+1)&
                  + cA(7,k,j,i)*x(kx  ,jx,ix-1) + cA(7,k  ,j  ,i+1)*x(kx  ,jx  ,ix+1)&
                  + cA(8,k,j,i)*x(kx-1,jx,ix-1) + cA(8,k+1,j  ,i+1)*x(kx+1,jx  ,ix+1)

                            endif
                            y(ky,jy,iy) = yy
                         enddo
                      enddo
                   enddo

                   ! step 3: APPLY THE RESTRICTION = average the 8 y's
                   !z1 = y(1,1,1)+y(1,1,2)+y(1,2,1)+y(1,2,2)
                   !z2 = y(2,1,1)+y(2,1,2)+y(2,2,1)+y(2,2,2)
                   i=2*i2-1
                   j=2*j2-1
                   k=2*k2-1
                   sm = msk(j,i)+msk(j,i+1)+msk(j+1,i)+msk(j+1,i+1)
                   cff = c3(sm)*0.5_8
          
                   z1 = y(1,1,1)*msk(j  ,i)+y(1,1,2)*msk(j  ,i+1)&
                       +y(1,2,1)*msk(j+1,i)+y(1,2,2)*msk(j+1,i+1)
                   z2 = y(2,1,1)*msk(j  ,i)+y(2,1,2)*msk(j  ,i+1)&
                       +y(2,2,1)*msk(j+1,i)+y(2,2,2)*msk(j+1,i+1)

                   if(k2<nz2)then
                      cA2(l,k2,j2,i2) = ( z1+z2 )*cff
                   else
                      cA2(l,k2,j2,i2) = ( z1+z2/2 )*cff
                   endif

                else ! the 'l' neighbour is outside the domain
                   ! the value of this coefficient should be
                   ! insensitive
                   cA2(l,k2,j2,i2) = 0. ! or better: NAN
                endif
             enddo ! end of l-loop: coarse point/direction
          enddo
!          else
!             cA2(:,:,j2,i2) = 0.
!          endif
       enddo
    enddo

  end subroutine coarsen_matrix_3D_linear_avg_8

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_3D_linear_avg_27(lev,cA,cA2,nx2,ny2,nz2) ! from lev to lev+1

    integer(kind=ip),intent(in):: lev
    integer(kind=ip):: nx2, ny2, nz2! on lev+1
    real(kind=rp), dimension(:,:,:,:), pointer :: cA,cA2

    !local
    integer(kind=ip):: nx, ny, nz
    integer(kind=ip):: m, l
    integer(kind=ip):: k2, j2, i2
    integer(kind=ip):: k1, j1, i1
    integer(kind=ip):: k, j, i
    integer(kind=ip):: kx, jx, ix
    integer(kind=ip):: ky, jy, iy
    integer(kind=ip):: dk, dj, di
    integer(kind=ip):: dk2, dj2, di2
    real(kind=rp)   :: cx, cy, cz, z1, z2
    real(kind=rp), dimension(4,4,4) :: x
    real(kind=rp), dimension(2,2,2) :: y
    real(kind=rp), dimension(-1:1,1:4) :: cI
    real(kind=rp) :: yy
    real(kind=rp),dimension(0:4)::c3
    real(kind=rp) :: cff
    integer(kind=ip):: sm

    integer(kind=1),dimension(:,:),pointer::msk

    nx = nx2*2
    ny = ny2*2
    nz = nz2*2

    data c3/0.,1.,0.5,0.3333333333333333333333333333,0.25/
      
    msk => grid(lev)%rmask


    ! weights for linear interpolation, cI has 2 entries 
    !
    ! cI(di,ix), di={-1,0,1}, ix={1,2,3,4}
    !
    ! cI(di,ix) is the weight to put on point x(ix) for coarse grid point
    ! located at di
    !
    cI(:,:)=0._8

    cI(-1,1)=0.75_8
    cI(-1,2)=0.25_8

    cI(0,1)=0.25_8
    cI(0,2)=0.75_8
    cI(0,3)=0.75_8
    cI(0,4)=0.25_8

    cI(+1,3)=0.25_8
    cI(+1,4)=0.75_8


    ! i1,j1,k1 are index on the fine grid
    ! i2,j2,k2 are index on the coarse grid

    ! l direction

    ! di2,dj2,dk2 index shift on the coarse grid (corresponding to 'l')

    ! map relating 'l' to directions
    ! l=1 is the central point (di=dj=dk=0)
    !
    !
    ! map for di = 0
    ! +---+---+---+---+---+---+
    ! | o | o | o | o | o | o |
    ! +--10---+---2---+---5---+
    ! | o | x | x | x | x | o |   ^
    ! +---+---+---+---+---+---+   |
    ! | o | x | y | y | x | o |   |
    ! +--18---+---1---+--13---+  k axis
    ! | o | x | y | y | x | o |   |
    ! +---+---+---+---+---+---+   |
    ! | o | x | x | x | x | o |
    ! +--26---+---3---+--21---+
    ! | o | o | o | o | o | o |
    ! +---+---+---+---+---+---+
    !
    !      -- j axis --> 

    !    cA(01)(k,j,i) -> p(k+0,j+0,i+0)  ***
    !    cA(02)(k,j,i) -> p(k+1,j+0,i+0)  ***
    !    cA(03)(k,j,i) -> p(k-1,j+0,i+0)  *** ! cA3(k+1,j,i) == cA2(k,j,i) ?

    !    cA(04)(k,j,i) -> p(k+1,j+1,i+1)
    !    cA(05)(k,j,i) -> p(k+1,j+1,i+0)
    !    cA(06)(k,j,i) -> p(k+1,j+1,i-1)
    !    cA(07)(k,j,i) -> p(k+1,j+0,i+1)
    !    cA(08)(k,j,i) -> p(k+1,j+0,i-1)
    !    cA(09)(k,j,i) -> p(k+1,j-1,i+1)
    !    cA(10)(k,j,i) -> p(k+1,j-1,i+0)
    !    cA(11)(k,j,i) -> p(k+1,j-1,i-1)

    !    cA(12)(k,j,i) -> p(k+0,j+1,i+1)
    !    cA(13)(k,j,i) -> p(k+0,j+1,i+0)
    !    cA(14)(k,j,i) -> p(k+0,j+1,i-1)
    !    cA(15)(k,j,i) -> p(k+0,j+0,i+1)
    !    cA(16)(k,j,i) -> p(k+0,j+0,i-1)
    !    cA(17)(k,j,i) -> p(k+0,j-1,i+1)
    !    cA(18)(k,j,i) -> p(k+0,j-1,i+0)
    !    cA(19)(k,j,i) -> p(k+0,j-1,i-1)

    !    cA(20)(k,j,i) -> p(k-1,j+1,i+1)
    !    cA(21)(k,j,i) -> p(k-1,j+1,i+0)
    !    cA(22)(k,j,i) -> p(k-1,j+1,i-1)
    !    cA(23)(k,j,i) -> p(k-1,j+0,i+1)
    !    cA(24)(k,j,i) -> p(k-1,j+0,i-1)
    !    cA(25)(k,j,i) -> p(k-1,j-1,i+1)
    !    cA(26)(k,j,i) -> p(k-1,j-1,i+0)
    !    cA(27)(k,j,i) -> p(k-1,j-1,i-1)



    do i2 = 1,nx2
       do j2 = 1,ny2
!          if(grid(lev+1)%rmask(j2,i2)==1)then
          do k2 = 1,nz2
             !
             do l=1,27! loop over coarse points=each connection
                ! loc relates the index l to each direction
                ! direction takes value in [-1,0,+1]
                di2 = loc(l,3)
                dj2 = loc(l,2)
                dk2 = loc(l,1)

                ! compute cA2 only if this is a valid direction in z
                if ((k2+dk2>=1).and.(k2+dk2<=nz2)) then


                   ! step 1: INTERPOLATE COARSE POINT 'l' on x's
                   !
                   ! x is a 4x4x4 array of fine points surrounding (k2,j2,i2)
                   !
                   do ix=1,4                    
                      cx = cI(di2,ix)
                      i = (i2+di2)*2+ix-3
                      do jx=1,4
                         cy = cI(dj2,jx)    
                         j = (j2+dj2)*2+jx-3
                         do kx=1,4
                            cz = cI(dk2,kx)
                            if(k2==1)then
                               ! forces x to 0 because x(k) is below the bottom boundary
                               if ((dk2==0).and.(kx==1)) cz = 0._8 
                               if ((dk2==0).and.(kx==2)) cz = 1._8 ! nearest interpolation
                            endif
                            if(k2==nz2)then
                               ! forces x to 0 because x(k) is above the top boundary
                               if ((dk2==0).and.(kx==4)) cz = 0._8 
                               if ((dk2==0).and.(kx==3)) cz = 0.5_8 ! nearest interpolation
                            endif
                            x(kx,jx,ix) = cx*cy*cz*grid(lev)%rmask(j,i)
                         enddo
                      enddo
                   enddo

                   ! step 2: COMPUTE FINE GRID LAPLACIAN on y's
                   !
                   ! y is a 2x2x2 array of fine points surrounding (k2,j2,i2)
                   !
                   do iy=1,2
                      i1 = i2*2+iy-2 ! absolute position on the fine grid
                      ix = iy+1 ! ix,jx,kx are coordinate of the "x" / the 4x4x4 cube
                      do jy=1,2
                         j1 = j2*2+jy-2
                         jx=jy+1
                         do ky=1,2
                            k1 = k2*2+ky-2
                            kx=ky+1
                            yy = 0._8
                            do m=1,27 ! assume that cA is already not symmetric
                               di=loc(m,3)
                               dj=loc(m,2)
                               dk=loc(m,1)
                               yy = yy + cA(m,k1,j1,i1)*x(kx+dk,jx+dj,ix+di)
                            enddo
                            y(ky,jy,iy) = yy
                         enddo
                      enddo
                   enddo

                   ! step 3: APPLY THE RESTRICTION = average the 8 y's
                   i=2*i2-1
                   j=2*j2-1
                   k=2*k2-1
                   sm = msk(j,i)+msk(j,i+1)+msk(j+1,i)+msk(j+1,i+1)
                   cff = c3(sm)*0.5_8
          
                   z1 = y(1,1,1)*msk(j  ,i)+y(1,1,2)*msk(j  ,i+1)&
                       +y(1,2,1)*msk(j+1,i)+y(1,2,2)*msk(j+1,i+1)
                   z2 = y(2,1,1)*msk(j  ,i)+y(2,1,2)*msk(j  ,i+1)&
                       +y(2,2,1)*msk(j+1,i)+y(2,2,2)*msk(j+1,i+1)
                   if(k2<nz2)then
                      cA2(l,k2,j2,i2) = ( z1+z2 )*cff
                   else
                      cA2(l,k2,j2,i2) = ( z1+z2/2 )*cff
                   endif

                else ! the 'l' neighbour is outside the domain
                   ! the value of this coefficient should be
                   ! insensitive
                   cA2(l,k2,j2,i2) = 0. ! or better: NAN
                endif
             enddo ! end of l-loop: coarse point/direction
          enddo
!          else
!             cA2(:,:,j2,i2) = 0.
!          endif
       enddo
    enddo


  end subroutine coarsen_matrix_3D_linear_avg_27

end module mg_define_matrix
