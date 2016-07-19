module mg_relax

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_define_matrix

  implicit none

    real(kind=rp), dimension(200) :: zrhs, zd, zud, zgam, zxc

contains

  !----------------------------------------
  subroutine relax(lev,nsweeps)
    integer(kind=ip), intent(in):: lev
    integer(kind=ip), intent(in):: nsweeps

!    real(kind=rp),dimension(:,:,:), pointer:: p
!    real(kind=rp),dimension(:,:,:), pointer:: b
!    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nd

!    p  => grid(lev)%p
!    b  => grid(lev)%b
!    cA => grid(lev)%cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
!    nd = size(cA(:,:,:,:),dim=1)

    if (grid(lev)%nz == 1) then

       call relax_2D_5(lev,grid(lev)%p,grid(lev)%b,grid(lev)%cA,nsweeps,nx,ny)

    else

       select case(trim(relax_method))

       case('Gauss-Seidel','GS')
          call relax_3D_8_GS(lev,grid(lev)%p,grid(lev)%b,grid(lev)%cA,nsweeps,nx,ny,nz)

       case('Red-Black','RB')
          call relax_3D_8_RB(lev,grid(lev)%p,grid(lev)%b,grid(lev)%cA,nsweeps,nx,ny,nz)

       case('Four-Color','FC')
          call relax_3D_8_FC(lev,grid(lev)%p,grid(lev)%b,grid(lev)%cA,nsweeps,nx,ny,nz)

       end select

    end if

  end subroutine relax

  !----------------------------------------
  subroutine relax_2D_5(lev,p,b,cA,nsweeps,nx,ny)

    integer(kind=ip)                         , intent(in)   :: lev
    real(kind=rp),dimension(1,0:ny+1,0:nx+1)  , intent(inout):: p
    real(kind=rp),dimension(1,0:ny+1,0:nx+1)  , intent(in)   :: b
    real(kind=rp),dimension(5,1,0:ny+1,0:nx+1), intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip)           :: i,j,k, it,rb
    integer(kind=ip)            :: ib,ie,jb,je,rbb,rbe,rbi
    real(kind=rp) :: z,gamma,g1,g2

    gamma = 1._8
    g1 = gamma
    g2 = 1._8 - gamma

    k=1

    do it = 1,nsweeps
       if (mod(it,1) == 0) then
          ib = 1 
          ie = nx
          jb = 1
          je = ny
       else
          ib = 0
          ie = nx+1
          jb = 0
          je = ny+1
       endif

       if (red_black) then
          rbb = 1
          rbe = 2
          rbi = 2
       else
          rbb = 0
          rbe = 0
          rbi = 1
       endif

       do rb = rbb,rbe
          do i = ib,ie
             do j = jb+mod(i+rb,rbi),je,rbi

                z =    b(k,j,i)                                           &
                     - cA(2,k,j,i)*p(k,j-1,i  ) - cA(2,k,j+1,i  )*p(k,j+1,  i)&
                     - cA(3,k,j,i)*p(k,j  ,i-1) - cA(3,k,j  ,i+1)*p(k,j  ,i+1)&
                     - cA(4,k,j,i)*p(k,j-1,i-1) - cA(4,k,j+1,i+1)*p(k,j+1,i+1)&
                     - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1)

                p(k,j,i) = z / cA(1,k,j,i)

             enddo
          enddo
       enddo

       call fill_halo(lev,p,nx,ny,1)

    enddo

  end subroutine relax_2D_5

 !----------------------------------------
  subroutine relax_3D_8_GS(lev,p,b,cA,nsweeps,nx,ny,nz)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(inout):: p
    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(in)   :: b
    real(kind=rp),dimension(8,nz,0:ny+1,0:nx+1), intent(in)   :: cA

    integer(kind=ip)            :: i,j,it

    call tic(lev,'relax_3D_8_GS')

    ! add a loop on smoothing
    do it = 1,nsweeps

       do i = 1, nx
          do j = 1, ny

             call relax_3D_8_heart(p,b,cA,i,j,nx,ny,nz)

          enddo ! j
       enddo    ! i

       call fill_halo(lev,p,nx,ny,nz)

    enddo  !it

    call toc(lev,'relax_3D_8_GS')

  end subroutine relax_3D_8_GS

 !----------------------------------------
  subroutine relax_3D_8_RB(lev,p,b,cA,nsweeps,nx,ny,nz)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(inout):: p
    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(in)   :: b
    real(kind=rp),dimension(8,nz,0:ny+1,0:nx+1), intent(in)   :: cA
!    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
!    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
!    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip) :: i,j,it
    integer(kind=ip) :: rb

    call tic(lev,'relax_3D_8_RB')

    ! add a loop on smoothing
    do it = 1,nsweeps

       do rb = 1, 2 ! Red black loop
          do i = 1, nx
             do j = 1+mod(i+rb,2),ny,2

                call relax_3D_8_heart(p,b,cA,i,j,nx,ny,nz)

             enddo ! j
          enddo    ! i

          call fill_halo(lev,p,nx,ny,nz)

       enddo       ! red-black

    enddo  !it

    call toc(lev,'relax_3D_8_RB')

  end subroutine relax_3D_8_RB

 !----------------------------------------
  subroutine relax_3D_8_FC(lev,p,b,cA,nsweeps,nx,ny,nz)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

!    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
!    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
!    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(inout):: p
    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(in)   :: b
    real(kind=rp),dimension(8,nz,0:ny+1,0:nx+1), intent(in)   :: cA

    integer(kind=ip)            :: i,j,it
    integer(kind=ip)            :: fc1,fc2

    call tic(lev,'relax_3D_8_FC')

    ! add a loop on smoothing
    do it = 1,nsweeps

       do fc1 = 1, 2 ! 
          do fc2 = 1, 2 ! 
             do i = 1 + mod(fc1-1,2), nx, 2
                do j = 1 + mod(fc2-1,2), ny, 2

                   call relax_3D_8_heart(p,b,cA,i,j,nx,ny,nz)

                enddo ! j
             enddo    ! i

             call fill_halo(lev,p,nx,ny,nz)

          enddo  ! fc2
       enddo  ! fc1

    enddo  !it

    call toc(lev,'relax_3D_8_FC')

  end subroutine relax_3D_8_FC

  !----------------------------------------
  subroutine relax_3D_8_heart(p,b,cA,i,j,nx,ny,nz)

    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(inout):: p
    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(in)   :: b
    real(kind=rp),dimension(8,nz,0:ny+1,0:nx+1), intent(in)   :: cA
!    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
!    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
!    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip)                          , intent(in)  :: i
    integer(kind=ip)                          , intent(in)  :: j
    integer(kind=ip)                          , intent(in)  :: nx,ny,nz

    !- Local -!
    integer(kind=ip) :: k
    real*8 :: bet
!    real(kind=rp), dimension(nz) :: zrhs, zd, zud, zgam

    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    k=1 !lower level
    zrhs(k) = b(k,j,i)                                              &
         - cA(3,k,j,i)*p(k+1,j-1,i)                                &
         - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
         - cA(5,k+1,j+1,i)*p(k+1,j+1,i) &
         - cA(6,k,j,i)*p(k+1,j,i-1)                                &
         - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
         - cA(8,k+1,j,i+1)*p(k+1,j,i+1) 

    if (cmatrix == 'real') then
       !- Exception for the redefinition of the coef for the bottom level
       zrhs(k) = zrhs(k) &
            - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1) &
            - cA(8,k,j,i)*p(k,j-1,i-1) - cA(8,k,j+1,i+1)*p(k,j+1,i+1)
    endif

    zd(k)   = cA(1,k,j,i)
    zud(k)  = cA(2,k+1,j,i)

    do k = 2,nz-1 !interior levels
       zrhs(k) = b(k,j,i) &
            - cA(3,k,j,i)*p(k+1,j-1,i) - cA(3,k-1,j+1,i)*p(k-1,j+1,i) &
            - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
            - cA(5,k,j,i)*p(k-1,j-1,i) - cA(5,k+1,j+1,i)*p(k+1,j+1,i) &
            - cA(6,k,j,i)*p(k+1,j,i-1) - cA(6,k-1,j,i+1)*p(k-1,j,i+1) &
            - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
            - cA(8,k,j,i)*p(k-1,j,i-1) - cA(8,k+1,j,i+1)*p(k+1,j,i+1) 
       zd(k)   = cA(1,k,j,i)
       zud(k)  = cA(2,k+1,j,i)
    enddo

    k=nz !upper level
    zrhs(k) = b(k,j,i)                                              &
         - cA(3,k-1,j+1,i)*p(k-1,j+1,i) &
         - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
         - cA(5,k,j,i)*p(k-1,j-1,i)                                &
         - cA(6,k-1,j,i+1)*p(k-1,j,i+1) &
         - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
         - cA(8,k,j,i)*p(k-1,j,i-1) 
    zd(k)   = cA(1,k,j,i)

    call tridiag(nz,zd,zud,zrhs,p(:,j,i)) !solve for vertical_coeff_matrix.p1d=rhs    

  end subroutine relax_3D_8_heart

  !----------------------------------------
  subroutine tridiag(l,d,dd,b,xc)
    !     Axc = b
    !     Solve tridiagonal system
    implicit none
    !     IMPORT/EXPORT
    integer                   ,intent(in)  :: l
    real(kind=rp),dimension(l),intent(in)  :: d,b
    real(kind=rp),dimension(l),intent(in)  :: dd
    real(kind=rp),dimension(l),intent(out) :: xc
    !     LOCAL
    integer                  :: k
    real(kind=rp),dimension(l):: gam
    real(kind=rp)             :: bet

    bet   = 1._8/d(1)
    xc(1) = b(1)*bet
    do k=2,l
       gam(k)= dd(k-1)*bet
       bet     = 1._8/(d(k)-dd(k-1)*gam(k))
       xc(k) = (b(k)-dd(k-1)*xc(k-1))*bet
    enddo
    do k=l-1,1,-1
       xc(k) = xc(k)-gam(k+1)*xc(k+1)
    enddo
    !    endif
  end subroutine tridiag

  !----------------------------------------
  subroutine compute_residual(lev,res)
    integer(kind=ip), intent(in) :: lev
    real(kind=rp)   , intent(out):: res

!    real(kind=rp),dimension(grid(lev)%nz,0:grid(lev)%ny+1,0:grid(lev)%nx+1)  :: p
!    real(kind=rp),dimension(:,:,:)  , pointer:: p
!    real(kind=rp),dimension(:,:,:)  , pointer:: b
!    real(kind=rp),dimension(:,:,:)  , pointer:: r
!    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nd
    real(kind=rp) ::resloc

!    p  => grid(lev)%p
!    b  => grid(lev)%b
!    r  => grid(lev)%r
!    cA => grid(lev)%cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
!    nd = size(cA(:,:,:,:),dim=1)

    if (grid(lev)%nz == 1) then

!       call compute_residual_2D_5(res,p,b,r,cA,nx,ny)

    else

!       call compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)
       call tic(lev,'residual_3D_8')
    
       call compute_residual_3D_8(res,grid(lev)%p,grid(lev)%b,grid(lev)%r,grid(lev)%cA,nx,ny,nz)
       call toc(lev,'residual_3D_8')
    end if

    call fill_halo(lev,grid(lev)%r,nx,ny,nz)

    if (lev >-1) then
       resloc=res
       call global_sum(lev,resloc,res)
       res = sqrt(res)
    else
       res = -999._8
    endif

  end subroutine compute_residual

  !----------------------------------------
  subroutine compute_residual_2D_5(res,p,b,r,cA,nx,ny)
    !TODO
    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip) :: i,j,k
    real(kind=rp)  :: z

    res = 0._8

    k=1

    do i = 1,nx
       do j = 1,ny

          z = b(k,j,i) - cA(1,k,j,i)*p(k,j,i)                           &
               - cA(2,k,j,i)*p(k,j-1,i  ) - cA(2,k,j+1,i  )*p(k,j+1,  i)&
               - cA(3,k,j,i)*p(k,j  ,i-1) - cA(3,k,j  ,i+1)*p(k,j  ,i+1)&
               - cA(4,k,j,i)*p(k,j-1,i-1) - cA(4,k,j+1,i+1)*p(k,j+1,i+1)&
               - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1)

          r(k,j,i) = z
          !          res = max(res,abs(r(k,j,i)))
          res = res+z*z

       enddo
    enddo

  end subroutine compute_residual_2D_5

  !----------------------------------------
  subroutine compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)
!  subroutine compute_residual_3D_8(res,lev)

    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(inout):: p
    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(in)   :: b
    real(kind=rp),dimension(nz,0:ny+1,0:nx+1)  , intent(inout)   :: r
    real(kind=rp),dimension(8,nz,0:ny+1,0:nx+1), intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz
!!$    integer(kind=ip)                        , intent(in)   :: lev

    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    integer(kind=ip)           :: i,j,k
!    integer(kind=ip)           :: nx, ny, nz


!    real(kind=rp),dimension(:,:,:)  , pointer:: p
!    real(kind=rp),dimension(:,:,:)  , pointer:: b
!    real(kind=rp),dimension(:,:,:)  , pointer:: r
!    real(kind=rp),dimension(:,:,:,:), pointer:: cA
!    nx = grid(lev)%nx
!    ny = grid(lev)%ny
!    nz = grid(lev)%nz

!    p  => grid(lev)%p
!    b  => grid(lev)%b
!    r  => grid(lev)%r
!    cA => grid(lev)%cA

    res = 0._8

    do i = 1,nx
       do j = 1,ny

             k=1 !lower level
             r(k,j,i) = b(k,j,i)                                           &
                  - cA(1,k,j,i)*p(k,j,i)                                   &
                  - cA(2,k+1,j,i)*p(k+1,j,i)                               &
                  - cA(3,k,j,i)*p(k+1,j-1,i)                               &
                  - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
                  - cA(5,k+1,j+1,i)*p(k+1,j+1,i)                           &
                  - cA(6,k,j,i)*p(k+1,j,i-1)                               &
                  - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
                  - cA(8,k+1,j,i+1)*p(k+1,j,i+1)

             if (cmatrix == 'real') then
             !- Exception for the redefinition of the coef for the bottom level
             r(k,j,i) = r(k,j,i) &
                  - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1) &
                  - cA(8,k,j,i)*p(k,j-1,i-1) - cA(8,k,j+1,i+1)*p(k,j+1,i+1)
             endif

             res = res+r(k,j,i)*r(k,j,i)

             do k = 2,nz-1 !interior levels
                r(k,j,i) = b(k,j,i)                                           &
                     - cA(1,k,j,i)*p(k,j,i)                                   &
                     - cA(2,k,j,i)*p(k-1,j,i)   - cA(2,k+1,j,i)*p(k+1,j,i)    &
                     - cA(3,k,j,i)*p(k+1,j-1,i) - cA(3,k-1,j+1,i)*p(k-1,j+1,i)&
                     - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
                     - cA(5,k,j,i)*p(k-1,j-1,i) - cA(5,k+1,j+1,i)*p(k+1,j+1,i)&
                     - cA(6,k,j,i)*p(k+1,j,i-1) - cA(6,k-1,j,i+1)*p(k-1,j,i+1)&
                     - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
                     - cA(8,k,j,i)*p(k-1,j,i-1) - cA(8,k+1,j,i+1)*p(k+1,j,i+1)

                res = res+r(k,j,i)*r(k,j,i)
             enddo

             k=nz !upper level
             r(k,j,i) = b(k,j,i)                                           &
                  - cA(1,k,j,i)*p(k,j,i)                                   &
                  - cA(2,k,j,i)*p(k-1,j,i)                                 &
                  - cA(3,k-1,j+1,i)*p(k-1,j+1,i)&
                  - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
                  - cA(5,k,j,i)*p(k-1,j-1,i)                               &
                  - cA(6,k-1,j,i+1)*p(k-1,j,i+1)&
                  - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
                  - cA(8,k,j,i)*p(k-1,j,i-1)

             res = res+r(k,j,i)*r(k,j,i)


       enddo
    enddo

  end subroutine compute_residual_3D_8

end module mg_relax
